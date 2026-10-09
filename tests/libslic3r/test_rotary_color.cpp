#include <catch2/catch_all.hpp>

#include "libslic3r/RotaryColor.hpp"
#include "libslic3r/PresetBundle.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Format/bbs_3mf.hpp"
#include "libslic3r/Semver.hpp"
#include "test_utils.hpp"

#include <limits>

using namespace Slic3r;
using Catch::Matchers::WithinAbs;

static DynamicPrintConfig project_config()
{
    auto config = DynamicPrintConfig::full_print_config();
    config.set_num_filaments(3);
    config.set_deserialize_strict({
        {"rotary_color_enabled", "1"}, {"rotary_color_filament", "2"},
        {"rotary_color_speed", "23.5"}, {"filament_rotary_angle", "-450,0,47.5"},
        {"filament_colour", "#FF0000;#00FF00;#0000FF"},
        {"filament_diameter", "1.7,1.75,1.8"}, {"filament_map", "1,3,2"},
        {"nozzle_diameter", "0.4,0.4,0.4"},
        {"nozzle_temperature", "190,215,230"}, {"filament_flow_ratio", "0.9,0.98,1.1"}
    });
    return config;
}

TEST_CASE("Rotary appearances inherit one material and tool while retaining their colors", "[RotaryColor]")
{
    auto config = project_config();
    const auto colors = config.option<ConfigOptionStrings>("filament_colour")->values;
    apply_rotary_color_config(config);
    for (double diameter : config.option<ConfigOptionFloats>("filament_diameter")->values)
        REQUIRE_THAT(diameter, WithinAbs(1.75, 1e-9));
    REQUIRE(config.option<ConfigOptionInts>("nozzle_temperature")->values == std::vector<int>{215, 215, 215});
    REQUIRE(config.option<ConfigOptionInts>("filament_map")->values == std::vector<int>{3, 3, 3});
    REQUIRE(config.option<ConfigOptionStrings>("filament_colour")->values == colors);
    REQUIRE_FALSE(config.opt_bool("enable_prime_tower"));
    const auto once = config;
    apply_rotary_color_config(config);
    REQUIRE(config == once);
}

TEST_CASE("Rotary settings reject missing angles and nonfinite motion parameters", "[RotaryColor]")
{
    PrintConfig config;
    config.apply(project_config(), true);
    REQUIRE(validate_rotary_color_config(config).empty());
    SECTION("Missing parent") {
        config.rotary_color_filament.value = 4;
        REQUIRE_FALSE(validate_rotary_color_config(config).empty());
    }
    SECTION("Missing angle") {
        config.filament_rotary_angle.values.pop_back();
        REQUIRE_FALSE(validate_rotary_color_config(config).empty());
    }
    SECTION("Invalid speed") {
        config.rotary_color_speed.value = GENERATE(0., -1., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN());
        REQUIRE_FALSE(validate_rotary_color_config(config).empty());
    }
    SECTION("Invalid angle") {
        config.filament_rotary_angle.values[1] = std::numeric_limits<double>::quiet_NaN();
        REQUIRE_FALSE(validate_rotary_color_config(config).empty());
    }
}

TEST_CASE("Removing a rotary color keeps other angles attached to their slots", "[RotaryColor]")
{
    PresetBundle bundle;
    bundle.set_num_filaments(3);
    bundle.project_config.apply_only(project_config(), rotary_color_config_keys());
    bundle.update_num_filaments(0);
    REQUIRE(bundle.project_config.opt_int("rotary_color_filament") == 1);
    REQUIRE_THAT(bundle.project_config.option<ConfigOptionFloats>("filament_rotary_angle")->values[1], WithinAbs(47.5, 1e-9));
    REQUIRE(bundle.project_config.opt_bool("rotary_color_enabled"));
    bundle.update_num_filaments(0);
    REQUIRE_FALSE(bundle.project_config.opt_bool("rotary_color_enabled"));
}

TEST_CASE("Rotary project settings survive a native 3MF round trip", "[RotaryColor][3mf]")
{
    auto config = project_config();
    Model model;
    auto *object = model.add_object();
    object->add_volume(make_cube(2., 2., 2.));
    model.add_default_instances();
    ScopedTemporaryDir backup("orca_rotary");
    model.set_backup_path(backup.string());
    ScopedTemporaryFile file(".3mf");
    const auto path = file.string();
    PlateData plate;
    plate.plate_index = 0;
    StoreParams params;
    params.path = path.c_str();
    params.model = &model;
    params.config = &config;
    params.plate_data_list.push_back(&plate);
    params.strategy = SaveStrategy::Zip64 | SaveStrategy::Silence;
    REQUIRE(store_bbs_3mf(params));
    Model loaded_model;
    ScopedTemporaryDir loaded_backup("orca_rotary_loaded");
    loaded_model.set_backup_path(loaded_backup.string());
    DynamicPrintConfig loaded;
    ConfigSubstitutionContext substitutions{ForwardCompatibilitySubstitutionRule::Enable};
    PlateDataPtrs plates;
    std::vector<Preset *> presets;
    bool bbl = false, orca = false;
    Semver version;
    REQUIRE(load_bbs_3mf(path.c_str(), &loaded, &substitutions, &loaded_model, &plates, &presets,
                         &bbl, &orca, &version, nullptr, LoadStrategy::LoadModel | LoadStrategy::LoadConfig));
    for (auto *item : plates)
        delete item;
    for (auto *preset : presets)
        delete preset;
    for (const auto &key : rotary_color_config_keys()) {
        REQUIRE(loaded.has(key));
        REQUIRE(*loaded.option(key) == *config.option(key));
    }
}
