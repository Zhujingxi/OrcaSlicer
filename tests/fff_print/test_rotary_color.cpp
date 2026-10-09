#include <catch2/catch_all.hpp>

#include "libslic3r/GCodeWriter.hpp"
#include "libslic3r/GCode.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/RotaryColor.hpp"
#include "test_helpers.hpp"
#include "test_utils.hpp"

#include <fstream>
#include <set>

using namespace Slic3r;
using namespace Slic3r::Test;
using Catch::Matchers::WithinAbs;

static DynamicPrintConfig rotary_config()
{
    return multifilament_config(3, {
        {"rotary_color_enabled", "1"}, {"rotary_color_filament", "1"},
        {"filament_rotary_angle", "0,47.5,-450"}, {"rotary_color_speed", "23"},
        {"filament_map", "1,1,1"}, {"gcode_flavor", "klipper"},
        {"nozzle_hrc", "55"},
        {"layer_height", "0.2"}, {"initial_layer_print_height", "0.2"},
        {"enable_prime_tower", "1"}, {"machine_start_gcode", "G90\n"},
        {"machine_end_gcode", ""}, {"change_filament_gcode", ";PHYSICAL_SWAP\n"},
        {"filament_start_gcode", ";SPOOL_START\n"}, {"filament_end_gcode", ";SPOOL_END\n"},
        {"enable_arc_fitting", "0"}, {"brim_type", "no_brim"}, {"skirt_loops", "0"}
    });
}

TEST_CASE("Virtual angle switches preserve physical extrusion and retraction state", "[RotaryColor]")
{
    const bool relative = GENERATE(false, true);
    const bool shared = GENERATE(false, true);
    auto config = rotary_config();
    config.set_key_value("use_relative_e_distances", new ConfigOptionBool(relative));
    config.set_key_value("single_extruder_multi_material", new ConfigOptionBool(shared));
    apply_rotary_color_config(config);
    PrintConfig print_config;
    print_config.apply(config, true);
    GCodeWriter writer;
    writer.apply_print_config(print_config);
    writer.set_extruders({0, 1, 2});
    writer.set_extruder(0);
    writer.travel_to_xyz({12., 34., 2.});
    writer.extrude_to_xyz({13., 34., 2.}, 5.);
    writer.filament()->retract(0.8, 0.15);
    const double e_before = writer.filament()->E();
    const auto position = writer.get_position();
    const std::string change = writer.set_extruder(1);
    REQUIRE(change.find("G1 A47.5") != std::string::npos);
    REQUIRE(change.find("F1380") != std::string::npos);
    REQUIRE(change.find("\nT") == std::string::npos);
    REQUIRE(change.find("G92") == std::string::npos);
    REQUIRE_THAT(writer.filament()->E(), WithinAbs(e_before, 1e-9));
    REQUIRE_THAT(writer.filament()->unretract(), WithinAbs(0.95, 1e-9));
    REQUIRE_THAT((writer.get_position() - position).norm(), WithinAbs(0., 1e-9));
    REQUIRE(writer.set_extruder(1).empty());
    REQUIRE(writer.set_extruder(2).find("G1 A-450") != std::string::npos);
    REQUIRE(writer.set_extruder(0).find("G1 A0") != std::string::npos);
}

TEST_CASE("Rotary moves preserve canonical XY and restore the linear feed", "[RotaryColor]")
{
    auto config = rotary_config();
    apply_rotary_color_config(config);
    PrintConfig print_config;
    print_config.apply(config, true);
    GCodeWriter writer;
    writer.apply_print_config(print_config);
    writer.set_extruders({0, 1});
    writer.set_extruder(0);
    writer.set_speed(1234.);
    const auto change = writer.set_extruder(1);
    REQUIRE(change.find("G1 F1234") != std::string::npos);
    REQUIRE(change.find(" X") == std::string::npos);
    REQUIRE(change.find(" Y") == std::string::npos);
    REQUIRE(writer.travel_to_xy({17.25, 63.75}).find("X17.25 Y63.75") != std::string::npos);
}

TEST_CASE("Equal angle colors change their preview identity without another rotation", "[RotaryColor]")
{
    auto config = rotary_config();
    config.set_deserialize_strict("filament_rotary_angle", "31,31,31");
    apply_rotary_color_config(config);
    PrintConfig print_config;
    print_config.apply(config, true);
    GCodeWriter writer;
    writer.apply_print_config(print_config);
    writer.set_extruders({0, 1});
    writer.set_extruder(0);
    const auto change = writer.set_extruder(1);
    REQUIRE(change.find(";ROTARY_COLOR:1") != std::string::npos);
    REQUIRE(change.find("G1 A") == std::string::npos);
    REQUIRE(change.find("T1") == std::string::npos);
}

TEST_CASE("Sliced virtual colors rotate without filament swaps or a purge tower", "[RotaryColor]")
{
    auto config = rotary_config();
    if (GENERATE(false, true))
        config.set_deserialize_strict("filament_colour", "#FF0000;#FF0000;#FF0000");
    const auto gcode = slice_with_object_overrides({cube(2.), cube(2.), cube(2.)}, config,
        {{{"extruder", "1"}}, {{"extruder", "2"}}, {{"extruder", "3"}}});
    REQUIRE(gcode.find("G1 A47.5") != std::string::npos);
    REQUIRE(gcode.find("G1 A-450") != std::string::npos);
    GCodeReader reader;
    int selections = 0;
    reader.parse_buffer(gcode, [&](GCodeReader &, const GCodeReader::GCodeLine &line) {
        if (line.cmd_is("T0"))
            ++selections;
        REQUIRE_FALSE(line.cmd_is("T1"));
        REQUIRE_FALSE(line.cmd_is("T2"));
        REQUIRE(line.comment() != "PHYSICAL_SWAP");
    });
    REQUIRE(selections == 1);
    REQUIRE(gcode.find(";TYPE:Prime tower") == std::string::npos);
}

TEST_CASE("Rotary preview colors do not count as physical filament changes", "[RotaryColor]")
{
    const bool full_export = GENERATE(false, true);
    auto dynamic = rotary_config();
    apply_rotary_color_config(dynamic);
    FullPrintConfig config;
    config.apply(dynamic, true);
    ScopedTemporaryFile file(".gcode");
    {
        std::ofstream out(file.string());
        if (full_export)
            out << slice_with_object_overrides({cube(2.), cube(2.), cube(2.)}, dynamic,
                {{{"extruder", "1"}}, {{"extruder", "2"}}, {{"extruder", "3"}}});
        else
            out << "G90\nM83\nT0\n;ROTARY_COLOR:0\nG1 X10 Y10 Z0.2 F600\nG1 X20 E1\n"
                   ";ROTARY_COLOR:1\nG1 A47.5 F1380\nG1 F600\nG1 X30 E1\n"
                   ";ROTARY_COLOR:2\nG1 A-450 F1380\nG1 F600\nG1 X40 E1\n";
    }
    GCodeProcessor processor;
    processor.apply_config(config);
    processor.process_file(file.string());
    const auto &result = processor.get_result();
    std::set<unsigned> colors;
    for (const auto &move : result.moves)
        if (move.type == EMoveType::Extrude)
            colors.insert(move.extruder_id);
    REQUIRE(colors == std::set<unsigned>{0, 1, 2});
    REQUIRE(result.print_statistics.total_filament_changes == 0);
}

TEST_CASE("Disabled rotary settings retain ordinary tool selection", "[RotaryColor][Regression]")
{
    auto config = rotary_config();
    config.set_key_value("rotary_color_enabled", new ConfigOptionBool(false));
    const auto original = config;
    apply_rotary_color_config(config);
    REQUIRE(config == original);
    PrintConfig print_config;
    print_config.apply(config, true);
    GCodeWriter writer;
    writer.apply_print_config(print_config);
    writer.set_extruders({0, 1});
    writer.set_extruder(0);
    const auto change = writer.set_extruder(1);
    REQUIRE(change.find("T1") != std::string::npos);
    REQUIRE(change.find("G1 A") == std::string::npos);
}

TEST_CASE("Cooling retains A moves whose angular feed equals the preceding feed", "[RotaryColor][Cooling]")
{
    auto config = rotary_config();
    apply_rotary_color_config(config);
    PrintConfig print_config;
    print_config.apply(config, true);
    GCode generator;
    generator.apply_print_config(print_config);
    generator.writer().set_extruders({0, 1, 2});
    CoolingBuffer cooling(generator);
    const std::string output = cooling.process_layer(
        "G1 X10 Y10 F1380\nG1 A47.5 F1380\nG1 F1380\nG1 A-450 F1380\n", 0, true);
    REQUIRE(output.find("G1 A47.5 F1380") != std::string::npos);
    REQUIRE(output.find("G1 A-450 F1380") != std::string::npos);
}

TEST_CASE("Rotary startup selects only the configured physical tool", "[RotaryColor]")
{
    auto config = rotary_config();
    config.set_deserialize_strict({
        {"nozzle_diameter", "0.4,0.4,0.4"}, {"filament_map", "3,3,3"},
        {"rotary_color_filament", "3"},
        {"machine_start_gcode", ";START_TOOL={initial_tool}\n;USED={is_extruder_used[0]},{is_extruder_used[1]},{is_extruder_used[2]}\n"}
    });
    const auto output = slice_with_object_overrides({cube(2.), cube(2.)}, config,
        {{{"extruder", "1"}}, {{"extruder", "2"}}});
    REQUIRE(output.find(";START_TOOL=2") != std::string::npos);
    INFO(output.substr(output.find(";START_TOOL="), 100));
    REQUIRE(output.find(";USED=false,false,true") != std::string::npos);
    GCodeReader reader;
    unsigned selections = 0;
    reader.parse_buffer(output, [&](GCodeReader &, const GCodeReader::GCodeLine &line) {
        REQUIRE_FALSE(line.cmd_is("T0"));
        REQUIRE_FALSE(line.cmd_is("T1"));
        if (line.cmd_is("T2"))
            ++selections;
    });
    REQUIRE(selections == 1);
}
