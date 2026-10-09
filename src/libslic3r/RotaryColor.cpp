#include "RotaryColor.hpp"
#include "Preset.hpp"
#include "Exception.hpp"

#include <algorithm>
#include <cmath>
#include <memory>

namespace Slic3r {

const t_config_option_keys &rotary_color_config_keys()
{
    static const t_config_option_keys keys{"rotary_color_enabled", "rotary_color_filament", "rotary_color_speed", "filament_rotary_angle"};
    return keys;
}

std::string validate_rotary_color_config(const PrintConfig &config)
{
    if (!config.rotary_color_enabled.value)
        return {};
    const size_t count = config.filament_diameter.size();
    if (config.rotary_color_filament.value < 1 || size_t(config.rotary_color_filament.value) > count)
        return "The physical filament for rotary colors is missing.";
    const int tool = config.filament_map.get_at(size_t(config.rotary_color_filament.value - 1));
    if (tool < 1 || size_t(tool) > config.nozzle_diameter.size())
        return "The physical tool for rotary colors is missing.";
    if (config.filament_rotary_angle.size() != count)
        return "Each rotary filament color must have an explicit A-axis angle.";
    if (!std::isfinite(config.rotary_color_speed.value) || config.rotary_color_speed.value <= 0. ||
        !std::isfinite(config.rotary_color_speed.value * 60.))
        return "Rotary speed must be finite and greater than zero.";
    for (double angle : config.filament_rotary_angle.values)
        if (!std::isfinite(angle))
            return "Rotary filament angles must be finite.";
    if (std::any_of(config.filament_is_mixed.values.begin(), config.filament_is_mixed.values.end(), [](auto v) { return v != 0; }))
        return "Rotary colors use one physical filament; remove mixed-material slots before enabling rotary colors.";
    return {};
}

void apply_rotary_color_config(DynamicPrintConfig &config)
{
    const auto *enabled = config.option<ConfigOptionBool>("rotary_color_enabled");
    if (!enabled || !enabled->value)
        return;
    const auto *colors = config.option<ConfigOptionStrings>("filament_colour");
    const auto *parent = config.option<ConfigOptionInt>("rotary_color_filament");
    if (!colors || !parent || parent->value < 1 || size_t(parent->value) > colors->size())
        throw SlicingError("The physical filament for rotary colors is missing.");
    const size_t count = colors->size(), source = size_t(parent->value - 1);
    const auto *mapping = config.option<ConfigOptionInts>("filament_map");
    const auto *nozzles = config.option<ConfigOptionFloats>("nozzle_diameter");
    if (!mapping || mapping->empty() || !nozzles || mapping->get_at(source) < 1 ||
        size_t(mapping->get_at(source)) > nozzles->size())
        throw SlicingError("The physical tool for rotary colors is missing.");
    const auto *angles = config.option<ConfigOptionFloats>("filament_rotary_angle");
    const auto *speed = config.option<ConfigOptionFloat>("rotary_color_speed");
    if (!angles || angles->size() != count ||
        std::any_of(angles->values.begin(), angles->values.end(), [](double angle) { return !std::isfinite(angle); }))
        throw SlicingError("Each rotary color must have an explicit finite A-axis angle.");
    if (!speed || !std::isfinite(speed->value * 60.) || speed->value <= 0. || std::round(speed->value * 60. * 1e9) == 0.)
        throw SlicingError("Rotary speed must be positive and representable in G-code.");
    // Unexpanded multi-nozzle presets concatenate variant columns. Keep every
    // column of the physical parent's preset, with fresh per-slot ownership.
    const auto *self_index = config.option<ConfigOptionInts>("filament_self_index");
    std::vector<size_t> parent_columns;
    if (self_index)
        for (size_t i = 0; i < self_index->size(); ++i)
            if (self_index->values[i] == parent->value)
                parent_columns.push_back(i);

    for (const auto &key : Preset::filament_options()) {
        auto *destination = dynamic_cast<ConfigOptionVectorBase *>(config.option(key));
        if (!destination || destination->empty())
            continue;
        std::unique_ptr<ConfigOption> original(destination->clone());
        const bool variants = filament_options_with_variant.count(key) != 0 && !parent_columns.empty() &&
                              self_index && destination->size() == self_index->size();
        const size_t width = variants ? parent_columns.size() : 1;
        destination->resize(count * width);
        for (size_t i = 0; i < count; ++i)
            for (size_t j = 0; j < width; ++j)
                destination->set_at(original.get(), i * width + j, variants ? parent_columns[j] : source);
    }
    if (self_index && !parent_columns.empty()) {
        auto &indices = config.option<ConfigOptionInts>("filament_self_index")->values;
        indices.clear();
        for (size_t i = 0; i < count; ++i)
            indices.insert(indices.end(), parent_columns.size(), int(i + 1));
    }
    for (const char *key : {"filament_map", "filament_nozzle_map", "filament_volume_map"}) {
        if (auto *mapping = config.option<ConfigOptionInts>(key)) {
            const int value = mapping->get_at(source);
            mapping->values.assign(count, value);
        }
    }
    // No material changes occur between these color slots.
    config.set_key_value("enable_prime_tower", new ConfigOptionBool(false));
    config.set_key_value("filament_map_mode", new ConfigOptionEnum<FilamentMapMode>(fmmManual));
    config.set_key_value("enable_filament_dynamic_map", new ConfigOptionBool(false));
    config.set_key_value("enable_pre_heating", new ConfigOptionBool(false));
    config.set_key_value("ooze_prevention", new ConfigOptionBool(false));
}

} // namespace Slic3r
