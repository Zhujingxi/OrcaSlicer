#pragma once

#include "PrintConfig.hpp"

namespace Slic3r {

const t_config_option_keys &rotary_color_config_keys();
// Only material/tool aliasing belongs here. Firmware owns the A-axis geometry.
std::string validate_rotary_color_config(const PrintConfig &config);
void apply_rotary_color_config(DynamicPrintConfig &config);

} // namespace Slic3r
