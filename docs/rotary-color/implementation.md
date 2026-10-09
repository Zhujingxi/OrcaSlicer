# Rotary virtual filament angles

Each virtual filament selects an absolute A-axis angle of one physical spool/tool.
The printer owns rotary geometry and XY compensation.
There is no slicer-side center, coordinate transform, optical solver or collision
planner in this implementation.

Target: OrcaSlicer 2.5.0-dev. Firmware inspected:
[SnapmakerU1-Extended-Firmware](https://github.com/Zhujingxi/SnapmakerU1-Extended-Firmware/tree/rotary-table), branch `rotary-table`, commit
`2e875805fbb8a0a2c8934be900eed5855bbe3997`, patch
`overlays/firmware-extended/11-patch-klipper/patches/home/lava/klipper/09_rotary_table.patch`.
This firmware accepts `G1 A` and transforms XY internally.

## Use

1. Configure and reference the printer's rotary A axis in firmware.
2. Add the desired color slots with the filament **+** button and set their
   preview swatches. They are appearances of the same loaded spool.
3. Click **Rotary color printing...** below the **Filament / Project Filaments** list. Enable the
   mode, choose the material's physical filament slot and hardware tool, enter
   rotation speed and an absolute angle for each color, then press OK.
4. Paint or assign parts with the ordinary color tools and slice. Preview uses
   those color slots; export emits `G1 A<angle> F<degrees-per-minute>` on changes.
   No extra tool change or material purge is emitted for a virtual color change.

Angles are absolute, unwrapped device values: `47.5`, `-90` and `450` are distinct
valid requests. There is no shortest-path or wrap selection. The initial tool
selection and ordinary startup/shutdown scripts still run. Initial startup must
establish the desired A reference; custom scripts must preserve that reference
and return to the slicer's normal absolute XYZ mode. Retraction and Z hop use
the selected filament's existing settings. Rotation restores the preceding F.

The mode defaults to off, angles to zero and requested speed to 30 degrees/s.
The speed is a user preference, not a hardware capability or maximum. Project
settings persist in 3MF; legacy projects without these keys keep rotary mode off.
Removing the selected physical filament disables the mode. Turning the mode off
restores ordinary material/tool interpretation of the palette.

## Source integration

* Project settings in `PrintConfig` and `PresetBundle`: enable, physical filament,
  rotary speed, and per-color absolute angles. Disabled is the legacy behavior.
* Native filament sidebar: edit rotary settings and angles for the existing color
  slots; reuse normal painting and 3MF storage.
* `PrintApply`: inherit the selected physical filament's material settings and
  tool mapping while retaining color slots for segmentation and preview.
* `GCode`/`GCodeWriter`: virtual switches preserve extrusion state, emit `G1 A…`,
  and skip physical filament-change scripts and purge towers. Exported XY stays
  in the ordinary slicer coordinate frame.
* `GCodeProcessor`: a rotary color annotation changes preview color without
  recording a physical tool/material change.
* `tests/fff_print`: command generation, E continuity, canonical XY, cooling,
  virtual colors, physical-tool startup and disabled regressions.
* `tests/libslic3r`: material inheritance, validation and project persistence.

The initial mode uses one physical spool/tool for all color slots. Existing
filament swatches are user-selected preview colors, not an optical prediction.
Firmware configuration, A-axis reference, travel limits and compensated motion
are printer responsibilities. Machine preview and time estimation do not model
the compensated XY motion or rotation duration.

## Verification

Linux Release builds of `OrcaSlicer`, `fff_print_tests` and `libslic3r_tests` pass.
The application also passes a `--help` startup smoke check.

* FFF tests: `[RotaryColor],[GCodeWriter],[Cooling],[MultiFilament]`, randomized
  order: 31 cases, 14,455 assertions passed.
* Core tests: `[RotaryColor],[3mf]`, randomized order: 12 cases, 116 assertions passed.

Coverage includes identical swatches with different angles, an unpainted physical
parent, a nonzero physical tool, exported-file preview reload, and 3MF persistence.
The GUI was launched with the Snapmaker U1 profile. The visible rotary button
and its settings dialog were opened and checked. Physical printing has not been
verified.
