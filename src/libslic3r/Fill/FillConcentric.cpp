#include "../ClipperUtils.hpp"
#include "../ExPolygon.hpp"
#include "../Surface.hpp"
#include "../VariableWidth.hpp"
#include "Arachne/WallToolPaths.hpp"

#include "FillConcentric.hpp"
#include <libslic3r/ShortestPath.hpp>
#include "../Geometry.hpp"

namespace Slic3r {

// Gap angle is now read from params.config->concentric_gap_angle

// Helper to find the closest intersection of a ray from center with the loop segments
// The ray starts at center and extends in the direction of the angle.
static bool find_intersection_aligned_with_angle(
    const Points& points, 
    bool is_closed_polygon, // true if implicit close (Polygon), false if explicit (Polyline with start==end)
    const Point& center, 
    double angle, 
    size_t& out_idx, 
    Point& out_point
) {
    Vec2d dir(std::cos(angle), std::sin(angle));
    Vec2d c = center.cast<double>();
    
    double min_dist2 = std::numeric_limits<double>::max();
    bool found = false;
    
    size_t count = points.size();
    // For closed polygons, there's an implicit segment from last to first.
    // For explicitly closed polylines (start==end), skip the last point.
    size_t segments = is_closed_polygon ? count : (count > 1 ? count - 1 : 0);
    
    for (size_t i = 0; i < segments; ++i) {
        Vec2d p1 = points[i].cast<double>();
        Vec2d p2 = points[(i + 1) % count].cast<double>();
        
        Vec2d seg_dir = p2 - p1;
        
        // Solve for intersection: c + t * dir = p1 + s * seg_dir
        // This is a 2x2 linear system.
        double denom = dir.x() * seg_dir.y() - dir.y() * seg_dir.x();
        if (std::abs(denom) < EPSILON) 
            continue; // Parallel or collinear
        
        Vec2d diff = p1 - c;
        double t = (diff.x() * seg_dir.y() - diff.y() * seg_dir.x()) / denom;
        double s = (diff.x() * dir.y() - diff.y() * dir.x()) / denom;
        
        // t >= 0: intersection is in the positive direction of the ray
        // 0 <= s <= 1: intersection is within the segment
        if (t >= 0 && s >= 0 && s <= 1.0) {
            Vec2d intersection = c + t * dir;
            double d2 = t * t; // Distance squared along ray (proportional)
            if (d2 < min_dist2) {
                min_dist2 = d2;
                out_idx = i;
                out_point = Point(coord_t(intersection.x()), coord_t(intersection.y()));
                found = true;
            }
        }
    }
    return found;
}


void FillConcentric::_fill_surface_single(
    const FillParams                &params, 
    unsigned int                     thickness_layers,
    const std::pair<float, Point>   &direction, 
    ExPolygon                        expolygon,
    Polylines                       &polylines_out)
{
    // no rotation is supported for this infill pattern
    BoundingBox bounding_box = expolygon.contour.bounding_box();
    
    coord_t min_spacing = scale_(this->spacing);
    coord_t distance = coord_t(min_spacing / params.density);
    
    if (params.density > 0.9999f && !params.dont_adjust) {
        distance = this->_adjust_solid_spacing(bounding_box.size()(0), distance);
        this->spacing = unscale<double>(distance);
    }

    Polygons   loops = to_polygons(expolygon);
    ExPolygons last { std::move(expolygon) };
    while (! last.empty()) {
        last = offset2_ex(last, -(distance + min_spacing/2), +min_spacing/2);
        append(loops, to_polygons(last));
    }

    // generate paths from the outermost to the innermost, to avoid
    // adhesion problems of the first central tiny loops
    loops = union_pt_chained_outside_in(loops);
    
    // Reverse loops to print from center (innermost) to outside (outermost)
    std::reverse(loops.begin(), loops.end());
    
    // Calculate center of bounding box for gap alignment
    Point center = bounding_box.center();
    
    // Convert config angle from degrees to radians and normalize to [-PI, PI]
    double gap_angle_degrees = params.config ? params.config->concentric_gap_angle.value : 90.0;
    double target_angle = gap_angle_degrees * M_PI / 180.0;
    while (target_angle > M_PI) target_angle -= 2.0 * M_PI;
    while (target_angle <= -M_PI) target_angle += 2.0 * M_PI;

    // Calculate the average centroid of all loops to find a stable "consensus center"
    Point average_center = bounding_box.center(); // Default to bbox center if no loops
    if (!loops.empty()) {
        Vec2d center_sum(0, 0);
        size_t count = 0;
        for (const Polygon &loop : loops) {
            center_sum += loop.centroid().cast<double>();
            count++;
        }
        if (count > 0) {
            Vec2d avg = center_sum / double(count);
            average_center = Point(avg.x(), avg.y());
        }
    }
    
    // Split paths at points aligned with the radial line from center at the specified angle
    // Split paths at points aligned with the radial line from center at the specified angle
    size_t iPathFirst = polylines_out.size();
    for (Polygon &loop : loops) {
        // Use the average center of ALL loops. 
        // This ensures straight gaps (fixed center) while being centrally balanced.
        Point center = average_center;
        
        size_t split_idx = 0;
        Point intersection;
        bool found = find_intersection_aligned_with_angle(loop.points, true, center, target_angle, split_idx, intersection);

        if (found) {
            // Insert the intersection point and split there
            loop.points.insert(loop.points.begin() + split_idx + 1, intersection);
            polylines_out.emplace_back(loop.split_at_index(split_idx + 1));
        } else {
            // Fallback to nearest vertex
            size_t best_idx = 0;
            double best_diff = std::numeric_limits<double>::max();
            
            for (size_t i = 0; i < loop.points.size(); ++i) {
                double dx = loop.points[i].x() - center.x();
                double dy = loop.points[i].y() - center.y();
                
                double angle = std::atan2(dy, dx);
                double diff = std::abs(angle - target_angle);
                if (diff > M_PI) diff = 2.0 * M_PI - diff;
                
                if (diff < best_diff) {
                    best_diff = diff;
                    best_idx = i;
                }
            }
            polylines_out.emplace_back(loop.split_at_index(best_idx));
        }
    }

    // clip the paths to prevent the extruder from getting exactly on the first point of the loop
    // Keep valid paths only.
    size_t j = iPathFirst;
    for (size_t i = iPathFirst; i < polylines_out.size(); ++ i) {
        polylines_out[i].clip_end(this->loop_clipping);
        if (polylines_out[i].is_valid()) {
            if (j < i)
                polylines_out[j] = std::move(polylines_out[i]);
            ++ j;
        }
    }
    if (j < polylines_out.size())
        polylines_out.erase(polylines_out.begin() + j, polylines_out.end());
    //TODO: return ExtrusionLoop objects to get better chained paths,
    // otherwise the outermost loop starts at the closest point to (0, 0).
    // We want the loops to be split inside the G-code generator to get optimum path planning.
}

void FillConcentric::_fill_surface_single(const FillParams& params,
    unsigned int                   thickness_layers,
    const std::pair<float, Point>& direction,
    ExPolygon                      expolygon,
    ThickPolylines& thick_polylines_out)
{
    assert(params.use_arachne);
    assert(this->print_config != nullptr && this->print_object_config != nullptr);

    // no rotation is supported for this infill pattern
    Point   bbox_size = expolygon.contour.bounding_box().size();
    coord_t min_spacing = scaled<coord_t>(this->spacing);

    if (params.density > 0.9999f && !params.dont_adjust) {
        coord_t                loops_count = std::max(bbox_size.x(), bbox_size.y()) / min_spacing + 1;
        Polygons               polygons = offset(expolygon, float(min_spacing) / 2.f);

        double min_nozzle_diameter = *std::min_element(print_config->nozzle_diameter.values.begin(), print_config->nozzle_diameter.values.end());
        Arachne::WallToolPathsParams input_params;
        input_params.min_bead_width = 0.85 * min_nozzle_diameter;
        input_params.min_feature_size = 0.25 * min_nozzle_diameter;
        input_params.wall_transition_length = 1.0 * min_nozzle_diameter;
        input_params.wall_transition_angle = 10;
        input_params.wall_transition_filter_deviation = 0.25 * min_nozzle_diameter;
        input_params.wall_distribution_count = 1;

        Arachne::WallToolPaths wallToolPaths(polygons, min_spacing, min_spacing, loops_count, 0, params.layer_height, input_params);

        std::vector<Arachne::VariableWidthLines>    loops = wallToolPaths.getToolPaths();
        std::vector<const Arachne::ExtrusionLine*> all_extrusions;
        for (Arachne::VariableWidthLines& loop : loops) {
            if (loop.empty())
                continue;
            for (const Arachne::ExtrusionLine& wall : loop)
                all_extrusions.emplace_back(&wall);
        }

        // Reverse extrusions to print from center (innermost) to outside (outermost)
        std::reverse(all_extrusions.begin(), all_extrusions.end());
        
        // Calculate center of bounding box for gap alignment
        Point center = expolygon.contour.bounding_box().center();
        
        // Convert config angle from degrees to radians and normalize to [-PI, PI]
        double gap_angle_degrees = params.config ? params.config->concentric_gap_angle.value : 90.0;
        double target_angle = gap_angle_degrees * M_PI / 180.0;
        while (target_angle > M_PI) target_angle -= 2.0 * M_PI;
        while (target_angle <= -M_PI) target_angle += 2.0 * M_PI;

        // Calculate the average centroid of all Arachne loops
        Point average_center = expolygon.contour.bounding_box().center(); // Default
        if (!all_extrusions.empty()) {
            Vec2d center_sum(0, 0);
            size_t count = 0;
            for (const Arachne::ExtrusionLine* extrusion : all_extrusions) {
                if (extrusion->empty() || !extrusion->is_closed) continue;
                // Convert to Polygon for stable centroid
                ThickPolyline tp = Arachne::to_thick_polyline(*extrusion);
                Polygon poly(tp.points);
                center_sum += poly.centroid().cast<double>();
                count++;
            }
            if (count > 0) {
                Vec2d avg = center_sum / double(count);
                average_center = Point(avg.x(), avg.y());
            }
        }
        
        // Split paths at points aligned with the radial line from center at the specified angle
        size_t firts_poly_idx = thick_polylines_out.size();
        for (const Arachne::ExtrusionLine* extrusion : all_extrusions) {
            if (extrusion->empty())
                continue;

            ThickPolyline thick_polyline = Arachne::to_thick_polyline(*extrusion);
            if (extrusion->is_closed) {
                // Use the average center of ALL loops to ensure straight gaps.
                Point center = average_center;

                size_t split_idx = 0;
                Point intersection;
                // Note: ThickPolyline points form a closed loop (start==end) if is_closed is true.
                // find_intersection expects false for is_closed_polygon if it's explicitly closed (Polyline).
                bool found = find_intersection_aligned_with_angle(thick_polyline.points, false, center, target_angle, split_idx, intersection);

                if (found) {
                    // Interpolate width
                    double w_start = thick_polyline.width[2 * split_idx];
                    double w_end   = thick_polyline.width[2 * split_idx + 1];
                    
                    const Point& p1 = thick_polyline.points[split_idx];
                    const Point& p2 = thick_polyline.points[split_idx + 1];
                    double total_len = (p2 - p1).cast<double>().norm();
                    double part_len  = (intersection - p1).cast<double>().norm();
                    double t = (total_len > EPSILON) ? part_len / total_len : 0.0;
                    
                    double w_interp = w_start + t * (w_end - w_start);
                    
                    // Modify widths: replace [2*i, 2*i+1] with [2*i, w_interp] and [w_interp, 2*i+1]
                    // Actually, we insert two new entries.
                    // Original segment i: w[2*i] -> w[2*i+1]
                    // New segment i: w[2*i] -> w_interp
                    // New segment i+1: w_interp -> w[2*i+1]
                    
                    thick_polyline.width[2 * split_idx + 1] = w_interp;
                    thick_polyline.width.insert(thick_polyline.width.begin() + 2 * split_idx + 2, w_interp);
                    thick_polyline.width.insert(thick_polyline.width.begin() + 2 * split_idx + 3, w_end); 
                    
                    // Insert point
                    thick_polyline.points.insert(thick_polyline.points.begin() + split_idx + 1, intersection);
                    
                    thick_polyline.start_at_index(split_idx + 1);
                } else {
                    // Fallback to nearest vertex
                    size_t best_idx = 0;
                    double best_diff = std::numeric_limits<double>::max();
                    
                    for (size_t i = 0; i < thick_polyline.points.size(); ++i) {
                        double dx = thick_polyline.points[i].x() - center.x();
                        double dy = thick_polyline.points[i].y() - center.y();
                        
                        double angle = std::atan2(dy, dx);
                        double diff = std::abs(angle - target_angle);
                        if (diff > M_PI) diff = 2.0 * M_PI - diff;
                        
                        if (diff < best_diff) {
                            best_diff = diff;
                            best_idx = i;
                        }
                    }
                    thick_polyline.start_at_index(best_idx);
                }
            }
            thick_polylines_out.emplace_back(std::move(thick_polyline));
        }

        // clip the paths to prevent the extruder from getting exactly on the first point of the loop
        // Keep valid paths only.
        size_t j = firts_poly_idx;
        for (size_t i = firts_poly_idx; i < thick_polylines_out.size(); ++i) {
            thick_polylines_out[i].clip_end(this->loop_clipping);
            if (thick_polylines_out[i].is_valid()) {
                if (j < i)
                    thick_polylines_out[j] = std::move(thick_polylines_out[i]);
                ++j;
            }
        }
        if (j < thick_polylines_out.size())
            thick_polylines_out.erase(thick_polylines_out.begin() + int(j), thick_polylines_out.end());

        // Note: removed reorder_by_shortest_traverse() to preserve center-to-outside order
    }
    else {
        Polylines polylines;
        this->_fill_surface_single(params, thickness_layers, direction, expolygon, polylines);
        append(thick_polylines_out, to_thick_polylines(std::move(polylines), min_spacing));
    }
}

} // namespace Slic3r
