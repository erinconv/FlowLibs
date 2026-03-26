#pragma once

#include "stdefx.h"

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include <richdem/common/Array2D.hpp>
#include <richdem/common/constants.hpp>

namespace stream_vectorization {

struct StreamVectorizationParams {
    std::string input_dem;
    std::string input_flowdirs;
    std::string input_accumulation;
    std::string input_streams;
    std::string input_stream_order;
    std::string input_watersheds;
    std::string input_subcatchments;
    std::string output_vector;
    bool use_tau_format = false;
};

enum class NodeType {
    Headwater,
    Junction,
    Outlet,
};

struct NodeInfo {
    int id = 0;
    NodeType type = NodeType::Headwater;
    int x = 0;
    int y = 0;
};

inline const char *node_type_name(const NodeType type) {
    switch (type) {
        case NodeType::Headwater:
            return "headwater";
        case NodeType::Junction:
            return "junction";
        case NodeType::Outlet:
            return "outlet";
    }

    return "unknown";
}

inline std::map<int, richdem::flowdir_t> build_direction_decoder(const bool use_tau_format) {
    constexpr std::array<int, 9> esri_dirs = {0, 16, 32, 64, 128, 1, 2, 4, 8};
    constexpr std::array<int, 9> tau_dirs = {0, 5, 4, 3, 2, 1, 8, 7, 6};
    constexpr std::array<richdem::flowdir_t, 9> internal_dirs = {0, 1, 2, 3, 4, 5, 6, 7, 8};

    std::map<int, richdem::flowdir_t> decoder;
    const auto &src_dirs = use_tau_format ? tau_dirs : esri_dirs;

    for (std::size_t index = 0; index < src_dirs.size(); ++index) {
        decoder.emplace(src_dirs[index], internal_dirs[index]);
    }

    return decoder;
}

inline void pixel_to_world(
    const std::vector<double> &geotransform,
    const int x,
    const int y,
    double &world_x,
    double &world_y) {
    const double pixel_x = static_cast<double>(x) + 0.5;
    const double pixel_y = static_cast<double>(y) + 0.5;

    world_x = geotransform[0] + pixel_x * geotransform[1] + pixel_y * geotransform[2];
    world_y = geotransform[3] + pixel_x * geotransform[4] + pixel_y * geotransform[5];
}

inline bool is_stream_cell(const richdem::Array2D<uint8_t> &stream_grid, const int x, const int y) {
    return stream_grid.inGrid(x, y) && !stream_grid.isNoData(x, y) && stream_grid(x, y) == 1;
}

inline bool flows_into(
    const richdem::Array2D<richdem::flowdir_t> &flowdirs,
    const int from_x,
    const int from_y,
    const int to_x,
    const int to_y) {
    if (!flowdirs.inGrid(from_x, from_y) || flowdirs.isNoData(from_x, from_y)) {
        return false;
    }

    const auto dir = flowdirs(from_x, from_y);
    if (dir == richdem::NO_FLOW) {
        return false;
    }

    return from_x + richdem::d8x[dir] == to_x && from_y + richdem::d8y[dir] == to_y;
}

inline int count_upstream_stream_cells(
    const richdem::Array2D<richdem::flowdir_t> &flowdirs,
    const richdem::Array2D<uint8_t> &stream_grid,
    const int x,
    const int y) {
    int upstream_count = 0;

    for (int dir = 1; dir <= 8; ++dir) {
        const int nx = x + richdem::d8x[dir];
        const int ny = y + richdem::d8y[dir];
        if (!is_stream_cell(stream_grid, nx, ny)) {
            continue;
        }
        if (flows_into(flowdirs, nx, ny, x, y)) {
            ++upstream_count;
        }
    }

    return upstream_count;
}

inline bool has_stream_downstream(
    const richdem::Array2D<richdem::flowdir_t> &flowdirs,
    const richdem::Array2D<uint8_t> &stream_grid,
    const int x,
    const int y,
    int &downstream_x,
    int &downstream_y) {
    downstream_x = x;
    downstream_y = y;

    if (flowdirs.isNoData(x, y)) {
        return false;
    }

    const auto dir = flowdirs(x, y);
    if (dir == richdem::NO_FLOW) {
        return false;
    }

    downstream_x = x + richdem::d8x[dir];
    downstream_y = y + richdem::d8y[dir];

    return is_stream_cell(stream_grid, downstream_x, downstream_y);
}

inline OGRSpatialReference build_spatial_reference(const std::string &projection_wkt) {
    OGRSpatialReference spatial_ref;
    if (!projection_wkt.empty()) {
        spatial_ref.importFromWkt(projection_wkt.c_str());
    }
    return spatial_ref;
}

inline void add_common_node_fields(OGRLayer *layer) {
    OGRFieldDefn field_node_id("node_id", OFTInteger);
    OGRFieldDefn field_node_type("node_type", OFTString);
    OGRFieldDefn field_grid_x("grid_x", OFTInteger);
    OGRFieldDefn field_grid_y("grid_y", OFTInteger);
    OGRFieldDefn field_fac("fac", OFTReal);
    OGRFieldDefn field_stream_order("stream_ord", OFTInteger);
    OGRFieldDefn field_watershed("watershed", OFTInteger);
    OGRFieldDefn field_subcatch("subcatch", OFTInteger);
    OGRFieldDefn field_elev("elev", OFTReal);

    layer->CreateField(&field_node_id);
    layer->CreateField(&field_node_type);
    layer->CreateField(&field_grid_x);
    layer->CreateField(&field_grid_y);
    layer->CreateField(&field_fac);
    layer->CreateField(&field_stream_order);
    layer->CreateField(&field_watershed);
    layer->CreateField(&field_subcatch);
    layer->CreateField(&field_elev);
}

inline void add_common_segment_fields(OGRLayer *layer) {
    OGRFieldDefn field_segment_id("segment_id", OFTInteger);
    OGRFieldDefn field_from_node("from_node", OFTInteger);
    OGRFieldDefn field_to_node("to_node", OFTInteger);
    OGRFieldDefn field_from_type("from_type", OFTString);
    OGRFieldDefn field_to_type("to_type", OFTString);
    OGRFieldDefn field_cell_count("cell_count", OFTInteger);
    OGRFieldDefn field_length("length_m", OFTReal);
    OGRFieldDefn field_drop("drop_z", OFTReal);
    OGRFieldDefn field_slope("slope", OFTReal);
    OGRFieldDefn field_fac_start("fac_start", OFTReal);
    OGRFieldDefn field_fac_end("fac_end", OFTReal);
    OGRFieldDefn field_fac_max("fac_max", OFTReal);
    OGRFieldDefn field_stream_order("stream_ord", OFTInteger);
    OGRFieldDefn field_watershed("watershed", OFTInteger);
    OGRFieldDefn field_subcatch("subcatch", OFTInteger);

    layer->CreateField(&field_segment_id);
    layer->CreateField(&field_from_node);
    layer->CreateField(&field_to_node);
    layer->CreateField(&field_from_type);
    layer->CreateField(&field_to_type);
    layer->CreateField(&field_cell_count);
    layer->CreateField(&field_length);
    layer->CreateField(&field_drop);
    layer->CreateField(&field_slope);
    layer->CreateField(&field_fac_start);
    layer->CreateField(&field_fac_end);
    layer->CreateField(&field_fac_max);
    layer->CreateField(&field_stream_order);
    layer->CreateField(&field_watershed);
    layer->CreateField(&field_subcatch);
}

inline int export_stream_network_vector(const StreamVectorizationParams &params) {
    GDALAllRegister();

    richdem::Array2D<float> dem(params.input_dem);
    richdem::Array2D<int> encoded_flowdirs(params.input_flowdirs);
    richdem::Array2D<double> accumulation(params.input_accumulation);
    richdem::Array2D<uint8_t> stream_grid(params.input_streams);
    richdem::Array2D<int> stream_order(params.input_stream_order);
    richdem::Array2D<int> watershed_labels(params.input_watersheds);
    richdem::Array2D<int> subcatchment_labels(params.input_subcatchments);

    richdem::Array2D<richdem::flowdir_t> flowdirs(encoded_flowdirs, richdem::NO_FLOW);
    const auto decoder = build_direction_decoder(params.use_tau_format);

    for (int y = 0; y < encoded_flowdirs.height(); ++y) {
        for (int x = 0; x < encoded_flowdirs.width(); ++x) {
            if (encoded_flowdirs.isNoData(x, y)) {
                continue;
            }

            const auto match = decoder.find(encoded_flowdirs(x, y));
            if (match == decoder.end()) {
                std::cerr << "ERROR: Invalid flow direction value " << encoded_flowdirs(x, y)
                          << " at (" << x << ", " << y << ")." << std::endl;
                return -1;
            }
            flowdirs(x, y) = match->second;
        }
    }

    std::map<std::size_t, NodeInfo> node_lookup;
    int next_node_id = 1;

    auto ensure_node = [&](const int x, const int y, const NodeType type) -> const NodeInfo & {
        const auto cell_id = flowdirs.xyToI(x, y);
        auto found = node_lookup.find(cell_id);
        if (found == node_lookup.end()) {
            NodeInfo node;
            node.id = next_node_id++;
            node.type = type;
            node.x = x;
            node.y = y;
            found = node_lookup.emplace(cell_id, node).first;
        } else if (type == NodeType::Outlet || found->second.type == NodeType::Headwater) {
            found->second.type = type;
        }

        return found->second;
    };

    for (int y = 0; y < stream_grid.height(); ++y) {
        for (int x = 0; x < stream_grid.width(); ++x) {
            if (!is_stream_cell(stream_grid, x, y)) {
                continue;
            }

            const int upstream_count = count_upstream_stream_cells(flowdirs, stream_grid, x, y);
            int downstream_x = x;
            int downstream_y = y;
            const bool has_downstream = has_stream_downstream(flowdirs, stream_grid, x, y, downstream_x, downstream_y);

            if (upstream_count == 0) {
                ensure_node(x, y, NodeType::Headwater);
            } else if (upstream_count > 1) {
                ensure_node(x, y, NodeType::Junction);
            }

            if (!has_downstream) {
                ensure_node(x, y, NodeType::Outlet);
            }
        }
    }

    GDALDriver *driver = GetGDALDriverManager()->GetDriverByName("GPKG");
    if (driver == nullptr) {
        std::cerr << "ERROR: GPKG driver is not available." << std::endl;
        return -1;
    }

    driver->Delete(params.output_vector.c_str());
    GDALDataset *dataset = driver->Create(params.output_vector.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (dataset == nullptr) {
        std::cerr << "ERROR: Failed to create output dataset: " << params.output_vector << std::endl;
        return -1;
    }

    auto spatial_ref = build_spatial_reference(dem.projection);
    OGRLayer *nodes_layer = dataset->CreateLayer("river_nodes", &spatial_ref, wkbPoint, nullptr);
    OGRLayer *segments_layer = dataset->CreateLayer("river_segments", &spatial_ref, wkbLineString, nullptr);
    if (nodes_layer == nullptr || segments_layer == nullptr) {
        GDALClose(dataset);
        std::cerr << "ERROR: Failed to create output layers." << std::endl;
        return -1;
    }

    add_common_node_fields(nodes_layer);
    add_common_segment_fields(segments_layer);

    for (const auto &[cell_id, node] : node_lookup) {
        static_cast<void>(cell_id);

        OGRFeature *feature = OGRFeature::CreateFeature(nodes_layer->GetLayerDefn());
        feature->SetField("node_id", node.id);
        feature->SetField("node_type", node_type_name(node.type));
        feature->SetField("grid_x", node.x);
        feature->SetField("grid_y", node.y);
        feature->SetField("fac", accumulation(node.x, node.y));
        feature->SetField("stream_ord", stream_order(node.x, node.y));
        feature->SetField("watershed", watershed_labels(node.x, node.y));
        feature->SetField("subcatch", subcatchment_labels(node.x, node.y));
        feature->SetField("elev", dem(node.x, node.y));

        double world_x = 0.0;
        double world_y = 0.0;
        pixel_to_world(dem.geotransform, node.x, node.y, world_x, world_y);
        OGRPoint point(world_x, world_y);
        feature->SetGeometry(&point);

        if (nodes_layer->CreateFeature(feature) != OGRERR_NONE) {
            OGRFeature::DestroyFeature(feature);
            GDALClose(dataset);
            std::cerr << "ERROR: Failed to create node feature." << std::endl;
            return -1;
        }
        OGRFeature::DestroyFeature(feature);
    }

    int next_segment_id = 1;

    for (const auto &[cell_id, start_node] : node_lookup) {
        static_cast<void>(cell_id);

        int downstream_x = start_node.x;
        int downstream_y = start_node.y;
        if (!has_stream_downstream(flowdirs, stream_grid, start_node.x, start_node.y, downstream_x, downstream_y)) {
            continue;
        }

        OGRLineString line;
        std::vector<std::pair<int, int>> cells;
        cells.emplace_back(start_node.x, start_node.y);

        double start_world_x = 0.0;
        double start_world_y = 0.0;
        pixel_to_world(dem.geotransform, start_node.x, start_node.y, start_world_x, start_world_y);
        line.addPoint(start_world_x, start_world_y);

        int current_x = start_node.x;
        int current_y = start_node.y;
        NodeInfo end_node = start_node;
        double fac_max = accumulation(current_x, current_y);

        while (true) {
            int next_x = current_x;
            int next_y = current_y;
            const bool continue_downstream = has_stream_downstream(flowdirs, stream_grid, current_x, current_y, next_x, next_y);
            if (!continue_downstream) {
                end_node = ensure_node(current_x, current_y, NodeType::Outlet);
                break;
            }

            current_x = next_x;
            current_y = next_y;
            cells.emplace_back(current_x, current_y);
            fac_max = std::max(fac_max, accumulation(current_x, current_y));

            double world_x = 0.0;
            double world_y = 0.0;
            pixel_to_world(dem.geotransform, current_x, current_y, world_x, world_y);
            line.addPoint(world_x, world_y);

            const int upstream_count = count_upstream_stream_cells(flowdirs, stream_grid, current_x, current_y);
            int downstream2_x = current_x;
            int downstream2_y = current_y;
            const bool has_downstream2 = has_stream_downstream(flowdirs, stream_grid, current_x, current_y, downstream2_x, downstream2_y);
            if (upstream_count != 1 || !has_downstream2) {
                end_node = ensure_node(current_x, current_y, !has_downstream2 ? NodeType::Outlet : NodeType::Junction);
                break;
            }
        }

        OGRFeature *feature = OGRFeature::CreateFeature(segments_layer->GetLayerDefn());
        feature->SetField("segment_id", next_segment_id++);
        feature->SetField("from_node", start_node.id);
        feature->SetField("to_node", end_node.id);
        feature->SetField("from_type", node_type_name(start_node.type));
        feature->SetField("to_type", node_type_name(end_node.type));
        feature->SetField("cell_count", static_cast<int>(cells.size()));
        feature->SetField("length_m", line.get_Length());

        const double drop_z = dem(start_node.x, start_node.y) - dem(end_node.x, end_node.y);
        feature->SetField("drop_z", drop_z);
        const double length_value = line.get_Length();
        feature->SetField("slope", length_value > 0.0 ? drop_z / length_value : 0.0);
        feature->SetField("fac_start", accumulation(start_node.x, start_node.y));
        feature->SetField("fac_end", accumulation(end_node.x, end_node.y));
        feature->SetField("fac_max", fac_max);
        // Reach order should describe the traced segment itself, so use the
        // upstream/start node order rather than the downstream junction cell.
        feature->SetField("stream_ord", stream_order(start_node.x, start_node.y));
        feature->SetField("watershed", watershed_labels(end_node.x, end_node.y));
        feature->SetField("subcatch", subcatchment_labels(end_node.x, end_node.y));
        feature->SetGeometry(&line);

        if (segments_layer->CreateFeature(feature) != OGRERR_NONE) {
            OGRFeature::DestroyFeature(feature);
            GDALClose(dataset);
            std::cerr << "ERROR: Failed to create segment feature." << std::endl;
            return -1;
        }
        OGRFeature::DestroyFeature(feature);
    }

    GDALClose(dataset);
    return 0;
}

}  // namespace stream_vectorization
