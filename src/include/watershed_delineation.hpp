#include "stdefx.h"

#include <richdem/common/Array2D.hpp>
#include <richdem/common/constants.hpp>


/**
 * @brief Parameters for watershed delineation operations
 */
struct WatershedParams {
    std::string input_dem;              // Carved DEM file path
    std::string input_flowdirs;         // Flow directions file path
    std::string input_accumulation;     // Flow accumulation file path
    std::string output_prefix;          // Output file prefix
    float stream_threshold;             // Accumulation threshold for stream definition
    bool use_tau_format = false;                // True if flow directions are TauDEM format
    
    // Default constructor
    WatershedParams() 
        : stream_threshold(1000.0f), use_tau_format(false) {}
    
    // Parameterized constructor
    WatershedParams(const std::string& dem, const std::string& flowdirs, 
                   const std::string& accu, const std::string& output,
                   float threshold = 1000.0f, bool tau = false)
        : input_dem(dem), input_flowdirs(flowdirs), input_accumulation(accu),
          output_prefix(output), stream_threshold(threshold), use_tau_format(tau) {}
};

/**
 * @brief Structure to represent a pour point (watershed outlet)
 */
struct PourPoint {
    double geo_x;       // Geographic X coordinate (longitude in decimal degrees)
    double geo_y;       // Geographic Y coordinate (latitude in decimal degrees)
    int pixel_x;        // Pixel X coordinate (column) - computed from geo coords
    int pixel_y;        // Pixel Y coordinate (row) - computed from geo coords
    int id;             // Unique identifier for this watershed
    std::string name;   // Optional name for the watershed
    bool is_geographic; // True if coordinates are geographic (need conversion)
    
    // Constructor for geographic coordinates (decimal degrees)
    PourPoint(double lon, double lat, int watershed_id, const std::string& ws_name = "")
        : geo_x(lon), geo_y(lat), pixel_x(-1), pixel_y(-1), 
          id(watershed_id), name(ws_name), is_geographic(true) {}
    
    // Constructor for pixel coordinates (for backward compatibility)
    PourPoint(int px, int py, int watershed_id, const std::string& ws_name = "", bool is_pixel = true)
        : geo_x(0), geo_y(0), pixel_x(px), pixel_y(py), 
          id(watershed_id), name(ws_name), is_geographic(false) {}
};

/**
 * @brief Structure to hold watershed statistics
 */
struct WatershedStats {
    int watershed_id;
    double area_cells;          // Area in number of cells
    double area_km2;            // Area in square kilometers
    double mean_elevation;      // Mean elevation in watershed
    double min_elevation;       // Minimum elevation
    double max_elevation;       // Maximum elevation
    int stream_cells;           // Number of stream cells
    double total_stream_length; // Total stream length in km
    
    WatershedStats() 
        : watershed_id(0), area_cells(0), area_km2(0), 
          mean_elevation(0), min_elevation(0), max_elevation(0),
          stream_cells(0), total_stream_length(0) {}
};

/**
 * @brief Class for watershed and subcatchment delineation
 * 
 * This class provides methods to:
 * - Delineate watersheds from pour points
 * - Extract stream networks based on accumulation thresholds
 * - Identify and label subcatchments
 * - Calculate watershed statistics
 */
class WatershedDelineation {
private:
    WatershedParams params;  // Configuration parameters
    
public:
    /**
     * @brief Constructor
     * @param params Configuration parameters for watershed delineation
     */
    WatershedDelineation(const WatershedParams& params) : params(params) {}
    
    /**
     * @brief Destructor
     */
    ~WatershedDelineation() = default;
    
    /**
     * @brief Convert geographic coordinates to pixel coordinates
     * 
     * Uses the DEM's geotransform to convert from geographic coordinates
     * (longitude, latitude in decimal degrees) to pixel coordinates (column, row).
     * 
     * @param geotransform GDAL geotransform vector [origin_x, pixel_width, 0, origin_y, 0, pixel_height]
     * @param geo_x Geographic X coordinate (longitude)
     * @param geo_y Geographic Y coordinate (latitude)
     * @param pixel_x Output pixel X coordinate (column)
     * @param pixel_y Output pixel Y coordinate (row)
     */
    void convert_geographic_to_pixel(const std::vector<double>& geotransform,
                                     double geo_x, double geo_y,
                                     int& pixel_x, int& pixel_y) {
        // Geotransform: [origin_x, pixel_width, 0, origin_y, 0, pixel_height]
        // Forward transform: geo_x = gt[0] + pixel_x * gt[1] + pixel_y * gt[2]
        //                    geo_y = gt[3] + pixel_x * gt[4] + pixel_y * gt[5]
        // Inverse transform (assuming gt[2]=0 and gt[4]=0):
        pixel_x = static_cast<int>((geo_x - geotransform[0]) / geotransform[1]);
        pixel_y = static_cast<int>((geo_y - geotransform[3]) / geotransform[5]);
    }
    
    /**
     * @brief Snap pour point to maximum flow accumulation within search radius
     * 
     * Finds the cell with maximum flow accumulation within a search radius around
     * the specified pour point. This ensures the outlet is located at the actual
     * stream position even if the user-specified coordinates are slightly off.
     * 
     * Algorithm:
     * 1. Search all cells within radius of pour point
     * 2. Find cell with maximum accumulation value
     * 3. If significantly higher than current location, snap to it
     * 4. Report adjustment distance
     * 
     * @param accumulation Flow accumulation grid
     * @param pour_point Pour point to adjust (modified in place)
     * @param search_radius Search radius in pixels (default: 10)
     * @param min_improvement Minimum accumulation ratio to trigger snap (default: 1.5)
     * @return True if pour point was adjusted, false otherwise
     */
    template<typename T>
    bool snap_pourpoint_to_stream(const richdem::Array2D<T>& accumulation,
                                  PourPoint& pour_point,
                                  int search_radius = 10,
                                  double min_improvement = 1.5);
    
    /**
     * @brief Convert pour points from geographic to pixel coordinates
     * 
     * @param pour_points Vector of pour points (may contain geographic coordinates)
     * @param geotransform GDAL geotransform from the DEM
     * @param width DEM width in pixels
     * @param height DEM height in pixels
     */
    void convert_pourpoints_to_pixels(std::vector<PourPoint>& pour_points,
                                      const std::vector<double>& geotransform,
                                      int width, int height) {
        std::cout << "Converting pour point coordinates..." << std::endl;
        
        for (auto& pp : pour_points) {
            if (pp.is_geographic) {
                convert_geographic_to_pixel(geotransform, pp.geo_x, pp.geo_y, 
                                          pp.pixel_x, pp.pixel_y);
                
                std::cout << "  Pour point " << pp.id << ": ";
                std::cout << "(" << pp.geo_x << "°, " << pp.geo_y << "°) → ";
                std::cout << "pixel (" << pp.pixel_x << ", " << pp.pixel_y << ")";
                
                // Validate coordinates are within bounds
                if (pp.pixel_x < 0 || pp.pixel_x >= width || 
                    pp.pixel_y < 0 || pp.pixel_y >= height) {
                    std::cout << " [OUT OF BOUNDS!]";
                }
                std::cout << std::endl;
            }
        }
    }
    
    /**
     * @brief Delineate a single watershed from a pour point
     * 
     * This function traces upstream from a pour point following the inverse
     * of the flow directions to identify all cells that drain to the outlet.
     * 
     * Algorithm:
     * 1. Start at the pour point
     * 2. For each cell, check all 8 neighbors
     * 3. If a neighbor flows into the current cell, add it to the watershed
     * 4. Continue until no more upstream cells are found
     * 
     * @param flowdirs Flow direction grid (D8 format)
     * @param pour_point The outlet point for watershed delineation
     * @param watershed_label Output grid with watershed labels
     * @param watershed_id The ID to assign to this watershed
     */
    template<typename T>
    void delineate_watershed(const richdem::Array2D<richdem::flowdir_t>& flowdirs,
                            const PourPoint& pour_point,
                            richdem::Array2D<T>& watershed_label,
                            T watershed_id);
    
    /**
     * @brief Delineate multiple watersheds from a list of pour points
     * 
     * @param flowdirs Flow direction grid (D8 format)
     * @param pour_points Vector of pour points
     * @param watershed_label Output grid with watershed labels
     * @return Number of watersheds successfully delineated
     */
    template<typename T>
    int delineate_multiple_watersheds(const richdem::Array2D<richdem::flowdir_t>& flowdirs,
                                      const std::vector<PourPoint>& pour_points,
                                      richdem::Array2D<T>& watershed_label);
    
    /**
     * @brief Extract stream network based on accumulation threshold
     * 
     * Cells with accumulation >= threshold are classified as stream cells.
     * 
     * @param accumulation Flow accumulation grid
     * @param stream_grid Output binary grid (1 = stream, 0 = non-stream)
     * @param threshold Accumulation threshold for stream definition
     */
    template<typename T>
    void extract_stream_network(const richdem::Array2D<T>& accumulation,
                               richdem::Array2D<uint8_t>& stream_grid,
                               T threshold);
    
    /**
     * @brief Calculate Horton-Strahler stream order
     * 
     * Assigns stream order to each stream cell following the Horton-Strahler system:
     * - Headwater streams (no upstream) = order 1
     * - When two streams of order n join = order n+1
     * - When streams of different orders join = higher order continues
     * 
     * Algorithm:
     * 1. Count upstream stream cells for each cell
     * 2. Process streams from headwaters downstream
     * 3. Assign orders based on upstream stream orders
     * 
     * @param flowdirs Flow direction grid
     * @param stream_grid Binary stream grid
     * @param stream_order Output grid with stream orders
     */
    void calculate_stream_order(
        const richdem::Array2D<richdem::flowdir_t>& flowdirs,
        const richdem::Array2D<uint8_t>& stream_grid,
        richdem::Array2D<int>& stream_order);
    
    /**
     * @brief Identify stream junctions and endpoints
     * 
     * Finds cells where streams merge (junctions) and where streams start (headwaters).
     * These become subcatchment outlets.
     * 
     * @param flowdirs Flow direction grid
     * @param stream_grid Binary stream grid
     * @param watershed_mask Mask defining the watershed extent
     * @return Vector of pour points for subcatchments
     */
    template<typename T>
    std::vector<PourPoint> identify_subcatchment_outlets(
        const richdem::Array2D<richdem::flowdir_t>& flowdirs,
        const richdem::Array2D<uint8_t>& stream_grid,
        const richdem::Array2D<T>& watershed_mask);
    
    /**
     * @brief Identify subcatchments for stream network
     * 
     * Divides a watershed into subcatchments by:
     * 1. Finding stream junctions and endpoints
     * 2. Delineating contributing areas to each junction/endpoint
     * 3. Assigning unique IDs to each subcatchment
     * 
     * @param flowdirs Flow direction grid
     * @param stream_grid Binary stream grid
     * @param watershed_mask Mask defining the watershed extent
     * @param subcatchment_label Output grid with subcatchment IDs
     * @return Number of subcatchments identified
     */
    template<typename T>
    int identify_subcatchments(const richdem::Array2D<richdem::flowdir_t>& flowdirs,
                              const richdem::Array2D<uint8_t>& stream_grid,
                              const richdem::Array2D<T>& watershed_mask,
                              richdem::Array2D<int>& subcatchment_label);
    
    /**
     * @brief Calculate statistics for a watershed
     * 
     * @param dem Digital elevation model
     * @param watershed_label Watershed label grid
     * @param stream_grid Stream network grid
     * @param watershed_id The watershed ID to calculate stats for
     * @return WatershedStats structure with calculated statistics
     */
    template<typename T>
    WatershedStats calculate_watershed_stats(const richdem::Array2D<T>& dem,
                                            const richdem::Array2D<int>& watershed_label,
                                            const richdem::Array2D<uint8_t>& stream_grid,
                                            int watershed_id);
    
    /**
     * @brief Main execution function for watershed delineation workflow
     * 
     * This orchestrates the complete workflow:
     * 1. Load input data (DEM, flow directions, accumulation)
     * 2. Extract stream network
     * 3. Delineate watersheds from pour points
     * 4. Identify subcatchments
     * 5. Calculate statistics
     * 6. Save results
     * 
     * @param pour_points Vector of pour points for watershed delineation
     * @return 0 on success, negative on error
     */
    int run_watershed_delineation(const std::vector<PourPoint>& pour_points);
};

// ============================================================================
// IMPLEMENTATION - Non-Template Functions
// ============================================================================

void WatershedDelineation::calculate_stream_order(
    const richdem::Array2D<richdem::flowdir_t>& flowdirs,
    const richdem::Array2D<uint8_t>& stream_grid,
    richdem::Array2D<int>& stream_order)
{
    std::cout << "Calculating Horton-Strahler stream order..." << std::endl;
    
    // Initialize stream order grid (0 = non-stream)
    stream_order.setAll(0);
    
    // Count number of upstream stream cells for each stream cell
    richdem::Array2D<int> upstream_count(stream_grid.width(), stream_grid.height(), 0);
    
    for (int y = 0; y < stream_grid.height(); y++) {
        for (int x = 0; x < stream_grid.width(); x++) {
            if (stream_grid(x, y) != 1) {
                continue;
            }
            
            // Check where this stream cell flows
            if (flowdirs(x, y) == richdem::NO_FLOW) {
                continue;
            }
            
            int next_x = x + richdem::d8x[flowdirs(x, y)];
            int next_y = y + richdem::d8y[flowdirs(x, y)];
            
            // If flows to another stream cell, increment its upstream count
            if (flowdirs.inGrid(next_x, next_y) && stream_grid(next_x, next_y) == 1) {
                upstream_count(next_x, next_y)++;
            }
        }
    }
    
    // Find headwater cells (no upstream stream cells) and assign order 1
    std::queue<std::pair<int, int>> to_process;
    int headwater_count = 0;
    
    for (int y = 0; y < stream_grid.height(); y++) {
        for (int x = 0; x < stream_grid.width(); x++) {
            if (stream_grid(x, y) == 1 && upstream_count(x, y) == 0) {
                stream_order(x, y) = 1;  // Headwater = order 1
                to_process.push({x, y});
                headwater_count++;
            }
        }
    }
    
    std::cout << "  Found " << headwater_count << " headwater streams (order 1)" << std::endl;
    
    // Process streams from headwaters downstream
    int max_order = 1;
    int cells_processed = 0;
    
    while (!to_process.empty()) {
        auto [cx, cy] = to_process.front();
        to_process.pop();
        cells_processed++;
        
        int current_order = stream_order(cx, cy);
        
        // Find where this stream cell flows
        if (flowdirs(cx, cy) == richdem::NO_FLOW) {
            continue;
        }
        
        int next_x = cx + richdem::d8x[flowdirs(cx, cy)];
        int next_y = cy + richdem::d8y[flowdirs(cx, cy)];
        
        // Check if flows to another stream cell
        if (!flowdirs.inGrid(next_x, next_y) || stream_grid(next_x, next_y) != 1) {
            continue;
        }
        
        // Decrement upstream count for downstream cell
        upstream_count(next_x, next_y)--;
        
        // If all upstream cells processed, calculate order for downstream cell
        if (upstream_count(next_x, next_y) == 0) {
            // Find all upstream stream cells and their orders
            std::vector<int> upstream_orders;
            
            for (int dir = 1; dir <= 8; dir++) {
                int ux = next_x + richdem::d8x[dir];
                int uy = next_y + richdem::d8y[dir];
                
                if (!flowdirs.inGrid(ux, uy) || stream_grid(ux, uy) != 1) {
                    continue;
                }
                
                // Check if this neighbor flows into next cell
                if (flowdirs(ux, uy) != richdem::NO_FLOW) {
                    int flow_to_x = ux + richdem::d8x[flowdirs(ux, uy)];
                    int flow_to_y = uy + richdem::d8y[flowdirs(ux, uy)];
                    
                    if (flow_to_x == next_x && flow_to_y == next_y) {
                        upstream_orders.push_back(stream_order(ux, uy));
                    }
                }
            }
            
            // Apply Horton-Strahler rules
            if (upstream_orders.empty()) {
                stream_order(next_x, next_y) = 1;  // Shouldn't happen, but default to 1
            } else if (upstream_orders.size() == 1) {
                stream_order(next_x, next_y) = upstream_orders[0];  // Single tributary
            } else {
                // Sort orders
                std::sort(upstream_orders.begin(), upstream_orders.end(), std::greater<int>());
                
                // If two or more streams have the same highest order, increment
                if (upstream_orders.size() >= 2 && upstream_orders[0] == upstream_orders[1]) {
                    stream_order(next_x, next_y) = upstream_orders[0] + 1;
                } else {
                    // Otherwise, keep the highest order
                    stream_order(next_x, next_y) = upstream_orders[0];
                }
            }
            
            max_order = std::max(max_order, stream_order(next_x, next_y));
            to_process.push({next_x, next_y});
        }
    }
    
    std::cout << "  Processed " << cells_processed << " stream cells" << std::endl;
    std::cout << "  Maximum stream order: " << max_order << std::endl;
    
    // Count cells in each order
    std::map<int, int> order_counts;
    for (int y = 0; y < stream_order.height(); y++) {
        for (int x = 0; x < stream_order.width(); x++) {
            if (stream_order(x, y) > 0) {
                order_counts[stream_order(x, y)]++;
            }
        }
    }
    
    // std::cout << "\n  Stream order distribution:" << std::endl;
    // for (const auto& [order, count] : order_counts) {
    //     std::cout << "    Order " << order << ": " << count << " cells" << std::endl;
    // }
}

// ============================================================================
// IMPLEMENTATION - Template Functions
// ============================================================================

template<typename T>
bool WatershedDelineation::snap_pourpoint_to_stream(
    const richdem::Array2D<T>& accumulation,
    PourPoint& pour_point,
    int search_radius,
    double min_improvement)
{
    // Get current accumulation value at pour point
    if (!accumulation.inGrid(pour_point.pixel_x, pour_point.pixel_y)) {
        std::cerr << "ERROR: Pour point (" << pour_point.pixel_x << ", " 
                  << pour_point.pixel_y << ") is outside grid!" << std::endl;
        return false;
    }
    
    T current_accu = accumulation(pour_point.pixel_x, pour_point.pixel_y);
    
    // Search for maximum accumulation within radius
    T max_accu = current_accu;
    int best_x = pour_point.pixel_x;
    int best_y = pour_point.pixel_y;
    
    // Search in square around pour point
    for (int dy = -search_radius; dy <= search_radius; dy++) {
        for (int dx = -search_radius; dx <= search_radius; dx++) {
            int search_x = pour_point.pixel_x + dx;
            int search_y = pour_point.pixel_y + dy;
            
            // Skip if outside grid
            if (!accumulation.inGrid(search_x, search_y)) {
                continue;
            }
            
            // Check if within circular radius (Euclidean distance)
            double distance = std::sqrt(dx*dx + dy*dy);
            if (distance > search_radius) {
                continue;
            }
            
            // Check if this cell has higher accumulation
            T accu_value = accumulation(search_x, search_y);
            if (accu_value > max_accu) {
                max_accu = accu_value;
                best_x = search_x;
                best_y = search_y;
            }
        }
    }
    
    // Calculate improvement ratio
    double improvement = (current_accu > 0) ? 
        static_cast<double>(max_accu) / static_cast<double>(current_accu) : 
        std::numeric_limits<double>::infinity();
    
    // Snap to new location if significantly better
    if (best_x != pour_point.pixel_x || best_y != pour_point.pixel_y) {
        if (improvement >= min_improvement || current_accu == 0) {
            // Calculate adjustment distance
            int dx = best_x - pour_point.pixel_x;
            int dy = best_y - pour_point.pixel_y;
            double distance = std::sqrt(dx*dx + dy*dy);
            
            std::cout << "  Pour point " << pour_point.id << " adjusted:" << std::endl;
            std::cout << "    Original:  (" << pour_point.pixel_x << ", " << pour_point.pixel_y 
                      << ") accumulation = " << current_accu << std::endl;
            std::cout << "    Snapped:   (" << best_x << ", " << best_y 
                      << ") accumulation = " << max_accu << std::endl;
            std::cout << "    Distance:  " << distance << " pixels" << std::endl;
            std::cout << "    Improvement: " << improvement << "x" << std::endl;
            
            // Update pour point location
            pour_point.pixel_x = best_x;
            pour_point.pixel_y = best_y;
            
            return true;
        }
    }
    
    return false;
}

template<typename T>
void WatershedDelineation::delineate_watershed(
    const richdem::Array2D<richdem::flowdir_t>& flowdirs,
    const PourPoint& pour_point,
    richdem::Array2D<T>& watershed_label,
    T watershed_id)
{
    // std::cout << "  Delineating watershed " << watershed_id 
    //           << " from pour point (pixel: " << pour_point.pixel_x << ", " << pour_point.pixel_y << ")" << std::endl;
    
    // Validate pour point is within grid
    if (!flowdirs.inGrid(pour_point.pixel_x, pour_point.pixel_y)) {
        std::cerr << "ERROR: Pour point (pixel: " << pour_point.pixel_x << ", " 
                  << pour_point.pixel_y << ") is outside grid bounds!" << std::endl;
        return;
    }
    
    // Queue for BFS traversal
    std::queue<std::pair<int, int>> to_process;
    to_process.push({pour_point.pixel_x, pour_point.pixel_y});
    watershed_label(pour_point.pixel_x, pour_point.pixel_y) = watershed_id;
    
    int cells_processed = 0;
    
    // BFS to find all cells that drain to this pour point
    while (!to_process.empty()) {
        auto [cx, cy] = to_process.front();
        to_process.pop();
        cells_processed++;
        
        // Check all 8 neighbors
        for (int dir = 1; dir <= 8; dir++) {
            int nx = cx + richdem::d8x[dir];
            int ny = cy + richdem::d8y[dir];
            
            // Skip if out of bounds
            if (!flowdirs.inGrid(nx, ny)) {
                continue;
            }
            
            // Skip if already labeled
            if (watershed_label(nx, ny) != 0) {
                continue;
            }
            
            // Skip if no flow
            if (flowdirs(nx, ny) == richdem::NO_FLOW) {
                continue;
            }
            
            // Check if neighbor flows into current cell
            // The neighbor flows into current cell if following its flow direction
            // leads back to the current cell
            int neighbor_flow_dir = flowdirs(nx, ny);
            int flow_to_x = nx + richdem::d8x[neighbor_flow_dir];
            int flow_to_y = ny + richdem::d8y[neighbor_flow_dir];
            
            if (flow_to_x == cx && flow_to_y == cy) {
                watershed_label(nx, ny) = watershed_id;
                to_process.push({nx, ny});
            }
        }
    }
    
    // std::cout << "    Watershed " << watershed_id << " contains " 
    //           << cells_processed << " cells" << std::endl;
}

template<typename T>
int WatershedDelineation::delineate_multiple_watersheds(
    const richdem::Array2D<richdem::flowdir_t>& flowdirs,
    const std::vector<PourPoint>& pour_points,
    richdem::Array2D<T>& watershed_label)
{
    std::cout << "Delineating " << pour_points.size() << " watersheds..." << std::endl;
    
    int successful = 0;
    for (const auto& pp : pour_points) {
        delineate_watershed(flowdirs, pp, watershed_label, static_cast<T>(pp.id));
        successful++;
    }
    
    std::cout << "Successfully delineated " << successful << " watersheds" << std::endl;
    return successful;
}

template<typename T>
void WatershedDelineation::extract_stream_network(
    const richdem::Array2D<T>& accumulation,
    richdem::Array2D<uint8_t>& stream_grid,
    T threshold)
{
    std::cout << "Extracting stream network (threshold = " << threshold << " cells)..." << std::endl;
    
    int stream_cells = 0;
    
    #pragma omp parallel for collapse(2) reduction(+:stream_cells)
    for (int y = 0; y < accumulation.height(); y++) {
        for (int x = 0; x < accumulation.width(); x++) {
            if (accumulation(x, y) >= threshold) {
                stream_grid(x, y) = 1;
                stream_cells++;
            } else {
                stream_grid(x, y) = 0;
            }
        }
    }
    
    std::cout << "  Stream network contains " << stream_cells << " cells" << std::endl;
}

template<typename T>
std::vector<PourPoint> WatershedDelineation::identify_subcatchment_outlets(
    const richdem::Array2D<richdem::flowdir_t>& flowdirs,
    const richdem::Array2D<uint8_t>& stream_grid,
    const richdem::Array2D<T>& watershed_mask)
{
    std::cout << "Identifying stream junctions and endpoints..." << std::endl;
    std::vector<PourPoint> outlets;
    
    // Count upstream stream cells for each stream cell
    richdem::Array2D<int> stream_inflow(stream_grid.width(), stream_grid.height(), 0);
    
    for (int y = 0; y < stream_grid.height(); y++) {
        for (int x = 0; x < stream_grid.width(); x++) {
            // Skip non-stream cells
            if (stream_grid(x, y) != 1 || watershed_mask(x, y) == 0) {
                continue;
            }
            
            // Check where this stream cell flows to
            if (flowdirs(x, y) == richdem::NO_FLOW) {
                continue;
            }
            
            int next_x = x + richdem::d8x[flowdirs(x, y)];
            int next_y = y + richdem::d8y[flowdirs(x, y)];
            
            // If it flows to another stream cell, increment that cell's inflow count
            if (flowdirs.inGrid(next_x, next_y) && stream_grid(next_x, next_y) == 1) {
                stream_inflow(next_x, next_y)++;
            }
        }
    }
    
    // Identify junction, endpoint, and outlet cells
    // Note: outlet_id starts at 1 because 0 is reserved for NoData/background
    int outlet_id = 1;  // Start at 1, not 0 (0 = NoData)
    int junctions = 0;
    int headwaters = 0;
    int watershed_outlets = 0;
    
    for (int y = 0; y < stream_grid.height(); y++) {
        for (int x = 0; x < stream_grid.width(); x++) {
            if (stream_grid(x, y) != 1 || watershed_mask(x, y) == 0) {
                continue;
            }
            
            // Check if this stream cell flows out of the watershed or has no flow
            bool is_watershed_outlet = false;
            if (flowdirs(x, y) == richdem::NO_FLOW) {
                is_watershed_outlet = true;
            } else {
                int next_x = x + richdem::d8x[flowdirs(x, y)];
                int next_y = y + richdem::d8y[flowdirs(x, y)];
                // Flows outside grid or outside watershed
                if (!flowdirs.inGrid(next_x, next_y) || watershed_mask(next_x, next_y) == 0) {
                    is_watershed_outlet = true;
                }
            }
            
            if (is_watershed_outlet) {
                // Watershed outlet: stream cell that flows out of watershed
                outlets.push_back(PourPoint(x, y, outlet_id, "Outlet_" + std::to_string(outlet_id), false));
                outlet_id++;
                watershed_outlets++;
            }
            else if (stream_inflow(x, y) >= 2) {
                // Junctions: cells with 2+ upstream stream cells
                outlets.push_back(PourPoint(x, y, outlet_id, "Junction_" + std::to_string(outlet_id), false));
                outlet_id++;
                junctions++;
            }
            else if (stream_inflow(x, y) == 0) {
                // Headwaters: stream cells with no upstream stream cells
                outlets.push_back(PourPoint(x, y, outlet_id, "Headwater_" + std::to_string(outlet_id), false));
                outlet_id++;
                headwaters++;
            }
        }
    }
    
    std::cout << "  Found " << watershed_outlets << " watershed outlets, " 
              << junctions << " junctions, and " 
              << headwaters << " headwater points" << std::endl;
    std::cout << "  Total subcatchment outlets: " << outlets.size() << std::endl;
    
    return outlets;
}

template<typename T>
int WatershedDelineation::identify_subcatchments(
    const richdem::Array2D<richdem::flowdir_t>& flowdirs,
    const richdem::Array2D<uint8_t>& stream_grid,
    const richdem::Array2D<T>& watershed_mask,
    richdem::Array2D<int>& subcatchment_label)
{
    std::cout << "\nIdentifying subcatchments..." << std::endl;
    
    // Initialize subcatchment labels to 0 (NoData/background)
    // Actual subcatchment IDs will be 1, 2, 3, etc.
    subcatchment_label.setAll(0);
    
    // Find stream junctions and endpoints
    auto outlets = identify_subcatchment_outlets(flowdirs, stream_grid, watershed_mask);
    
    if (outlets.empty()) {
        std::cout << "  No subcatchment outlets found" << std::endl;
        return 0;
    }
    
    // Delineate subcatchment for each outlet
    std::cout << "Delineating " << outlets.size() << " subcatchments..." << std::endl;
    
    for (const auto& outlet : outlets) {
        // Only delineate within the watershed mask
        richdem::Array2D<int> temp_label(subcatchment_label.width(), subcatchment_label.height(), 0);
        delineate_watershed(flowdirs, outlet, temp_label, 1);
        
        // Apply watershed mask and merge into subcatchment_label
        for (int y = 0; y < subcatchment_label.height(); y++) {
            for (int x = 0; x < subcatchment_label.width(); x++) {
                if (temp_label(x, y) == 1 && watershed_mask(x, y) != 0 && subcatchment_label(x, y) == 0) {
                    subcatchment_label(x, y) = outlet.id;
                }
            }
        }
    }
    
    // Count cells in each subcatchment and identify non-empty ones
    std::map<int, int> subcatchment_sizes;
    for (int y = 0; y < subcatchment_label.height(); y++) {
        for (int x = 0; x < subcatchment_label.width(); x++) {
            if (subcatchment_label(x, y) > 0) {
                subcatchment_sizes[subcatchment_label(x, y)]++;
            }
        }
    }
    
    // Report statistics
    int outlets_found = outlets.size();
    int subcatchments_with_cells = subcatchment_sizes.size();
    int empty_subcatchments = outlets_found - subcatchments_with_cells;
    
    std::cout << "\n=== Subcatchment Statistics ===" << std::endl;
    std::cout << "  Stream outlets identified: " << outlets_found << std::endl;
    std::cout << "  Subcatchments with cells:  " << subcatchments_with_cells << std::endl;
    // if (empty_subcatchments > 0) {
    //     std::cout << "  Empty subcatchments:       " << empty_subcatchments 
    //               << " (outlets claimed by downstream subcatchments)" << std::endl;
    // }
    
    // Return number of non-empty subcatchments
    return subcatchments_with_cells;
}

template<typename T>
WatershedStats WatershedDelineation::calculate_watershed_stats(
    const richdem::Array2D<T>& dem,
    const richdem::Array2D<int>& watershed_label,
    const richdem::Array2D<uint8_t>& stream_grid,
    int watershed_id)
{
    WatershedStats stats;
    stats.watershed_id = watershed_id;
    stats.min_elevation = std::numeric_limits<double>::max();
    stats.max_elevation = std::numeric_limits<double>::lowest();
    
    double sum_elevation = 0.0;
    int cell_count = 0;
    int stream_count = 0;
    
    // Calculate statistics
    for (int y = 0; y < watershed_label.height(); y++) {
        for (int x = 0; x < watershed_label.width(); x++) {
            if (watershed_label(x, y) == watershed_id) {
                cell_count++;
                double elev = static_cast<double>(dem(x, y));
                sum_elevation += elev;
                
                if (elev < stats.min_elevation) {
                    stats.min_elevation = elev;
                }
                if (elev > stats.max_elevation) {
                    stats.max_elevation = elev;
                }
                
                if (stream_grid(x, y) == 1) {
                    stream_count++;
                }
            }
        }
    }
    
    stats.area_cells = cell_count;
    stats.mean_elevation = (cell_count > 0) ? (sum_elevation / cell_count) : 0.0;
    stats.stream_cells = stream_count;
    
    // Calculate area in km2 (assumes geotransform is available)
    // This is a simplified calculation - actual implementation would use dem.geotransform
    double cell_size_x = 30.0;  // Default: 30m (will be updated from geotransform)
    double cell_size_y = 30.0;
    stats.area_km2 = (cell_count * cell_size_x * cell_size_y) / 1000000.0;
    
    // Estimate stream length (simplified)
    stats.total_stream_length = (stream_count * cell_size_x) / 1000.0;  // km
    
    return stats;
}

// ============================================================================
// NON-TEMPLATE IMPLEMENTATION
// ============================================================================

int WatershedDelineation::run_watershed_delineation(const std::vector<PourPoint>& pour_points_input)
{
    std::cout << "\n=== WATERSHED DELINEATION ===" << std::endl;
    std::cout << "Loading input data..." << std::endl;
    
    try {
        // Load DEM
        richdem::Array2D<float> dem(params.input_dem);
        std::cout << "  DEM: " << dem.width() << " x " << dem.height() << " cells" << std::endl;
        std::cout << "  Geotransform: [" << dem.geotransform[0] << ", " << dem.geotransform[1] 
                  << ", " << dem.geotransform[2] << ", " << dem.geotransform[3] 
                  << ", " << dem.geotransform[4] << ", " << dem.geotransform[5] << "]" << std::endl;
        
        // Load flow accumulation
        richdem::Array2D<int> accumulation(params.input_accumulation);
        std::cout << "  Accumulation: " << accumulation.width() << " x " 
                  << accumulation.height() << " cells" << std::endl;
        
        // Make a copy of pour points for conversion
        std::vector<PourPoint> pour_points = pour_points_input;
        
        // Convert geographic coordinates to pixel coordinates
        convert_pourpoints_to_pixels(pour_points, dem.geotransform, dem.width(), dem.height());
        
        // Snap pour points to maximum flow accumulation
        std::cout << "\nSnapping pour points to stream locations..." << std::endl;
        int snapped_count = 0;
        for (auto& pp : pour_points) {
            if (snap_pourpoint_to_stream(accumulation, pp, 10, 1.5)) {
                snapped_count++;
            }
        }
        if (snapped_count > 0) {
            std::cout << "Adjusted " << snapped_count << " pour point(s) to stream locations" << std::endl;
        } else {
            std::cout << "All pour points are already at optimal locations" << std::endl;
        }
        
        // Load and convert flow directions
        std::cout << "  Loading flow directions..." << std::endl;
        richdem::Array2D<int> flowdirs_orig(params.input_flowdirs);
        richdem::Array2D<richdem::flowdir_t> flowdirs(flowdirs_orig, richdem::NO_FLOW);
        
        // Convert flow directions from ESRI/TauDEM format to dephier format
        std::map<int, richdem::flowdir_t> dirs_map;
        constexpr std::array<int, 9> d8_esri = {0, 16, 32, 64, 128, 1, 2, 4, 8};
        constexpr std::array<richdem::flowdir_t, 9> d8_dephier = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        constexpr std::array<uint8_t, 9> d8_tau = {0, 5, 4, 3, 2, 1, 8, 7, 6};
        
        if (params.use_tau_format) {
            std::transform(d8_tau.begin(), d8_tau.end(), d8_dephier.begin(),
                         std::inserter(dirs_map, dirs_map.end()),
                         [](int const& x, richdem::flowdir_t i) {
                           return std::make_pair(x, i);
                         });
        } else {
            std::transform(d8_esri.begin(), d8_esri.end(), d8_dephier.begin(),
                         std::inserter(dirs_map, dirs_map.end()),
                         [](int const& x, richdem::flowdir_t i) {
                           return std::make_pair(x, i);
                         });
        }
        
        // Convert flow directions
        std::cout << "  Converting flow directions from " 
                  << (params.use_tau_format ? "TauDEM" : "ESRI/ArcGIS") 
                  << " format..." << std::endl;
        
        int converted_count = 0;
        int nodata_count = 0;
        for (int y = 0; y < flowdirs.height(); y++) {
            for (int x = 0; x < flowdirs.width(); x++) {
                if (!flowdirs_orig.isNoData(x, y)) {
                    try {
                        flowdirs(x, y) = dirs_map.at(flowdirs_orig(x, y));
                        converted_count++;
                    } catch (const std::out_of_range&) {
                        std::cerr << "WARNING: Invalid flow direction value " 
                                  << flowdirs_orig(x, y) << " at (" << x << ", " << y << ")" << std::endl;
                    }
                } else {
                    nodata_count++;
                }
            }
        }
        std::cout << "  Converted " << converted_count << " cells, " 
                  << nodata_count << " NoData cells" << std::endl;
        flowdirs_orig.clear();
        
        // Extract stream network
        std::cout << "\n=== STREAM NETWORK EXTRACTION ===" << std::endl;
        richdem::Array2D<uint8_t> stream_grid(dem.width(), dem.height(), 0);
        stream_grid.geotransform = dem.geotransform;
        stream_grid.projection = dem.projection;
        extract_stream_network(accumulation, stream_grid, 
                              static_cast<int>(params.stream_threshold));
        
        // Save stream network
        std::cout << "Saving stream network..." << std::endl;
        stream_grid.setNoData(0);  // Set NoData value to 0 (background)
        stream_grid.saveGDAL(params.output_prefix + "_streams.tif");
        
        // Calculate Horton-Strahler stream order
        std::cout << "\n=== STREAM ORDER CLASSIFICATION ===" << std::endl;
        richdem::Array2D<int> stream_order(dem.width(), dem.height(), 0);
        stream_order.geotransform = dem.geotransform;
        stream_order.projection = dem.projection;
        calculate_stream_order(flowdirs, stream_grid, stream_order);
        
        // Save stream order
        std::cout << "Saving stream order..." << std::endl;
        stream_order.setNoData(0);  // Set NoData value to 0 (non-stream)
        stream_order.saveGDAL(params.output_prefix + "_stream_order.tif");
        
        // Delineate watersheds
        std::cout << "\n=== WATERSHED DELINEATION ===" << std::endl;
        richdem::Array2D<int> watershed_label(dem.width(), dem.height(), 0);
        watershed_label.geotransform = dem.geotransform;
        watershed_label.projection = dem.projection;
        
        int num_watersheds = delineate_multiple_watersheds(flowdirs, pour_points, watershed_label);
        
        // Save watershed labels
        std::cout << "Saving watershed labels..." << std::endl;
        watershed_label.setNoData(0);  // Set NoData value to 0 (background)
        watershed_label.saveGDAL(params.output_prefix + "_watersheds.tif");
        
        // Identify and delineate subcatchments
        std::cout << "\n=== SUBCATCHMENT DELINEATION ===" << std::endl;
        richdem::Array2D<int> subcatchment_label(dem.width(), dem.height(), 0);
        subcatchment_label.geotransform = dem.geotransform;
        subcatchment_label.projection = dem.projection;
        
        int num_subcatchments = identify_subcatchments(flowdirs, stream_grid, watershed_label, subcatchment_label);
        
        if (num_subcatchments > 0) {
            std::cout << "Saving subcatchment labels..." << std::endl;
            subcatchment_label.setNoData(0);  // Set NoData value to 0 (background)
            subcatchment_label.saveGDAL(params.output_prefix + "_subcatchments.tif");
        }
        
        // Calculate and print statistics for each watershed
        std::cout << "\n=== WATERSHED STATISTICS ===" << std::endl;
        for (const auto& pp : pour_points) {
            auto stats = calculate_watershed_stats(dem, watershed_label, stream_grid, pp.id);
            
            // std::cout << "\nWatershed " << stats.watershed_id;
            // if (!pp.name.empty()) {
            //     std::cout << " (" << pp.name << ")";
            // }
            // std::cout << ":" << std::endl;
            // std::cout << "  Area: " << stats.area_cells << " cells (" 
            //          << stats.area_km2 << " km²)" << std::endl;
            // std::cout << "  Elevation: min=" << stats.min_elevation 
            //          << " m, max=" << stats.max_elevation 
            //          << " m, mean=" << stats.mean_elevation << " m" << std::endl;
            // std::cout << "  Stream cells: " << stats.stream_cells << std::endl;
            // std::cout << "  Stream length: ~" << stats.total_stream_length << " km" << std::endl;
        }
        
        std::cout << "\n=== WATERSHED DELINEATION COMPLETED ===" << std::endl;
        std::cout << "Output files created:" << std::endl;
        std::cout << "  " << params.output_prefix << "_streams.tif - Stream network (binary)" << std::endl;
        std::cout << "  " << params.output_prefix << "_stream_order.tif - Horton-Strahler stream order" << std::endl;
        std::cout << "  " << params.output_prefix << "_watersheds.tif - Watershed labels" << std::endl;
        if (num_subcatchments > 0) {
            std::cout << "  " << params.output_prefix << "_subcatchments.tif - Subcatchment labels (" 
                      << num_subcatchments << " non-empty subcatchments)" << std::endl;
        }
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return -1;
    }
}


