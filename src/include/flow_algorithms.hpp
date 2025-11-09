#include "stdefx.h"

#include <dephier/dephier.hpp>
#include <richdem/common/Array2D.hpp>
#include <richdem/common/constants.hpp>



struct Params {
    std::string input_file;
    std::string output_prefix;
    float ocean_level;
    bool write_extra_files = false;
    bool modify_mode = false;
    bool use_tau_format = false;
};

class FlowAlgorithms{

    // This class contains a series of custom flow
    // algorithms that uses the results from the Depression Hierarchy framework
    // to produce consistent flow accumulations and flow directions

    public:
        // Constructor and destructor with default implementations
        FlowAlgorithms(const Params &params) : params(params) {}
        ~FlowAlgorithms() = default;

        /**
         * @brief This functions performs pre-treatment of the input DEM by identifying and labeling border pixels as ocean cells
         * 
         * @tparam T: Type of the input DEM raster.
         * @tparam V: Type of the ocean level.
         * @param topo: Input DEM raster.
         * @param label: Depression hierarchy labels.
         * @param ocean_level: Ocean level.
         */
        template<typename A, typename B>
        void initialize_ocean_labels(const richdem::Array2D<A> &topo, richdem::Array2D<B> &label);

         /**
         * @brief Computes the carving of the input DEM raster based on the depression hierarchy
         * 
         * @tparam T: Type of the input DEM raster.
         * @param topo:  Input DEM raster. 
         * @param label: Depression hierarchy labels
         * @param flowdirs : Flow directions from the depression hierarchy
         */
        template<typename T>
        void compute_carving(richdem::Array2D<T> &topo,
            richdem::Array2D<richdem::dephier::dh_label_t> &label,
            richdem::Array2D<richdem::flowdir_t> &flowdirs);

        int compute_flow_accumulation(const std::string &input_flow_dirs);

    private:
        Params params;  // Configuration parameters for the algorithms
        /**
         * @brief Performs the carving
         * 
         */
        template<typename T>
        void perform_carving(richdem::Array2D<T> &topo, const richdem::Array2D<richdem::flowdir_t> &flowdirs,
            const T max_elevation, richdem::Array2D<T> &accu, richdem::Array2D<richdem::dephier::dh_label_t> &nips);

        template <typename T>
        void perform_modified_carving(const richdem::Array2D<T> &topo, richdem::Array2D<richdem::flowdir_t> &flowdirs,
            richdem::Array2D<richdem::dephier::dh_label_t> &label, const T max_elevation,
            richdem::Array2D<T> &accu, richdem::Array2D<richdem::dephier::dh_label_t> &nips);

        /**
          * @brief Computes the number of inflow pixels (NIPS) for carving propagation
          * 
          * This function prepares the NIPS grid needed for the carving algorithm to
          * properly handle flow confluences when propagating minimum pit elevations
          * downstream through the depression hierarchy.
          * 
          * @tparam T: Type of the input DEM raster.
          * @param topo: Input DEM raster.
          * @param flowdirs: Flow directions from the depression hierarchy.
          * @param max_elevation: Maximum elevation for initialization.
          * @param accu: Carving accumulation grid (not modified by this function).
          * @param nips: Output - Number of inflow pixels for each cell (modified).
          */
        template <typename T>
        void compute_carving_accumulation(const richdem::Array2D<T> &topo, const richdem::Array2D<richdem::flowdir_t> &flowdirs,
             const T max_elevation, richdem::Array2D<T> &accu, richdem::Array2D<richdem::dephier::dh_label_t> &nips);
         

        /**
         * @brief Compute outlets and pit elevations
         * 
         * @tparam T: Type of the input DEM raster.
         * @param topo: Input DEM raster.
         * @param label: Depression hierarchy labels.
         * @param flowdirs: Flow directions from the depression hierarchy.
         */
        template<class T>
        void compute_outlet_and_pit_elevations(const richdem::dephier::DepressionHierarchy<T> depression_hierarchy,
            richdem::Array2D<T> &topo, const T max_elevation, std::vector<bool> &check_value, std::vector<T> &pit_min);
        
        template<class T>
        void find_inlet_and_outlet_cells(const richdem::dephier::DepressionHierarchy<T> depression_hierarchy,
            richdem::Array2D<T> &topo, richdem::Array2D<richdem::dephier::dh_label_t> &label,
            richdem::Array2D<T> &outcell_pitelev, richdem::Array2D<T> &incell_pitelev,
            const T max_elevation, std::vector<bool> &check_value, std::vector<T> &pit_min);

};


// Template function implementations (must be in header for templates)
template<typename A, typename B>
void FlowAlgorithms::initialize_ocean_labels(const richdem::Array2D<A> &topo, 
    richdem::Array2D<B> &label){
    std::cout << "Initializing ocean labels..." << std::endl;

    // Start by labeling ocean cells and nodata cells
    #pragma omp parallel for
    for (size_t i = 0; i <label.size(); i++){
        if (topo.isNoData(i) || topo(i) == params.ocean_level){
            label(i) = richdem::dephier::OCEAN;
        }
    }

    // Label edge pixels as ocean
    #pragma omp parallel for
    for (int y = 0; y < topo.height(); y++){
        label(topo.xyToI(0, y)) = richdem::dephier::OCEAN;
        label(topo.xyToI(topo.width() - 1, y)) = richdem::dephier::OCEAN;
    }
    #pragma omp parallel for
    for (int x = 0; x < topo.width(); x++){
        label(topo.xyToI(x, 0)) = richdem::dephier::OCEAN;
        label(topo.xyToI(x, topo.height() - 1)) = richdem::dephier::OCEAN;
    }
}

int FlowAlgorithms::compute_flow_accumulation(const std::string &input_flowdirs){
    std::cout << "Running FLOW ACCUMULATION algorithm" << std::endl;
    
    // Build flow direction format conversion map
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

    // Load flow directions from file
    std::cout << "Loading flow directions..." << std::endl;
    richdem::Array2D<int> flowdirs_orig(input_flowdirs);
    richdem::Array2D<richdem::flowdir_t> flowdirs(flowdirs_orig, richdem::NO_FLOW);
    
    std::cout << "Data width  = " << flowdirs_orig.width() << std::endl;
    std::cout << "Data height = " << flowdirs_orig.height() << std::endl;
    std::cout << "Data cells  = " << flowdirs_orig.numDataCells() << std::endl;

    // Convert flow directions from input format to dephier format
    auto vals_are_ok = true;
    #pragma omp parallel for collapse(2)
    for (int y = 0; y < flowdirs.height(); y++) {
        for (int x = 0; x < flowdirs.width(); x++) {
            if (!flowdirs_orig.isNoData(x, y)) {
                auto dir_in = flowdirs_orig(x, y);
                try {
                    flowdirs(x, y) = dirs_map.at(dir_in);
                } catch (const std::out_of_range&) {
                    if (vals_are_ok) {
                        std::cout << "ERROR: direction " << dir_in << " is not valid. "
                                << "Check that inputs are " 
                                << (params.use_tau_format ? "TauDEM" : "ESRI") 
                                << " flow directions" << std::endl;
                        vals_are_ok = false;
                    }
                }
            }
        }
    }

    if (!vals_are_ok) {
        return -1;
    }
    
    dirs_map.clear();

    // Calculate number of inflowing pixels (NIPS) for each cell
    std::cout << "Calculating inflowing pixels..." << std::endl;
    richdem::Array2D<uint16_t> nips(flowdirs_orig, 0);
    flowdirs_orig.clear();  // Free memory

    for (int y = 0; y < nips.height(); y++) {
        for (int x = 0; x < nips.width(); x++) {
            if (flowdirs(x, y) == richdem::NO_FLOW) {
                continue;
            }
            int my_nx = x + richdem::d8x[flowdirs(x, y)];
            int my_ny = y + richdem::d8y[flowdirs(x, y)];
            if (flowdirs.inGrid(my_nx, my_ny)) {
                nips(my_nx, my_ny)++;
            }
        }
    }

    // Initialize accumulation grid (cell count)
    std::cout << "Starting standard accumulation (cell count)..." << std::endl;
    richdem::Array2D<int> accu(input_flowdirs);
    accu.setAll(1);

    // Perform flow accumulation
    for (int y = 0; y < flowdirs.height(); y++) {
        for (int x = 0; x < flowdirs.width(); x++) {
            if (flowdirs(x, y) == richdem::NO_FLOW) {
                continue;
            }
            if (nips(x, y) != 0) {
                continue;
            }
            
            auto current_cell_accu = accu(x, y);
            int next_x = x + richdem::d8x[flowdirs(x, y)];
            int next_y = y + richdem::d8y[flowdirs(x, y)];
            
            while (flowdirs.inGrid(next_x, next_y)) {
                accu(next_x, next_y) = current_cell_accu + accu(next_x, next_y);
                if (nips(next_x, next_y) > 1) {
                    nips(next_x, next_y)--;
                    break;
                }
                if (flowdirs(next_x, next_y) == 0) {
                    break;
                }
                current_cell_accu = accu(next_x, next_y);
                int dx = richdem::d8x[flowdirs(next_x, next_y)];
                int dy = richdem::d8y[flowdirs(next_x, next_y)];
                next_x += dx;
                next_y += dy;
            }
        }
    }

    // Save result
    std::cout << "Saving accumulation result..." << std::endl;
    accu.saveGDAL(params.output_prefix + "_accu.tif");
    
    nips.clear();
    
    std::cout << "ACCU algorithm completed successfully!" << std::endl;
    
    return 0;
}

template<typename T>
void FlowAlgorithms::compute_carving(richdem::Array2D<T> &topo,
    richdem::Array2D<richdem::dephier::dh_label_t> &label,
    richdem::Array2D<richdem::flowdir_t> &flowdirs){
    
    std::cout << "Generating depression hierarchy... " << std::endl;
    // Generate the depression hierarchy
    richdem::dephier::DepressionHierarchy<T> depression_hierarchy = richdem::dephier::GetDepressionHierarchy<T, richdem::Topology::D8>(topo, label, flowdirs);
    const T max_elevation = topo.max() + 100;
    
    // Allocate the elevation grids
    richdem::Array2D<T> outcell_pitelev(topo);
    richdem::Array2D<T> incell_pitelev(topo, max_elevation);

    // Initialize outcell_pitelev grids
    #pragma omp parallel for collapse(2)
    for (int y = 0; y < topo.height(); y++){
        for (int x = 0; x < topo.width(); x++){
            outcell_pitelev(x,y) = max_elevation;
        }
    }

    // Initialize tracking vectors:
    // - check_value: flags which depressions need to be re-examined in next iteration
    // - pit_min: stores the minimum pit elevation found for each depression (including upstream)
    std::vector<bool> check_value(depression_hierarchy.size(), true);
    std::vector<T> pit_min(depression_hierarchy.size(), max_elevation);
    // Compute outlets and pit elevations
    compute_outlet_and_pit_elevations(depression_hierarchy, topo, max_elevation, check_value, pit_min);
    // Find inlet and outlet cells
    find_inlet_and_outlet_cells(depression_hierarchy, topo, label, outcell_pitelev, incell_pitelev, max_elevation, check_value, pit_min);
 
    // ==============================================================================
    // COMBINE INLET AND OUTLET ELEVATIONS TO CREATE CARVING SURFACE
    // ==============================================================================
    // PURPOSE:
    //   Merge the inlet and outlet elevation grids to create a unified carving 
    //   surface (accu). This surface represents the target elevations for carving
    //   that ensure smooth water flow through the depression hierarchy.
    //
    // PROCESS:
    //   For each cell, take the MINIMUM of:
    //   - outcell_pitelev: elevation at outlet cells (where water exits depressions)
    //   - incell_pitelev: elevation at inlet cells (where water enters depressions)
    //
    //   This ensures the carved surface has no barriers to flow at either entry
    //   or exit points of depressions.
    // ==============================================================================
    
    std::cout << "Combining inlet and outlet elevations..." << std::endl;
    
    // Create an alias for outcell_pitelev called accu (accumulation/carving surface)
    // This is a reference, so any changes to accu will modify outcell_pitelev
    richdem::Array2D<T> &accu = outcell_pitelev;
 
    // Merge inlet elevations into the carving surface
    for (unsigned y = 0; y < label.height(); y++){
        for (unsigned x = 0; x < label.width(); x++){
            // If inlet cell has valid elevation and is lower than current value
            if (incell_pitelev(x, y) != max_elevation && incell_pitelev(x, y) < accu(x, y)){
                // Use the lower inlet elevation for smoother carving
                accu(x, y) = incell_pitelev(x, y);
            }
        }
    }
    
    // Release memory from incell_pitelev, since it is no longer needed
    incell_pitelev.clear();
    
    // Initialize grid to track number of inflow pixels (nips) for carving propagation
    richdem::Array2D<richdem::dephier::dh_label_t> nips(topo, 0);
    
    // Compute NIPS and prepare for carving accumulation
    compute_carving_accumulation(topo, flowdirs, max_elevation, accu, nips);
 
    // Perform the carving on the DEM
    perform_carving(topo, flowdirs, max_elevation, accu, nips);

    // Save the carved DEM
    if (!params.modify_mode) {
        // Standard mode: save carved DEM only
        std::cout << "Saving carved DEM..." << std::endl;
        accu.geotransform = topo.geotransform;
        if (params.write_extra_files) {
          accu.saveGDAL(params.output_prefix + "-carve-min.tif");
        }
        topo.saveGDAL(params.output_prefix + "-dem-carved.tif");
    } else{
        // Modified mode: regenerate flow directions and save both DEM and directions
        std::cout << "Generating flow directions on carved DEM..." << std::endl;
        perform_modified_carving(topo, flowdirs,label, max_elevation, accu, nips);
        
        // Save the carved DEM
        std::cout << "Saving carved DEM..." << std::endl;
        topo.geotransform = accu.geotransform;
        topo.saveGDAL(params.output_prefix + "-mod-dem-carved.tif");
        
        // Convert directions to ArcGIS format
        auto dirs_arcgis = label;
        #pragma omp parallel for collapse(2)
        for (int y = 0; y < dirs_arcgis.height(); y++) {
          for (int x = 0; x < dirs_arcgis.width(); x++) {
            dirs_arcgis(x, y) = richdem::d8_arcgis[flowdirs(x, y)];
          }
        }

        std::cout << "Saving flow directions..." << std::endl;
        dirs_arcgis.geotransform = topo.geotransform;
        dirs_arcgis.saveGDAL(params.output_prefix + "-mod-dirs.tif");
        dirs_arcgis.clear();

        // Save flow directions (convert to int16 for GDAL compatibility)
        // std::cout << "Saving flow directions..." << std::endl;
        // richdem::Array2D<int16_t> flowdirs_int(flowdirs.width(), flowdirs.height());
        // flowdirs_int.geotransform = topo.geotransform;
        // flowdirs_int.projection = topo.projection;
        // for (unsigned i = 0; i < flowdirs.size(); i++){
        //     flowdirs_int(i) = static_cast<int16_t>(flowdirs(i));
        // }
        // // Convert
        // flowdirs_int.saveGDAL(params.output_prefix + "-mod-flowdirs.tif");
        
        std::cout << "Modified carving complete!" << std::endl;
    }
}

template <typename T>
void FlowAlgorithms::perform_modified_carving(const richdem::Array2D<T> &topo, richdem::Array2D<richdem::flowdir_t> &flowdirs,
    richdem::Array2D<richdem::dephier::dh_label_t> &label, const T max_elevation,
    richdem::Array2D<T> &accu, richdem::Array2D<richdem::dephier::dh_label_t> &nips){
    
    // ==============================================================================
    // REGENERATE DEPRESSION HIERARCHY ON CARVED DEM
    // ==============================================================================
    // PURPOSE:
    //   After carving the DEM, regenerate the depression hierarchy and flow 
    //   directions on the carved surface to get accurate flow routing.
    //
    // NOTE:
    //   - topo has already been carved (modified in place by perform_carving)
    //   - accu contains the carved elevations (passed from compute_carving)
    //   - We need to regenerate flowdirs and labels on the carved terrain
    // ==============================================================================
    
    std::cout << "Regenerating depression hierarchy on carved DEM..." << std::endl;
    
    // Reset flow directions
    #pragma omp parallel for collapse(2)
    for (int y = 0; y < flowdirs.height(); y++) {
    for (int x = 0; x < flowdirs.width(); x++) {
        flowdirs(x, y) = richdem::NO_FLOW;
    }
    }

    // Reset labels to NO_DEP (required by GetDepressionHierarchy)
    std::cout << "Resetting labels..." << std::endl;
    label.setAll(richdem::dephier::NO_DEP);
    
    // Reinitialize ocean labels on the carved DEM
    initialize_ocean_labels(topo, label);
    
    // Regenerate depression hierarchy using the carved DEM
    richdem::dephier::GetDepressionHierarchy<T, richdem::Topology::D8>(topo, label, flowdirs);

    // Recompute NIPS for the new flow directions
    std::cout << "Recomputing NIPS for carved DEM..." << std::endl;
    nips.setAll(0);
    for (int y = 0; y < topo.height(); y++){
        for (int x = 0; x < topo.width(); x++){
            if (flowdirs(x, y) == richdem::NO_FLOW){
                continue;
            }
            int next_x = x + richdem::d8x[flowdirs(x, y)];
            int next_y = y + richdem::d8y[flowdirs(x, y)];
            if (topo.inGrid(next_x, next_y)){
                nips(next_x, next_y)++;
            }
        }
    }

    // Propagate minimum elevations using the new flow directions
    std::cout << "Propagating elevations with new flow directions..." << std::endl;
    for (int y = 0; y < topo.height(); y++){
        for (int x = 0; x < topo.width(); x++){
            if (flowdirs(x,y) == richdem::NO_FLOW){
                continue;
            }
            if (nips(x,y) != 0){
                continue;
            }

            T current_cell_min = accu(x,y);
            int next_x = x + richdem::d8x[flowdirs(x, y)];
            int next_y = y + richdem::d8y[flowdirs(x, y)];

            while (topo.inGrid(next_x, next_y)){
                accu(next_x, next_y) = std::min(current_cell_min, accu(next_x, next_y));
                if (nips(next_x, next_y) > 1){
                    nips(next_x, next_y)--;
                    break;
                }
                if (flowdirs(next_x, next_y) == 0){
                    break;
                }

                current_cell_min = accu(next_x, next_y);
                int dx = richdem::d8x[flowdirs(next_x, next_y)];
                int dy = richdem::d8y[flowdirs(next_x, next_y)];
                next_x += dx;
                next_y += dy;
            }
        }
    }
}


template <typename T>
void FlowAlgorithms::perform_carving(richdem::Array2D<T> &topo, const richdem::Array2D<richdem::flowdir_t> &flowdirs,
    const T max_elevation, richdem::Array2D<T> &accu, richdem::Array2D<richdem::dephier::dh_label_t> &nips){
    
    // ==============================================================================
    // PROPAGATE MINIMUM CARVING ELEVATIONS DOWNSTREAM
    // ==============================================================================
    // PURPOSE:
    //   Propagate the minimum pit elevations (stored in accu) downstream along
    //   flow paths. This ensures that each cell along a flow path knows the 
    //   minimum elevation it needs to be carved to for water to flow from pits
    //   through the depression hierarchy to the ocean.
    //
    // ALGORITHM:
    //   1. Start from flow sources (cells with nips=0 - hilltops/ridges)
    //   2. Trace downstream following flow directions
    //   3. At each step, propagate the minimum of current accu and upstream accu
    //   4. Handle confluences (cells with multiple inflows) by decrementing nips
    //   5. Stop when hitting cells that still need upstream contributions
    //
    // NIPS USAGE:
    //   - nips=0: Flow source, start tracing here
    //   - nips=1: Single inflow, continue downstream propagation
    //   - nips>1: Multiple inflows, wait for all contributors (decrement & stop)
    //
    // RESULT:
    //   After this step, accu[x,y] contains the minimum elevation needed at (x,y)
    //   to ensure continuous flow from all upstream pits to the ocean.
    // ==============================================================================
    
    std::cout << "Propagating carving elevations downstream..." << std::endl;
    
    // ----------------------------------------------------------------------
    // STEP 1: Trace flow paths from all source cells (nips=0)
    // ----------------------------------------------------------------------
    
    for (int y = 0; y < topo.height(); y++){
        for (int x = 0; x < topo.width(); x++){
            
            // Skip cells with no flow direction (sinks, edge cells, nodata)
            if (flowdirs(x, y) == richdem::NO_FLOW){
                continue;
            }
            
            // Skip cells that still need upstream contributions (nips > 0)
            // Only start tracing from flow sources (nips == 0)
            if (nips(x, y) != 0){
                continue;
            }
            
            // ------------------------------------------------------------------
            // STEP 2: Initialize downstream trace from this source cell
            // ------------------------------------------------------------------
            
            // Current minimum elevation to propagate downstream
            T current_cell_accu = accu(x, y);
            
            // Calculate first downstream neighbor
            int next_x = x + richdem::d8x[flowdirs(x, y)];
            int next_y = y + richdem::d8y[flowdirs(x, y)];
            
            // ------------------------------------------------------------------
            // STEP 3: Trace downstream until hitting a stopping condition
            // ------------------------------------------------------------------
            
            while (topo.inGrid(next_x, next_y)){
                
                // Update downstream cell with minimum elevation
                // (propagate the lowest pit elevation needed for flow)
                accu(next_x, next_y) = std::min(current_cell_accu, accu(next_x, next_y));
                
                // Check if downstream cell has multiple inflows (confluence)
                if (nips(next_x, next_y) > 1){
                    // Decrement inflow count - this path has contributed
                    nips(next_x, next_y)--;
                    // Stop here - cell still needs other upstream contributors
                    break;
                }
                
                // Check if downstream cell is a sink or has no flow direction
                if (flowdirs(next_x, next_y) == 0){
                    break;
                }
                
                // Update current minimum for next iteration
                current_cell_accu = accu(next_x, next_y);
                
                // Move to next downstream cell
                const int dx = richdem::d8x[flowdirs(next_x, next_y)];
                const int dy = richdem::d8y[flowdirs(next_x, next_y)];
                next_x += dx;
                next_y += dy;
            }
        }
    }
    
    // ==============================================================================
    // APPLY CARVING TO THE DEM
    // ==============================================================================
    // PURPOSE:
    //   Modify the DEM elevations by carving down to the computed minimum
    //   elevations (accu). This removes barriers to flow and creates continuous
    //   flow paths through the depression hierarchy.
    //
    // PROCESS:
    //   1. For each cell, take minimum of original DEM and carving surface
    //   2. Update accu grid to store final carved elevations (for output)
    //   3. Set cells at max_elevation to 0 (unprocessed/edge cells)
    // ==============================================================================
    
    std::cout << "Applying carving to DEM..." << std::endl;
    
    #pragma omp parallel for collapse(2)
    for (unsigned y = 0; y < accu.height(); y++){
        for (unsigned x = 0; x < accu.width(); x++){
            
            // Carve DEM: take minimum of original elevation and carving target
            topo(x, y) = std::min(topo(x, y), accu(x, y));
            
            // Update accu grid for output
            if (accu(x, y) == max_elevation){
                // Unprocessed cells (no carving needed) - set to 0
                accu(x, y) = 0;
            } else {
                // Store the final carved elevation
                accu(x, y) = topo(x, y);
            }
        }
    }
    
    std::cout << "Carving complete!" << std::endl;
}

template <typename T>
void FlowAlgorithms::compute_carving_accumulation(const richdem::Array2D<T> &topo, const richdem::Array2D<richdem::flowdir_t> &flowdirs,
    const T max_elevation, richdem::Array2D<T> &accu, richdem::Array2D<richdem::dephier::dh_label_t> &nips){
    
    // ==============================================================================
    // CALCULATE NUMBER OF INFLOW PIXELS (NIPS) FOR CARVING PROPAGATION
    // ==============================================================================
    // PURPOSE:
    //   Compute how many upstream cells flow INTO each cell. This is essential for
    //   the carving propagation algorithm, which must trace flow paths while 
    //   handling confluences (stream junctions) correctly.
    //
    // CONTEXT:
    //   This function specifically prepares the NIPS grid for the carving algorithm
    //   (perform_carving), NOT for general flow accumulation calculations. The NIPS
    //   grid tells us which cells are flow sources vs. confluences, enabling proper
    //   downstream propagation of minimum pit elevations.
    //
    // ALGORITHM:
    //   1. Traverse all cells in the DEM
    //   2. For each cell, follow its flow direction to find downstream neighbor
    //   3. Increment the neighbor's inflow count (nips)
    //   4. Result: nips[x,y] = number of cells that flow into (x,y)
    //
    // USAGE IN CARVING:
    //   - nips=0: Flow sources (hilltops/ridges) - start carving propagation here
    //   - nips=1: Single inflow - continue propagating minimum elevations
    //   - nips>1: Confluences - must wait for all upstream paths before continuing
    //
    // NOTE:
    //   This is adapted from D8 flow accumulation methodology but used specifically
    //   for coordinating the downstream propagation of carving elevations.
    // ==============================================================================
    
    std::cout << "Calculating number of inflow pixels..." << std::endl;
    
    // Traverse all cells to count inflows
    for (unsigned y = 0; y < topo.height(); y++){
        for (unsigned x = 0; x < topo.width(); x++){
            
            // Skip cells with no flow direction (sinks, nodata, etc.)
            if (flowdirs(x, y) == richdem::NO_FLOW){
                continue;
            }
            
            // Calculate downstream neighbor coordinates using D8 flow direction
            // flowdirs(x,y) is a value 1-8 indicating one of 8 cardinal/diagonal directions
            const int neighbor_x = x + richdem::d8x[flowdirs(x, y)];
            const int neighbor_y = y + richdem::d8y[flowdirs(x, y)];
            
            // Increment inflow count for downstream neighbor (if in grid bounds)
            if (topo.inGrid(neighbor_x, neighbor_y)){
                nips(neighbor_x, neighbor_y)++;
            }
        }
    }
    
    std::cout << "Number of inflow pixels calculated." << std::endl;
}

template<typename T>
void FlowAlgorithms::find_inlet_and_outlet_cells(const richdem::dephier::DepressionHierarchy<T> depression_hierarchy,
    richdem::Array2D<T> &topo, richdem::Array2D<richdem::dephier::dh_label_t> &label,
    richdem::Array2D<T> &outcell_pitelev, richdem::Array2D<T> &incell_pitelev,
    const T max_elevation, std::vector<bool> &check_value, std::vector<T> &pit_min){
    
    // ==============================================================================
    // IDENTIFY INLET AND OUTLET CELLS FOR CARVING
    // ==============================================================================
    // PURPOSE:
    //   For each depression, identify two critical cell types:
    //   1. OUTLET CELL: Where water exits the depression (spill point)
    //   2. INLET CELL: Where water enters from downstream child depressions
    //   
    //   Both cells are assigned the minimum pit elevation to create smooth carving
    //   paths through the DEM that respect the depression hierarchy.
    //
    // CONTEXT:
    //   After this function, the grids outcell_pitelev and incell_pitelev will be
    //   combined to create the final carving surface, ensuring water flows 
    //   continuously from pits → inlets → outlets → ocean.
    // ==============================================================================
    
    std::cout << "Finding inlet and outlet cells..." << std::endl;
    
    // Process each depression in the hierarchy (skip ocean at index 0)
    for (unsigned x = 1; x < depression_hierarchy.size(); x++){
        
        // ----------------------------------------------------------------------
        // STEP 1: Mark the outlet cell with minimum pit elevation
        // ----------------------------------------------------------------------
        // The outlet (spill point) is where this depression overflows into 
        // its parent depression. We mark it with the minimum pit elevation
        // found among all depressions that drain through this point.
        
        // Get outlet cell coordinates for this depression
        int out_x, out_y;
        richdem::dephier::flat_c_idx outlet_now = depression_hierarchy[x].out_cell;
        label.iToxy(outlet_now, out_x, out_y);
        
        // Assign minimum pit elevation to outlet cell (if lower than current)
        if (outcell_pitelev(out_x, out_y) > pit_min[x]) {
            outcell_pitelev(out_x, out_y) = pit_min[x];
        }
        
        // ----------------------------------------------------------------------
        // STEP 2: Determine which depression to search for inlet cell
        // ----------------------------------------------------------------------
        // We need to find the inlet cell - the cell where water enters this
        // depression from downstream. This is a neighbor of the outlet that
        // belongs to a child depression.
        //
        // Special case: Geolinks represent lateral merges between depressions.
        // If the outlet cell IS the geolink, we set target_label = -1 to use
        // a different search strategy.
        
        int target_label;
        if (label(out_x, out_y) == depression_hierarchy[x].geolink) {
            // Outlet is at geolink - use alternative search (current depression)
            target_label = -1;
        } else {
            // Normal case - search for geolinked depression
            target_label = depression_hierarchy[x].geolink;
        }
        
        // Track the best inlet candidate (lowest elevation neighbor)
        T min_elevation = max_elevation + 1;
        int keepx = -1;  // Inlet cell x coordinate
        int keepy = -1;  // Inlet cell y coordinate
        
        // ----------------------------------------------------------------------
        // STEP 3: Search 8 neighbors of outlet to find inlet cell
        // ----------------------------------------------------------------------
        // The inlet must be a neighbor of the outlet that belongs to a 
        // depression draining into the current one.
        
        for (int n = 1; n <= 8; n++) {
            int target_label_now = target_label;
            
            // Get neighbor coordinates using D8 directions
            const int my_nx = out_x + richdem::d8x[n];
            const int my_ny = out_y + richdem::d8y[n];
            
            // Skip if neighbor is outside grid bounds
            if (!topo.inGrid(my_nx, my_ny)) {
                continue;
            }
            
            // ------------------------------------------------------------------
            // STEP 3a: Trace neighbor's depression hierarchy
            // ------------------------------------------------------------------
            // Check if this neighbor's depression eventually drains to the
            // current depression by traversing up the parent chain.
            
            richdem::dephier::dh_label_t neighbour_label = label(my_nx, my_ny);
            
            // Traverse parent chain to see if it connects to current depression
            while (neighbour_label != richdem::dephier::OCEAN && target_label_now == -1) {
                if (neighbour_label == x) {
                    // Found it! This neighbor drains into current depression
                    target_label_now = neighbour_label;
                } else {
                    // Check if we've reached ocean parent (stop condition)
                    if (depression_hierarchy[neighbour_label].ocean_parent) {
                        break;
                    }
                    // Move up to parent depression
                    neighbour_label = depression_hierarchy[neighbour_label].parent;
                }
            }
            
            // ------------------------------------------------------------------
            // STEP 3b: Keep neighbor with lowest elevation that connects
            // ------------------------------------------------------------------
            // If this neighbor's depression chain connects to current depression,
            // and it has the lowest elevation so far, mark it as inlet candidate.
            
            if (neighbour_label == target_label_now) {
                if (topo(my_nx, my_ny) < min_elevation) {
                    keepx = my_nx;
                    keepy = my_ny;
                    min_elevation = topo(my_nx, my_ny);
                }
            }
        }
        
        // ----------------------------------------------------------------------
        // STEP 4: Assign minimum pit elevation to inlet cell
        // ----------------------------------------------------------------------
        // If we found a valid inlet cell, mark it with the minimum pit elevation.
        // This ensures smooth transition from child depression into parent.
        
        if (keepx >= 0) {
            if (incell_pitelev(keepx, keepy) > pit_min[x]) {
                incell_pitelev(keepx, keepy) = pit_min[x];
            }
        }
    }
}

template<typename T>
void FlowAlgorithms::compute_outlet_and_pit_elevations(const richdem::dephier::DepressionHierarchy<T> depression_hierarchy,
    richdem::Array2D<T> &topo,const T max_elevation, std::vector<bool> &check_value, std::vector<T> &pit_min){

    // ==============================================================================
    // COMPUTE MINIMUM PIT ELEVATIONS FOR ALL DEPRESSIONS
    // ==============================================================================
    // PURPOSE: 
    //   For each depression in the hierarchy, find the LOWEST pit elevation of any
    //   depression that eventually drains into it (including itself and all upstream
    //   child depressions). This minimum pit elevation will be used later to set 
    //   the carving target elevation for that depression.
    //
    // ALGORITHM:
    //   Iteratively propagate minimum pit values upstream through the depression 
    //   hierarchy (from children to parents) until convergence. The algorithm 
    //   handles both parent-child relationships (vertical hierarchy) and geolinks
    //   (lateral connections between depressions at similar levels).
    //
    // INPUT:
    //   - depression_hierarchy: The computed depression hierarchy structure
    //   - topo: DEM raster with elevation values
    //   - max_elevation: Maximum elevation for initialization
    //   - check_value: Flag array (modified by this function)
    //   - pit_min: Output array to store minimum pit elevations (modified)
    //
    // OUTPUT:
    //   - pit_min[i] = lowest pit elevation among depression i and all its children
    // ==============================================================================
    
    std::cout << "Computing minimum pit elevations..." << std::endl;
    
    int round = 1;
    while (true){
        // ----------------------------------------------------------------------
        // Setup for this iteration round
        // ----------------------------------------------------------------------
        // Save which depressions were flagged from previous round (geolinks)
        std::vector<bool> check_value_from_geolink = check_value;
        // Reset flags - only depressions that get updated will be flagged
        std::fill(check_value.begin(), check_value.end(), false);

        // Progress reporting for long-running computations
        // if ((round % 100 == 0) || (round % 10 == 0 && round < 100) || round < 3) {
        //     std::cout << "Processing round " << round << std::endl;
        // }
        
        // Track if any updates occurred this round (convergence check)
        bool done = true;
        
        // ----------------------------------------------------------------------
        // Process each depression in the hierarchy
        // ----------------------------------------------------------------------
        for (unsigned x = 1; x < depression_hierarchy.size(); x++){
            // Skip depressions that don't need checking this round (optimization)
            if (round != 1 && !check_value_from_geolink[x]) {
                continue;
            }
            done = false;  // At least one depression is being processed
            
            // Get the pit cell location and elevation for this depression
            richdem::dephier::flat_c_idx dep_pit = depression_hierarchy[x].pit_cell;
            int dx, dy;
            topo.iToxy(dep_pit, dx, dy);
            T pit_elev_now = topo(dx, dy);
            
            // ----------------------------------------------------------------------
            // Traverse up the hierarchy propagating minimum pit elevation
            // ----------------------------------------------------------------------
            // Start at current depression and work upstream through parent chain
            richdem::dephier::dh_label_t current_label = x;
            while (current_label != richdem::dephier::OCEAN){
                // Mark this depression as checked (prevent redundant work)
                if (check_value[current_label]){
                    check_value[current_label] = false;
                }
                
                // Update minimum pit elevation for this depression
                if (pit_min[current_label] > pit_elev_now){
                    // Found a lower pit elevation - update it
                    pit_min[current_label] = pit_elev_now;
                } else {
                    // This depression already has the minimum - use it
                    pit_elev_now = pit_min[current_label];
                }

                // ----------------------------------------------------------------------
                // Propagate minimum to geolinked depression (lateral connection)
                // ----------------------------------------------------------------------
                // Geolinks connect depressions at the same level that eventually merge
                richdem::dephier::dh_label_t geolink = depression_hierarchy[current_label].geolink;
                if (pit_min[geolink] > pit_elev_now){
                    // Geolinked depression has a higher pit min - update it
                    pit_min[geolink] = pit_elev_now;
                    // Flag geolink for re-checking since it was updated
                    check_value[geolink] = true;
                }
                // --------------------------------------------------------------
                // Check if this depression flows to ocean (stop condition)
                // --------------------------------------------------------------
                // If the current depression flows to the ocean, stop the iteration
                // --------------------------------------------------------------
                if (depression_hierarchy[current_label].ocean_parent){
                    break;
                }
                // --------------------------------------------------------------
                // Move to parent depression and continue propagating upstream
                // --------------------------------------------------------------
                current_label = depression_hierarchy[current_label].parent;
            }
        }
        round++;
        if (done){
            break;
        }
    }
}