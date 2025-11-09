
 #include "include/watershed_delineation.hpp"
 #include "include/stdefx.h"
 
 /**
  * @brief Print usage information
  */
 void PrintUsage(const char* program_name) {
     std::cout << "\n=== Watershed Delineation - Standalone Executable ===" << std::endl;
     std::cout << "\nUsage:" << std::endl;
     std::cout << "  " << program_name 
               << " <carved_dem> <flow_dirs> <accumulation> <output_prefix> <pourpoints_file> <stream_threshold>"
               << std::endl;
     
     std::cout << "\nArguments:" << std::endl;
     std::cout << "  carved_dem       Path to carved DEM file (from carve command)" << std::endl;
     std::cout << "  flow_dirs        Path to flow directions file (from carve with mod flag)" << std::endl;
     std::cout << "  accumulation     Path to flow accumulation file (from accu command)" << std::endl;
     std::cout << "  output_prefix    Prefix for output files" << std::endl;
     std::cout << "  pourpoints_file  Path to file with pour point coordinates" << std::endl;
     std::cout << "  stream_threshold Minimum accumulation to define streams (e.g., 10000)" << std::endl;
     
     std::cout << "\nPour Points File Format:" << std::endl;
     std::cout << "  longitude latitude id name" << std::endl;
     std::cout << "  (one pour point per line, coordinates in decimal degrees)" << std::endl;
     std::cout << "\nExample Pour Points File:" << std::endl;
     std::cout << "  -69.9408235 48.2644502 1 Main_Watershed" << std::endl;
     std::cout << "  -69.8523456 48.1234567 2 North_Tributary" << std::endl;
     
     std::cout << "\nExample Usage:" << std::endl;
     std::cout << "  " << program_name 
               << " dem_carved.tif flowdirs.tif accu.tif output pourpoints.txt 10000"
               << std::endl;
     
     std::cout << "\nOutputs:" << std::endl;
     std::cout << "  <output_prefix>_streams.tif        Binary stream network (0/1)" << std::endl;
     std::cout << "  <output_prefix>_stream_order.tif   Horton-Strahler stream order" << std::endl;
     std::cout << "  <output_prefix>_watersheds.tif     Watershed labels" << std::endl;
     std::cout << "  <output_prefix>_subcatchments.tif  Subcatchment labels" << std::endl;
     
     std::cout << "\nNotes:" << std::endl;
     std::cout << "  - Pour points are automatically snapped to nearest stream (within 10 pixels)" << std::endl;
     std::cout << "  - Subcatchments are automatically delineated based on stream junctions" << std::endl;
     std::cout << "  - Stream order is calculated using Horton-Strahler classification" << std::endl;
     std::cout << "  - All output rasters use NoData = 0 for background values" << std::endl;
     std::cout << std::endl;
 }
 
 /**
  * @brief Load pour points from file
  * 
  * File format: longitude latitude id name
  * Lines starting with # are comments
  * 
  * @param filename Path to pour points file
  * @return Vector of pour points with geographic coordinates
  */
 std::vector<PourPoint> LoadPourPoints(const std::string& filename) {
     std::vector<PourPoint> pour_points;
     std::ifstream file(filename);
     
     if (!file.is_open()) {
         std::cerr << "ERROR: Could not open pour points file: " << filename << std::endl;
         return pour_points;
     }
     
     std::cout << "Loading pour points from: " << filename << std::endl;
     
     std::string line;
     int line_number = 0;
     
     while (std::getline(file, line)) {
         line_number++;
         
         // Skip empty lines and comments
         if (line.empty() || line[0] == '#') {
             continue;
         }
         
         std::istringstream iss(line);
         double lon, lat;
         int id;
         std::string name;
         
         // Parse line: longitude latitude id name
         if (!(iss >> lon >> lat >> id)) {
             std::cerr << "WARNING: Invalid format at line " << line_number 
                       << " (skipping): " << line << std::endl;
             continue;
         }
         
         // Name is optional (rest of the line)
         std::getline(iss, name);
         // Trim leading whitespace from name
         size_t start = name.find_first_not_of(" \t");
         if (start != std::string::npos) {
             name = name.substr(start);
         } else {
             name = "Watershed_" + std::to_string(id);
         }
         
         // Create pour point with geographic coordinates
         PourPoint pp(0, 0, id, name);  // Initialize with pixel coords (0,0) temporarily
         pp.geo_x = lon;
         pp.geo_y = lat;
         pp.is_geographic = true;
         
         pour_points.push_back(pp);
         
         std::cout << "  Pour point " << id << ": (" << lon << ", " << lat 
                   << ") - " << name << std::endl;
     }
     
     file.close();
     
     std::cout << "Loaded " << pour_points.size() << " pour point(s)" << std::endl;
     
     return pour_points;
 }
 
 /**
  * @brief Main function
  */
 int main(int argc, char** argv) {
     // Check arguments
     if (argc < 7) {
         std::cerr << "ERROR: Not enough arguments specified" << std::endl;
         PrintUsage(argv[0]);
         return -1;
     }
     
     // Parse command line arguments
     std::string carved_dem = argv[1];
     std::string flow_dirs = argv[2];
     std::string accumulation = argv[3];
     std::string output_prefix = argv[4];
     std::string pourpoints_file = argv[5];
     float stream_threshold = std::stof(argv[6]);
     
     // Display configuration
     std::cout << "Configuration:" << std::endl;
     std::cout << "  Carved DEM:       " << carved_dem << std::endl;
     std::cout << "  Flow Directions:  " << flow_dirs << std::endl;
     std::cout << "  Accumulation:     " << accumulation << std::endl;
     std::cout << "  Output Prefix:    " << output_prefix << std::endl;
     std::cout << "  Pour Points File: " << pourpoints_file << std::endl;
     std::cout << "  Stream Threshold: " << stream_threshold << " cells" << std::endl;
     std::cout << std::endl;
     
     // Load pour points
     std::vector<PourPoint> pour_points = LoadPourPoints(pourpoints_file);
     
     if (pour_points.empty()) {
         std::cerr << "ERROR: No valid pour points loaded from file" << std::endl;
         return -1;
     }
     
     // Create watershed parameters
     WatershedParams params(
         carved_dem,
         flow_dirs,
         accumulation,
         output_prefix,
         stream_threshold,
         false  // Use ESRI/ArcGIS flow direction format (default)
     );
     
     // Create watershed delineation object
     WatershedDelineation watershed_delineator(params);
     
     // Run watershed delineation
     std::cout << "\nStarting watershed delineation..." << std::endl;
     int result = watershed_delineator.run_watershed_delineation(pour_points);
     
     // Report results
     if (result == 0) {
         std::cout << "\n========================================" << std::endl;
         std::cout << "  SUCCESS!" << std::endl;
         std::cout << "========================================" << std::endl;
         std::cout << "\nOutput files created:" << std::endl;
         std::cout << "  " << output_prefix << "_streams.tif        - Stream network (binary)" << std::endl;
         std::cout << "  " << output_prefix << "_stream_order.tif   - Horton-Strahler stream order" << std::endl;
         std::cout << "  " << output_prefix << "_watersheds.tif     - Watershed labels" << std::endl;
         std::cout << "  " << output_prefix << "_subcatchments.tif  - Subcatchment labels" << std::endl;
         std::cout << std::endl;
     } else {
         std::cout << "\n========================================" << std::endl;
         std::cout << "  FAILED!" << std::endl;
         std::cout << "========================================" << std::endl;
         std::cerr << "ERROR: Watershed delineation failed with error code: " << result << std::endl;
         std::cout << std::endl;
     }
     
     return result;
 }
 
 