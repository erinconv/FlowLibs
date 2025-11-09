
#include "include/flow_algorithms.hpp"
#include <iostream>
#include <string>
#include <cstdlib>

// Print usage information
void PrintUsage(const char* prog_name) {
  std::cout << "\n=== Carving Algorithm ===" << std::endl;
  std::cout << "\nUsage:" << std::endl;
  std::cout << "  " << prog_name << "<Input DEM> <Output Prefix> <Ocean Level> [mod|modify]" << std::endl;
  std::cout << "\nExamples:" << std::endl;
  std::cout << "  " << prog_name << "dem.tif output 0 mod" << std::endl;
  std::cout << std::endl;
}

// Parse command line arguments and dispatch to appropriate algorithm
int main(int argc, char** argv) {
  // Check minimum arguments
  if (argc < 2) {
    std::cerr << "ERROR: Not enough arguments specified" << std::endl;
    PrintUsage(argv[0]);
    return -1;
  }
  

  // Setup parameters
  Params params;
  params.input_file = argv[1];
  params.output_prefix = argv[2];
  params.ocean_level = std::stof(argv[3]);
  // modify_mode has default value = false in struct definition

  // Load the DEM
  richdem::Array2D<float> topo(params.input_file);
  // Initialize the label and flow dirs
  richdem::Array2D<richdem::dephier::dh_label_t> label(topo, richdem::dephier::NO_DEP);
  richdem::Array2D<richdem::flowdir_t> flowdirs(topo, richdem::NO_FLOW);

  // Parse optional modify mode argument
  params.modify_mode = true;
  
  // Create FlowAlgorithms object AFTER setting all params
  FlowAlgorithms flow_algorithms(params);
    
  // Compute the pre-treatment, where ocean labels are initialized
  flow_algorithms.initialize_ocean_labels(topo, label);
  
  // Run carving algorithm
  flow_algorithms.compute_carving(topo, label, flowdirs);
  
  // Compute the flow accumulation
  if (params.modify_mode = true){
    // std::string flow_dirs_file = params.output_prefix + "-mod-flowdirs.tif";
    std::string flow_dirs_file = params.output_prefix + "-mod-dirs.tif";
    int a = flow_algorithms.compute_flow_accumulation(flow_dirs_file);
    if (a != 0){
      std::cerr << "ERROR: Flow accumulation failed" << std::endl;
      return a;
    }
  } 
  
  std::cout << "Carve algorithm completed successfully!" << std::endl;
  
  return 0;
}

