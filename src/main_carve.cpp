#include "include/stdefx.h"

#include <richdem/common/Array2D.hpp>
#include <richdem/common/version.hpp>
#include <richdem/common/constants.hpp>
#include <richdem/depressions/Lindsay2016.hpp>
#include <richdem/flats/flat_resolution.hpp>
#include <richdem/methods/d8_methods.hpp>

namespace {

void PrintUsage(const char *program_name) {
    std::cout << "\n=== Complete Breaching Workflow ===" << std::endl;
    std::cout << "\nUsage:" << std::endl;
    std::cout << "  " << program_name << " <input_dem> <output_prefix> <ocean_level>" << std::endl;
    std::cout << "\nArguments:" << std::endl;
    std::cout << "  input_dem      Path to the input DEM raster" << std::endl;
    std::cout << "  output_prefix  Prefix used to write outputs" << std::endl;
    std::cout << "  ocean_level    Reserved for future workflow variants" << std::endl;
    std::cout << "\nOutputs:" << std::endl;
    std::cout << "  <prefix>_DEM_breached.tif" << std::endl;
    std::cout << "  <prefix>_DIR.tif" << std::endl;
    std::cout << "  <prefix>_FAC.tif" << std::endl;
    std::cout << std::endl;
}

template <typename ElevationType>
int RunWorkflow(const std::string &input_dem, const std::string &output_prefix, const float ocean_level) {
    static_cast<void>(ocean_level);

    std::cout << "Loading DEM from: " << input_dem << std::endl;
    richdem::Array2D<ElevationType> elevations(input_dem);

    std::cout << "Data width  = " << elevations.width() << std::endl;
    std::cout << "Data height = " << elevations.height() << std::endl;
    std::cout << "Data cells  = " << elevations.numDataCells() << std::endl;

    std::cout << "Running complete breaching..." << std::endl;
    richdem::CompleteBreaching_Lindsay2016<richdem::Topology::D8>(elevations);

    const auto breached_output = output_prefix + "_DEM_breached.tif";
    std::cout << "Saving breached DEM to: " << breached_output << std::endl;
    elevations.saveGDAL(breached_output);

    std::cout << "Computing D8 flow directions..." << std::endl;
    richdem::Array2D<uint8_t> flowdirs;
    richdem::barnes_flat_resolution_d8(elevations, flowdirs, false);

    richdem::Array2D<int> flowdirs_esri(flowdirs, 0);
    #pragma omp parallel for collapse(2)
    for (int y = 0; y < flowdirs.height(); ++y) {
        for (int x = 0; x < flowdirs.width(); ++x) {
            if (flowdirs.isNoData(x, y)) {
                flowdirs_esri(x, y) = flowdirs_esri.noData();
            } else {
                flowdirs_esri(x, y) = richdem::d8_arcgis[flowdirs(x, y)];
            }
        }
    }

    const auto flowdirs_output = output_prefix + "_DIR.tif";
    std::cout << "Saving ESRI D8 flow directions to: " << flowdirs_output << std::endl;
    flowdirs_esri.saveGDAL(flowdirs_output);

    std::cout << "Computing D8 flow accumulation..." << std::endl;
    richdem::Array2D<double> accumulation(flowdirs, 0.0);
    richdem::d8_flow_accum(flowdirs, accumulation);

    const auto accumulation_output = output_prefix + "_FAC.tif";
    std::cout << "Saving D8 flow accumulation to: " << accumulation_output << std::endl;
    accumulation.saveGDAL(accumulation_output);

    std::cout << "Complete breaching workflow completed successfully!" << std::endl;
    return 0;
}

}  // namespace

int main(int argc, char **argv) {
    richdem::PrintRichdemHeader(argc, argv);

    if (argc < 4) {
        std::cerr << "ERROR: Not enough arguments specified" << std::endl;
        PrintUsage(argv[0]);
        return -1;
    }

    return RunWorkflow<float>(argv[1], argv[2], std::stof(argv[3]));
}
