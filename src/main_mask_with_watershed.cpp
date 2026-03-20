#include "include/watershed_masking.hpp"

#include "include/stdefx.h"

#include <richdem/common/version.hpp>

namespace {

void PrintUsage(const char *program_name) {
    std::cout << "\n=== Mask Raster With Watershed ===" << std::endl;
    std::cout << "\nUsage:" << std::endl;
    std::cout << "  " << program_name << " <watershed_raster> <input_raster> <output_raster>" << std::endl;
    std::cout << std::endl;
}

}  // namespace

int main(int argc, char **argv) {
    richdem::PrintRichdemHeader(argc, argv);

    if (argc != 4) {
        PrintUsage(argv[0]);
        return -1;
    }

    return watershed_masking::dispatch_mask_raster_with_watershed(argv[1], argv[2], argv[3]);
}
