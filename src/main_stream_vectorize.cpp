#include "include/stream_vectorization.hpp"

#include <richdem/common/version.hpp>

namespace {

void PrintUsage(const char *program_name) {
    std::cout << "\n=== Stream Network Vectorization ===" << std::endl;
    std::cout << "\nUsage:" << std::endl;
    std::cout << "  " << program_name
              << " <dem> <flow_dirs> <accumulation> <streams> <stream_order> <watersheds> <subcatchments> <output_gpkg>"
              << std::endl;
    std::cout << std::endl;
}

}  // namespace

int main(int argc, char **argv) {
    richdem::PrintRichdemHeader(argc, argv);

    if (argc != 9) {
        PrintUsage(argv[0]);
        return -1;
    }

    stream_vectorization::StreamVectorizationParams params;
    params.input_dem = argv[1];
    params.input_flowdirs = argv[2];
    params.input_accumulation = argv[3];
    params.input_streams = argv[4];
    params.input_stream_order = argv[5];
    params.input_watersheds = argv[6];
    params.input_subcatchments = argv[7];
    params.output_vector = argv[8];

    return stream_vectorization::export_stream_network_vector(params);
}
