#pragma once

#include "stdefx.h"

#include <gdal_priv.h>
#include <richdem/common/Array2D.hpp>
#include <richdem/common/gdal.hpp>

namespace watershed_masking {

template <typename MaskType, typename RasterType>
void mask_array_with_watershed(
    const richdem::Array2D<MaskType> &watershed,
    const richdem::Array2D<RasterType> &raster,
    richdem::Array2D<RasterType> &masked) {
    if (watershed.width() != raster.width() || watershed.height() != raster.height()) {
        throw std::runtime_error("Watershed raster dimensions do not match input raster dimensions.");
    }

    masked = richdem::Array2D<RasterType>(raster);
    masked.setNoData(raster.noData());

    for (int y = 0; y < raster.height(); ++y) {
        for (int x = 0; x < raster.width(); ++x) {
            const bool outside_watershed =
                watershed.isNoData(x, y) || watershed(x, y) == 0;

            masked(x, y) = outside_watershed ? raster.noData() : raster(x, y);
        }
    }
}

template <typename RasterType>
int mask_raster_with_watershed(
    const std::string &watershed_raster,
    const std::string &input_raster,
    const std::string &output_raster) {
    richdem::Array2D<int> watershed(watershed_raster);
    richdem::Array2D<RasterType> raster(input_raster);

    if (watershed.width() != raster.width() || watershed.height() != raster.height()) {
        std::cerr << "ERROR: Watershed raster dimensions do not match input raster dimensions." << std::endl;
        return -1;
    }

    richdem::Array2D<RasterType> masked;
    mask_array_with_watershed(watershed, raster, masked);

    masked.saveGDAL(output_raster);
    return 0;
}

inline int dispatch_mask_raster_with_watershed(
    const std::string &watershed_raster,
    const std::string &input_raster,
    const std::string &output_raster) {
    switch (richdem::peekGDALType(input_raster)) {
        case GDT_Byte:
            return mask_raster_with_watershed<uint8_t>(watershed_raster, input_raster, output_raster);
        case GDT_UInt16:
            return mask_raster_with_watershed<uint16_t>(watershed_raster, input_raster, output_raster);
        case GDT_Int16:
            return mask_raster_with_watershed<int16_t>(watershed_raster, input_raster, output_raster);
        case GDT_UInt32:
            return mask_raster_with_watershed<uint32_t>(watershed_raster, input_raster, output_raster);
        case GDT_Int32:
            return mask_raster_with_watershed<int32_t>(watershed_raster, input_raster, output_raster);
        case GDT_Float32:
            return mask_raster_with_watershed<float>(watershed_raster, input_raster, output_raster);
        case GDT_Float64:
            return mask_raster_with_watershed<double>(watershed_raster, input_raster, output_raster);
        default:
            std::cerr << "ERROR: Unsupported raster type." << std::endl;
            return -1;
    }
}

}  // namespace watershed_masking
