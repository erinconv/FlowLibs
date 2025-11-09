# FlowLibs

This library is a collection of C++ implementations of different hydrological algorithms for terrain analysis. It performs basic watershed and subwatershed delineation, computes flow direction, and calculates flow accumulation. The main purpose of this library is to serve as a support tool for the [pycequeau](https://github.com/erinconv/pycequeau) library, whose main objective is to perform physiographic analysis for obtaining input data for the CEQUEAU model.

## Can it be used as standalone software?

Although its main purpose is to serve as a supporting library for Pycequeau, this library is designed to function as standalone software. It can be easily downloaded and compiled for more general purposes related to hydrological analysis and physiographic basin description.

## How to use it?

The easiest way to use this library is through Docker. You can build it locally by running the following command:

```bash
  docker build -t erinconv/flowlibs:latest -f docker/Dockerfile .
```

Alternatively, you can pull the latest build from Docker Hub:

```bash
  docker pull erinconv/flowlibs:latest
```

Once the algorithms are either built or pulled, this repository provides different Python wrappers to facilitate the execution of the algorithms.

### Running the carving algorithm

```bash
python python/run_carving.py path/to/your/DEM.tif
```

This command will produce the flow directions (`run1-mod-dirs.tif`), flow accumulation (`run1_accu.tif`), and the carved DEM (`run1-mod-dem-carved.tif`).

### Running the watershed delineation

```bash
python python/run_watersheds.py path/to/your/carving/outputs --pour_points X Y id basin_name
```
Where `X` and `Y` correspond to the point coordinate, `id` is the label that will be assigned to the delineated watershed, and `basin_name` is the name of the basin. This command will delineate watersheds based on the provided pour point coordinates.