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

## Documentation site

This repository now includes a MkDocs Material documentation site in [`docs/`](docs/).

Install the documentation dependencies and run it locally with:

```bash
pip install -r docs/requirements.txt
mkdocs serve
```

The most detailed page for the new vector export is:

- `docs/vector-outputs.md`

Once the algorithms are either built or pulled, this repository provides different Python wrappers to facilitate the execution of the algorithms.

### Running the carving algorithm

```bash
python python/run_carving.py path/to/your/DEM.tif
```

This command now runs a complete-breaching workflow. It produces:

- `run1_DEM_breached.tif`
- `run1_DIR.tif`
- `run1_FAC.tif`

The `DIR` raster is written in ESRI/ArcGIS D8 format.

### Running the watershed delineation

```bash
python python/run_watersheds.py path/to/your/carving/outputs --pour_points X Y id basin_name
```
Where `X` and `Y` correspond to the point coordinate, `id` is the label that will be assigned to the delineated watershed, and `basin_name` is the name of the basin. This command will delineate watersheds based on the provided pour point coordinates.

### Exporting the river network to vector format

After running the carving workflow and watershed delineation, you can export the stream network to a GeoPackage:

```bash
python python/run_stream_vectorize.py path/to/your/carving/outputs
```

By default, this command expects the following rasters in the folder:

- `run1_DEM_breached.tif`
- `run1_DIR.tif`
- `run1_FAC.tif`
- `run1_streams.tif`
- `run1_stream_order.tif`
- `run1_watersheds.tif`
- `run1_subcatchments.tif`

It writes:

- `run1_river_network.gpkg`

The GeoPackage contains two layers:

- `river_segments`
- `river_nodes`

The exported vector network includes hydrological attributes such as stream order, contributing accumulation, watershed and subcatchment IDs, segment length, elevation drop, and slope.
