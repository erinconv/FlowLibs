# FlowLibs

This library is a collection of C++ implementations of different hydrological algorithms for terrain analysis. It performs basic watershed and subwatershed delineation, flow direction computation, and flow accumulation. The main purpose of this library is to serve as a support tool for the [pycequeau](https://github.com/erinconv/pycequeau) library, whose main purpose is to perform physiographic analysis for obtaining input data for the CEQUEAU model.

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