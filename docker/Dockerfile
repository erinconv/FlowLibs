# Build and run Depression Hierarchies with GDAL in Docker
FROM ubuntu:22.04 AS build

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      build-essential \
      cmake \
      ninja-build \
      git \
      libgdal-dev \
      libomp-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy the repository into the image
COPY . /app

# Ensure submodules are present
RUN git submodule update --init --recursive || true

# Configure and build with Make, Release, GDAL enabled
RUN mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DUSE_GDAL=On .. && \
    make -j $(nproc)

# Runtime stage (keeps image smaller while retaining GDAL runtime)
FROM ubuntu:22.04 AS runtime

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      gdal-bin \
      libgdal30 \
      libgomp1 \
      libstdc++6 \
      libomp5 && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /work

# Copy the built binaries (both main executable and standalone watershed tool)
COPY --from=build /app/build/main_carve.exe /usr/local/bin/main_carve
COPY --from=build /app/build/watershed_delineation.exe /usr/local/bin/watershed_delineation.exe

# Add all binaries to PATH for easy access
ENV PATH="/usr/local/bin:${PATH}"

# Default entrypoint is the main dephier binary; args are CLI args
# To use main_carve.exe or watershed_delineation.exe, override entrypoint: docker run --entrypoint main_carve.exe ... or docker run --entrypoint watershed_delineation.exe ...
# ENTRYPOINT ["/usr/local/bin/FlowLibs"]
CMD ["--help"]