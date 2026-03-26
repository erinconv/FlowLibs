#!/usr/bin/env python3
"""
Export the stream network to a vector dataset.
"""

import argparse
import os
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Export stream network rasters to a GeoPackage vector dataset"
    )
    parser.add_argument(
        "input_path",
        type=str,
        help="Path to the folder containing the workflow and watershed outputs",
    )
    parser.add_argument("--dem-name", type=str, default="run1_DEM_breached.tif")
    parser.add_argument("--flow-dirs-name", type=str, default="run1_DIR_masked.tif")
    parser.add_argument("--accumulation-name", type=str, default="run1_FAC_masked.tif")
    parser.add_argument("--streams-name", type=str, default="run1_streams_masked.tif")
    parser.add_argument("--stream-order-name", type=str, default="run1_stream_order.tif")
    parser.add_argument("--watersheds-name", type=str, default="run1_watersheds.tif")
    parser.add_argument("--subcatchments-name", type=str, default="run1_subcatchments.tif")
    parser.add_argument("--output-name", type=str, default="run1_river_network_masked.gpkg")
    parser.add_argument("--image", type=str, default="erinconv/flowlibs:latest")

    args = parser.parse_args()

    input_path = os.path.abspath(args.input_path)
    if not os.path.isdir(input_path):
        print(f"ERROR: Input folder not found: {input_path}", file=sys.stderr)
        return 1

    required_files = [
        args.dem_name,
        args.flow_dirs_name,
        args.accumulation_name,
        args.streams_name,
        args.stream_order_name,
        args.watersheds_name,
        args.subcatchments_name,
    ]

    for file_name in required_files:
        file_path = os.path.join(input_path, file_name)
        if not os.path.isfile(file_path):
            print(f"ERROR: Missing required file: {file_path}", file=sys.stderr)
            return 1

    print("Using masked hydrology rasters for vectorization:")
    print(f"  flow directions: {args.flow_dirs_name}")
    print(f"  accumulation:    {args.accumulation_name}")
    print(f"  streams:         {args.streams_name}")
    print(f"  output:          {args.output_name}")
    print()

    cmd = [
        "docker",
        "run",
        "--rm",
        "--entrypoint",
        "stream_vectorize.exe",
        "-v",
        f"{input_path}:/out",
        args.image,
        f"/out/{args.dem_name}",
        f"/out/{args.flow_dirs_name}",
        f"/out/{args.accumulation_name}",
        f"/out/{args.streams_name}",
        f"/out/{args.stream_order_name}",
        f"/out/{args.watersheds_name}",
        f"/out/{args.subcatchments_name}",
        f"/out/{args.output_name}",
    ]

    print("Running command:")
    print(" ".join(cmd))
    print()

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"ERROR: Stream vectorization failed with code {exc.returncode}", file=sys.stderr)
        return exc.returncode

    print(f"Vector network: {os.path.join(input_path, args.output_name)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
