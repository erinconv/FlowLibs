#!/usr/bin/env python3
"""
Mask all output rasters in a folder using a watershed raster.
"""

import argparse
import os
import sys

import docker_run


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Mask all output rasters in a folder using a watershed raster"
    )
    parser.add_argument(
        "input_path",
        type=str,
        help="Path to the folder containing watershed outputs",
    )
    parser.add_argument(
        "--watershed-name",
        type=str,
        default="run1_WAT.tif",
        help="Watershed raster name inside the folder",
    )
    parser.add_argument(
        "--suffix",
        type=str,
        default="_masked",
        help="Suffix added before the .tif extension for masked rasters",
    )
    parser.add_argument(
        "--image",
        type=str,
        default="erinconv/flowlibs:latest",
        help="Docker image to run",
    )

    args = parser.parse_args()

    input_path = os.path.abspath(args.input_path)
    if not os.path.isdir(input_path):
        print(f"ERROR: Input folder not found: {input_path}", file=sys.stderr)
        return 1

    watershed_path = os.path.join(input_path, args.watershed_name)
    if not os.path.isfile(watershed_path):
        print(f"ERROR: Watershed raster not found: {watershed_path}", file=sys.stderr)
        return 1

    rasters_to_mask = []
    for filename in sorted(os.listdir(input_path)):
        if not filename.lower().endswith(".tif"):
            continue
        if filename == args.watershed_name:
            continue
        if filename.endswith(f"{args.suffix}.tif"):
            continue
        rasters_to_mask.append(filename)

    if not rasters_to_mask:
        print("No rasters found to mask.")
        return 0

    print(f"Watershed raster: {args.watershed_name}")
    print("Rasters to mask:")
    for raster_name in rasters_to_mask:
        print(f"  {raster_name}")

    for raster_name in rasters_to_mask:
        stem, ext = os.path.splitext(raster_name)
        output_name = f"{stem}{args.suffix}{ext}"

        cmd = [
            "docker",
            "run",
            "--rm",
            "--entrypoint",
            "mask_with_watershed.exe",
            "-v",
            f"{input_path}:/out",
            args.image,
            f"/out/{args.watershed_name}",
            f"/out/{raster_name}",
            f"/out/{output_name}",
        ]

        print("\nRunning command:")
        print(" ".join(cmd))

        rc = docker_run.run_docker(cmd)
        if rc != 0:
            print(f"ERROR: Failed masking {raster_name} with code {rc}", file=sys.stderr)
            return rc

    print("\nMasked rasters created successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
