#!/usr/bin/env python3
"""
Watershed Delineation Script
Runs watershed delineation using the dephier Docker container
"""

import os
import subprocess
import sys
import argparse

def create_pourpoints_file(output_file, pour_points):
    """
    Create a pour points file from a list of tuples
    
    Args:
        output_file: Path to save the pour points file
        pour_points: List of tuples (longitude, latitude, id, name)
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("# Pour Points File\n")
        f.write("# Format: longitude latitude id name\n")
        f.write("# Coordinates are in decimal degrees (e.g., -69.9408235 48.2644502)\n")
        f.write("# Longitude: negative = West, positive = East\n")
        f.write("# Latitude: negative = South, positive = North\n\n")

        for pp in pour_points:
            if len(pp) == 4:
                lon, lat, id, name = pp
                f.write(f"{lon} {lat} {id} {name}\n")
            elif len(pp) == 3:
                lon, lat, id = pp
                f.write(f"{lon} {lat} {id}\n")

    print(f"Created pour points file: {output_file}")


def main():
    """_summary_

    Returns:
        _type_: _description_
    """
    parser = argparse.ArgumentParser("Python wrapper for executing the watershed delineation algorithm from FlowLibs")

    parser.add_argument("input_path",
                        type = str,
                        help = "Path to the Folder containing the input files")
    parser.add_argument("--dem_name",
                        type = str,
                        help = "Name of the input DEM. Preferable the CARVED DEM for more accurate results.",
                        default = "run1-mod-dem-carved.tif")
    parser.add_argument("--flow_dirs_name",
                        type = str,
                        help = "NAme of the flow directions file.",
                        default = "run1-mod-dirs.tif")
    parser.add_argument("--accumulation_name",
                        type = str,
                        help = "Name of the accumulation file.",
                        default = "run1_accu.tif")
    parser.add_argument("--output_prefix",
                        type = str,
                        default = "run1",
                        help = "Output prefix that will be used to name output files.")
    parser.add_argument("--pour_points",
                        action="append",
                        nargs = 4,
                        metavar = ("X", "Y", "ID", "Basin name"),
                        help = "Specify a pour point as x y name. Can be used multiple times.")
    parser.add_argument("--stream_threshold",
                        type = int,
                        default=10000,
                        help = "Thereshold to define stream and subcatcment delineation. DEFAULT = 10000 cells")

    args = parser.parse_args()

    if args.pour_points:
        print("Input pour points:")
        for ppoint in args.pour_points:
            # if isinstance(ppoint, tuple):
            # 69.8961205°W 48.4220364°N 
            #     pass
            # else:
            #     raise TypeError(f"Invalid input for pour point: {ppoint}")
            x, y, idx, name = ppoint
            try:
                x = float(x)
                y = float(y)
                idx = int(idx)
                print(f"  Point Name: {name}, Coordinates: ({x}, {y})")
            except ValueError as e:
                raise ValueError(f"  Error with point {name}: X or Y were not valid numbers. {e}") from e
    else:
        raise ValueError("No pour point specified")

    input_path = args.input_path
    pour_points = args.pour_points
    stream_threshold = args.stream_threshold
    watershed_output_prefix = args.output_prefix
    
    # Input files (from previous carve and accu commands)
    carved_dem = os.path.join(input_path, args.dem_name)
    flow_dirs = os.path.join(input_path, args.flow_dirs_name)
    accumulation = os.path.join(input_path, args.accumulation_name)

    # Docker image name
    image = "erinconv/flowlibs:latest"

    # ========== DEFINE POUR POINTS ==========
    # Define pour points using geographic coordinates (decimal degrees)
    # Format: (longitude, latitude, id, name)
    # Longitude: negative = West, positive = East
    # Latitude: negative = South, positive = North
    
    # Example: 69.9408235°W 48.2644502°N
    # Convert W to negative: -69.9408235
    # 69.9439113°W 48.2642379°N 
    # pour_points = [
    #     (-69.9439113, 48.2642379, 1, "Main_Watershed"),
    #     # (-69.9472252, 48.2633507, 1, "Main_Watershed"),
    #     # Add more pour points as needed:
    #     # (-69.5, 48.5, 2, "Tributary_1"),
    #     # (-70.0, 48.0, 3, "Tributary_2"),
    # ]

    # Create pour points file
    pourpoints_file = os.path.join(input_path, "pourpoints.txt")
    create_pourpoints_file(pourpoints_file, pour_points)

    # ========== CHECK INPUT FILES ==========

    if not os.path.isfile(carved_dem):
        print(f"ERROR: Carved DEM not found: {carved_dem}", file=sys.stderr)
        print("Run 'carve' command first with 'mod' flag", file=sys.stderr)
        return 1

    if not os.path.isfile(flow_dirs):
        print(f"ERROR: Flow directions not found: {flow_dirs}", file=sys.stderr)
        print("Run 'carve' command first with 'mod' flag", file=sys.stderr)
        return 1

    if not os.path.isfile(accumulation):
        print(f"ERROR: Accumulation not found: {accumulation}", file=sys.stderr)
        print("Run 'accu' command first", file=sys.stderr)
        return 1

    # ========== RUN WATERSHED DELINEATION ==========

    print("\n" + "="*60)
    print("WATERSHED DELINEATION")
    print("="*60)
    print(f"Carved DEM:      {carved_dem}")
    print(f"Flow Directions: {flow_dirs}")
    print(f"Accumulation:    {accumulation}")
    print(f"Output Prefix:   {watershed_output_prefix}")
    print(f"Stream Threshold: {stream_threshold} cells")
    print(f"Pour Points:     {len(pour_points)}")
    print("="*60 + "\n")

    # Build Docker command
    cmd = [
        "docker",
        "run",
        "--rm",
        "--entrypoint",
        "watershed_delineation.exe",
        "-v",
        f"{input_path}:/out",
        image,
        f"/out/{os.path.basename(carved_dem)}",
        f"/out/{os.path.basename(flow_dirs)}",
        f"/out/{os.path.basename(accumulation)}",
        f"/out/{watershed_output_prefix}",
        f"/out/{os.path.basename(pourpoints_file)}",
        str(stream_threshold),
    ]

    print("Running command:")
    print(" ".join(cmd))
    print()

    try:
        result = subprocess.run(cmd, check=True)

        if result.returncode == 0:
            print("\n" + "="*60)
            print("SUCCESS!")
            print("="*60)
            print("\nOutput files:")
            print(f"  {watershed_output_prefix}_streams.tif        - Stream network (binary)")
            print(f"  {watershed_output_prefix}_stream_order.tif   - Horton-Strahler stream order")
            print(f"  {watershed_output_prefix}_watersheds.tif     - Watershed labels")
            print(f"  {watershed_output_prefix}_subcatchments.tif  - Subcatchment labels")
            print(f"\nLocation: {input_path}")

        return result.returncode

    except subprocess.CalledProcessError as e:
        print(f"\nERROR: Watershed delineation failed with code {e.returncode}", file=sys.stderr)
        return e.returncode


if __name__ == "__main__":
    raise SystemExit(main())

