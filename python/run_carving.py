import os
import argparse
import subprocess
import sys


def main() -> int:
    """Python wrapper for executing the complete-breaching workflow from FlowLibs

    Raises:
        FileNotFoundError: If the inut file is not found

    Returns:
        int: Exit code (0 for success)
    """
    parser = argparse.ArgumentParser("Python wrapper for executing the carving workflow from FlowLibs")

    parser.add_argument("input_path",
                        type =str,
                        help = "Path to the input DEM File")
    parser.add_argument("--output-prefix",
                        type = str,
                        default = "run1",
                        help = "Output prefix used to name output rasters")
    parser.add_argument("--ocean-level",
                        type = float,
                        default = 0.0,
                        help = "Reserved for future workflow variants")

    args = parser.parse_args()

    # Assign parsed arguments to the variablesd
    input_path = args.input_path
    if not os.path.isfile(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}", file=sys.stderr)

    mount_dir = os.path.dirname(input_path)
    input_name = os.path.basename(input_path)
    image = "erinconv/flowlibs:latest"

    # Create an output directory alongside the input
    output_dir = os.path.join(mount_dir, "dephier_outputs")
    os.makedirs(output_dir, exist_ok=True)

    # Output prefix (inside container) and ocean level
    output_prefix_in_container = f"/out/{args.output_prefix}"
    ocean_level = str(args.ocean_level)

    # Build docker command:
    # - Mount input directory read-only at /data
    # - Mount output directory read-write at /out
    # - Invoke the main complete-breaching workflow
    cmd = [
        "docker",
        "run",
        "--rm",
        "--entrypoint", "main_carve",  # Note: no .exe extension in container
        "-v",
        f"{mount_dir}:/data:ro",
        "-v",
        f"{output_dir}:/out",
        image,
        f"/data/{input_name}",
        output_prefix_in_container,
        ocean_level,
    ]

    print("Running command:")
    print(" ".join(cmd))
    print()

    try:
        completed = subprocess.run(cmd, check=True)
        if completed.returncode != 0:
            return completed.returncode
    except subprocess.CalledProcessError as e:
        return e.returncode

    print(f"Breached DEM: {os.path.join(output_dir, args.output_prefix + '_DEM_breached.tif')}")
    print(f"Flow directions: {os.path.join(output_dir, args.output_prefix + '_DIR.tif')}")
    print(f"Accumulation: {os.path.join(output_dir, args.output_prefix + '_FAC.tif')}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
