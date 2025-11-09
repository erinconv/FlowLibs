import os
import argparse
import subprocess
import sys


def main() -> int:
    """Python wrapper for executing the carving algorithm from FlowLibs

    Raises:
        FileNotFoundError: If the inut file is not found

    Returns:
        int: Exit code (0 for success)
    """
    parser = argparse.ArgumentParser("Python wrapper for executing the carving algorithm from FlowLibs")

    parser.add_argument("input_path",
                        type =str,
                        help = "Path to the input DEM File")

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
    output_prefix_in_container = "/out/run1"
    ocean_level = "0"

    # Build docker command:
    # - Mount input directory read-only at /data
    # - Mount output directory read-write at /out
    # - Invoke: dephier <command> /data/<input> /out/run1 <ocean_level> [additional args]
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

    # if command == "2":
    #     # Add 'mod' to enable saving flow directions
    cmd.append("mod")
    print(" ".join(cmd))

    try:
        completed = subprocess.run(cmd, check=True)
        if completed.returncode != 0:
            return completed.returncode
    except subprocess.CalledProcessError as e:
        return e.returncode

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


