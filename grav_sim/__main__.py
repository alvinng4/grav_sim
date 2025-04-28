"""
Print compilation information for the grav_sim C library.

Usage:
    python -m grav_sim [-p PATH]

This script searches for the compiled C library in the parent directory of the current file
(or a user-specified path), loads it using ctypes, and calls the `print_compilation_info` 
function defined in the C library. Users can run this script to verify the installation 
of the grav_sim library.
"""

import argparse
import ctypes
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Print compilation information for the grav_sim C library.")
    parser.add_argument(
        "-p", "--path", type=str, default=None,
        help="Path to search for the grav_sim C library. Defaults to the parent directory of the current file."
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    search_path = Path(args.path).parent if args.path else Path(__file__).parent.parent
    c_lib_files = [str(p) for p in search_path.rglob("*libgrav_sim*")]
    if len(c_lib_files) == 0:
        raise FileNotFoundError(f"C library not found from path: {search_path}")

    c_lib_path = c_lib_files[0]
    c_lib: ctypes.CDLL = ctypes.cdll.LoadLibrary(c_lib_path)
    c_lib.print_compilation_info()
    print(f"C library location: {c_lib_path}")


if __name__ == "__main__":
    main()
