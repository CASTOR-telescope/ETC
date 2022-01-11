"""
constants.py

A file to track filepaths and telescope parameters.

Isaac Cheng - January 2022
"""

from pathlib import Path

# ------------------------- FILEPATHS (N.B. the trailing slash) ------------------------ #

# The location of the ETC repo
BASEPATH = str(Path(__file__).parent.parent) + "/"  # e.g., "/arc/home/IsaacCheng/ETC/"

# The directory to store output files. This directory must already exist and be
# discoverable by the container (e.g., via the "-v" bind mount for a local Docker build)
OUTPATH = "/arc/home/IsaacCheng/ETC_plots/"

# The directory containing the passband response files
PASSBAND_PATH = BASEPATH + "data/passbands/"


# -------------------------------------------------------------------------------------- #

# -------------------------------- TELESCOPE PARAMETERS -------------------------------- #
# -------------------------------------------------------------------------------------- #
