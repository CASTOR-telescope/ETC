"""
filepaths.py

Contains filepaths for CASTOR.

Isaac Cheng - 2022
"""

from os.path import dirname, abspath

# ------------------------- FILEPATHS (N.B. the trailing slash) ------------------------ #

# The location of the castor_etc package
__BASEPATH = dirname(abspath(__file__)) + "/"

# The directory to store output files. This directory must already exist and be
# discoverable by the container (e.g., via the "-v" bind mount for a local Docker build)
OUTPATH = "/arc/home/IsaacCheng/CASTOR/ETC_plots/"  # ! CHANGE ME !

# The directory containing the data files (e.g., passbands, sky background, etc.)
DATAPATH = __BASEPATH + "data/"

# -------------------------------------------------------------------------------------- #
