"""
Contains filepaths for CASTOR.
"""

from os.path import dirname, abspath, join

# ------------------------- FILEPATHS (N.B. no trailing slash) ------------------------- #

# The location of the castor_etc package
__BASEPATH = dirname(abspath(__file__))  # + "/"

# The directory to store output files. This directory must already exist and be
# discoverable by the container (e.g., via the "-v" bind mount for a local Docker build)
OUTPATH = "/arc/home/IsaacCheng/CASTOR/ETC_plots"  # ! CHANGE ME !

# The directory containing the data files (e.g., passbands, sky background, etc.)
DATAPATH = join(__BASEPATH, "data")

# -------------------------------------------------------------------------------------- #
