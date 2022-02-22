"""
CASTOR Exposure Time Calculator (ETC)
=====================================

A Python package for easy analysis of CASTOR performance and modifications to CASTOR
parameters in Python scripts. See the `ETC_frontend` GitHub repository for a graphical
user interface to complement this package.

Includes:
  1. Astronomical source generation and background noise estimation
  2. Telescope imaging chain simulation (both photometry and spectroscopy)
  3. Convenience functions for converting between useful quantities (e.g., flux to
     electron/s to AB magnitude)

Copyright 2022, CASTOR Mission Team
Author: Isaac Cheng
Contact: isaac.cheng.ca@gmail.com

---

GNU General Public License v3 (GNU GPLv3)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
"""

# REVIEW: Untested!

from setuptools import setup

long_description = __doc__.strip()  # Remove leading and trailing newlines

setup(
    name="castor_etc",
    version="0.1dev",
    description="CASTOR Exposure Time Calculator (ETC)",
    long_description=long_description,
    url="https://github.com/CASTOR-telescope/ETC",
    author="Isaac Cheng",
    author_email="isaac.cheng.ca@gmail.com",
    # maintainer="Isaac Cheng",
    # maintainer_email="isaac.cheng.ca@gmail.com"
    packages=["castor_etc"],
    install_requires=["numpy", "scipy", "astropy", "pandas", "photutils"],
    license="GPLv3",
    python_requires=">=3.9",
    platforms=["Linux"],  # only tested on Ubuntu. MacOS may be okay, Windows probably not
)
