"""
CASTOR Exposure Time Calculator (ETC)
=====================================

A modular, user-friendly Python package for easy estimation and analysis of CASTOR imaging
performance. See the [`ETC_frontend`](https://github.com/CASTOR-telescope/ETC_frontend)
GitHub repository for a graphical user interface to complement this package.

Includes:
  1. Astronomical source generation and background noise estimation
  2. Telescope imaging chain simulation, featuring a pixel-based photometry approach
  3. Convenience functions for converting between useful quantities (e.g., flux to
     electron/s to AB magnitude)

Author: Isaac Cheng
Contact: isaac.cheng.ca@gmail.com

---

        GNU General Public License v3 (GNU GPLv3)

(c) 2022.                            (c) 2022.
Government of Canada                 Gouvernement du Canada
National Research Council            Conseil national de recherches
Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
All rights reserved                  Tous droits réservés

NRC disclaims any warranties,        Le CNRC dénie toute garantie
expressed, implied, or               énoncée, implicite ou légale,
statutory, of any kind with          de quelque nature que ce
respect to the software,             soit, concernant le logiciel,
including without limitation         y compris sans restriction
any warranty of merchantability      toute garantie de valeur
or fitness for a particular          marchande ou de pertinence
purpose. NRC shall not be            pour un usage particulier.
liable in any event for any          Le CNRC ne pourra en aucun cas
damages, whether direct or           être tenu responsable de tout
indirect, special or general,        dommage, direct ou indirect,
consequential or incidental,         particulier ou général,
arising from the use of the          accessoire ou fortuit, résultant
software. Neither the name           de l'utilisation du logiciel. Ni
of the National Research             le nom du Conseil National de
Council of Canada nor the            Recherches du Canada ni les noms
names of its contributors may        de ses  participants ne peuvent
be used to endorse or promote        être utilisés pour approuver ou
products derived from this           promouvoir les produits dérivés
software without specific prior      de ce logiciel sans autorisation
written permission.                  préalable et particulière
                                     par écrit.

This file is part of the             Ce fichier fait partie du projet
FORECASTOR ETC project.              FORECASTOR ETC.

FORECASTOR ETC is free software:     FORECASTOR ETC est un logiciel
you can redistribute it and/or       libre ; vous pouvez le redistribuer
modify it under the terms of         ou le modifier suivant les termes de
the GNU General Public               la "GNU General Public
License as published by the          License" telle que publiée
Free Software Foundation,            par la Free Software Foundation :
either version 3 of the              soit la version 3 de cette
License, or (at your option)         licence, soit (à votre gré)
any later version.                   toute version ultérieure.

FORECASTOR ETC is distributed        FORECASTOR ETC est distribué
in the hope that it will be          dans l'espoir qu'il vous
useful, but WITHOUT ANY WARRANTY;    sera utile, mais SANS AUCUNE
without even the implied warranty    GARANTIE : sans même la garantie
of MERCHANTABILITY or FITNESS FOR    implicite de COMMERCIALISABILITÉ
A PARTICULAR PURPOSE. See the        ni d'ADÉQUATION À UN OBJECTIF
GNU General Public License for       PARTICULIER. Consultez la Licence
more details.                        Générale Publique GNU pour plus
                                     de détails.

You should have received             Vous devriez avoir reçu une
a copy of the GNU General            copie de la Licence Générale
Public License along with            Publique GNU avec FORECASTOR ETC ;
FORECASTOR ETC. If not, see          si ce n'est pas le cas, consultez :
<http://www.gnu.org/licenses/>.      <http://www.gnu.org/licenses/>.
"""

from setuptools import setup

long_description = __doc__.strip()  # Remove leading and trailing newlines

setup(
    name="castor_etc",
    version="1.3.2",  # see semantic versioning (<https://semver.org/spec/v2.0.0.html>)
    description="CASTOR Exposure Time Calculator (ETC)",
    long_description=long_description,
    url="https://github.com/CASTOR-telescope/ETC",
    author="FORECASTOR Team",
    # author_email="isaac.cheng.ca@gmail.com",
    # maintainer="Isaac Cheng",
    # maintainer_email="isaac.cheng.ca@gmail.com"
    packages=[
        "castor_etc",
        "castor_etc.data",
        "castor_etc.data.UVMOS_data",
        "castor_etc.data.galaxy_spectra",
        "castor_etc.data.passbands",
        "castor_etc.data.pickles_spectra",
        "castor_etc.data.UVMOS_data",
        "castor_etc.data.transit_data",
        "castor_etc.data.grism_data",
        "castor_etc.data.psfs",
        "castor_etc.data.sky_background",
    ],
    package_data={
        "castor_etc.data.galaxy_spectra": ["*.txt"],
        "castor_etc.data.grism_data": [
            "*_dispersion_u.txt",
            "*_dispersion_uv.txt",
            "*_efficiency_.u.txt",
            "*_efficiency_.uv.txt",
            "*_profile_uv.txt",
        ],
        "castor_etc.data.passbands": ["*.uv", "*.u", "*.g"],
        "castor_etc.data.pickles_spectra": ["dat/*.dat"],  # must use forward slash
        "castor_etc.data.psfs": ["*.fits"],
        "castor_etc.data.sky_background": ["*.fits", "*.txt"],
        "castor_etc.data.transit_data": [
            "*.txt",
            "instrument_data/*.csv",
            "instrument_data/transmission_functions/*.dat",
            "stellar_models/*.txt",
        ],
        "castor_etc.data.UVMOS_data": ["*.dat", "*.txt"],
    },
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "astropy",
        "pandas",
        "photutils",
        "tqdm",
        "scikit-image",
        "astroquery",
        "pytransit",
        "arviz",
        "celerite",
        "emcee",
        "corner",
        "spectres",
    ],  # Packages listed after the 'pytransit' package and before the 'spectres' package are pre-requisites to run the 'pytransit' package.
    license="GPLv3",
    python_requires=">=3.9",
    platforms=["Linux"],  # only tested on Ubuntu. MacOS and Windows likely okay.
)
