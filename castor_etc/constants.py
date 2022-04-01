"""
Useful constants. (N.B. use cgs units when possible)

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

import astropy.units as u

# ----------------------------------- SOME CONSTANTS ----------------------------------- #

# Speed of light in vacuum
LIGHTSPEED = 2.99792458e10 << (u.cm / u.s)  # cm/s

# Planck's constant
PLANCK_H = 6.62607015e-27 << (u.erg * u.s)  # erg.s

# Boltzmann's constant
K_B = 1.380649e-16 << (u.erg / u.K)  # erg/K

# The radius of the Sun
SUN_RADIUS = 6.96340e10 << u.cm  # cm

# Solar luminosity
SUN_LUMINOSITY = 3.828e33 << (u.erg / u.s)  # erg/s

# Coefficient for the conversion from flam (erg/cm^2/s/A) to photlam (photon/cm^2/s/A).
# Reference: <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf>
FLAM_TO_PHOTLAM = 5.0341165675427094e7  # photon/erg/A

# Coefficient for the conversion from fnu (erg/cm^2/s/Hz) to photlam (photon/cm^2/s/A).
# Reference: <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf>
FNU_TO_PHOTLAM = 1.5091901796421519e26  # photon.Hz/erg

# Coefficient for the conversion from fnu (erg/cm^2/s/Hz) to flam (erg/cm^2/s/A).
# Reference: <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf>
FNU_TO_FLAM = 2.99792458e18  # Hz/A

# Coefficient for the conversion from flam (erg/cm^2/s/A) to fnu (erg/cm^2/s/Hz).
# Reference: <https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf>
FLAM_TO_FNU = 3.3356409519815204e-19  # A/Hz

# 1 steradian to 1 arcsec^2
SR_TO_SQARCSEC = 4.254517029615221e10  # arcsec^2/sr

# 1 parsec
PC = 3.0856775815e18 << u.cm  # cm

# -------------------------------------------------------------------------------------- #
