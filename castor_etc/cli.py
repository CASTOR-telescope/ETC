#         GNU General Public License v3 (GNU GPLv3)
#
# (c) 2022.                            (c) 2022.
# Government of Canada                 Gouvernement du Canada
# National Research Council            Conseil national de recherches
# Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
# All rights reserved                  Tous droits réservés
#
# NRC disclaims any warranties,        Le CNRC dénie toute garantie
# expressed, implied, or               énoncée, implicite ou légale,
# statutory, of any kind with          de quelque nature que ce
# respect to the software,             soit, concernant le logiciel,
# including without limitation         y compris sans restriction
# any warranty of merchantability      toute garantie de valeur
# or fitness for a particular          marchande ou de pertinence
# purpose. NRC shall not be            pour un usage particulier.
# liable in any event for any          Le CNRC ne pourra en aucun cas
# damages, whether direct or           être tenu responsable de tout
# indirect, special or general,        dommage, direct ou indirect,
# consequential or incidental,         particulier ou général,
# arising from the use of the          accessoire ou fortuit, résultant
# software. Neither the name           de l'utilisation du logiciel. Ni
# of the National Research             le nom du Conseil National de
# Council of Canada nor the            Recherches du Canada ni les noms
# names of its contributors may        de ses  participants ne peuvent
# be used to endorse or promote        être utilisés pour approuver ou
# products derived from this           promouvoir les produits dérivés
# software without specific prior      de ce logiciel sans autorisation
# written permission.                  préalable et particulière
#                                      par écrit.
#
# This file is part of the             Ce fichier fait partie du projet
# FORECASTOR ETC project.              FORECASTOR ETC.
#
# FORECASTOR ETC is free software:     FORECASTOR ETC est un logiciel
# you can redistribute it and/or       libre ; vous pouvez le redistribuer
# modify it under the terms of         ou le modifier suivant les termes de
# the GNU General Public               la "GNU General Public
# License as published by the          License" telle que publiée
# Free Software Foundation,            par la Free Software Foundation :
# either version 3 of the              soit la version 3 de cette
# License, or (at your option)         licence, soit (à votre gré)
# any later version.                   toute version ultérieure.
#
# FORECASTOR ETC is distributed        FORECASTOR ETC est distribué
# in the hope that it will be          dans l'espoir qu'il vous
# useful, but WITHOUT ANY WARRANTY;    sera utile, mais SANS AUCUNE
# without even the implied warranty    GARANTIE : sans même la garantie
# of MERCHANTABILITY or FITNESS FOR    implicite de COMMERCIALISABILITÉ
# A PARTICULAR PURPOSE. See the        ni d'ADÉQUATION À UN OBJECTIF
# GNU General Public License for       PARTICULIER. Consultez la Licence
# more details.                        Générale Publique GNU pour plus
#                                      de détails.
#
# You should have received             Vous devriez avoir reçu une
# a copy of the GNU General            copie de la Licence Générale
# Public License along with            Publique GNU avec FORECASTOR ETC ;
# FORECASTOR ETC. If not, see          si ce n'est pas le cas, consultez :
# <http://www.gnu.org/licenses/>.      <http://www.gnu.org/licenses/>.


"""
cli.py

FORECASTOR ETC Command Line tool - built using the click module
"""

import click
import os
import sys
import platform
import numpy as np
import astropy
import scipy
import pandas
import castor_etc

@click.group()
@click.version_option(version=castor_etc.__version__)
def cli():
    """
    FORECASTOR ETC Command Line Interface.

    The CASTOR Exposure Time Calculator (ETC) helps users estimate the SNR 
    and performance of the CASTOR telescope across its UV, U, and G bands.
    This command line tool validates your installation, check local data assets 
    (passbands, PSFs), and gather information for bug reports.
    """
    pass

@cli.command()
def info():
    """Prints installation and system information for bug reports"""
    click.secho("--- Environment Info ---", fg="cyan", bold=True)
    click.echo(f"OS/Platform:    {platform.platform()}")
    click.echo(f"Python:         {sys.version.split()[0]}")
    click.echo(f"Executable:     {sys.executable}")
    
    click.secho("\n--- Package Info ---", fg="cyan", bold=True)
    click.echo(f"Version:        {castor_etc.__version__}")
    click.echo(f"Install Path:   {os.path.dirname(castor_etc.__file__)}")
    
    click.secho("\n--- Core Dependencies ---", fg="cyan", bold=True)
    click.echo(f"NumPy:          {np.__version__}")
    click.echo(f"SciPy:          {scipy.__version__}")
    click.echo(f"Astropy:        {astropy.__version__}")
    click.echo(f"Pandas:         {pandas.__version__}")

@cli.command()
def validate():
    """Verify data folder content to ensure data folders exist and have content

    NOTE: This only validates that local data assets exist in the package 
    directory. It does NOT verify if the package is correctly installed 
    in your Python environment. Use 'info' for installation details.
    """
    from castor_etc import verify_data_installation, DATAPATH
    
    click.secho(f"Validating data content in: {DATAPATH}", fg="yellow")
    
    results = verify_data_installation()
    
    if not results:
        click.secho("\n[!] Data directory is empty/not found.", fg="red", bold=True)
        click.echo(f"Expected path: {DATAPATH}")
        sys.exit(1)

    all_ok = True
    click.echo("--- Content Check ---")
    for folder, has_content in sorted(results.items()):
        if has_content:
            status = click.style("FOUND", fg="green")
            click.echo(f"  [{status}] {folder}")
        else:
            status = click.style("EMPTY", fg="red")
            click.echo(f"  [{status}] {folder}")
            all_ok = False
            
    if all_ok:
        click.secho("\nSuccess: local data assets found.", fg="green", bold=True)
    else:
        click.secho("\nWarning: missing local data assets", 
                    fg="red", bold=True)
        sys.exit(1)

if __name__ == "__main__":
    cli()