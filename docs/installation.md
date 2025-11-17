# Installation

## Pre-Requisites
- Python 3.8 or higher
- pip package manager

## Install directly from GitHub

Currently, the fastest way to install the package is done via GitHub.

```Bash
pip install git+https://github.com/CASTOR-telescope/ETC.git
```

## Install from source

You can obtain the full pacakge source code either through [cloning from GitHub](https://github.com/CASTOR-telescope/ETC) or [downloading the zipped source code from releases](https://github.com/CASTOR-telescope/ETC/releases).

Once downloaded, navigate to the `/ETC` folder and install the package via the following command:

```Bash
pip install .
```

You can also attempt to build the wheel files locally by running the following commands:

```Bash
pip install --upgrade build
python -m build
```

### Install in develop mode
If you want to install in
[develop mode](https://pip-python3.readthedocs.io/en/latest/reference/pip_install.html#install-editable),
use:

```bash
pip install -e .
```

In develop mode, the installation of the package simply links back to the
[`castor_etc`](castor_etc/) folder itself, meaning any changes made to this package will
be reflected in your environment.

## Additional installs

### Transit models 

If you are doing transit simulations, then after installing the Python package, you will need to
download these [stellar models](https://kona.ubishops.ca/jsikora/poet_stellar_models.tar.gz) into a
directory. When you are using `getStarData()` or `use_gaia_spectrum()` in
[spectrum.py](castor_etc/spectrum.py) and when you are using the `Observation` class in
[transit.py](castor_etc/transit.py), you will need to set the `stellar_model_dir` parameter in
these functions/class to the directory containing the stellar models on your local machine.
