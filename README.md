# CASTOR Exposure Time Calculator (ETC)

Isaac Cheng - 2022

**This package is under active development. It should be stable but please let me know if
there are any bugs or unusual behaviour!**

This Python package is an exposure time calculator for
[CASTOR](https://www.castormission.org/). Since CASTOR is currently in Phase 0, this tool
contains tunable parameters that are not normally exposed in most other exposure time
calculators; the user is highly encouraged to read the docstrings and view the example
notebooks prior to using this software. Likewise, there may be frequent updates to this
package as the mission matures (see the [changelog](CHANGELOG.md) for more details).

## Table of Contents

1. [Installation](#installation)
2. [Building the Docker image for CANFAR/local machine](docker/README.md)
3. [Viewing CANFAR Jupyter notebook session logs and Other Terminal
   Commands](docker/how_to_view_session_logs.md)
4. Getting started and specific use-case examples: See the
   [ETC_notebooks](https://github.com/CASTOR-telescope/ETC_notebooks) repository
   containing examples in Jupyter notebooks.
5. The browser-based graphical user interface to complement this ETC is located in the
   [ETC_frontend](https://github.com/CASTOR-telescope/ETC_frontend) repository.
6. [Known issues](#known-issues)
7. [Roadmap](#roadmap)
8. [Questions & Issues](#questions--issues)

## Installation

If you do not wish to clone the repository, the easiest way to install the `castor_etc`
package is via:

```bash
pip install git+https://github.com/CASTOR-telescope/ETC.git
```

Note that executing this command in a conda environment (i.e., by building from the
[environment file](docker/castor_etc_env.yml)) is recommended.

Alternatively, after cloning this repository via:

```bash
git clone https://github.com/CASTOR-telescope/ETC.git
```

You can then install the ETC package using one of the following commands, which should be
executed while inside the repository folder (i.e., `ETC/`).

To install the `castor_etc` package "normally", use either:

```bash
pip install .
```

or

```bash
python setup.py install
```

If you want to install in [develop
mode](https://pip-python3.readthedocs.io/en/latest/reference/pip_install.html#install-editable),
use:

```bash
pip install -e .
```

In develop mode, the installation of the package simply links back to the
[`castor_etc`](castor_etc/) folder itself, meaning any changes made to this package will
be reflected in your environment.

## Known Issues

- Rectangular aperture sometimes produces the incorrect number of pixels in the photometry
  aperture. I believe this is actually a bug in `astropy` since the elliptical aperture
  does not suffer from this problem and the aperture weights are clearly wrong (see the
  photometry source code for an example that highlights this flaw).
  - The package will raise a warning when the number of aperture pixels is incorrect
    (i.e., when it deviates by more than 0.1 pixels from the true value). When this
    happens, either slightly change the aperture dimensions or avoid using the rectangular
    aperture for now.

## Roadmap

Here are some future plans for the ETC:

- Use a sampled PSF to create a pixel-based kernel to convolve with the simulated images.
  Additionally, point sources should use the sampled PSF to calculate its encircled energy
  (we currently assume the point source PSF is a Gaussian).
- Add slitless grism and multi-object spectroscopic capabilities to the ETC

## Questions & Issues

Please reach out if you have any questions, either through email
([isaac.cheng.ca@gmail.com](mailto:isaac.cheng.ca@gmail.com)) or the [discussions
page](https://github.com/CASTOR-telescope/ETC/discussions). You can also ping me on Slack
or even set up an online video/audio call! Larger issues or feature requests can be posted
and tracked via the [issues page](https://github.com/CASTOR-telescope/ETC/issues).
