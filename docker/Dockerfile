FROM quay.io/jupyter/scipy-notebook:hub-5.2.1

LABEL maintainer="Charles Lee <leec16@myumanitoba.ca>"
LABEL description="FORECASTOR analysis environment"

# Install system dependencies as root
USER root

# Build argument for wheel file
ARG VERSION="2.0.0"
ARG WHEEL_URL=https://github.com/CASTOR-telescope/ETC/releases/download/v2.0.0/castor_etc-2.0.0-py3-none-any.whl

# Set to Notebook user
USER ${NB_UID}

# Install required Python packages
RUN pip install astropy ipyevents pylint black[jupyter] jupyterlab-git notebook==6.5.7
RUN jupyter nbextension enable --py --sys-prefix ipyevents

# Install from wheel file
RUN pip install --no-cache-dir "${WHEEL_URL}"