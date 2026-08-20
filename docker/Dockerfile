FROM images.canfar.net/skaha/astroml:latest

LABEL maintainer="Charles Lee <leec16@myumanitoba.ca>"
LABEL description="FORECASTOR analysis environment"

# Install system dependencies as root
USER root

# Build argument for wheel file
ARG WHEEL_URL=https://github.com/CASTOR-telescope/ETC/releases/download/v2.0.0/castor_etc-2.0.0-py3-none-any.whl

# Install from wheel file
RUN pip install --no-cache-dir "${WHEEL_URL}"