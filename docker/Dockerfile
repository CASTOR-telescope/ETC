FROM images.canfar.net/skaha/astroml:latest

# Container metadata
LABEL maintainer="Charles Lee <leec16@myumanitoba.ca"
LABEL description="FORECASTOR Exposure Time Calculator"
LABEL version="1.0.0"

# Install system dependencies from root
USER root

# Update system packages
RUN apt-get update && apt-get install -y system-packages

# Install UV
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

# Install ETC from source
RUN git clone https://github.com/CASTOR-telescope/ETC.git /tmp/tools && \
    cd /tmp/tools && \
    uv pip install --no-cache-dir . && \
    rm -rf /tmp/tools

