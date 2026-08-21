#!/bin/bash

# FORECASTOR local container build tools
#
# This is build script for any users hoping to build a container locally. Please note that this is NOT meant to be used for CANFAR builds.

## User Parameters
VERSION=$(date +%y.%m.%d)       # container version
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # get current working directory for this build script (https://stackoverflow.com/a/246128)
REPO_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"         # https://stackoverflow.com/a/8426110
# Alternatively, we can also run:
# git rev-parse --show-toplevel

## Stellar Model files

# This defines where the Docker container will store the stellar models ; by default, CASTOR ETC searchs for it under 3 different places:
#           Variation 1: stellar_model_dir = "/arc/projects/CASTOR/stellar_models" --> Default path (Working in the CANFAR server).
#   Variation 2: stellar_model_dir = <path to local stellar models directory>
#   Variation 3: stellar_model_dir = join(DATAPATH,"transit_data/stellar_models") --> This path should be used when building docker container locally.
DOCKER_STELLAR_MODEL_DIR="/opt/FORECASTOR/stellar_models"

## Custom Environment Parameters
CUSTOMIZE_ENV=no        # yes if custom; no otherwise
NB_USER=IsaacCheng      # Notebook username
NOTEBOOK_DIR=/arc/home/IsaacCheng/CASTOR/ETC        # Notebook directory
STELLAR_MODEL_DIR=""        # Change this to match where the stellar model files are installed locally
JUPYTER_ENABLE_LAB=yes
JUPYTER_TOKEN=""
GRANT_SUDO=yes
CHOWN_HOME=yes
CHOWN_HOME_OPTS="-R"

## Docker builds
cd ${REPO_DIR}  # necessary so Docker can access other folders within the repo
if [[ "$CUSTOMIZE_ENV" = "yes" ]]
then
    echo "Building with custom JupyterLab environment"
    docker build --build-arg NOTEBOOK_DIR=${NOTEBOOK_DIR} \
                 --build-arg NB_USER=${NB_USER} \
                 -t castor_etc:${VERSION} \
                 --build-arg CACHEBUST=$(date +%s) \
                 -f docker/local/Dockerfile.yesCustomEnv .
elif [[ "$CUSTOMIZE_ENV" = "no" ]]
then
    echo "Building with default JupyterLab environment"
    docker build -t castor_etc:${VERSION} \
                 --build-arg CACHEBUST=$(date +%s) \
                 -f docker/local/Dockerfile.noCustomEnv .
else
    echo "ERROR: CUSTOMIZE_ENV is must be yes or no"
    exit 1
fi
echo "Finishing building castor_etc:${VERSION}"

# Run the project
echo "Now running castor_etc_v${VERSION}..."

docker run --interactive \
        --rm \
        --tty \
        --env DISPLAY=host.docker.internal:0 \
        -p 8888:8888 \
        -v ${REPO_DIR}:${NOTEBOOK_DIR} \
        -v ${STELLAR_MODEL_DIR}:${DOCKER_STELLAR_MODEL_DIR} \
        --env JUPYTER_ENABLE_LAB=${JUPYTER_ENABLE_LAB} \
        --env JUPYTER_TOKEN=${JUPYTER_TOKEN} \
        --env NB_USER=${NB_USER} \
        --env CHOWN_HOME=${CHOWN_HOME} \
        --env CHOWN_HOME_OPTS=${CHOWN_HOME_OPTS} \
        --env GRANT_SUDO=${GRANT_SUDO} \
        --workdir ${NOTEBOOK_DIR} \
        --user root \
        --name castor_etc_v${VERSION} \
        -d castor_etc:${VERSION}

#
# Print the JupyterLab URL
#
# (wait until URL is generated. Also route stderr to stdout)
while ! docker logs castor_etc_v${VERSION} 2>&1 | grep -q "or http*"
do
    # ! FIXME: this does not catch any errors!
    docker logs castor_etc_v${VERSION}
    if docker logs castor_etc_v${VERSION} 2>&1 | grep -q "*Error*"
    then
        echo "Error: JupyterLab failed to start"
        exit 1
    fi
    sleep 1
done
# (output log containing URL)
docker logs castor_etc_v${VERSION}
#
echo "DONE! Use the URL above to access the JupyterLab instance for castor_etc_v${VERSION}."
