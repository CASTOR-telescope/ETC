#!/bin/bash
# Build script for the CASTOR exposure time calculator. Not meant to be used on CANFAR
echo "Building CASTOR exposure time calculator..."
#
# Set some parameters
#
VERSION=$(date +%y.%m.%d.%H%M)
# (following line from <https://stackoverflow.com/a/246128>)
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#
# Load custom parameters
#
source Docker_env
#
# Build the project
#
docker build --build-arg NOTEBOOK_DIR=${NOTEBOOK_DIR} \
             --build-arg NB_USER=${NB_USER} \
             -t castor_etc:${VERSION} \
             -f ${SCRIPT_DIR}/Dockerfile.${CUSTOMIZE_ENV}CustomEnv .
#
echo "Finishing building castor_etc:${VERSION}"
echo "Now running castor_etc_v${VERSION}..."
#
# Run the project
#
docker run --interactive \
        --ip 0.0.0.0 \
        --rm \
        --env DISPLAY=host.docker.internal:0 \
        -p 8888:8888 \
        -v ${SCRIPT_DIR}:${NOTEBOOK_DIR} \
        -v /arc/home/IsaacCheng/ETC_plots:/arc/home/IsaacCheng/ETC_plots \
        -v ${SCRIPT_DIR}/.vscode/extensions:/home/IsaacCheng/.vscode-server/extensions \
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
