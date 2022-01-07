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
# Build the project
#
docker build -t castor_etc:${VERSION} -f ${SCRIPT_DIR}/Dockerfile .
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
        -v ${SCRIPT_DIR}:/arc/projects/CASTOR/ETC \
        -v ${SCRIPT_DIR}/.vscode/extensions:/home/jovyan/.vscode-server/extensions \
        --env JUPYTER_ENABLE_LAB=yes \
        --env JUPYTER_TOKEN="" \
        --env GRANT_SUDO=yes \
        --user root \
        --name castor_etc_v${VERSION} \
        -d castor_etc:${VERSION}
#
# Print the JupyterLab URL
#
# (wait until URL is generated. Also route stderr to stdout)
while ! docker logs castor_etc_v${VERSION} 2>&1 | grep -q "or http*";
do
    sleep 1
done
# (output log containing URL)
docker logs castor_etc_v${VERSION}
#
echo "DONE! Use the URL above to access the JupyterLab instance for castor_etc_v${VERSION}."
