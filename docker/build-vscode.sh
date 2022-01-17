#!/bin/bash
# Build script for the CASTOR exposure time calculator. Not meant to be used on CANFAR
echo "Building CASTOR exposure time calculator..."
#
# Set some parameters
#
VERSION=$(date +%y.%m.%d.%H%M)
# (following line from <https://stackoverflow.com/a/246128>)
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# (following line from <https://stackoverflow.com/a/8426110>)
REPO_DIR="$(dirname "$SCRIPT_DIR")"
#
# Load custom parameters
#
source Docker_env
#
# Build the project
#
cd ${REPO_DIR}  # necessary so Docker can access other folders within the repo
if [[ "$CUSTOMIZE_ENV" = "yes" ]]
then
    echo "Building with custom JupyterLab environment"
    docker build --build-arg NOTEBOOK_DIR=${NOTEBOOK_DIR} \
                 --build-arg NB_USER=${NB_USER} \
                 -t castor_etc:${VERSION} \
                 -f docker/Dockerfile.yesCustomEnv .
elif [[ "$CUSTOMIZE_ENV" == "no" ]]
then
    echo "Building with default JupyterLab environment"
    docker build -t castor_etc:${VERSION} \
                 -f docker/Dockerfile.noCustomEnv .
else
    echo "ERROR: CUSTOMIZE_ENV is must be 'yes' or 'no'"
    exit 1
fi
#
echo "Finishing building castor_etc:${VERSION}"
echo "Now running castor_etc_v${VERSION}..."
#
# Run the project
#
docker run --interactive \
        --ip 0.0.0.0 \
        --rm \
        --tty \
        --env DISPLAY=host.docker.internal:0 \
        -p 8888:8888 \
        -v ${REPO_DIR}:${NOTEBOOK_DIR} \
        -v /arc/home/${NB_USER}/ETC_plots:/arc/home/${NB_USER}/ETC_plots \
        -v ${REPO_DIR}/.vscode/extensions:/home/${NB_USER}/.vscode-server/extensions \
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
