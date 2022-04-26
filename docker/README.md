# Docker Setup

Brief walkthrough on how to prepare the `castor_etc` package in a Jupyter lab Docker
container.

## CANFAR Build Instructions

1. Download the git repo:

   ```bash
   git clone https://github.com/CASTOR-telescope/ETC.git
   ```

2. Ensure you have [Docker](https://docs.docker.com/get-started/) installed.

3. To build an image ready to be deployed on [CANFAR](https://www.canfar.net/en/), simply
   run the following command _from this repository's top-level directory (i.e., `ETC/`)_.
   The requirement to run from the top-level directory is so that the Dockerfile can
   access files throughout this repo and is not limited to those in the [docker](./)
   folder.

   ```bash
   docker build -t castor_etc:<VERSION> \
                --build-arg CACHEBUST=$(date +%s) \
                -f docker/Dockerfile.noCustomEnv .
   ```

   where `<VERSION>` is the version you would like to tag the image with. The `CACHEBUST`
   argument ensures that the newest version of the `castor_etc` package is always being
   installed.

4. Then follow the instructions detailed on the skaha GitHub for [session
   containers](https://github.com/opencadc/skaha/tree/master/containers#publishing-skaha-containers).
   Remember to tag the pushed image as `notebook` on [Harbor](https://images.canfar.net)
   to be able to access it via the Science Portal drop-down menu!

## Local Build Instructions

To build the project _locally_ (i.e., not on [CANFAR](https://www.canfar.net/en/)), do the
following:

1. Download the git repo:

   ```bash
   git clone https://github.com/CASTOR-telescope/ETC.git
   ```

2. Ensure you have [Docker](https://docs.docker.com/get-started/) installed.

3. Modify the variables in [Docker_env](./Docker_env) to your desired values. Below
   is a brief explanation of the parameters:

   ```bash
   # The first 2 variables + the `VERSION` parameter affects the Docker image (i.e., they
   # are relevant to `docker build`) and the Docker containers.
   # The other settings only affect the individual containers (i.e., they are relevant to
   # `docker run`).

   NB_USER=IsaacCheng  # the username for the notebook. I recommend setting is equal to
                       # your CANFAR username. Ignored if CUSTOMIZE_ENV=no
   NOTEBOOK_DIR="/arc/home/IsaacCheng/CASTOR/ETC"  # the repo bind mount destination. I
                                                   # recommend setting this path to be the
                                                   # same path as if you ran this repo on
                                                   # CANFAR. Ignored if CUSTOMIZE_ENV=no

   CUSTOMIZE_ENV=yes  # yes or no. If yes, apply my custom notebook configuration

   JUPYTER_ENABLE_LAB=yes  # yes or no. If yes, use JupyterLab
   JUPYTER_TOKEN=""  # the token (similar to a password) for the JupterLab instance
   GRANT_SUDO=yes  # yes or no. If yes, give sudo privileges to notebook user
   CHOWN_HOME=yes  # yes or no. Must be yes if NB_USER is not "jovyan"
   CHOWN_HOME_OPTS="-R"  # must be "-R" if CHOWN_HOME=yes

   # If the following parameters exist, will overwrite the default values in the build script
   VERSION=0.1  # version tag for the Docker image
   ```

4. Run the [build.sh](./build.sh) script to build the image (i.e., run `./build.sh`).
   It should automatically determine the version and the locations of the repository and
   script (i.e., you can execute the script from any directory, not just the
   repository-level directory).

   The first time building this may take a while since the reference image I am using
   (`jupyter/scipy-notebook:hub-2.0.1`) creates a Docker image that is about 3 GB! Note
   that the actual size of the pushed image (i.e., on [Harbor](https://images.canfar.net))
   is much smaller (only around 1 GB).

5. Copy & paste the URL from Docker (e.g., something like
   `http://127.0.0.1:8888/lab?token=<RANDOM STRING>`)
   into a web browser to access JupyterLab. Note that any changes you make to the `ETC`
   directory in the Docker container _will_ be reflected in your local `ETC` directory!
   You also have sudo privileges and root access if you build using the script above.

6. Once you are done, shutdown JupyterLab (i.e., using the `File`>`Shutdown` buttons near
   the top left corner of the menu bar).

7. When JupyterLab has successfully shutdown, run the following command (TIP: use bash
   autocompletion to fill in the version number):

   ```bash
   docker image rm castor_etc:<VERSION>
   ```

   where `<VERSION>` is some string like `22.01.05.1234`. You can also check the version
   number by looking at the last line of the output from `./build.sh`, which should say
   something like "DONE! Use the URL above to access the JupyterLab instance for
   castor_etc_v22.01.05.1234".

Please reach out if you have any questions about this, either through email
([isaac.cheng.ca@gmail.com](mailto:isaac.cheng.ca@gmail.com)) or the [discussions
page](https://github.com/CASTOR-telescope/ETC/discussions). You can also ping me on Slack
or even set up an online video/audio call! Larger issues or feature requests can be posted
and tracked via the [issues page](https://github.com/CASTOR-telescope/ETC/issues).

### Tip: Copying Files Out of a Docker Container

Also remember you can copy files and folders within a container to the local
filesystem using the following command (TIP: use bash autocompletion to fill in the
version number):

```bash
docker cp castor_etc_v<VERSION>:<DOCKER_PATH_TO_COPY> <LOCAL_DESTINATION_PATH>
```

## Viewing the Docker Logs from CANFAR and Other Terminal Commands

Please see this [howto](./how_to_view_session_logs.md) document for instructions.
