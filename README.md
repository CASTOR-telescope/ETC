# CASTOR Exposure Time Calculator (ETC)

Isaac Cheng - January 2022

WARNING: Everything is super messy at the moment because I am stuck in debugging hell...

TODO: Figure out a way to ensure consistent absolute paths for project on CANFAR & local
systems...

## Local Build Instructions

To build the project _locally_ (i.e., not on [CANFAR](https://www.canfar.net/en/)), do the
following:

1. Download the git repo:

   ```bash
    git clone https://github.com/CASTOR-telescope/ETC.git
   ```

2. Ensure you have [Docker](https://docs.docker.com/get-started/) installed.

3. Modify the variables in [Docker_env](Docker_env) to your desired values. Below is a
   brief explanation of the parameters:

   ```bash
   # The first 3 variables + the `VERSION` parameter affects the Docker image (i.e., they
   # are relevant to `docker build`) and the Docker containers.
   # The other settings only affect the individual containers (i.e., they are relevant to
   # `docker run`).

   NB_USER=IsaacCheng  # the username for the notebook. I recommend setting is equal to
                       # your CANFAR username
   NOTEBOOK_DIR="/arc/home/IsaacCheng/ETC"  # the repo bind mount destination. I recommend
                                            # setting this path to be the same path as if
                                            # you ran this repo on CANFAR

   CUSTOMIZE_ENV=yes  # yes or no. If yes, apply my custom notebook configuration

   JUPYTER_ENABLE_LAB=yes  # yes or no. If yes, use JupyterLab
   JUPYTER_TOKEN=""  # the token (similar to a password) for the JupterLab instance
   GRANT_SUDO=yes  # yes or no. If yes, give sudo privileges to notebook user
   CHOWN_HOME=yes  # yes or no. Must be yes if NB_USER is not "jovyan"
   CHOWN_HOME_OPTS="-R"  # must be "-R" if CHOWN_HOME=yes

   # If the following parameters exist, will overwrite the default values in the build script
   VERSION=0.1  # version tag for the Docker image
   SCRIPT_DIR="/arc/home/IsaacCheng/ETC"  # the absolute path of this repo on local machine
   ```

   - You may also wish to save outputs and plots to a separate directory (i.e., not a
     subfolder in this repo), in which case you should add a bind mount in
     [build.sh](build.sh) and modify the `OUTPATH` variable in
     [`constants.py`](src/constants.py) to the proper mounted path. For example, add:

     ```bash
     -v /arc/local_directory/ETC_plots:/container_directory/ETC_plots
     ```

     to your [build.sh](build.sh) file and change `OUTPATH` in
     [`constants.py`](src/constants.py) to "`/container_directory/ETC_plots/`".

<!-- 3. Open the [Dockerfile](Dockerfile) and modify the `WORKDIR` value to be whichever path
   you would like this repo to be contained in. I recommend setting this path to be the
   same path as if you ran this repo on CANFAR.

   In the same [Dockerfile](Dockerfile), modify the `USERNAME` variable to equal your
   CANFAR username. -->

<!-- 4. Open [build.sh](build.sh) and modify the line containing "`/arc/home/IsaacCheng/ETC`"
   to be the same value as the `WORKDIR` variable from the [Dockerfile](Dockerfile). This
   line creates a [bind mount](https://docs.docker.com/storage/bind-mounts/) between this
   repository and the "virtual repository" in the Docker container so that any changes you
   make to the repository files in the Docker container will be reflected outside the
   container.

   - You may also wish to save outputs and plots to a separate directory (i.e., not a
     subfolder in this repo), in which case you should add a bind mount in
     [build.sh](build.sh) and modify the `OUTPATH` variable in
     [`constants.py`](src/constants.py) to the proper mounted path. For example, add:

     ```bash
     -v /arc/local_directory/ETC_plots:/container_directory/ETC_plots
     ```

     to your [build.sh](build.sh) file and change `OUTPATH` in
     [`constants.py`](src/constants.py) to "`/container_directory/ETC_plots/`". -->

4. Run the following command to build the image. It should automatically mount this
   git directory in the Docker container.

   ```bash
    ./build.sh
   ```

   The first time building this may take a while since the reference image I am using
   (`jupyter/scipy-notebook:hub-2.0.1`) creates a Docker container that is almost 3 GB!

   (In case you're wondering, the `build-vscode.sh` file is for my own use with VS Code
   that contains some paths + settings specific to my local machine.)

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

Please reach out if you have any questions about this; my email is
[isaac.cheng.ca@gmail.com](mailto:isaac.cheng.ca@gmail.com). You can also ping me on Slack
or even set up an online video/audio call! Larger issues or feature requests can be posted
and tracked via the [issues page](https://github.com/CASTOR-telescope/ETC/issues).

## Inside the Docker Container

Please modify the `OUTPATH` variable in [`constants.py`](src/constants.py) to the proper
directory. Also remember you can copy files and folders within a container to the local
filesystem using the following command (TIP: use bash autocompletion to fill in the
version number):

```bash
docker cp castor_etc_v<VERSION>:<DOCKER_PATH_TO_COPY> <LOCAL_DESTINATION_PATH>
```
