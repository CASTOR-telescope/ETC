# CASTOR Exposure Time Calculator (ETC)

Isaac Cheng - January 2022

## Local Build Instructions

To build the project _locally_ (i.e., not on [CANFAR](https://www.canfar.net/en/)), do the
following:

1. Download the git repo:

   ```bash
    git clone https://github.com/CASTOR-telescope/ETC.git
   ```

2. Ensure you have [Docker](https://docs.docker.com/get-started/) installed.

3. Run the following command to build the image. It should automatically mount this
   git directory in the resulting container.

   ```bash
    ./build.sh
   ```

   The first time building this may take a while since the reference image I am using
   (`jupyter/scipy-notebook:hub-2.0.1`) creates a Docker container that is almost 3 GB!

   (In case you're wondering, the `build-vscode.sh` file is for my own use with VS Code
   that contains some paths + settings specific to my local machine.)

4. Copy & paste the URL from Docker (e.g., something like
   `http://127.0.0.1:8888/lab?token=<RANDOM STRING>`)
   into a web browser to access JupyterLab. Note that any changes you make to the `ETC`
   directory in the Docker container _will_ be reflected in your local `ETC` directory!
   You also have sudo privileges and root access if you build using the script above.

5. Once you are done, shutdown JupyterLab (i.e., using the `File`>`Shutdown` buttons near
   the top left corner of the menu bar).

6. When JupyterLab has successfully shutdown, run the following command (TIP: use bash
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
