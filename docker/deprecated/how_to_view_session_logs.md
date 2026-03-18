# How to view notebook session logs on the CANFAR Science Platform (and other `skaha` terminal commands)

Viewing the session logs for your Docker notebook instance is invaluable for
troubleshooting and debugging any unexpected behaviour. This walkthrough is for accessing
the logs of an image pushed & launched on CANFAR. Note that accessing local Docker logs is
simply done via the command [`docker
logs`](https://docs.docker.com/engine/reference/commandline/logs/).

The steps below are summarized from: <https://github.com/opencadc/skaha/tree/master/doc>
and <https://github.com/opencadc/vostools/tree/master/vos>.

0. Ensure you have a [CANFAR](https://www.canfar.net/) account with access to the [Science
   Platform](https://github.com/opencadc/skaha/tree/master/doc#introduction-and-access).

1. Install [`vos`](https://github.com/opencadc/vostools/tree/master/vos):

   ```bash
   pip install vos
   ```

2. Generate a [proxy
   certificate](https://github.com/opencadc/skaha/tree/master/doc#proxy-certificates) for
   the `skaha` API:

   ```bash
   cadc-get-cert -u <CADC_username> --days-valid 30
   ```

   You will be required to input your CADC password.

3. To view the logs for a session, input the following command:

   ```bash
   curl -E ~/.ssl/cadcproxy.pem https://ws-uv.canfar.net/skaha/session/<sessionID>?view=logs
   ```

   where `<sessionID>` is the 8-character string you see in the URL after launching a
   notebook session. For example, if the URL is

   ```text
   https://ws-uv.canfar.net/notebook/abcde1f2/lab/tree/arc/home/IsaacCheng
   ```

   then the `<sessionID>` is `abcde1f2`.

   - Alternatively, you can get the `<sessionID>` from the terminal using:

     ```bash
     curl -E ~/.ssl/cadcproxy.pem https://ws-uv.canfar.net/skaha/session
     ```

     which will return a JSON list of parameters. The `<sessionID>` is the value
     corresponding to the `"id"` key.

     - The `<sessionID>` is also the
       [`?token=`](https://code.visualstudio.com/docs/datascience/jupyter-notebooks#_connect-to-a-remote-jupyter-server)
       parameter used to connect to remote Jupyter sessions. In fact, you can connect to a
       running notebook session without ever interacting with the Science Platform GUI
       simply through the `"connectURL"` field. This `"connectURL"` value is also the
       [URI](https://code.visualstudio.com/docs/datascience/jupyter-notebooks#_connect-to-a-remote-jupyter-server)
       VS Code uses to connect to a remote Jupyter server.

   - It is often useful to redirect the output (containing the session logs) to a file.
     This can be done via:

     ```bash
     curl -E ~/.ssl/cadcproxy.pem https://ws-uv.canfar.net/skaha/session/<sessionID>?view=logs > logs.txt
     ```

---

Following is a full list of commands taken from the
[`skaha` documentation](https://github.com/opencadc/skaha/tree/master/doc):

> To view all sessions and jobs:
>
> ```bash
> curl -E ~/.ssl/cadcproxy.pem https://ws-uv.canfar.net/skaha/session
> ```
>
> To view a single session or job:
>
> ```bash
> curl -E ~/.ssl/cadcproxy.pem https://ws-uv.canfar.net/skaha/session/<sessionID>
> ```
>
> To view logs for a session (this shows the complete output (stdout and stderr) for the
> image for the job):
>
> ```bash
> curl -E ~/.ssl/cadcproxy.pem https://ws-uv.canfar.net/skaha/session/<sessionID>?view=logs
> ```
>
>
>
> To view scheduling events for session (scheduling events will only be seen when there
> are issues scheduling the job on a node):
>
> ```bash
> curl -E ~/.ssl/cadcproxy.pem https://ws-uv.canfar.net/skaha/session/<sessionID>?view=events
> ```

---

## Starting and Stopping from the Terminal

- To start a session:

  ```bash
  curl -E ~/.ssl/cadcproxy.pem https://ws-uv.canfar.net/skaha/session \
       -d "name=<mySessionName>" \
       -d "image=images.canfar.net/<folderName>/<imageName>:<version>"
  ```

- To delete a session:

   ```bash
   curl -E ~/.ssl/cadcproxy.pem https://ws-uv.canfar.net/skaha/session/<sessionId> \
        -X DELETE
   ```
