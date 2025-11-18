CANFAR Science Portal
=====================

What is CANFAR?
----------------
`Canadian Advanced Network for Astronomical Research (CANFAR) <https://www.canfar.net/en/>`_ is a platform for Canadian researchers in astronomy and operated by the Canadian Astronomy Data Centre (CADC) and the Digital Research Alliance of Canada (DRAC).


It is highly recommended for any CASTOR developer affiliated with Canadian astronomy to utilize CANFAR for developing and testing of CASTOR ETC.

Launching on CANFAR
-------------------

#. Ensure you have a Canadian Astronomy Data Centre account (or
    `request one <https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/auth/request.html>`_ if you
    do not have one yet).
#. Go to `CANFAR <https://www.canfar.net/en/>`_ and sign in to the
    `Science Portal <https://www.canfar.net/science-portal/>`_. If you cannot access this,
    then you must send an email to `support@canfar.net <mailto:support@canfar.net>`_
    requesting access to the Science Portal.
#. Inside the `Science Portal <https://www.canfar.net/science-portal/>`_, click the "``+``"
    icon to launch a new session. Under "``type``", select "``notebook``". If multiple
    ``castor_etc`` versions are available, you can select the specific version you would like
    to use under the "``container image``" field. The version number is denoted by the string
    following the colon (e.g., ``images.canfar.net/castor/castor_etc:1.0.0`` means version
    ``1.0.0`` of the ``castor_etc`` notebook image).
#. Assign a name to your Jupyter notebook container and choose the maximum amount of
    memory (RAM) and maximum number of CPU cores you would like to have available for your
    notebook container. Note that RAM and CPU are shared resources. A reasonable amount to
    use, for example, would be 16 GB of RAM and 4 CPU cores. Larger RAM and CPU requests
    are not uncommon, but please be mindful of others using the Science Portal.
#. Click the blue "``Launch``" button to start your new ``castor_etc`` notebook session. It
    can take up to a minute to launch the session depending on which computing node you are
    assigned to and the last time the image was launched. Additionally, only 1 session of
    each type is allowed at a time and they automatically shut down after 14 consecutive
    days.
#. The `JupyterLab <https://jupyter.org/>`_ environment is set up to work "out of the box"
    with the ``castor_etc`` Python package. You do not need to install the ``castor_etc``
    package separately; everything is running inside a
    ```conda`` <https://docs.conda.io/en/latest/>`_ environment. In addition to ``LaTeX``, other
    convenience tools like ``git`` are preloaded as well. Feel free to suggest any changes to
    the default configuration.
