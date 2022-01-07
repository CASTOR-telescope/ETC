#
# Dockerfile for the CASTOR exposure time calculator. Template taken from
# <https://github.com/arfon/astropy-jupyterlab-docker>.
#
FROM jupyter/scipy-notebook:hub-2.0.1
LABEL maintainer="Isaac Cheng <isaac.cheng.ca@gmail.com>"
ENV NOTEBOOK_DIR="/arc/projects/CASTOR/ETC"
WORKDIR ${NOTEBOOK_DIR}
#
# Install astropy for science, ipyevents for interactive plots, and PyQt5 for GUI
#
RUN pip install astropy ipyevents PyQt5
RUN jupyter nbextension enable --py --sys-prefix ipyevents
#
# Custom JupyterLab configuration
#
# (or is it ${NB_USER}?)
USER ${NB_UID}
COPY misc/overrides.json /opt/conda/share/jupyter/lab/settings/
COPY misc/matplotlibrc ${HOME}/.config/matplotlib/matplotlibrc
# #
# # VS Code configuration
# # (see <https://code.visualstudio.com/remote/advancedcontainers/avoid-extension-reinstalls)
# #
# RUN mkdir -p ${HOME}/.vscode-server/extensions \
#     && chown -R ${NB_USER} ${HOME}/.vscode-server
