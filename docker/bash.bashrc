# From <https://github.com/opencadc/skaha/blob/master/containers/session-containers/astroml-notebook/config/bash.bashrc>

export PS1="${debian_chroot:+($debian_chroot)}\u \W \$ "
. /opt/conda/etc/profile.d/mamba.sh

source /opt/conda/etc/profile.d/conda.sh
conda activate base
eval "$(command conda shell.bash hook 2> /dev/null)"
