#!/bin/bash

# package installer for apt based Dockerfiles
# 1. Install apt packages silently
# 2. Update metadata before install
# 3. Remove unneeded files after install
# From <https://github.com/opencadc/skaha/blob/master/containers/session-containers/astroml-notebook/scripts/apt-install.sh>

[[ $# == 0 ]] && echo "usage: apt-install.sh [packages...]" >&2 && exit 0

if [[ -f $1 ]]; then
    packages="$(cat $1)"
else
    packages="$*"
fi
set -eu

apt-get update --yes -qq
apt-get update --yes --fix-missing
# gcc and g++ compilers are needed to install photutils and celerite python packages.
apt-get -y install gcc g++
DEBIAN_FRONTEND=noninteractive apt-get install --yes ${packages}
apt-get autoremove --purge -y
apt-get clean --yes
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*