#!/bin/bash
# Download the newest VADR model packages into versioned directories
# (vadr-models-<name>-<version>), the layout dfv/vadr2mss_config.py expects.
#
# Unlike download_vadr_models.sh, this does NOT fetch the versions the pipeline
# is pinned to (VADR_*_MODELS_VERSION in the Dockerfile): it fetches newer ones
# side by side, so they can be evaluated before the pins are moved. Existing
# directories are left untouched.
#
# usage: download_vadr_models_latest.sh [-d DIR] [-n] [MODEL ...]
#   -d DIR   destination (default: $VADRMODELDIR, else the current directory)
#   -n       dry run: print what would be downloaded and exit
#   MODEL    any of: sarscov2 mpxv rsv cox1 corona flavi calici flu hmpv
#            (default: all of them)
#
# Versions: VADR 1.7.3's vadr-install.sh pins, except mpxv (1.7-2 is on the FTP
# site but vadr-install.sh still pins 1.4.2-1). cox1 is not in vadr-install.sh;
# 1.2-2 is its newest release. Checked 2026-10-19; as of then the pins equal
# these versions, so this script only matters once newer packages appear.

set -euo pipefail

BASEURL="https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models"
ALLMODELS="sarscov2 mpxv rsv cox1 corona flavi calici flu hmpv"

model_version () {
    case "$1" in
        sarscov2) echo "1.6.3-1" ;;
        mpxv)     echo "1.7-2" ;;
        rsv)      echo "1.7-1" ;;
        cox1)     echo "1.2-2" ;;
        corona)   echo "1.3-3" ;;
        flavi)    echo "1.7-2" ;;
        calici)   echo "1.7-2" ;;
        flu)      echo "1.7-2" ;;
        hmpv)     echo "1.01" ;;  # GitHub release tag v1.01 (greninger-lab)
    esac
}

# FTP directory: three libraries are under <name>viridae
model_urldir () {
    case "$1" in
        corona|flavi|calici) echo "${1}viridae" ;;
        *)                   echo "$1" ;;
    esac
}

# Dockerfile ENV that selects the version
model_env () {
    case "$1" in
        sarscov2) echo "VADR_SCOV2_MODELS_VERSION" ;;
        mpxv)     echo "VADR_MPXV_MODELS_VERSION" ;;
        rsv)      echo "VADR_RSV_MODELS_VERSION" ;;
        cox1)     echo "VADR_COX1_MODELS_VERSION" ;;
        corona)   echo "VADR_CORONA_MODELS_VERSION" ;;
        flavi)    echo "VADR_FLAVI_MODELS_VERSION" ;;
        calici)   echo "VADR_CALCI_MODELS_VERSION" ;;
        flu)      echo "VADR_FLU_MODELS_VERSION" ;;
        hmpv)     echo "VADR_HMPV_MODELS_VERSION" ;;
    esac
}

usage () {
    sed -n '/^# usage:/,/^#            (default/p' "$0" | sed 's/^# \{0,1\}//'
    exit "$1"
}

DEST="${VADRMODELDIR:-$PWD}"
DRYRUN=0
while getopts "d:nh" opt; do
    case "$opt" in
        d) DEST="$OPTARG" ;;
        n) DRYRUN=1 ;;
        h) usage 0 ;;
        *) usage 1 >&2 ;;
    esac
done
shift $((OPTIND - 1))

SELECTED="${*:-$ALLMODELS}"
for m in $SELECTED; do
    case " $ALLMODELS " in
        *" $m "*) ;;
        *) echo "ERROR: unknown model '$m' (choose from: $ALLMODELS)" >&2; exit 1 ;;
    esac
done

if [ ! -d "$DEST" ]; then
    echo "ERROR: destination directory '$DEST' does not exist" >&2
    exit 1
fi
echo "Destination: $DEST"

TMPDIR_DL=$(mktemp -d "$DEST/.vadr-models-dl.XXXXXX")
trap 'rm -rf "$TMPDIR_DL"' EXIT

for m in $SELECTED; do
    v=$(model_version "$m")
    name="vadr-models-$m-$v"
    url="$BASEURL/$(model_urldir "$m")/$v/$name.tar.gz"
    # hMPV comes from GitHub; its tag archive unpacks to $name/ with the model
    # files one level down in hmpv/
    [ "$m" = hmpv ] && url="https://github.com/greninger-lab/vadr-models-hmpv/archive/refs/tags/v$v.tar.gz"
    if [ -d "$DEST/$name" ]; then
        echo "[skip] $name (already present)"
        continue
    fi
    if [ "$DRYRUN" = 1 ]; then
        echo "[would download] $url"
        continue
    fi
    echo "[download] $url"
    curl -fL -# --retry 3 -o "$TMPDIR_DL/$name.tar.gz" "$url"
    # extract into the temp dir first so an interrupted run never leaves a
    # half-extracted vadr-models-* directory behind
    tar -xzf "$TMPDIR_DL/$name.tar.gz" -C "$TMPDIR_DL"
    if [ "$m" = hmpv ]; then
        mv "$TMPDIR_DL/$name" "$TMPDIR_DL/$name.src" && mv "$TMPDIR_DL/$name.src/hmpv" "$TMPDIR_DL/$name"
    fi
    if [ ! -d "$TMPDIR_DL/$name" ]; then
        echo "ERROR: $name.tar.gz did not contain the expected top directory '$name/'" >&2
        exit 1
    fi
    mv "$TMPDIR_DL/$name" "$DEST/$name"
    rm -f "$TMPDIR_DL/$name.tar.gz"
done

[ "$DRYRUN" = 1 ] && exit 0

echo ""
echo "Done. To make the pipeline use these versions, set (Dockerfile ENV and"
echo "Singularity.def %environment):"
for m in $SELECTED; do
    echo "  $(model_env "$m")=\"$(model_version "$m")\""
done
