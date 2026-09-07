#!/usr/bin/env bash
###################################################################################
# start_server_noninteractive.sh
#
# Non-interactive variant of start_server.sh, for use from SLURM batch jobs
# (run_selenzymerf.sh) which cannot answer the `read` prompts in start_server.sh.
# Behavior is controlled by environment variables instead:
#
#   REMAKE_SELENZYME2  1|0 (default 1)  Rebuild the selenzyme2/ mount directory
#                                        (code + data) from scratch.
#   DELETE_UPLOADS     1|0 (default 0)  Only consulted when REMAKE_SELENZYME2=0;
#                                        clears selenzyme2/selenzyPro/uploads/*.
#   DATA_SOURCE_DIR    path (default "" = unset)
#                                        If set, its contents are copied into
#                                        selenzyme2/selenzyPro/data/ instead of
#                                        unzipping compressed_data/{data_2023,seqs}.zip.
#                                        Point this at a filtered DB directory
#                                        produced by generate_selenzyme_db.py
#                                        (plus seqs.fasta alongside it -- pass
#                                        --seqs_fasta to that script) to run the
#                                        server against a seed's leakage-free
#                                        reference DB. Only consulted when
#                                        REMAKE_SELENZYME2=1.
#
# Runs the already-built Apptainer image (build/selenzyme.sif, built from selenzyme.def) instead
# of Docker -- this cluster's SLURM nodes have no sudo/Docker, but do have the Apptainer module
# (confirmed: `module avail apptainer` lists 1.1.8/1.2.4/1.3.5/1.4.5). Apptainer shares the host's
# network namespace by default (no Docker-style -p remapping), so the Flask server started inside
# the container by selenzyme.def's %runscript (flaskform.py, hardcoded `app.run(port=5001)`) is
# reached directly at localhost:5001 on the node running this script, NOT the port 32784 the old
# Docker -p 32784:5001 mapping used -- callers must point --server_url at :5001.
#
# `sudo rm -rf`/`docker` calls from the Docker version are gone: rebuilding selenzyme2/ is a plain
# `rm -rf` (no sudo needed, it's a normal directory under $DIR, not a Docker-managed volume), and
# there is nothing to `docker build`/`docker rmi` since the image is prebuilt.
#
# NOTE: written to be logically correct and ready to run, but not itself executed (no apptainer/
# sbatch invoked while authoring it) -- the first real SLURM submission is the actual test of the
# networking assumption above; watch its log for "server container started" / a reachable port.
###################################################################################
set -euo pipefail

REMAKE_SELENZYME2="${REMAKE_SELENZYME2:-1}"
DELETE_UPLOADS="${DELETE_UPLOADS:-0}"
DATA_SOURCE_DIR="${DATA_SOURCE_DIR:-}"
SELENZYME2_DIR="${SELENZYME2_DIR:-selenzyme2}"  # separate mount dir per concurrent job (e.g. per seed) so parallel SLURM jobs sharing this filesystem path do not collide
SIF_PATH="${SIF_PATH:-build/selenzyme.sif}"

DIR=$(cd "$(dirname "$0")"; pwd)
cd "$DIR"

#### Building the directory to mount

if [ "$REMAKE_SELENZYME2" == "1" ]; then
    rm -rf "$SELENZYME2_DIR"
    mkdir "$SELENZYME2_DIR"

    echo "     COPYING selenzyme2 code"
    cp -rp gitcode2023/* "$SELENZYME2_DIR"
    mkdir -p "$SELENZYME2_DIR/selenzyPro"

    echo "     COPYING selenzyme2 data"
    mkdir -p "$SELENZYME2_DIR/selenzyPro/data"

    if [ -n "$DATA_SOURCE_DIR" ]; then
        if [ ! -d "$DATA_SOURCE_DIR" ]; then
            echo "ERROR: DATA_SOURCE_DIR=$DATA_SOURCE_DIR does not exist" >&2
            exit 1
        fi
        echo "     copying pre-built DB from $DATA_SOURCE_DIR"
        # -L dereferences symlinks so the mounted dir is self-contained even if
        # generate_selenzyme_db.py was run with --link_mode symlink (the default).
        cp -rpL "$DATA_SOURCE_DIR"/. "$SELENZYME2_DIR/selenzyPro/data/"
        if [ ! -f "$SELENZYME2_DIR/selenzyPro/data/seqs.fasta" ]; then
            echo "WARNING: no seqs.fasta found in DATA_SOURCE_DIR ($DATA_SOURCE_DIR)." \
                 "Pass --seqs_fasta to generate_selenzyme_db.py so the filtered DB" \
                 "directory is self-contained; falling back to unzipping seqs.zip." >&2
            unzip -q compressed_data/seqs.zip -d "$SELENZYME2_DIR/selenzyPro/data"
        fi
    else
        echo "unzipping data_2023.zip !!!"
        unzip -q compressed_data/data_2023.zip -d "$SELENZYME2_DIR/selenzyPro/data/"
        mv "$SELENZYME2_DIR/selenzyPro/data/data_2023"/* "$SELENZYME2_DIR/selenzyPro/data/"
        rmdir "$SELENZYME2_DIR/selenzyPro/data/data_2023"

        echo "unzipping seqs.zip"
        unzip -q compressed_data/seqs.zip -d "$SELENZYME2_DIR/selenzyPro/data"
    fi

    echo "     MAKING additional folders"
    mkdir -p "$SELENZYME2_DIR/selenzyPro/log"
    mkdir -p "$SELENZYME2_DIR/selenzyPro/uploads"
else
    if [ "$DELETE_UPLOADS" == "1" ]; then
        rm -rf "$SELENZYME2_DIR/selenzyPro/uploads"/*
    fi
fi


### The Apptainer bit

if [ ! -f "$SIF_PATH" ]; then
    echo "ERROR: $SIF_PATH not found. Build it first with:" >&2
    echo "  module load Apptainer && apptainer build $SIF_PATH selenzyme.def" >&2
    echo "(may require --fakeroot depending on cluster policy)" >&2
    exit 1
fi

module load Apptainer 2>/dev/null || true

# No image build/remove step needed (unlike the Docker version) -- $SIF_PATH is prebuilt.
# --bind mounts selenzyme2/ at the same /selenzyme2 path selenzyme.def's %runscript expects
# (matching the Dockerfile's -v "$DIR/selenzyme2:/selenzyme2"). No -p port mapping: Apptainer
# shares the host network namespace by default, so the %runscript's flaskform.py is reachable
# directly at localhost:5001 (see header note).
echo "     STARTING selenzyme2023 via Apptainer ($SIF_PATH)"
apptainer run --bind "$DIR/$SELENZYME2_DIR:/selenzyme2" "$SIF_PATH"

echo ""
echo "server exited"
