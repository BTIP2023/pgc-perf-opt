# Docker Guide

This directory contains the separate .Dockerfiles used to build isolated
environments for each research objective, including extended objectives.

Subdirectories will contain auxiliary files which were copied directly from 
the **rocker-org (n.d.)** repository, and also files that were made specifically
to build the images of this project (e.g. `scripts/{install_pgc_base.sh, install_pgc_aux.R}`).

These containers will resolve all dependency issues brought up in all the `README.md`.

**ATTN:** Please run all commands from the project root.

## Dockerfiles and their Images

### 1. Base Image (`pgc-perf-opt.Dockerfile`)

This image serves as the base image of all research objectives. This attempts to ensure full
compatibility with the project's entire codebase, from `code/` to `presentations/`. Notable features include utilities like `tree`, `curl` and `nano`, pre-installation of all required packages (`tidyverse` *et al.*) for the current project, and out-of-the-box support for Conda and Reticulate for the default user `rstudio` (accessible through RStudio Server, and may be setup for `root` terminal by running `Rscript reticulate::use_miniconda('r-reticulate')`).

To pull the latest version:

```bash
docker pull yshebron/pgc-perf-opt:latest
```

**Important:** This image *natively supports* the CPU vulnerability performance evaluation and optimization research objective.

Examples of resolved issues:

- `knitr` issue that prevented compiling of presentations.
- `save_image` issue with reticulate, plotly, and kaleido that caused the pipeline to hang.
- package namespace issues.
- tidyverse installation in Ubuntu.

Please see `pgc-perf-opt.Dockerfile` and `scripts/{install_pgc_base.sh, install_pgc_aux.R}` for further information.

Should you encounter these or any other issue, please open an Issue in our main
[repository](https://github.com/PGCInternship2023/pgc-perf-opt).

#### Usage

To build this image:

```bash
docker build -t <yourdockerusername>/pgc-perf-opt:<tag> . -f ./docker/pgc-perf-opt.Dockerfile
```

To start a container using this image:

```bash
docker compose -f ./docker/compose.yml up
```

You are **encouraged to modify** `compose.yml` for your purposes. In particular, you may
set the host-side port as follows: default setting is `127.0.0.1:0:8787`, where `127.0.0.1` makes the RStudio Server
only accessible from your computer (omit it to make it accessible from other computers on the network), and `0` makes
the port a random value (change it to a specific number if you wish).
You may also set the following environment variables:

- `PASSWORD=string` (password for RStudio server, note that username is always `rstudio`)
- `ROOT=bool` (whether to start the container as ROOT user)
- `DISABLE_AUTH=bool` (whether to authenticate RStudio Server users)
- `USERID=int` (default non-root user ID)
- `GROUPID=int` (default non-root group ID)

Note that there should always be a `DOCKER_RUNNING=true` environment variable. DO NOT MODIFY this.

To access the RStudio Server:

1. Open Docker Desktop and locate the `docker > docker-pgc-perf-opt-1` container.
2. Open the container details and click the host port and container port mapping (e.g. 8080:8787)
3. Login to RStudio Server (if DISABLE_AUTH is `FALSE`). Use username `rstudio` and your `compose.yml` password.

Running from the Docker CLI:

```bash
# Set Password:
docker run --rm -ti -e PASSWORD=yourpassword -p 8787:8787 pgc-perf-opt/cpu

# Set Root:
docker run --rm -ti -e ROOT=true -p 8787:8787 pgc-perf-opt/cpu

# Set RStudio Server Authentication:
docker run --rm -ti -e DISABLE_AUTH=true -p 127.0.0.1:8787:8787 pgc-perf-opt/cpu

# Set User ID and Group ID:
docker run --rm -ti -e USERID=1001 -e GROUPID=1001 -p 8787:8787 pgc-perf-opt/cpu

# Select Port to Expose (use `127.0.0.1` for `localhost` only):
-p 127.0.0.1:8787:8787
```

### 2. CPU benchmark (`cpu.Dockerfile`)

This section will contain information about the CPU benchmark container.

### 3. GPU benchmark (`gpu.Dockerfile`)

This section will contain information about the GPU benchmark container.

## Notes

Under consideration is the Conda platform for managing dev environments.

## Reference

rocker-org (n.d.). rocker-org/rocker-versioned2: Run current & prior versions of R using docker.
Retrieved July 18, 2023, from [https://github.com/rocker-org/rocker-versioned2](https://github.com/rocker-org/rocker-versioned2)
