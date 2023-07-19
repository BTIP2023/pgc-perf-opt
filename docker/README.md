# Docker Guide

This directory contains the separate .Dockerfiles used to build isolated
environments for each research objective, including extended objectives.

Subdirectories will contain auxiliary files which were copied directly from 
the **rocker-org (n.d.)** repository, and also files that were made specifically
to build the images of this project (e.g. `scripts/{install_pgc_base.sh, install_pgc_aux.R}`).

These containers will resolve all dependency issues brought up in all the `README.md`.

**ATTN:** Please run all commands from the project root.

## Dockerfiles and their Images

### Base Image (`pgc-perf-opt.Dockerfile`)

This image serves as the base image of all research objectives. This attempts to ensure full
compatibility with the project's entire codebase, from `code/` to `presentations/`.

To pull the latest version:

```bash
docker pull yshebron/pgc-perf-opt:latest
```

**Important:** This image *natively supports* the CPU vulnerability performance evaluation and optimization objective.

Examples of resolved issues:

- `knitr` issue that prevented compiling of presentations.
- `save_image` issue with reticulate, plotly, and kaleido that caused the pipeline to hang.
- package namespace issues.
- tidyverse installation in Ubuntu.

Should you encounter these or any other issue, please open an Issue in our main
[repository](https://github.com/PGCInternship2023/pgc-perf-opt).

#### Usage

To build this image:

```bash
docker build -t <yourdockerusername>/pgc-perf-opt:<tag> . -f ./docker/pgc-perf-opt.Dockerfile
```

To start a container using this image:

```bash
docker-compose -f ./docker/compose.yml up
```

You are **encouraged to modify** `./docker/compose.yml` for your purposes.

Other commands:

```bash
# Password:
docker run --rm -ti -e PASSWORD=yourpassword -p 8787:8787 pgc-perf-opt/cpu
# Root:
docker run --rm -ti -e ROOT=true -p 8787:8787 pgc-perf-opt/cpu
# Disable Authentication:
docker run --rm -ti -e DISABLE_AUTH=true -p 127.0.0.1:8787:8787 pgc-perf-opt/cpu
# User ID and Group ID:
docker run --rm -ti -e USERID=1001 -e GROUPID=1001 -p 8787:8787 pgc-perf-opt/cpu
# To not expose:
-p 127.0.0.1:8787:8787
```

### CPU benchmark (`cpu.Dockerfile`)

This section will contain information about the CPU benchmark container.

### GPU benchmark (`gpu.Dockerfile`)

This section will contain information about the GPU benchmark container.

## Notes

Under consideration is the Conda platform for managing dev environments.

## Reference

rocker-org (n.d.). rocker-org/rocker-versioned2: Run current & prior versions of R using docker.
Retrieved July 18, 2023, from [https://github.com/rocker-org/rocker-versioned2](https://github.com/rocker-org/rocker-versioned2)
