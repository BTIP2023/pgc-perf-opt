
FROM yshebron/pgc-perf-opt:latest

# Flexiblas installation
# https://github.com/mpimd-csc/flexiblas
# https://www.mpi-magdeburg.mpg.de/2427308/flexiblas_install

# Environment variables to skip debian prompts
ENV DEBIAN_FRONTEND noninteractive
# ENV DEBCONF_NOWARNINGS yes

ARG NCPUS=16

# Install dependencies
RUN --mount=type=cache,target=/var/cache/apt apt-get update
RUN apt-get upgrade -y
# For Flexiblas
RUN apt-get install -y debhelper cmake gfortran gcc libatlas-base-dev libopenblas-dev build-essential fakeroot \
	libopenblas-openmp-dev libopenblas-pthread-dev libopenblas-serial-dev \
	libopenblas64-openmp-dev libopenblas64-pthread-dev \
	libopenblas64-serial-dev libblis-openmp-dev libblis-pthread-dev \
	libblis-serial-dev libblis64-openmp-dev libblis64-pthread-dev \
	libblis64-serial-dev libblas-dev libblas64-dev liblapack-dev \
	liblapack64-dev libmkl-dev fakeroot

# RUN apt-get install -y --no-install-recommends \
#     debhelper \
#     cmake \
#     gfortran \
#     gcc \
#     libatlas-base-dev \
#     libopenblas-dev \
#     build-essential \
#     fakeroot

# More lib dependencies
# RUN apt-get install -y --no-install-recommends \
#     libopenblas-openmp-dev \
#     libopenblas-serial-dev \
#     libopenblas64-openmp-dev \
#     libopenblas64-pthread-dev \
#     libopenblas64-serial-dev \
#     libblis-openmp-dev \
#     libblis-pthread-dev \
#     libblis-serial-dev \
#     libblis64-openmp-dev \
#     libblis64-pthread-dev \
#     libblis64-serial-dev \
#     libblas64-dev \
#     liblapack64-dev \
#     libmkl-dev
    
RUN apt-get upgrade -y

# Download and extract source
# Source will be in /home/rstudio/pgc-perf-opt/flexiblas/flexiblas-3.3.0
RUN mkdir /home/rstudio/pgc-perf-opt/flexiblas
WORKDIR /home/rstudio/pgc-perf-opt/flexiblas
RUN curl https://csc.mpi-magdeburg.mpg.de/mpcsc/software/flexiblas/flexiblas-3.2.1.tar.gz | tar -xz
WORKDIR /home/rstudio/pgc-perf-opt/flexiblas/flexiblas-3.3.0
RUN fakeroot dpkg-buildpackage -us -uc -j${NCPUS}
RUN dpkg -i ../libflexiblas-*.deb
RUN debian/rules clean

WORKDIR /home/rstudio/pgc-perf-opt