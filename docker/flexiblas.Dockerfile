
FROM yshebron/pgc-perf-opt:latest

COPY . /home/rstudio/pgc-perf-opt

# Flexiblas installation
# Environment variables to skip debian prompts
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NOWARNINGS yes

# Install dependencies
RUN ./docker/scripts/install_flexiblas.sh

# Download and extract source
# Source will be in /home/rstudio/pgc-perf-opt/flexiblas/flexiblas-3.3.0
RUN mkdir /home/rstudio/pgc-perf-opt/flexiblas
WORKDIR /home/rstudio/pgc-perf-opt/flexiblas
RUN curl https://csc.mpi-magdeburg.mpg.de/mpcsc/software/flexiblas/flexiblas-3.3.0.tar.gz | tar -xz
WORKDIR /home/rstudio/pgc-perf-opt/flexiblas/flexiblas-3.3.0
RUN fakeroot dpkg-buildpackage -us -uc
RUN dpkg -i ../libflexiblas-*.deb
RUN debian/rules clean

WORKDIR /home/rstudio/pgc-perf-opt