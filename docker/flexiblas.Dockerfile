
FROM yshebron/pgc-perf-opt:latest

COPY . /home/rstudio/pgc-perf-opt
WORKDIR /home/rstudio/pgc-perf-opt

# Flexiblas installation
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NOWARNINGS yes
RUN ./docker/scripts/install_flexiblas.sh
RUN curl https://csc.mpi-magdeburg.mpg.de/mpcsc/software/flexiblas/flexiblas-3.3.0.tar.gz | tar -xz
WORKDIR /home/rstudio/pgc-perf-opt/flexiblas-3.3.0
RUN fakeroot dpkg-buildpackage -us -uc
RUN dpkg -i ../libflexiblas-*.deb
RUN debian/rules cleana
WORKDIR /home/rstudio/pgc-perf-opt