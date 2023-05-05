# Here is some help for using the Dockerfile. For more information refer to
# docker.com
#
# Build the container:
# docker build --tag ubuntu:mesh2hrtf .
#
# Start the container with out an external volume:
# docker run -dit --name mesh2hrtf ubuntu:mesh2hrtf /bin/bash
#
# Start the container with an external volume mounted at /home/data:
# docker run -dit --name mesh2hrtf -v '/local/folder:/home/data' ubuntu:mesh2hrtf /bin/bash

FROM ubuntu:22.04

# All commands are run from this path
WORKDIR /home

# Configure tzdata and timezone if problems appear during `apt-get install cmake`
#ENV TZ=Europe/Berlin
#RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install required dependencies (remove specific versions if desired)
RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get upgrade -y
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y cmake \
                   make \
                   build-essential\
                   libx11-dev \
                   libxrandr-dev \
                   libxinerama-dev \
                   libxcursor-dev \
                   libxi-dev \
                   libgl1-mesa-dev \
                   libglu1-mesa-dev \
                   libeigen3-dev

# copy NumCalc from Mesh2HRTF git repo to docker container
COPY NumCalc /home/NumCalc

# build NumCalc
RUN cd /home/NumCalc/src && make

# add symbolic linc for NumCalc
RUN ln -s /home/NumCalc/bin/NumCalc /usr/local/bin

# copy pmp-library from Mesh2HRTF git repo to docker container
COPY Mesh2Input/Meshes/GradingHybrid/pmp-library /home/pmp-library

# build the pmp-library
RUN cd /home/pmp-library && mkdir build && cd build && cmake .. && make && make install

# add libpmp to the search path
RUN echo /usr/local/lib > /etc/ld.so.conf.d/local.conf
RUN ldconfig

# add symbolic link for hrtf-mesh-grading
RUN ln -s /home/pmp-library/build/hrtf_mesh_grading /usr/local/bin
