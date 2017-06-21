##############################################
# Dockerfile for 
# https://github.com/AntonioJBT/py27_packaging
##############################################


############
# Base image
############

# FROM python:3-onbuild 
# FROM ubuntu:17.04

FROM jfloff/alpine-python:2.7
# https://github.com/jfloff/alpine-python
# This is a minimal Python 3 image that can start from python or bash

# Or simply run:
# docker run --rm -ti jfloff/alpine-python bash
# docker run --rm -ti jfloff/alpine-python python hello.py

#########
# Contact
#########
MAINTAINER Antonio Berlanga-Taylor <a.berlanga@imperial.ac.uk>


#########################
# Update/install packages
#########################
# Install system dependencies
# For Alpine see:
# https://wiki.alpinelinux.org/wiki/Alpine_Linux_package_management
RUN apk update && apk upgrade \
    && apk add \
    vim \
    sudo

#    wget \
#    bzip2 \
#    unzip \
#    git \ # Already in Alpine Python


#########################
# Install Python packages
#########################

RUN pip install --upgrade pip future

#########################
# Install package to test 
#########################

RUN cd home \
    && git clone https://github.com/AntonioJBT/py27_packaging.git \
    && cd py27_packaging \
    && python setup.py install \
    && cd ..

###############################
# Install external dependencies
###############################


############################
# Default action to start in
############################
# Only one CMD is read (if several only the last one is executed)
#ENTRYPOINT ['/xxx']
#CMD echo "Hello world"
CMD ["/bin/sh"]

# Create a shared folder between docker container and host
#VOLUME ["/shared/data"]
