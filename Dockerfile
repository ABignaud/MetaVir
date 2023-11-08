# syntax=docker/dockerfile:1

# Start with metabat2 docker
FROM metabat/metabat:latest

LABEL Name=metavir Version=1.1.0

WORKDIR /app
COPY ./ /app

# install python and python-pip
RUN apt-get update && apt-get install -y python3-pip

# Install python dependencies
RUN pip3 install -Ur requirements.txt

# Build checkV database
RUN checkv download_database ./

# Install metator
RUN pip3 install .

ENTRYPOINT ["metavir"]
