FROM ubuntu:24.04
RUN apt -qq update && \
    apt install -y \
    git \
    cmake \
    python3 \
    python3-pip \
    libmpfr-dev \
    libgmp-dev \
    libboost-dev \
    libboost-filesystem-dev \
    libboost-timer-dev \
    libglu1-mesa-dev \
    libxrender1 \
    libxcursor1 \
    libxft2 \
    libxinerama1
CMD ["/bin/bash"]
ENV MAIN_DIR=/root/moving_heat_source
WORKDIR ${MAIN_DIR}
RUN git clone --recurse-submodules \
    https://github.com/ordinary-slim/moving_heat_source . && \
    python3 -m pip install --break-system-packages \
    -r python_requirements.txt && \
    python3 setup.py develop
