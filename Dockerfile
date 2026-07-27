FROM ubuntu:stonking-20260612
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    wget \
    && rm -rf /var/lib/apt/lists/*
ENV PIP_BREAK_SYSTEM_PACKAGES=1
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/miniconda && \
    rm miniconda.sh
ENV PATH="/opt/miniconda/bin:$PATH"    
ENV EMU_DATABASE_DIR="/opt/emu_db"
RUN pip install osfclient && \
    mkdir -p ${EMU_DATABASE_DIR} && \
    cd ${EMU_DATABASE_DIR} && \
    osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar && \
    tar -xvf emu.tar && \
    rm -f ._*
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
RUN conda create -y -n emu_env python=3.7 emu
ENV PATH="/opt/miniconda/envs/emu_env/bin:$PATH"
WORKDIR /app
CMD ["emu"]