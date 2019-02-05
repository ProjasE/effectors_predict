FROM rocker/rstudio:3.4.0
LABEL Name=rstudio_ss
EXPOSE 8787

ENV USER=rstudio
ENV PASSWORD=temp123

RUN apt-get update -y && \
    apt-get install -y \
    vim \
    git \
    build-essential \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    libssh2-1-dev \
    zlib1g-dev \
    qawk

ENV ROOT=TRUE

COPY Effectors /home/Effectors
WORKDIR /home/Effectors/

RUN bash secret_santa_installer_fixed.sh
ENV PATH=\
"/home/Effectors/SS_dep/WoLFPSort/bin:/home/Effectors/SS_dep/targetp-1.1:/home/Effectors/SS_dep/signalp-4.1:/home/Effectors/SS_dep/signalp-2.0:/home/Effectors/SS_dep/tmhmm-2.0c/bin:/home/Effectors/SS_dep/signalp-3.0:$PATH"

RUN Rscript install_SecretSanta.R