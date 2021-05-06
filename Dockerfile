FROM ubuntu

RUN apt-get update \
    && apt-get install -y git python3 python3-pip gcc make wget vim\
    && pip3 install numpy scipy matplotlib seaborn

WORKDIR /root
RUN wget https://boostorg.jfrog.io/artifactory/main/release/1.67.0/source/boost_1_67_0.tar.bz2 \
    && tar -jxvf boost_1_67_0.tar.bz2

RUN git clone https://github.com/qkirikigaku/SigTracer
WORKDIR /root/SigTracer
RUN  make compile_docker

RUN echo 'alias python="python3"' >> ~/.bashrc
