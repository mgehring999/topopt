FROM ubuntu:latest AS base

RUN apt-get update && apt-get install -y locales wget unzip && rm -rf /var/lib/apt/lists/* \
    && localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8
ENV LANG en_US.utf8

RUN cd /root/ && wget https://ampl.com/dl/open/ipopt/ipopt-linux64.zip && unzip ipopt-linux64.zip
RUN mv /root/ipopt /bin
RUN . /root/.bashrc && ipopt

