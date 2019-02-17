FROM fpco/stack-build:lts-5.11

ADD . /src
WORKDIR /src
RUN stack install

FROM ubuntu:18.04
COPY --from=0 /root/.local/bin/phybin /usr/bin/phybin
