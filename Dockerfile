FROM ubuntu:bionic as builder
WORKDIR /
RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential ca-certificates git zlib1g cmake libgsl-dev libxerces-c-dev xsdcxx \
  && rm -rf /var/lib/apt/lists/* \
  && git clone --depth 1 --branch refactor https://github.com/SwissTPH/openmalaria.git \
  && cd openmalaria \
  && ./build.sh -c -b=master -r -a=openmalariaRelease

FROM ubuntu:bionic as openmalaria
WORKDIR /om/
RUN apt-get update && apt-get install -y --no-install-recommends \
  libgsl23 libxerces-c3.2 \
  && rm -rf /var/lib/apt/lists/*
COPY --from=builder /openmalaria/openmalariaRelease/* /om/
ENTRYPOINT ["/om/openMalaria"]  

# docker build -t openmalaria .
# docker run -it --rm -v /home/acavelan/git/openmalaria/nhh:/root/nhh openmalaria -p nhh -s scenarioNonHumanHosts.xml
