FROM ubuntu:bionic as builder
WORKDIR /
RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential ca-certificates git zlib1g cmake libgsl-dev libxerces-c-dev xsdcxx libboost-dev \
  && rm -rf /var/lib/apt/lists/* \
  && git clone --depth 1 --branch schema-43.0 https://github.com/SwissTPH/openmalaria.git \
  && cd openmalaria \
  # && echo "schema-30.3" > version.txt \
  && mkdir build \
  && cd build \
  && cmake -DCMAKE_BUILD_TYPE=Release .. \
  && make -j8 \
  && cd .. \
  && mkdir openmalariaRelease \
  && cp build/openMalaria openmalariaRelease/ \
  && cp build/schema/scenario_current.xsd openmalariaRelease/scenario_43.xsd \
  && cp build/schema/scenario_current.xsd openmalariaRelease/ \
  && cp test/*.csv openmalariaRelease/

FROM ubuntu:bionic as openmalaria
WORKDIR /om/
RUN apt-get update && apt-get install -y --no-install-recommends \
  libgsl23 libxerces-c3.2 \
  && rm -rf /var/lib/apt/lists/*
COPY --from=builder /openmalaria/openmalariaRelease/* ./
ENTRYPOINT ["/om/openMalaria"]  

# docker build -t openmalaria .
# docker run -it --rm -v /home/acavelan/git/openmalaria/nhh:/root/nhh openmalaria -p nhh -s scenarioNonHumanHosts.xml
