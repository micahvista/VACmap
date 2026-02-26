FROM ubuntu:22.04

ARG VACMAP_REF=main
ARG MICROMAMBA_VERSION=1.5.10-0

ENV DEBIAN_FRONTEND=noninteractive
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV PATH=/opt/micromamba/envs/vacmap_env/bin:${PATH}
ENV MPLCONFIGDIR=/tmp/matplotlib

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    bzip2 \
    ca-certificates \
    curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN curl -L "https://conda.anaconda.org/conda-forge/linux-64/micromamba-${MICROMAMBA_VERSION}.tar.bz2" -o /tmp/micromamba.tar.bz2 \
    && tar -xjf /tmp/micromamba.tar.bz2 -C /tmp bin/micromamba \
    && install -m 755 /tmp/bin/micromamba /usr/local/bin/micromamba \
    && rm -rf /tmp/micromamba.tar.bz2 /tmp/bin

RUN git clone --depth 1 --branch "${VACMAP_REF}" --single-branch https://github.com/micahvista/VACmap.git /opt/VACmap \
    && micromamba create -y -n vacmap_env -f /opt/VACmap/VACmap_environment.yml \
    && micromamba run -n vacmap_env python -m pip install --no-cache-dir --no-build-isolation /opt/VACmap \
    && if ! head -n 1 /opt/micromamba/envs/vacmap_env/bin/vacmap | grep -q '^#!'; then \
         sed -i '1i#!/opt/micromamba/envs/vacmap_env/bin/python' /opt/micromamba/envs/vacmap_env/bin/vacmap; \
         chmod +x /opt/micromamba/envs/vacmap_env/bin/vacmap; \
       fi \
    && micromamba clean --all --yes \
    && rm -rf /opt/VACmap \
    && apt-get purge -y --auto-remove build-essential git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p "${MPLCONFIGDIR}" \
    && vacmap --help >/dev/null

ENTRYPOINT ["vacmap"]
CMD ["--help"]
