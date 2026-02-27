FROM ubuntu:22.04

ARG VACMAP_REF=c486fe256d3df478485166975ec5b718422a242c
ARG MICROMAMBA_VERSION=1.5.10-0
ARG MICROMAMBA_SHA256=e69fe45b21f37fb807d0a2c1ce110f4bd2cabe7aae79bfe92a758125dcbd54f2

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ENV DEBIAN_FRONTEND=noninteractive
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV PATH=/opt/micromamba/envs/vacmap_env/bin:${PATH}
ENV MPLCONFIGDIR=/tmp/matplotlib

# hadolint ignore=DL3008
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    bzip2 \
    ca-certificates \
    curl \
    && curl -fsSL "https://conda.anaconda.org/conda-forge/linux-64/micromamba-${MICROMAMBA_VERSION}.tar.bz2" -o /tmp/micromamba.tar.bz2 \
    && echo "${MICROMAMBA_SHA256}  /tmp/micromamba.tar.bz2" | sha256sum -c - \
    && tar -xjf /tmp/micromamba.tar.bz2 -C /tmp bin/micromamba \
    && install -m 755 /tmp/bin/micromamba /usr/local/bin/micromamba \
    && rm -rf /tmp/micromamba.tar.bz2 /tmp/bin \
    && git clone https://github.com/micahvista/VACmap.git /opt/VACmap \
    && git -C /opt/VACmap checkout --detach "${VACMAP_REF}" \
    && micromamba create -y -n vacmap_env -f /opt/VACmap/VACmap_environment.yml \
    && micromamba run -n vacmap_env python -m pip install --no-cache-dir --no-build-isolation /opt/VACmap \
    && if ! head -n 1 /opt/micromamba/envs/vacmap_env/bin/vacmap | grep -q '^#!'; then \
         sed -i '1i#!/opt/micromamba/envs/vacmap_env/bin/python' /opt/micromamba/envs/vacmap_env/bin/vacmap; \
         chmod +x /opt/micromamba/envs/vacmap_env/bin/vacmap; \
       fi \
    && micromamba clean --all --yes \
    && rm -rf /opt/VACmap \
    && mkdir -p "${MPLCONFIGDIR}" \
    && useradd --create-home --uid 10001 --shell /usr/sbin/nologin vacmap \
    && chown -R vacmap:vacmap "${MPLCONFIGDIR}" /home/vacmap \
    && vacmap --help >/dev/null \
    && apt-get purge -y --auto-remove build-essential git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

USER 10001:10001
ENTRYPOINT ["vacmap"]
CMD ["--help"]

LABEL version=1.0.0
