# ==============================================================================
# SMART — Somatic Mutation Annotation and Reporting Tool
# ==============================================================================
# Multi-stage build:
#   Stage 1: VEP installation with plugins
#   Stage 2: Final image with all tools combined
# ==============================================================================

ARG VEP_VERSION=114.0
FROM ensemblorg/ensembl-vep:release_${VEP_VERSION} AS vep-base

# ==============================================================================
# Final stage: combine everything into one image
# ==============================================================================
FROM ubuntu:22.04

ARG VEP_VERSION=114.0
LABEL maintainer="SMART Pipeline Team"
LABEL description="SMART — Somatic Mutation Annotation and Reporting Tool"
LABEL smart.vep.version="${VEP_VERSION}"

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# ---------- System dependencies ----------
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    wget \
    curl \
    gzip \
    bzip2 \
    tabix \
    bgzip \
    bcftools \
    samtools \
    openjdk-17-jre-headless \
    python3 \
    python3-pip \
    perl \
    cpanminus \
    libdbi-perl \
    libdbd-mysql-perl \
    libwww-perl \
    libjson-perl \
    libmodule-build-perl \
    libarchive-zip-perl \
    libbio-perl-perl \
    libtry-tiny-perl \
    libset-intervaltree-perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libhts-dev \
    git \
    unzip \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# ---------- GATK ----------
ARG GATK_VERSION=4.6.0.0
ENV GATK_VERSION=${GATK_VERSION}
LABEL smart.gatk.version="${GATK_VERSION}"
RUN mkdir -p /opt/gatk && \
    wget -q "https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip" \
         -O /tmp/gatk.zip && \
    unzip -q /tmp/gatk.zip -d /opt/ && \
    mv /opt/gatk-${GATK_VERSION}/* /opt/gatk/ && \
    rm -rf /tmp/gatk.zip /opt/gatk-${GATK_VERSION} && \
    chmod +x /opt/gatk/gatk
ENV PATH="/opt/gatk:${PATH}"

# ---------- VEP (copy from vep-base stage) ----------
COPY --from=vep-base /opt/vep /opt/vep
ENV PATH="/opt/vep/src/ensembl-vep:${PATH}"
ENV PERL5LIB="/opt/vep/src/ensembl-vep:/opt/vep/src/ensembl-vep/modules:${PERL5LIB:-}"

# ---------- Python dependencies ----------
COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt

# ---------- OncoKB Annotator ----------
RUN git clone https://github.com/oncokb/oncokb-annotator.git /opt/oncokb-annotator && \
    pip3 install --no-cache-dir -r /opt/oncokb-annotator/requirements/common.txt || true

# ---------- Pipeline scripts ----------
COPY scripts/ /opt/smart/scripts/
RUN chmod +x /opt/smart/scripts/*.sh /opt/smart/scripts/*.py 2>/dev/null || true

# ---------- Working directory ----------
WORKDIR /data

# ---------- Entrypoint ----------
COPY entrypoint.sh /opt/smart/entrypoint.sh
RUN chmod +x /opt/smart/entrypoint.sh

ENTRYPOINT ["/opt/smart/entrypoint.sh"]
CMD ["--help"]
