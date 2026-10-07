# Vibrion Sentinel - Multi-stage Docker Build
# Optimized for genomic surveillance pipelines

# Stage 1: Build environment with all dependencies
FROM mambaorg/micromamba:1.5-jammy AS builder

USER root
WORKDIR /app

# Copy environment specification
COPY environment.yml /app/environment.yml

# Create conda environment
RUN micromamba create -y -n vibrion -f /app/environment.yml && \
    micromamba clean --all --yes

# Stage 1b: Bake databases with versions recorded.
# The per-sample manifest (liability guardrails §5) requires tool+DB
# versions, and "latest" is not a version — so both DBs are fetched here,
# at image build time, and their versions are written to /app/db_versions.txt.
ENV PATH="/opt/conda/envs/vibrion/bin:$PATH"
RUN mkdir -p /app/db/amrfinderplus && \
    amrfinder_update --force_update --database /app/db/amrfinderplus && \
    { echo "# tool+database versions baked at image build time"; \
      echo "built: $(date -u +%Y-%m-%dT%H:%M:%SZ)"; \
      amrfinder --version; \
      medaka --version; \
      echo "medaka_models_prefetched:"; } > /app/db_versions.txt && \
    for model in r1041_e82_400bps_fast_g615 r1041_e82_400bps_hac_g615 r1041_e82_400bps_sup_g615; do \
      medaka tools download_models --models "$model" >> /app/db_versions.txt 2>&1 || \
        { echo "FATAL: medaka model prefetch failed for $model -- confirm the"; \
          echo "exact model names against the pinned medaka (medaka tools list)"; \
          echo "and the prefetch invocation, then rebuild."; exit 1; }; \
      echo "  - $model" >> /app/db_versions.txt; \
    done && \
    cat /app/db_versions.txt
# NOTE: the `medaka tools download_models` invocation above must be verified
# against the pinned medaka at build time; the loop fails the build loudly
# (not silently) if the command or a model name is wrong.

# Stage 2: Runtime image
FROM mambaorg/micromamba:1.5-jammy AS runtime

USER root
WORKDIR /app

# Copy conda environment from builder
COPY --from=builder /opt/conda/envs/vibrion /opt/conda/envs/vibrion

# Copy baked databases + version record from the DB stage
COPY --from=builder /app/db /app/db
COPY --from=builder /app/db_versions.txt /app/db_versions.txt
ENV AMRFINDERPLUS_DB=/app/db/amrfinderplus

# Set environment variables
ENV PATH="/opt/conda/envs/vibrion/bin:$PATH"
ENV CONDA_DEFAULT_ENV=vibrion
ENV VIBRION_HOME=/app
ENV VIBRION_DATA=/data
ENV VIBRION_OUTPUT=/output

# Create directories
RUN mkdir -p /data /output /app/logs

# Copy application code
COPY workflow/ /app/workflow/
COPY backend/ /app/backend/
COPY scripts/ /app/scripts/
COPY data/references/ /app/data/references/
COPY data/metadata/ /app/data/metadata/

# Copy configuration
COPY workflow/config/config.yaml /app/workflow/config/config.yaml

# Expose ports
EXPOSE 8000 8888

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import snakemake; print('OK')" || exit 1

# Default command: run fast triage
ENTRYPOINT ["python", "/app/scripts/fast_triage.py"]
CMD ["--help"]

# Stage 3: ingest service (Phase 1). Extends the pipeline runtime image so the
# worker has the full bioinformatics toolchain (snakemake, flye, medaka, ...)
# and the API shares the exact same codebase. Default target stays `runtime`;
# build this stage explicitly for the service.
FROM runtime AS service

USER root

# WeasyPrint (Phase 3 PDF reports) renders via Pango/Cairo — system
# libraries the conda env does not provide. Installed here, before the pip
# deps, so the image can render /jobs/{id}/report.pdf.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libpango-1.0-0 \
        libpangoft2-1.0-0 \
        libharfbuzz-subset0 \
        libjpeg-dev \
        libopenjp2-7-dev \
        libffi-dev \
        fonts-dejavu-core \
    && rm -rf /var/lib/apt/lists/*

# Service Python deps (FastAPI app + worker). The vibrion conda env's pip is
# used so there is one Python per image.
COPY service/requirements.txt /app/service/requirements.txt
RUN /opt/conda/envs/vibrion/bin/pip install --no-cache-dir \
        -r /app/service/requirements.txt

# Service code (app, worker, key CLI).
COPY service/ /app/service/

# The worker writes pipeline configs referencing this Snakefile; already in
# the image via the runtime stage (/app/workflow). REPO_DIR defaults to /app.
ENV REPO_DIR=/app
ENV WORKER_WORKDIR=/tmp/sentinel-work

EXPOSE 8000
