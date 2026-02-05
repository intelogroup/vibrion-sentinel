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

# Stage 2: Runtime image
FROM mambaorg/micromamba:1.5-jammy AS runtime

USER root
WORKDIR /app

# Copy conda environment from builder
COPY --from=builder /opt/conda/envs/vibrion /opt/conda/envs/vibrion

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
