FROM python:3.10-slim-bookworm

LABEL org.opencontainers.image.title="Molecular Transport Audit Python Core"
LABEL org.opencontainers.image.description="Python scientific environment for Molecular Transport Audit"
LABEL org.opencontainers.image.source="https://github.com/Helen2813/molecular-transport-audit"

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PIP_NO_CACHE_DIR=1

WORKDIR /opt/mta

COPY pyproject.toml README.md ./
COPY core ./core

RUN python -m pip install .

RUN python -c "import numpy, pandas, yaml, transport_audit; print('transport_audit', transport_audit.__version__)"

RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

CMD ["python3", "--version"]
