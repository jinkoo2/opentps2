FROM python:3.11-slim

ENV POETRY_HOME=/opt/poetry \
    PATH="/opt/poetry/bin:$PATH" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

# Install system deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    libgl1 \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# Install Poetry
RUN python -m venv $POETRY_HOME \
    && $POETRY_HOME/bin/pip install --upgrade pip setuptools wheel \
    && $POETRY_HOME/bin/pip install poetry

WORKDIR /opt/project

# Copy only manifest files first (better Docker caching)
COPY pyproject.toml poetry.lock ./opentps/

# Install dependencies
RUN cd opentps && poetry install --no-root

# Copy full source
COPY . .

# Set PYTHONPATH so tests can import both submodules
ENV PYTHONPATH="/opt/project/opentps/opentps_core:/opt/project/opentps/opentps_gui"

CMD ["/bin/bash"]
