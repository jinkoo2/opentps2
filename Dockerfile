FROM python:3.11-slim

ENV POETRY_HOME=/opt/poetry \
    PATH="/opt/poetry/bin:$PATH" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

# Install system deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libglib2.0-0 \
    libgl1-mesa-glx \
    libxcb-icccm4 \
    libxcb-image0 \
    libxcb-keysyms1 \
    libxcb-randr0 \
    libxcb-render-util0 \
    libxcb-xinerama0 \
    libxcb-xkb1 \
    libxkbcommon-x11-0 \
    libxcb-xfixes0 \
    libxcb-shape0 \
    libxcb-shm0 \
    libx11-xcb1 \
    libxcb-cursor0 \
    libxcb-util1 \
    libxcb-glx0 \
    libxcb-dri2-0 \
    libxcb-dri3-0 \
    libxcb-present0 \
    git \
    xvfb \
    --no-install-recommends && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


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
