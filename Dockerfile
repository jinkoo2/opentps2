FROM python:3.11-slim

MAINTAINER OpenTPS Developers
LABEL maintainer="OpenTPS Developers"
LABEL description="Docker image for OpenTPS"

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libglib2.0-0 \
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

RUN pip install --upgrade pip setuptools wheel poetry

WORKDIR /opt/project

COPY . .

RUN poetry config virtualenvs.create false \
    && poetry install --no-root \
    && pip install -e .

CMD ["/bin/bash"]