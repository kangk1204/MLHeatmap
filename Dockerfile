FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

RUN apt-get update && \
    apt-get install -y --no-install-recommends libgomp1 && \
    rm -rf /var/lib/apt/lists/*

COPY pyproject.toml README.md LICENSE requirements-lock.txt ./
COPY src ./src

# Install the pinned, tested versions first so the image reproduces the
# manuscript results exactly, then install the package (deps already satisfied).
RUN python -m pip install --upgrade pip && \
    python -m pip install -r requirements-lock.txt && \
    python -m pip install ".[full]"

EXPOSE 8765

HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=5 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8765/api/v1/capabilities', timeout=3).read()"

CMD ["mlheatmap", "--host", "0.0.0.0", "--no-browser"]
