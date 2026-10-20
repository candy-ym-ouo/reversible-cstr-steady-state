ARG BASE_IMAGE=python:3.12-slim
FROM ${BASE_IMAGE}

WORKDIR /app
COPY . /app

RUN mkdir -p /logs/verifier \
    && chown -R 10001:10001 /app /logs

USER 10001:10001
ENV PYTHONUNBUFFERED=1

CMD ["python3", "reactor_report.py", "--case", "cases/base.json"]
