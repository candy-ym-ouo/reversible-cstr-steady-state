# Validation record

Run date: 2026-10-20, Asia/Shanghai.

## Container build

The image builds successfully with Python 3.12 and starts as UID/GID `10001`.

```text
docker build \
  --build-arg BASE_IMAGE=docker.m.daocloud.io/library/python:3.12-slim \
  -t reversible-cstr-steady-state:local environment
```

The base-image override was needed only because this host could not reach Docker Hub. The Dockerfile default remains `python:3.12-slim`.

Inspection of `/app` found only the Dockerfile, README, case file, and the imperfect program. It did not contain `tests/` or `solution/`.

## Harbor runs

```text
harbor run -p . -a nop -e docker -y -q -o work/harbor-jobs
harbor run -p . -a oracle -e docker -y -q -o work/harbor-jobs
```

| Run | Trials | Exceptions | Mean reward | Record |
|---|---:|---:|---:|---|
| `nop` | 1 | 0 | 0.000 | `validation/harbor-nop-result.json` |
| `oracle` | 1 | 0 | 1.000 | `validation/harbor-oracle-result.json` |

The direct container runs produced the same rewards:

| Candidate | Reward | Record |
|---|---:|---|
| Fresh environment | 0.000 | `validation/container-nop-report.json` |
| Reference implementation | 1.000 | `validation/container-oracle-report.json` |
| Missing/crashing program | 0.000 | `validation/container-crash-report.json` |

## Scoring parts

- `0.30`: temperature-dependent reversible kinetics and fixed-temperature mass balance.
- `0.30`: coupled sensible heat and utility heat balance.
- `0.20`: discovery of all three steady states in the multiplicity case.
- `0.20`: local stability labels from the concentration/temperature Jacobian.

## Target-model trial

An isolated copy containing only `instruction.md` and `environment/` was given to `deepseek-v4-flash` through the local Codex CLI. The model did not receive `tests/` or `solution/`.

| Candidate | Reward | Result | Record |
|---|---:|---|---|
| `deepseek-v4-flash` | 0.300 | found 2 of 3 steady states in the close-root case; the third branch was missed | `validation/model-deepseek-v4-flash-report.json` |

The model implementation did not pass the all-roots or stability checks, so it does not receive full credit. The missing `net_rate_kmol_m3_s` field was added to the visible output contract after that trial; the root-count failure is independent of that field.
