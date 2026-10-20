#!/usr/bin/env bash
set -u

reward_file="/logs/verifier/reward.txt"
mkdir -p /logs/verifier
rm -f /logs/verifier/reward.json
printf '0\n' > "$reward_file"

on_signal() {
  printf '0\n' > "$reward_file"
  exit 0
}
trap on_signal INT TERM HUP

python3 /tests/test_outputs.py
status=$?

if [ "$status" -ne 0 ]; then
  printf '0\n' > "$reward_file"
fi

exit 0
