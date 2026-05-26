#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOLS_DIR="${ROOT_DIR}/tools"
mkdir -p "${TOOLS_DIR}"

# Install piawka locally from source if no binary is on PATH.
if ! command -v piawka >/dev/null 2>&1; then
  if [[ ! -d "${TOOLS_DIR}/piawka/.git" ]]; then
    git clone --depth=1 https://github.com/novikovalab/piawka.git "${TOOLS_DIR}/piawka"
  else
    git -C "${TOOLS_DIR}/piawka" pull --ff-only
  fi
fi

# Create a Python conda env for pixy + aggregation stack.
if ! conda env list | awk '{print $1}' | grep -qx 'piawka-paper-py'; then
  conda create -n piawka-paper-py python=3.11 pip -y
fi
conda run -n piawka-paper-py python -m pip install --upgrade pip
conda run -n piawka-paper-py python -m pip install \
  "git+https://github.com/ksamuk/pixy.git@master" pandas matplotlib typing_extensions

echo "Bootstrap complete."
echo "piawka_bin=${TOOLS_DIR}/piawka/piawka"
echo "pixy_bin=/Users/ntikhomirov/mambaforge/envs/piawka-paper-py/bin/pixy"
