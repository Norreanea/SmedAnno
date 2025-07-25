#!/usr/bin/env bash
set -euo pipefail

# Clean any old/broken entries
sudo rm -f /etc/apt/sources.list.d/{nvidia-container-toolkit,libnvidia-container,nvidia-docker,nvidia-container-runtime}.list 2>/dev/null

# Key + repo (generic deb repo)
ARCH=$(dpkg --print-architecture)
sudo mkdir -p /etc/apt/keyrings
sudo wget -qO /etc/apt/keyrings/nvidia-container-toolkit.asc \
  https://nvidia.github.io/libnvidia-container/gpgkey

cat <<EOF | sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
deb [arch=${ARCH} signed-by=/etc/apt/keyrings/nvidia-container-toolkit.asc] \
https://nvidia.github.io/libnvidia-container/stable/deb/${ARCH} /
# deb [arch=${ARCH} signed-by=/etc/apt/keyrings/nvidia-container-toolkit.asc] \
# https://nvidia.github.io/libnvidia-container/experimental/deb/${ARCH} /
EOF

# Install
sudo apt-get update
sudo apt-get install -y nvidia-container-toolkit

# Wire Docker to NVIDIA runtime
sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker

# Sanity test (inside a CUDA container)
docker run --rm --gpus all ubuntu nvidia-smi
