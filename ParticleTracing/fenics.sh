#!/bin/bash
echo "[!] Begining FEniCS Simulation..."
DATE=`date +%Y-%m-%d`
date
echo "  [+] Configuring DOCKER data volume"
docker volume create --name instant-cache > /dev/null 2>&1
echo "  [+] Mounting persistent DOCKER one-shot container"
docker run --rm -v instant-cache:/home/fenics/.instant -v $(pwd):/home/fenics/shared -w /home/fenics/shared quay.io/fenicsproject/stable "$@"
echo "[!] Computation Complete"
DATE=`date +%Y-%m-%d`
date
