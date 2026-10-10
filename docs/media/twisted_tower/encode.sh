#!/bin/bash
# usage: encode.sh <frames dir> <out.mp4>
set -euo pipefail
FF=/tmp/claude-0/-home-user-SparLab/d56879da-faaf-5259-b306-3c582599460d/scratchpad/linkedin/venv/lib/python3.11/site-packages/imageio_ffmpeg/binaries/ffmpeg-linux-x86_64-v7.0.2
"$FF" -hide_banner -loglevel error -y -framerate 30 -i "$1/%04d.png" \
  -c:v libx264 -preset slow -crf 18 -profile:v high -pix_fmt yuv420p \
  -movflags +faststart -r 30 "$2"
