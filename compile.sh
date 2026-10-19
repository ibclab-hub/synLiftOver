#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
"${CXX:-g++}" -std=c++11 -O2 src/step2/synLiftOver_new.cpp -o src/step2/synLiftOver.out
perl -c src/step2/synLiftOver.pl
