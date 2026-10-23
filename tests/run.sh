#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "$0")/.." && pwd)"
build_dir="$(mktemp -d)"
trap 'rm -rf "$build_dir"' EXIT
cd "$build_dir"
"${FC:-gfortran}" -g -fcheck=all -ffpe-trap=invalid,zero,overflow \
    "$repo_dir/src/units.f90" "$repo_dir/src/forces.f90" "$repo_dir/src/energies.f90" \
    "$repo_dir/src/integration.f90" "$repo_dir/tests/test_physics.f90" -o test_physics
./test_physics
"${FC:-gfortran}" -g -fcheck=all -ffpe-trap=invalid,zero,overflow \
    "$repo_dir/src/units.f90" "$repo_dir/src/forces.f90" "$repo_dir/src/energies.f90" \
    "$repo_dir/src/integration.f90" "$repo_dir/tests/test_units.f90" -o test_units
./test_units
