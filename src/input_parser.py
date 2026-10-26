#!/usr/bin/env python3
"""Parse and validate the simulation input before allocating Fortran arrays."""
import argparse
import math


def parse_config(file_name):
    config = {}
    with open(file_name) as source:
        for number, line in enumerate(source, 1):
            line = line.split('#', 1)[0].strip()
            if not line:
                continue
            if ':' not in line:
                raise ValueError(f'Line {number}: expected key: value')
            key, value = line.split(':', 1)
            config[key.strip()] = value.strip()
    return config


def config_values(config):
    n_atoms = int(config.get('n_atoms', '30'))
    n_steps = int(config.get('n_steps', '1000'))
    dt = float(config.get('dt', '0.0001'))  # ps (0.1 fs)
    box_length = float(config.get('box_length', '10.0'))  # Angstrom
    stride = int(config.get('stride', '1'))
    if n_atoms <= 0 or n_atoms % 3:
        raise ValueError('n_atoms must be a positive multiple of 3')
    if n_steps <= 0 or stride <= 0:
        raise ValueError('n_steps and stride must be positive integers')
    if not math.isfinite(dt) or dt <= 0:
        raise ValueError('dt must be finite and positive (ps)')
    if not math.isfinite(box_length) or box_length <= 0:
        raise ValueError('box_length must be finite and positive (Angstrom)')
    return n_atoms, n_steps, dt, box_length, stride


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('file_name', help='Path to the input file')
    args = parser.parse_args()
    try:
        values = config_values(parse_config(args.file_name))
    except (OSError, ValueError) as error:
        parser.error(str(error))
    with open('parsed_config.tmp', 'w') as destination:
        destination.write('\n'.join(map(str, values)) + '\n')


if __name__ == '__main__':
    main()
