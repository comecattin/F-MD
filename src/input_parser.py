#!/usr/bin/env python3
"""Parsing of the input file."""
import argparse

def parse_config(file_name):
    
    config = {}

    with open(file_name, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            key, value = line.split(':', 1)
            value = value.split('#', 1)[0].strip()
            config[key.strip()] = value.strip()
        
    return config



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parsing of the input file")
    parser.add_argument("file_name", help="Path to the input file")
    args = parser.parse_args()
    config = parse_config(args.file_name)

    dt = config.get('dt', '0.001')
    n_steps = config.get('n_steps', '1000')
    n_atoms = config.get('n_atoms', '30')
    box_length = config.get('box_length', '10.0')
    stride = config.get('stride', '1')

    with open('parsed_config.tmp', 'w') as f:
        f.write(f'{n_atoms}\n')
        f.write(f'{n_steps}\n')
        f.write(f'{dt}\n')
        f.write(f'{box_length}\n')
        f.write(f'{stride}\n')
    