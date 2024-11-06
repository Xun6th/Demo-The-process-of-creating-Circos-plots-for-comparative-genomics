#!/usr/bin/env python3
import argparse

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile:
        numbers = [int(line.strip()) for line in infile]

    if not numbers:
        return

    groups = []
    current_group = [numbers[0]]

    for i in range(1, len(numbers)):
        if numbers[i] == numbers[i - 1] + 1:
            current_group.append(numbers[i])
        else:
            groups.append(current_group)
            current_group = [numbers[i]]

    groups.append(current_group)  # Append the last group

    with open(output_file, 'w') as outfile:
        for group in groups:
            if len(group) > 1:
                outfile.write(f"{group[0]}\t{group[-1]}\n")
            else:
                outfile.write(f"{group[0]}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a file to merge consecutive numbers.")
    parser.add_argument("input_file", help="Input file path")
    parser.add_argument("output_file", help="Output file path")

    args = parser.parse_args()
    process_file(args.input_file, args.output_file)

