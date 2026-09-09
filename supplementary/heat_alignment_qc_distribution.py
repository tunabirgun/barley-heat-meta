import argparse
from pathlib import Path
import runpy
import sys

from bx.bitset_builders import MAX, binned_bitsets_from_file
from qcmodule import BED

parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True)
parser.add_argument('-r', required=True)
args = parser.parse_args()
chromosomes = set()
maximum = 0
with open(args.r) as stream:
    for line in stream:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.split()
        chromosomes.add(fields[0])
        maximum = max(maximum, int(fields[2]))
# Native upstream/downstream windows extend at most 10 kb.
capacity = max(MAX, maximum + 10001)
lengths = dict.fromkeys(chromosomes, capacity)

def bitsets(intervals):
    lines = ('\t'.join(map(str, interval)) + '\n' for interval in intervals)
    return binned_bitsets_from_file(lines, lens=lengths)

BED.binned_bitsets_from_list = bitsets
sys.argv = [str(Path(sys.executable).with_name('read_distribution.py')), '-i', args.i, '-r', args.r]
runpy.run_path(sys.argv[0], run_name='__main__')
