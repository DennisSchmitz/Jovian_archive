#!/usr/bin/env python3
"""Generate sample sheet for Jovian pipeline.

Usage:
  generate_sample_sheet.py <source_dir>

<source_dir> is a directory containing input fastq files with typical
filenames as used in the legacy (non-automated) process. Output will
be a sample sheet in YAML, for example:

  sample1_id:
    R1:
      path/to/sample_R1.fq.gz
    R2:
      path/to/sample_R2.fq.gz
  sample2_id:
  ...
"""

import argparse
import pathlib
import re
import yaml


fq_pattern = re.compile("(.*)_R?(1|2).*\.f(ast)?q(\.gz)?")


def main(args):
    assert args.dir.is_dir(), "Argument must be a directory."
    
    samples = {}
    for file_ in args.dir.iterdir():
        if file_.is_dir():
            continue
        match = fq_pattern.fullmatch(file_.name)
        sample = samples.setdefault(match.group(1), {})
        sample["R{}".format(match.group(2))] = str(file_)

    print(yaml.dump(samples, default_flow_style=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=pathlib.Path, 
                       help="Directory where input files are located")
    main(parser.parse_args())
