#!/usr/bin/env python3
"""
add_parent_dir_prefix.py

Read a text file containing filenames (one per line) or md5sum lines
and prepend the parent/sample directory to each filename.

Examples:
  01001CM_122022_S1_L005_I1_001.fastq.gz -> 01001CM_122022/01001CM_122022_S1_L005_I1_001.fastq.gz
  <md5>  01001CM_122022_S1_L005_I1_001.fastq.gz -> <md5>  01001CM_122022/01001CM_122022_S1_L005_I1_001.fastq.gz

Usage:
  python3 scripts/add_parent_dir_prefix.py input.txt -o output.txt

If -o / --output is omitted, output is written to stdout.
"""
from __future__ import annotations

import argparse
import re
import sys
from typing import Optional


MD5_RE = re.compile(r"^[0-9a-fA-F]{32}$")
SAMPLE_RE = re.compile(r"^(.*?)_S\d+_")


def sample_from_filename(fn: str) -> str:
    """Extract sample prefix from filename.

    Strategy:
    - If filename matches r'^(.*?)_S\d+_' take the captured group (works for 10x names)
    - Else, fallback to joining the first two underscore-separated parts if available
      (e.g. '01001CM_122022_...' -> '01001CM_122022')
    - If only one part, return that part.
    """
    m = SAMPLE_RE.match(fn)
    if m:
        return m.group(1)
    parts = fn.split("_")
    if len(parts) >= 2:
        return parts[0] + "_" + parts[1]
    return parts[0]


def transform_line(line: str) -> Optional[str]:
    """Transform a single line. Returns new line or None for blank lines.

    Preserves an md5 followed by whitespace if present.
    """
    line = line.rstrip("\n")
    if not line.strip():
        return None

    # Split into tokens; filename is usually the last token
    tokens = line.split()
    if not tokens:
        return None

    if MD5_RE.match(tokens[0]) and len(tokens) >= 2:
        # md5 present, filename should be last token
        md5 = tokens[0]
        filename = tokens[-1]
        # If filename already contains a '/', assume already prefixed
        if "/" in filename:
            newfn = filename
        else:
            sample = sample_from_filename(filename)
            newfn = f"{sample}/{filename}"
        # Preserve two spaces between md5 and filename (common md5sum format)
        return f"{md5}  {newfn}"

    # otherwise assume whole line is a filename
    filename = tokens[-1]
    if "/" in filename:
        return filename
    sample = sample_from_filename(filename)
    return f"{sample}/{filename}"


def main(argv=None):
    parser = argparse.ArgumentParser(description="Prefix filenames with sample parent directory")
    parser.add_argument("input", help="Input text file with filenames (one per line) or md5sum lines")
    parser.add_argument("-o", "--output", help="Write output to this file (defaults to stdout)")
    args = parser.parse_args(argv)

    out_stream = None
    try:
        with open(args.input, "r") as fh:
            lines = fh.readlines()
    except Exception as e:
        print(f"Error reading input file '{args.input}': {e}", file=sys.stderr)
        sys.exit(2)

    if args.output:
        out_stream = open(args.output, "w")
    else:
        out_stream = sys.stdout

    try:
        for line in lines:
            new = transform_line(line)
            if new is None:
                print(file=out_stream)
            else:
                print(new, file=out_stream)
    finally:
        if args.output and out_stream is not None:
            out_stream.close()


if __name__ == "__main__":
    main()
