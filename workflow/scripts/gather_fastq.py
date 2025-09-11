from pathlib import Path
import re
import sys
from typing import List, Tuple, Dict

## FASTQ DETECTION DEFINITIONS ################################################

"""
Objective: given an input directory of FASTQ files, return lists of R1 and R2
files, paired by mate. This should be robust to commonly used file extensions:
.fastq / .FASTQ / .fq, etc. and compression formats. If any files are unmatched
they are returned separately.
"""

# --- Configuration ---

COMP_EXTS = (".gz", ".bz2", ".zst", ".xz")
BASE_EXTS = (".fastq", ".fq")

# Try explicit tokens first; then bare _1/_2 (no lookbehind needed).

READ_TOKEN_PATTERNS = [
    re.compile(r"([._-])R(?P<read>[12])([._-]|$)", re.IGNORECASE),     # _R1_
    re.compile(r"([._-])read(?P<read>[12])([._-]|$)", re.IGNORECASE),  # _read1_
    re.compile(r"([._-])(?P<read>[12])([._-]|$)"),                     # _1_
]

STRIP_DECOS = [
    re.compile(r"([._-])S\d+([._-]|$)", re.IGNORECASE),   # _S1_, _S23_
    re.compile(r"([._-])L\d{3}([._-]|$)", re.IGNORECASE), # _L001_
    re.compile(r"([._-])\d{3}([._-]|$)"),                 # trailing _001
]

def strip_fastq_exts(name: str) -> str:
    """Returns the file stripped of common FASTQ extensions (case insensitive).
    
    Strips common compression format extensions first, then FASTQ extensions.
    """

    low = name.lower()
    for cext in COMP_EXTS:
        if low.endswith(cext):
            name = name[: -len(cext)]
            low = name.lower()
            break
    for bext in BASE_EXTS:
        if low.endswith(bext):
            name = name[: -len(bext)]
            break
    return name

def detect_read_and_key(stem: str) -> Tuple[str | None, str | None]:
    """Determines if read 1 or 2 and returns shared identifier for FASTQ mates
    
    Tries to look for *_R1_*, *_read1_*, and *_1_* identifiers in turn, returns a
    shared ID stripped of lane IDs and other common  
    """
    base = stem + "."
    for rgx in READ_TOKEN_PATTERNS:
        m = rgx.search(base)
        if not m:
            continue
        read = m.group("read")
        # replace only the matched token with {READ}
        key = rgx.sub(lambda mm: f"{mm.group(1)}{{READ}}{mm.group(3)}", base, count=1)
        for deco in STRIP_DECOS:
            key = deco.sub(r"\1\2", key)
        key = re.sub(r"[._-]+", "_", key).strip("_")
        return read, key
    return None, None

def lane_sort_key(p: Path):
    s = p.name
    lane = re.search(r"_L(\d{3})", s)
    sidx = re.search(r"_S(\d+)", s)
    chunk = re.search(r"_([0-9]{3})(?=\.(?:fastq|fq)(?:\.(?:gz|bz2|zst|xz))?$)", s, re.IGNORECASE)
    return (
        int(lane.group(1)) if lane else 0,
        int(sidx.group(1)) if sidx else 0,
        int(chunk.group(1)) if chunk else 0,
        s,
    )

def find_fastq_pairs(input_dir: str, recursive: bool = True
                    ) -> Tuple[List[str], List[str], Dict[str, List[str]]]:
    root = Path(input_dir)
    files = root.rglob("*") if recursive else root.glob("*")

    def is_fastq(p: Path) -> bool:
        low = p.name.lower()
        return any(low.endswith(be) for be in BASE_EXTS) or any(
            low.endswith(be + ce) for be in BASE_EXTS for ce in COMP_EXTS
        )

    buckets: Dict[str, Dict[str, List[Path]]] = {}
    for p in files:
        if not (p.is_file() and is_fastq(p)):
            continue
        stem = strip_fastq_exts(p.name)
        read, key = detect_read_and_key(stem)
        if not (read and key):
            continue
        buckets.setdefault(key, {"1": [], "2": []})[read].append(p)

    r1_paths: List[str] = []
    r2_paths: List[str] = []
    unmatched_R1: List[str] = []
    unmatched_R2: List[str] = []

    for key in sorted(buckets):
        g = buckets[key]
        g["1"].sort(key=lane_sort_key)
        g["2"].sort(key=lane_sort_key)
        n = min(len(g["1"]), len(g["2"]))
        r1_paths.extend(str(p) for p in g["1"][:n])
        r2_paths.extend(str(p) for p in g["2"][:n])
        unmatched_R1.extend(str(p) for p in g["1"][n:])
        unmatched_R2.extend(str(p) for p in g["2"][n:])

    return r1_paths, r2_paths, {"R1": unmatched_R1, "R2": unmatched_R2}

## END ########################################################################


# --- Example ---
# if __name__ == "__main__":
#     r1, r2, unmatched = find_fastq_pairs("/path/to/fastqs")
#     print("Pairs:")
#     for a, b in zip(r1, r2):
#         print(a, "<->", b)
#     print("Unmatched:", unmatched)


if __name__ == "__main__":
    r1, r2, unmatched = find_fastq_pairs(sys.argv[1], recursive=True)

    print("Pairs:")
    for a, b in zip(r1, r2):
        print(a, "<->", b)
    print("Unmatched:", unmatched)

    with open(sys.argv[2], 'w') as f:
        f.write('\n'.join(r1))

    with open(sys.argv[3], 'w') as f:
        f.write('\n'.join(r2))
