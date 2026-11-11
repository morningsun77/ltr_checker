#!/usr/bin/env python
# coding: utf-8
"""Emit a candidate FASTA + GFF3 from a genome-coordinate .scn and the genome FASTA.

This is the candidate-set deliverable that lets ltr_checker replace running LTR_Finder /
LTRharvest / LtrDetector separately: the merged candidate elements as sequences (FASTA)
and as annotations (GFF3), in addition to the LTR_retriever-ready candidates.scn.

seq-nr in the .scn is the 0-based index into chroms.txt (the FASTA record order). Element
and LTR coordinates are 1-based genome coordinates (LTRharvest -tabout columns):
    s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim seq-nr

Sequence is pulled by random access through the genome's .fai index (built by
`samtools faidx`), so a 14 Gb genome is handled without loading it into memory. The GFF3
is always written (coordinates only); the FASTA is skipped with --no-fasta.
"""
import argparse
import os
import sys


class Faidx:
    """Minimal random-access FASTA reader backed by a samtools .fai index."""

    def __init__(self, genome):
        self.genome = genome
        fai = genome + ".fai"
        if not os.path.exists(fai):
            raise FileNotFoundError(
                f"{fai} not found; build it with `samtools faidx {genome}`")
        self.index = {}
        for line in open(fai):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                # name -> (length, byte-offset of first base, bases/line, bytes/line)
                self.index[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
        self._fh = open(genome, "rb")

    def fetch(self, name, start, end):
        """Return the 1-based inclusive [start, end] subsequence of record `name` (str)."""
        if name not in self.index:
            return ""
        length, offset, lbases, lwidth = self.index[name]
        start = max(1, start)
        end = min(length, end)
        if end < start:
            return ""
        s0, e0 = start - 1, end - 1
        byte_start = offset + (s0 // lbases) * lwidth + (s0 % lbases)
        byte_end = offset + (e0 // lbases) * lwidth + (e0 % lbases) + 1
        self._fh.seek(byte_start)
        raw = self._fh.read(byte_end - byte_start)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii", "replace")

    def close(self):
        self._fh.close()


def load_chrom_order(chroms_txt):
    return [l.strip() for l in open(chroms_txt) if l.strip()]


def convert(scn, genome, chroms_txt, out_fasta, out_gff, write_fasta=True, source="ltr_checker"):
    order = load_chrom_order(chroms_txt)
    fa = Faidx(genome) if write_fasta else None
    fh_fa = open(out_fasta, "w") if write_fasta else None
    n = 0
    with open(out_gff, "w") as gff:
        gff.write("##gff-version 3\n")
        for line in open(scn):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            try:
                s, e = int(f[0]), int(f[1])
                lL_s, lL_e = int(f[3]), int(f[4])
                rL_s, rL_e = int(f[6]), int(f[7])
                sim = f[9]
                seqnr = int(f[-1])
            except (IndexError, ValueError):
                continue
            if e < s:
                s, e = e, s
            if seqnr < 0 or seqnr >= len(order):
                continue
            chrom = order[seqnr]
            n += 1
            eid = f"{chrom}:{s}-{e}"
            # GFF3: the element, then its two LTRs. Strand is '.' (unknown): the .scn
            # (LTRharvest -tabout form) carries no strand column, so it is not available
            # at this stage. Emitting '.' is correct for unknown strand; a hardcoded '+'
            # would be wrong for the ~half of elements that are on the minus strand.
            attrs = f"ID=cand{n};Name={eid};ltr_identity={sim}"
            gff.write(f"{chrom}\t{source}\tLTR_retrotransposon\t{s}\t{e}\t{sim}\t"
                      f".\t.\t{attrs}\n")
            gff.write(f"{chrom}\t{source}\tlong_terminal_repeat\t{lL_s}\t{lL_e}\t.\t"
                      f".\t.\tParent=cand{n};side=5prime\n")
            gff.write(f"{chrom}\t{source}\tlong_terminal_repeat\t{rL_s}\t{rL_e}\t.\t"
                      f".\t.\tParent=cand{n};side=3prime\n")
            if write_fasta:
                seq = fa.fetch(chrom, s, e)
                fh_fa.write(f">cand{n}::{eid}\n")
                for i in range(0, len(seq), 60):
                    fh_fa.write(seq[i:i + 60] + "\n")
    if fa is not None:
        fa.close()
    if fh_fa is not None:
        fh_fa.close()
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scn", required=True, help="genome-coordinate candidates .scn")
    ap.add_argument("--genome", required=True, help="genome FASTA (needs a .fai index)")
    ap.add_argument("--chrom-order", required=True, help="chroms.txt (FASTA record order)")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-gff", required=True)
    ap.add_argument("--no-fasta", action="store_true", help="write only the GFF3")
    a = ap.parse_args()
    n = convert(a.scn, a.genome, a.chrom_order, a.out_fasta, a.out_gff,
                write_fasta=not a.no_fasta)
    print(f"[scn_to_fasta_gff] {n} candidate elements -> {a.out_gff}"
          + ("" if a.no_fasta else f", {a.out_fasta}"), file=sys.stderr)


if __name__ == "__main__":
    main()
