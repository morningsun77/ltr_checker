"""Minimal streaming FASTA reader (no deps). Yields (name, seq) with name = header up to first whitespace."""
def read_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)

def write_fasta(fh, name, seq, width=0):
    fh.write(">" + name + "\n")
    if width and width > 0:
        for i in range(0, len(seq), width):
            fh.write(seq[i:i+width] + "\n")
    else:
        fh.write(seq + "\n")
