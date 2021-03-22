import logging

logger = logging.getLogger("Hydraslayer")


def get_consensusbase(bases, mincov=3):
    """

    :param mincov:
    :type bases: list
    """

    bases = "".join(bases)
    a = bases.count('A')
    t = bases.count('T')
    c = bases.count('C')
    g = bases.count('G')
    n = bases.count("N") + bases.count('-')

    counts = [(a, 'A'), (t, 'T'), (c, 'C'), (g, 'G')]
    s_dic = sorted(counts, key=lambda x: x[0], reverse=True)

    max = s_dic[0]

    if max[0] < mincov:
         return "N"
    else:
        return max[1]

def get_gap_pos_alignment(record):
    sequence = str(record.seq)
    N_pos = [x for x, nuc in enumerate(sequence) if nuc == "N"]

    return N_pos


def extract_positions(seq, positions):
    return "".join([seq[idx] for idx in positions])


def pretty_print_alignment(fasta_sequences):
    alignment_len = len(fasta_sequences[0])

    for x in range(alignment_len):
        row = []
        for alignment in fasta_sequences:
            row.append(alignment[x])

        print(" ".join(row))


def get_consensus(fasta_seqs, mincov):
    """Get the per-position consensus sequence, excluding gaps.
        All read sequences (not the assembly sequence) are merged into a consensus sequence.
    """
    consensus = []
    alignment_len = len(fasta_seqs[0])

    for x in range(alignment_len):

        bases = [fasta[x] for fasta in fasta_seqs]
        consensus.append(get_consensusbase(bases, mincov))

   # logger.debug(f'Consensus {"".join(consensus)}')
    return "".join(consensus)
