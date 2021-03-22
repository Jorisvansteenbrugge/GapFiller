#!/usr/bin/python3

from Bio import SeqIO
import re
import logging
import Hydraslayer.Utility as Utility

logger = logging.getLogger("Hydraslayer")



class Gap:

    def __init__(self, gap_seq, contig, start, length, bound):
        logger.debug(f"Gap Identified in scaffold {contig.id} of {length}bp, starting at {start}")
        self.gap_seq = gap_seq
        self.contig  = contig
        self.start   = start
        self.length  = length

        self.extended_gap_seq = self._extend_gap_seq(bound)

    def __str__(self):
        return ">{0}\n{1}\n".format(self.contig.id, self.extended_gap_seq)

    def __repr__(self):
        return self.__str__()

    def _extend_gap_seq(self, bound=1000):
        if self.start - bound < 0:
            l_start = 0
        else:
            l_start = self.start - bound

        end = self.start + self.length + bound

        gapseq = str(self.contig.seq)[l_start: end + 1]
        return gapseq


class Assembly:

    def __init__(self, assembly_file):
        self.assembly_file = assembly_file
        self.contigs = self._get_contigs(assembly_file)
        self.gaps = []

    def _get_contigs(self, assembly_file):
        return tuple(SeqIO.parse(assembly_file, 'fasta'))

    def get_size(self):
        return sum([len(contig.seq) for contig in self.contigs])

    def resolve_gap(self, gap, gap_filling):

        contigname = gap.contig.id

        contig_idx = [idx for idx, val in enumerate(self.contigs) if val.id == contigname][0]

        contig_seq_old = str(self.contigs[contig_idx].seq)

        assert len(gap.gap_seq) == len(gap_filling), "someting went wrong (<0.o>)"

        contig_seq = contig_seq_old.replace(gap.gap_seq, gap_filling)

        self.contigs[contig_idx].seq = contig_seq

    def _identify_gaps(self, contig, minlen, bound):

        seq = str(contig.seq)
        if "N" * minlen not in seq:
            return None

        s = re.search(r'N+', seq).group(0)
        self.gaps.append(Gap(s, contig, seq.index(s), len(s), bound))

    def get_total_gap_length(self):
        N_c = 0

        for record in SeqIO.parse(self.assembly_file, 'fasta'):
            seq = str(record.seq)
            N_c += seq.count("N")

        return N_c

    def get_gaps(self, minlen, bound):

        for contig in self.contigs:
            self._identify_gaps(contig, minlen, bound)

        return self.gaps

    def write_to_output(self, outfile):
        with open(outfile, 'w') as out_file:
            for contig in self.contigs:
                out_file.write(f">{contig.id}\n{contig.seq}\n")

