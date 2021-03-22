#!/usr/bin/env python

import logging
from Bio import SeqIO
from io import StringIO
from random import sample
from subprocess import Popen, PIPE
import Hydraslayer.Utility as Utility

logger = logging.getLogger("Hydraslayer")


def revcomp(seq):
    maptable = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join([maptable[nuc] for nuc in seq[::-1]])


def derive_gap_position(gap, bound):
    start = gap.start-bound
    end   = gap.start+gap.length+bound

    contig_name = gap.contig.id.replace("|", "\\|")
    return f"{contig_name}:{start}-{end}"


def _get_matching_regions(bamfile, gap_pos):


    cmd = f"samtools view {bamfile} {gap_pos}"
    logger.debug(f"Samtools cmd: {cmd}")
    p = Popen(cmd, stdout=PIPE, shell=True)
    out = p.communicate()[0].decode().split("\n")

    return out


def _samrecords_to_fasta(sam_records, random_sample=None):
    fasta_fmt = ""
    c = 1

    if random_sample:
        if random_sample > len(sam_records):
            random_sample = len(sam_records)
        sam_records = sample(sam_records, random_sample)

    for record in sam_records:

        try:
            seq = record.split()[9]
            fasta_fmt += f">read{c}\n{seq}\n"

        except IndexError:
            pass

        c += 1

    return fasta_fmt


def _msa_kalign(fasta_reads):

    cmd = "kalign > /dev/stdout 2> /dev/null "
    #logger.debug(f"kalign2: {cmd}")
    p = Popen(cmd, shell=True, stdout=PIPE, stdin=PIPE)
    out = p.communicate(input=fasta_reads.encode())[0].decode()
    return out


def extract_gap_region(alignments_fasta):
    alignments = StringIO(alignments_fasta)

    sliced_alignments = []
    # get the position of the gap
    for record in SeqIO.parse(alignments, 'fasta'):
        if "read" not in record.id:
            N_positions = Utility.get_gap_pos_alignment(record)
            sliced_alignments.append(Utility.extract_positions(str(record.seq), N_positions))

    alignments = StringIO(alignments_fasta)
    for record in SeqIO.parse(alignments, 'fasta'):
          if "read" in record.id:
            sliced_alignments.append(Utility.extract_positions(str(record.seq), N_positions))

    return sliced_alignments


def determine_improvement(consensus, gap):
    if gap.gap_seq.count("N") > consensus.count("N"):  # then the new seq is better
        return True
    else:
        return False


class Fixer:

    def __init__(self, bamfile, arguments):
        self.bamfile = bamfile
        self.arguments = arguments

    def fix_gap(self, gap):
        gap_pos_fmt = derive_gap_position(gap, self.arguments.bounds)
        matching_reads = _get_matching_regions(self.bamfile, gap_pos_fmt)
        fasta_seqs = _samrecords_to_fasta(matching_reads, random_sample=self.arguments.align)
        fasta_seqs += f">gapseq\n{gap.extended_gap_seq}\n"

        alignments_fasta = _msa_kalign(fasta_seqs)

        sliced_alignments = extract_gap_region(alignments_fasta)
        assembly_seq = sliced_alignments.pop(0)

        # all read sequences are merged
        consensus = Utility.get_consensus(sliced_alignments, self.arguments.mincov)

        if determine_improvement(consensus, gap):
            return (True, consensus)
        else:
            return (False, None)