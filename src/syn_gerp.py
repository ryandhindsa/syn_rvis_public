#!/nfs/goldstein/software/python3.6.1-x86_64_shared/bin/python
"""
Purpose: calculate synGEPR
Author: Ryan Dhindsa
"""

# Imports --------------------------------------------------------------------
import re
import argparse
import gzip
import tabix
import logging
import sys
from Bio import SeqIO
from itertools import zip_longest
from annotate_gerp import query_gerp
from annotate_rscu import read_rscu_f
from filter_snvs import exclude_regions
from utils import read_hgnc_ensembl_map
from get_fh import get_fh
import pyfaidx
import numpy


# Globals --------------------------------------------------------------------
CHROMS = [str(c) for c in range(1, 23)] + ["X", "Y"]
CHROMS_SET = set(CHROMS)

FOURFOLD_DEGEN = set(["CTT", "CTC", "CTA", "CTG",
                      "GTT", "GTC", "GTA", "GTG",
                      "TCT", "TCC", "TCA", "TCG",
                      "CCT", "CCC", "CCA", "CCG",
                      "ACT", "ACC", "ACA", "ACG",
                      "GCT", "GCC", "GCA", "GCG",
                      "CGT", "CGC", "CGA", "CGG",
                      "GGT", "GGC", "GGA", "GGG"])

GTF_HEADER = ("chromosome", "source", "feature", "start", "end", "score",
              "strand", "frame", "attribute")

INTRONIC_BP_TO_EXCLUDE = 6
EXONIC_BP_TO_EXCLUDE = 2 


# Logger ---------------------------------------------------------------------
logging_level = logging.DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(logging_level)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging_level)
handler.setFormatter(logging.Formatter(
     "%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
logger.addHandler(handler)


# Classes --------------------------------------------------------------------
class Gene:
    def __init__(self, gene, chrom, strand):
        """
        Initializes class with gene name, chromosome, and strand
        :param gene: gene name
        :param chrom: chromosome that gene is on
        :param strand: strand that gene is on
        """
        self.gene = gene
        self.strand = strand
        self.chrom = str(chrom)
        self.cds = [] 
        self.bed = []

    def add_cds(self, start, end):
        """
        Adds start and end position for each exon in gene to self.cds
        """
        new_coords = range(start, end + 1)
        self.cds += new_coords

    def get_codon(self, genome, coords):
        """
        Retrieves codon by concatenating each nucleotide in the coords tuple
        :param genome: genome faidx object
        :param coords: tuple w/ three positions
        :return codon: the codon sequence
        """
        codon = ""
        for i in coords:
            nt = genome[self.chrom][i-1]

            if self.strand == "-":
                nt = nt.complement.seq
            else:
                nt = nt.seq

            codon += str(nt)

        return codon

    def calc_syn_gerp(self, genome, gerp_header, gerp_tb, rscu_d, 
                      bases_to_exclude):
        """
        Calculates mean GERP for all 4-fold degenerate sites in gene
        Iterates through all exons in self.cds and filters out codons that
        are not found in FOURFOLD_DEGEN. Then calculates mean GERP for
        these sites
        :param genome: pyfaidx genome object
        :param gerp_header: header of gerp file
        :param gerp_tb: gerp tabix object
        """
        if self.strand == "+":
            self.cds.sort()

        elif self.strand == "-":
            self.cds.sort(reverse=True)

        else:
            raise ValueError("Unrecognized strand: {}".format(self.strand))

        fourfold_gerp_scores = []

        i = 1
        for codon_coords in zip_longest(*[iter(self.cds)]*3):
            if None in codon_coords:
                continue

            codon = self.get_codon(genome, codon_coords)
            syn_pos = codon_coords[2]
            
            gerp = float(query_gerp(gerp_header, gerp_tb, 
                                    self.chrom, syn_pos))
            
            if codon in FOURFOLD_DEGEN and syn_pos not in bases_to_exclude:
                fourfold_gerp_scores.append(float(gerp))

            # get rscu score of codon for bed file
            if self.gene in rscu_d.keys():
                if codon not in rscu_d[self.gene].keys():
                    rscu = "NA"
                    
                rscu = rscu_d[self.gene][codon]
            else:
                rscu = "NA"

            bed_line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".\
                    format(self.chrom, syn_pos, self.strand,
                           self.gene, i, codon, rscu, gerp)

            i += 1
            self.bed.append(bed_line)

        self.syn_gerp = numpy.mean(fourfold_gerp_scores)
    

# Functions ------------------------------------------------------------------
def get_attribute_dict(gtf_attribute):
    """
    Converts values in "attribute" GTF column to a dictionary
    :param gtf_attribute: string of GTF attributes
    """
    d = {}
    for x in gtf_attribute.split(";"):
        if len(x) < 2:
            continue

        k, v = x.split()
        d[k] = v
    return(d)


def read_ensembl_gtf(ensembl_gtf_fh, hgnc_map, canonical_transcripts):
    """
    Reads in ENSEMBL GTF and stores CDS coordinates for each gene in 
    Gene object
    :param ensembl_gtf_fh: ensembl gtf
    :param hgnc_ensembl: map of ensembl --> HGNC conversions
    :returns all_genes: a dict of Gene objects
    """

    all_genes = {}
    
    with gzip.open(ensembl_gtf_fh, 'rt') as ensembl_gtf:
        for line in ensembl_gtf:
            if line.startswith("#"):
                continue

            gtf_dict = dict(zip(GTF_HEADER, line.strip().split("\t")))

            chrom = gtf_dict["chromosome"]

            if chrom not in CHROMS_SET:
                continue

            if gtf_dict["feature"] != "CDS":
                continue

            attr_dict = get_attribute_dict(gtf_dict["attribute"])
            transcript = attr_dict["transcript_id"].strip('"')

            if transcript in canonical_transcripts:
                start, end = [int(p) for p in (gtf_dict["start"], 
                                               gtf_dict["end"])]

                strand = gtf_dict["strand"]
                gene = hgnc_map[transcript]

                if gene in all_genes:
                    all_genes[gene].add_cds(start, end)

                else:
                    all_genes[gene] = Gene(gene, chrom, strand)
                    all_genes[gene].add_cds(start, end)
        
    return all_genes


def calc_all(all_genes, bases_to_exclude,
             rscu_fh, gerp_fp, genome_fa, syn_gerp_out, bed_out):
    """
    Calculates mean gerp score for all Gene objects contained in list
    all_genes and writes values to outfile
    :param all_genes: dict of Gene objects
    :param gerp_fp: path to tabix-indexed gerp file
    :param genome_fa: path to reference genome fasta that has been indexed
        via samtools faidx
    :param outfile: path to output file
    """
    with gzip.open(gerp_fp, 'rt') as gerp_f:
        gerp_header = gerp_f.readline()
        gerp_header = gerp_header.strip().split("\t")

    gerp_tb = tabix.open(gerp_fp)
    genome = pyfaidx.Fasta(genome_fa)
    
    rscu = read_rscu_f(rscu_fh)
    
    syn_gerp_out.write("#GENE\tSYN_GERP\n")

    bed_out.write("#CHROM\tPOS\tSTRAND\tGENE\tCDS_POS\tCODON\tRSCU\tGERP\n")

    for gene_obj in all_genes.values():
        gene_obj.calc_syn_gerp(genome, gerp_header, gerp_tb, rscu, 
                               bases_to_exclude)
        syn_gerp_out.write("{}\t{}\n".format(gene_obj.gene, \
                                             gene_obj.syn_gerp)) 
        for line in gene_obj.bed:
            bed_out.write(line)


# Main -----------------------------------------------------------------------
if __name__ == "__main__":
    # Argparse
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-e","--ensembl_gtf", 
                        help="Input file")
    
    parser.add_argument("--rscu", 
                        help="RSCU file", 
                        default="data/per_gene_rscu.txt",
                        type=get_fh)

    parser.add_argument("-he", "--hgnc-ensembl",
                        help="Path to HGNC ENSEMBL ID map",
                        default="data/ensembl/hgnc_canonical_ens75.txt",
                        type=argparse.FileType("r"))

    parser.add_argument("-r", "--ref-genome", 
                        help="Path to indexed reference genome fasta")

    parser.add_argument("-g", "--gerp", 
                        help="Path to GERP file",
                        default="data/gerp/gerp_merged.txt.gz")


    parser.add_argument("-o", "--output",
                        help="Path to output file",
                        type=argparse.FileType("w"))

    parser.add_argument("-bo", "--bed-out",
                        help="Path to output file which will contain GERP "
                        "and RSCU scores for each fourfold degenerate site" ,
                        type=argparse.FileType("w"))

    args = parser.parse_args()
    
    hgnc_map = read_hgnc_ensembl_map(args.hgnc_ensembl)
    canonical_transcripts = set(hgnc_map.keys())

    # Set canonical_only to FALSE, b/c even though we are only calculating
    # synGERP for canonical transcripts, want to eliminate any signal
    # that may be conserved due to splicing in alternative transcripts
    with gzip.open(args.ensembl_gtf, "rt") as ensembl_gtf:
        bases_to_exclude = exclude_regions(ensembl_gtf, 
                                           canonical_transcripts,
                                           EXONIC_BP_TO_EXCLUDE, 
                                           INTRONIC_BP_TO_EXCLUDE,
                                           canonical_only=False)

    all_genes = read_ensembl_gtf(args.ensembl_gtf, hgnc_map,
                                 canonical_transcripts)
    
    calc_all(all_genes, bases_to_exclude, 
             args.rscu, args.gerp, args.ref_genome, args.output,
             args.bed_out)
