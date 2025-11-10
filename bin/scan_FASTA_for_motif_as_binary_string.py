import argparse
import re
import pysam
import logging

def getParams():
    '''Parse parameters from the command line'''
    parser = argparse.ArgumentParser(description='''
    This script scans a FASTA file for a motif and reports two binary CDT matrices (0/1) indicating the presence of the motif.
    It reports the position of the last base in the motif match for both sense and antisense strands.
    
    Example: python scan_FASTA_for_motif_as_binary_matrix.py -i input.fasta -m AG -o output
    ''')
    parser.add_argument('-i', '--input', metavar='input', required=True, help='Input FASTA file')
    parser.add_argument('-m', '--motif', metavar='iupac_motif', required=True, help='Motif to search for (IUPAC notation supported)')
    parser.add_argument('-o', '--output', metavar='output', required=True, help='Base name for output CDT files (e.g., output will create output_sense.cdt and output_anti.cdt)')

    args = parser.parse_args()
    return args

# IUPAC nucleotide codes
iupac_dict = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
}

def iupac_to_regex(motif):
    '''Converts an IUPAC motif to a regex pattern.'''
    return ''.join(iupac_dict.get(base, base) for base in motif)

def reverse_complement(sequence):
    '''Returns the reverse complement of a DNA sequence.'''
    complement = str.maketrans("ACGT", "TGCA")
    return sequence.translate(complement)[::-1]

def scan_fasta_for_motif(fasta_file, motif, output_base):
    '''Scans both strands of a FASTA file and writes two separate CDT files for sense and antisense strand matches.'''
    
    regex_pattern = iupac_to_regex(motif)
    motif_regex = re.compile(regex_pattern)
    
    try:
        fasta = pysam.FastaFile(fasta_file)
    except FileNotFoundError:
        print(f"Error: FASTA file {fasta_file} not found.")
        return

    with open(f'{output_base}_sense.cdt', 'w') as sense_file, open(f'{output_base}_anti.cdt', 'w') as anti_file:
        header_written = False

        for seq_id in fasta.references:
            sequence = fasta.fetch(seq_id)
            
            # Get sequence length
            n = len(sequence)
            
            # Create vectors for both strands
            sense_matches = ['0'] * n
            anti_matches = ['0'] * n
            
            # Scan sense strand
            for match in motif_regex.finditer(sequence):
                match_end = match.start() + len(motif) - 1
                sense_matches[match_end] = '1'
            
            # Scan antisense strand
            rev_comp_sequence = reverse_complement(sequence)
            for match in motif_regex.finditer(rev_comp_sequence):
                match_end = match.start() + len(motif) - 1
                anti_matches[len(sequence) - match_end - 1] = '1'
            
            # Write header once
            if not header_written:
                # Calculate the symmetric range from -(n-2)/2 to n/2
                start_pos = -(n - 2) // 2
                end_pos = n // 2
                header = 'YORF\tNAME\t' + '\t'.join(map(str, range(start_pos, end_pos + 1))) + '\n'
                sense_file.write(header)
                anti_file.write(header)
                header_written = True
            
            # Write results to sense and antisense CDT files
            sense_file.write(f'{seq_id}\t{seq_id}\t' + '\t'.join(sense_matches) + '\n')
            anti_file.write(f'{seq_id}\t{seq_id}\t' + '\t'.join(anti_matches) + '\n')

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    args = getParams()
    logging.info(f"Starting motif scan on {args.input} for motif {args.motif}")
    scan_fasta_for_motif(args.input, args.motif, args.output)
