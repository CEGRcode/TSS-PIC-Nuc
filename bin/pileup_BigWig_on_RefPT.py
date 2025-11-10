import sys
import argparse
import pyBigWig

def getParams():
    '''Parse parameters from the command line.'''
    parser = argparse.ArgumentParser(description='''
This script will pile up BigWig scores.

Example: python pileup_BigWig_on_RefPT.py -i INPUT.bw -r REF.bed -o OUTPUT.bed''')
    parser.add_argument('-i', '--input', metavar='bigwig_fn', required=True, help='A BigWig file of signal to pile up.')
    parser.add_argument('-r', '--reference', metavar='bed_fn', required=True, help='A BED file with coordinates to calculate distances from.')
    parser.add_argument('-o', '--output', metavar='bed_fn', required=True, help='A BED file to output the distance scores.')

    return parser.parse_args()

def loadBedGraph(bg_fn):
    '''Load BedGraph data into a dictionary.'''
    coord2score = {}
    with open(bg_fn, 'r') as reader:
        for line in reader:
            if line.startswith("chrom"):
                continue
            tokens = line.strip().split('\t')
            coord2score.setdefault(tokens[0], {})
            for i in range(int(tokens[1]), int(tokens[2])):
                coord2score[tokens[0]][i] = float(tokens[3])
    return coord2score

if __name__ == "__main__":
    '''Gets per-base BigWig signal within BED coordinate regions to build a pileup.'''
    
    args = getParams()
    
    # Initialize header indicator
    headerWritten = False

    # Open BigWig file
    bw = pyBigWig.open(args.input)

    # Open output file for writing
    with open(args.output, 'w') as writer:
        # Parse reference BED file
        with open(args.reference, 'r') as reader:
            for line in reader:
                tokens = line.strip().split('\t')

                # Extract reference information
                chrstr = tokens[0]
                start = int(tokens[1])
                stop = int(tokens[2])
                my_id = tokens[3]
                strand = tokens[5]

                # Write header if not yet written
                if not headerWritten:
                    writer.write("YORF\tNAME\t" + "\t".join([str(i) for i in range(stop - start)]) + "\n")
                    headerWritten = True

                # Fetch values from BigWig
                scoreDict = {}
                fetched_intervals = bw.intervals(chrstr, start, stop)
                if fetched_intervals is not None:
                    for interval in fetched_intervals:
                        for index in range(interval[0], interval[1]):
                            scoreDict[index] = interval[2]

                # Prepare values for output
                values = [str(scoreDict.get(i, 0)) for i in range(start, stop)]
                if strand == "-":
                    values.reverse()  # Reverse values for negative strand

                # Write results to output
                writer.write(f"{my_id}\t{my_id}\t" + "\t".join(values) + "\n")

    # Close the BigWig file
    bw.close()
