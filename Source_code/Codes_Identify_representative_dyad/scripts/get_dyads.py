import gzip
import tempfile

import click
import pandas as pd
import io

def get_closer_maxima(fin, fout):
    """
    Select which of the putative dyads will be selected as the most representative one.
    We will select the dyad that has more reads and is close to the maximum. In case of a tie, we will select the one closer to the peak.
    :param fin:
    :param fout:
    :return:
    """
    df = pd.read_csv(fin, sep='\t',
                     names=['chr', 'start', 'end', 'reads', 'chr2', 'start2', 'end2', 'ID', 'score_wav', 'overlapp'])

    # get the distance of the midpoint to the wavelet peak
    df['distance_to_center_wavelet'] = (df['start2'] - df['start'] + 30).abs()  # need to remove the 30 bp we have added previously

    # sort them so that we will get those with higher reads, in case of tie we get the distance to the peak
    df.sort_values(by=['reads', 'distance_to_center_wavelet', 'start'], ascending=[False, True, True], inplace=True)

    # pandas is using quicksort algorithm,
    # which might not be a stable sort (meaning that the relative order of equal sort items is not preserved)
    # For that, we are including the start column in the sort

    intersected_file_sorted = tempfile.NamedTemporaryFile()

    # this is to make it faster
    df.to_csv(intersected_file_sorted.name, sep='\t', header=False, index=False, compression='gzip')

    done = set()
    keep_lines = []
    with gzip.open(intersected_file_sorted.name, 'rt') as infile:
        for line in infile:
            line_spl = line.rstrip().split('\t')

            # if the nucleosome-wavelet is not done
            if line_spl[7] not in done:

                # write it in the file
                out = '{}\t{}\t{}\t{}\t{}\t{}\n'.format(line_spl[0], line_spl[1], line_spl[2], line_spl[3], line_spl[8],
                                                        line_spl[10])
                keep_lines.append(out)

                # add it to the black list to avoid two dyads within one dyad.
                done.add(line_spl[7])

    # remove around 1% putative overlapping nucleosomes. Be aware this is a fix introduced after the Tag Paper was created
    selected = pd.read_csv(io.StringIO('\n'.join(keep_lines)), delim_whitespace=True,
                           names = ["chr","p1","p2","reads","smoothing_score","distance"])

    set_data = set()
    for chrom, data in selected.groupby(by='chr', sort=False):
        forbidden = set()
        for i, row in data.iterrows():
            p2 = row['p2']
            if p2 not in forbidden:
                set_data.add(i)

                # set of forbidden positions around the dyad
                for p in range(p2 - 147, p2 + 148):
                    forbidden.add(p)

    non_overlapp = selected.loc[set_data]
    non_overlapp.to_csv(fout, sep='\t', header=False, index=False,
                        compression='gzip')



@click.command()
@click.argument('input', metavar='INPUT', type=click.Path(exists=True))
@click.argument('output', metavar='OUTPUT', type=click.Path())
def cli(input, output):
    """
    Get the position of the dyads by sorting them by number of read and distance to 'expected' center
    """
    get_closer_maxima(input, output)


if __name__ == '__main__':
    cli()
