import functools
import gzip
from itertools import islice

import click
import pandas as pd


OUTPUT_NAME = None


def read_file(f):
    """
    Read the input file
    :param f: file
    :return: dataframe
    """
    df = pd.read_csv(f, sep='\t', names=['Chr', 'Pos-1', 'Position', 'Counts'], usecols=['Chr', 'Position', 'Counts'])
    return df


@functools.lru_cache(100)
def kernel_triweight(n):
    """
    Formula of the triweighted kernel
    :param n:
    :return:
    """
    
    x = (n-15)/15
    if abs(x) > 1:
        k = 0
    else:
        k = (35/32)*(1-x**2)**3  # we are including the hash the kernel
        
    return k


def window(seq, n=2):
    """
    Returns a sliding window (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...
    """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def kernel_density_reads(data, output_file):
    """
    Apply a kernel density smoothing to all the mapped reads
    :param data:
    :param output_file:
    :return:
    """

    chrom = data['Chr'].values[0]
    
    # get a dictionary where each position has a count
    dic = dict(zip(data.Position, data.Counts))
    
    # get a list where each position has the dyad count
    list_counts = []
    for i in range(max(data['Position']) + 30):
        list_counts.append(dic.get(i, 0))

    # define outfile
    with gzip.open(output_file, 'wt') as outfile:
        
        # this will return a window of N nucleotides from the big list
        for ix, win in enumerate(window(list_counts, 31)):

            val = 0
            # position and readcount of each window
            for n, c in enumerate(win):
                if c > 0:
                    # apply function
                    k = kernel_triweight(n)
                    # add value
                    val += k*c
            # write the value to the position in the middle of the window
            out = '{}\t{}\t{}\t{}\n'.format(chrom, ix+14, ix+15, val)
            outfile.write(out)


@click.command()
@click.argument('input', metavar='<BED FILE>', type=click.Path(exists=True))
@click.argument('output', metavar='<OUT FILE>', type=click.Path())
def cli(input, output):
    """
    Compute the kernel-smoothed dyad count.

    The input file must contain data of a single chromosome
    """
    kernel_density_reads(read_file(input), output)


if __name__ == '__main__':
    cli()
