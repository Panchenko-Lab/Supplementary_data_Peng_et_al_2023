import gzip

import click
import pandas as pd


def filter_(nucleosomes_file, sizes, output_file, l=58):
    df = pd.read_csv(sizes, sep='\t', names=['chr', 'size'])
    chrom_sizes = dict(zip(df['chr'], df['size']))

    with gzip.open(nucleosomes_file, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
        for line in infile:
            line_spl = line.rstrip().split('\t')
            chrom = line_spl[0]
            start = int(line_spl[1])
            if (start > l) & (start < chrom_sizes[chrom]-l):
                outfile.write(line)


@click.command()
@click.argument('infile', metavar='NUCLEOSOMES', type=click.Path(exists=True))
@click.argument('sizes', metavar='<CHR SIZES>', type=click.Path(exists=True))
@click.argument('outfile', metavar='OUTPUT', type=click.Path())
def cli(infile, sizes, outfile):
    """
    Remove nucleosomes failling in the first or last 58 bp
    """
    filter_(infile, sizes, outfile, l=58)


if __name__ == '__main__':
    cli()
