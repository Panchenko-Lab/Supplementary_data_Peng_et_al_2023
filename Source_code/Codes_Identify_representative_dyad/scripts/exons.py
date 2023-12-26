import click
import numpy as np
import pandas as pd


def get_exons(infile, CGS, outfile):

    df_cosmic = pd.read_csv(CGS)
    cancer_genes = set(df_cosmic['Gene Symbol'].tolist())

    df = pd.read_csv(infile,
                     sep='\t', names=['chr', 'type', 'class', 'start', 'end', 'd1', 'strand', 'd2', 'info'])
    df['transcript_id'] = df['info'].apply(lambda x: x.split(';')[1].split(' ')[2])
    df['exon_number'] = df['info'].apply(lambda x: int(x.split(';')[8].split(' ')[2]))
    df['gene_id'] = df['info'].apply(lambda x: x.split(';')[4].split(' ')[2].replace('"', ''))

    # remove cancer genes
    df = df[~df['gene_id'].isin(cancer_genes)]

    # remove unwanted exons (the first and the last)
    all_coordinates = []
    for gene, data in df.groupby(by='transcript_id'):
        not_wanted = [np.max(data['exon_number'].tolist()), np.min(data['exon_number'].tolist())]

        # this will remove the first and the last exon
        coord = data[~data['exon_number'].isin(not_wanted)][['chr', 'start', 'end']]
        if len(coord) > 0:
            all_coordinates.append(coord)

    df_wanted = pd.concat(all_coordinates)
    df_wanted.sort_values(by=['chr', 'start'], inplace=True)
    df_wanted.to_csv(outfile, columns=['chr', 'start', 'end'], sep='\t', header=False, index=False,
                     compression='gzip')


@click.command()
@click.argument('input', metavar='<REGIONS FILE>')
@click.argument('cgs', metavar='<CANCER GENE SENSUS>')
@click.argument('output', metavar='<OUT FILE>')
def cli(input, cgs, output):
    """
    Get the coordinates of the exons removing genes in the Cancer Gene Sensus
    and remove first and least exons of the rest
    """
    get_exons(input, cgs, output)


if __name__ == '__main__':
    cli()
