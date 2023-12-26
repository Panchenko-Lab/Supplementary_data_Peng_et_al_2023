import click
import pandas as pd


def normalize_kernel(file_in, file_out):
    """
    This will apply a rolling window to normalize the score of each position
    :param file_in:
    :param file_out:
    :return:
    """

    df = pd.read_csv(file_in, sep='\t', names=['chr', 'pos1', 'pos2', 'score'])

    # do the rolling window
    df['centered'] = df['score'].rolling(151, center=True).sum()
    
    df.dropna(inplace=True)
    wanted = df[df['score'] > 0].copy()
    
    wanted['smoothed_score'] = wanted['score'] / wanted['centered']
    wanted[['chr', 'pos1', 'pos2', 'smoothed_score']].to_csv(file_out, sep='\t', header=False, index=False,
                                                             compression='gzip')


@click.command()
@click.argument('input', metavar='<IN FILE>', type=click.Path(exists=True))
@click.argument('output', metavar='<OUT FILE>', type=click.Path())
def cli(input, output):
    """
    Correct the kernel-smoothed dyad count to compute the
    stringency metric.

    Input and output parameters are names that will construct file names
    as <name>.chrX.gz (where X represents the chromosome)
    """
    normalize_kernel(input, output)


if __name__ == '__main__':
    cli()
