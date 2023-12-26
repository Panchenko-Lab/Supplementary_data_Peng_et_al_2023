import gzip

import click


def remove(f_in, f_out):

    bad_set = set()
    done = set()

    with gzip.open(f_in, 'rt') as infile:
        for line in infile:

            line_spl = line.split('\t')
            distance = line_spl[10]
            # if the closest one is within 150 bp

            if int(distance) <= 147:
                id1 = '{}_{}'.format(line_spl[0], line_spl[1])
                id2 = '{}_{}'.format(line_spl[5], line_spl[6])
                bad_set.add(id1)
                bad_set.add(id2)

    with gzip.open(f_in, 'rt') as infile, gzip.open(f_out, 'wt') as outfile:
        for line in infile:
            line_spl = line.split('\t')
            id1 = '{}_{}'.format(line_spl[0], line_spl[1])
            id_complete = '\t'.join(line_spl[0:5])
            if (id1 not in bad_set) and (id_complete not in done):
                outfile.write('{}\n'.format(id_complete))
                done.add(id_complete)


@click.command()
@click.argument('input', metavar='<IN FILE>')
@click.argument('output', metavar='<OUT FILE>')
def cli(input, output):
    """
    Input and output parameters are names that will construct file names
    as <name>.chrX.gz (where X represents the chromosome)
    """
    remove(input, output)


if __name__ == '__main__':
    cli()
