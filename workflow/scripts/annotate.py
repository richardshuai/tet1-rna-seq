import sys
import pyensembl
import pandas as pd
import os


def main(tsv_path):
    """
    annotate.py

    Manually annotate a file with gene symbols from Ensemble ID.
    :param tsv_path: Path for TSV to annotate
    """
    data = pyensembl.Genome(
        reference_name='GRCm38',
        annotation_name='mus_musculus',
        gtf_path_or_url='resources/genome.gtf')
    data.index()

    df = annotate(tsv_path, data)

    parent_dir = os.path.dirname(tsv_path)
    name, ext = os.path.splitext(os.path.basename(tsv_path))
    out_dir = os.path.join(parent_dir, "{}_anno{}".format(name, ext))

    df.to_csv(out_dir, sep='\t')


def annotate(tsv_path, data):
    df = pd.read_table(tsv_path, sep='\t', index_col=0)
    df['gene_symbol'] = df.index
    df['gene_symbol'] = df['gene_symbol'].apply(lambda id: data.gene_by_id(id).gene_name)
    df.index = df['gene_symbol']
    df = df.drop('gene_symbol', axis=1)
    return df


if __name__ == '__main__':
    tsv_path = sys.argv[1]
    main(tsv_path)
