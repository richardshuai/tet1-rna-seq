import sys
import pyensembl
import pandas as pd

def main():
    sys.stderr = open(snakemake.log[0], "w")
    data = pyensembl.Genome(
        reference_name='GRCm38',
        annotation_name='mus_musculus',
        gtf_path_or_url=snakemake.params.gtf)
    data.index()
    df = annotate(snakemake.input.dge, data)
    df.to_csv(snakemake.output.table, sep='\t')


def annotate(tsv_path, data):
    df = pd.read_table(tsv_path, sep='\t', index_col=0)
    df['gene_symbol'] = df.index
    df['gene_symbol'] = df['gene_symbol'].apply(lambda id: data.gene_by_id(id).gene_name)
    df.index = df['gene_symbol']
    df = df.drop('gene_symbol', axis=1)
    return df


if __name__ == '__main__':
    main()