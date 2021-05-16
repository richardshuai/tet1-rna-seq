rule count_matrix:
    input:
        get_star_output,
    output:
        "results/counts/all.tsv",
    log:
        "logs/count-matrix.log",
    params:
        samples=units["sample_name"].tolist(),
        strand=get_strandedness(units),
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"


rule deseq2_init:
    input:
        counts="results/counts/all.tsv",
    output:
        "results/deseq2/all.rds",
    params:
        samples=config["samples"],
        model=config["diffexp"]["model"],
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule deseq2:
    input:
        "results/deseq2/all.rds",
    output:
        table=report(
            "results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"
        ),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst"),
    params:
        contrast=get_contrast,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log",
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"

rule pca:
    input:
        "results/deseq2/all.rds",
    output:
        report("results/pca.svg", "../report/pca.rst"),
    params:
        pca_labels=config["pca"]["labels"],
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/pca.log",
    script:
        "../scripts/plot-pca.R"

rule tmm_norm:
    input:
        counts="results/counts/all.tsv"
    output:
        table="results/counts/TMM_normalized.tsv"
    params:
        samples=config["samples"],
    conda:
        "../envs/edgeR.yaml"
    log:
        "logs/tmm_norm.log"
    script:
        "../scripts/normalize_TMM.R"

rule make_table:
    input:
        counts="results/counts/TMM_normalized.tsv",
        dge="results/diffexp/{contrast}.diffexp.tsv"
    output:
        table="results/tables/{contrast}.table.tsv",
    params:
        gtf="resources/genome.gtf",
    conda:
        "../envs/pyensembl.yaml"
    log:
        "logs/{contrast}.make_table.log"
    script:
        "../scripts/make_table.py"