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


rule tcc:
    input:
        counts="results/counts/all.tsv",
    output:
        dge="results/tcc/diffexp.tsv"
    params:
        samples=config["samples"],
    conda:
        "../envs/tcc.yaml"
    log:
        "logs/deseq2/init.log",
    threads: get_tcc_threads()
    script:
        "../scripts/tcc.R"


rule annotate_table:
    input:
        dge="results/tcc/diffexp.tsv"
    output:
        table="results/tables/table.tsv",
    params:
        gtf="resources/genome.gtf",
    conda:
        "../envs/pyensembl.yaml"
    log:
        "logs/annotate_table.log"
    script:
        "../scripts/annotate_table.py"


rule heatmap:
    input:
        dge="results/tables/table.tsv",
    params:
        samples=config["samples"]
    output:
        heatmap="results/heatmap.jpeg"
    conda:
        "../envs/pheatmap.yaml"
    log:
        "logs/heatmap.log"
    script:
        "../scripts/plot-heatmap.R"