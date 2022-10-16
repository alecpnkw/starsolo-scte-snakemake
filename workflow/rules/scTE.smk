rule scte_install:
    output: 
        directory("resources/scTE")
    conda:
        "../envs/star-scte.yaml"
    shell:
        """
        git clone https://github.com/alecpnkw/scTE.git "resources/scTE"
        cd resources/scTE
        python setup.py install
        """

rule scte_build:
    input:
        "resources/scTE"
    output:
        "resources/scte_build/{genome}.{mode}.idx"
    params:
        prefix = "resources/scte_build/{genome}"
    conda:
        "../envs/star-scte.yaml"
    threads: 4
    resources:
        mem_mb = 5000,
        disk_mb = 20000
    shell:
        """
        scTE_build \
            -g {wildcards.genome} \
            -m {wildcards.mode} \
            -o {params.prefix}
        """

rule clean_bam:
    input:
        "results/starsolo/{dataset}_{genome}/Aligned.sortedByCoord.out.bam"
    conda:
        "../envs/star-scte.yaml"
    output:
        "results/clean_bam/{dataset}_{genome}/Aligned.sortedByCoord.out.bam"
    shell:
        """
        samtools view {input} -h | awk '/^@/ || /CB:/' | samtools view -h -b > {output}
        """

rule scte_quant:
    input:
        bam = "results/clean_bam/{dataset}_{genome}/Aligned.sortedByCoord.out.bam",
        index = "resources/scte_build/{genome}.{mode}.idx",
    output:
        "results/scte_quant/{dataset}.{genome}.{mode}.csv"
    params:
        expect = config["scte_expect_cells"],
        min_counts = config["scte_min_counts"],
        prefix = "results/scte_quant/{dataset}.{genome}.{mode}"
    resources:
        mem_mb = 13000,
        disk_mb = 55000
    conda:
        "../envs/star-scte.yaml"
    threads: 24
    shell:
        """
        scTE -i {input.bam} \
            -o {params.prefix} \
            -x {input.index} \
            -p 96 \
            --keeptmp True \
            --expect-cells {params.expect} \
            --min_counts {params.min_counts}
        """
    