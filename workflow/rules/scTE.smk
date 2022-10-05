## need to init scte... 
rule scte_install:
    output: 
        directory("resources/scTE")
    conda:
        "envs/star-scte.yaml"
    shell:
        """
        git clone https://github.com/JiekaiLab/scTE.git "resources/scTE"
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
        "envs/star-scte.yaml"
    shell:
        """
        scTE_build \
            -g {wildcards.genome} \
            -m {wildcards.mode} \
            -o {params.prefix}
        """

rule scte_quant:
    input:
        bam = "results/starsolo/{dataset}_{genome}/Aligned.sortedByCoord.out.bam",
        index = "resources/scte_build/{genome}.{mode}.idx",
    output:
        directory("results/scte_quant/{dataset}.{genome}.{mode}")
    params:
        genome = config["genome"]["name"],
        expect = 30000,
        min_counts = 3000
    conda:
        "envs/star-scte.yaml"
    threads: 48
    shell:
        """"
        scTE -i {input.bam} \
            -o {output} \
            -g {params.genome} \
            -x {input.index} \
            -p {threads} \
            --keeptmp True \
            --expect-cells {params.expect} \
            --min_counts {params.min_counts}
        """
    