rule star_genome_generate:
    input:
        fasta = config["genome"]["fasta"],
        gtf = config["genome"]["gtf"]
    output:
        directory("resources/star_genome_generate/{genome}")
    threads: 4
    conda: "../envs/star-scte.yaml"
    resources:
        mem_mb = 40000,
        disk_mb = 50000,
        walltime = 360
    shell:
        """
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --genomeSAsparseD 3
        """

def STARSOLO_INPUT(wildcards):
    """
    Returns list of R1 fastqs, R2 fastqs, and genome build given a path to FASTQ directory in the .samples.csv 
    file.
    """
    fastq_dir = samples.loc[wildcards.dataset, "fastqs"]
    r1, r2, unmatched = scripts.gather_fastq.find_fastq_pairs(fastq_dir)
    return {
        "R1_fastqs": r1,
        "R2_fastqs": r2,
        "genome": "resources/star_genome_generate/{genome}"
    }

rule starsolo:
    input: unpack(STARSOLO_INPUT)
    output:
        # see STAR manual for additional output files
        "results/starsolo/{dataset}_{genome}/Aligned.sortedByCoord.out.bam",
        "results/starsolo/{dataset}_{genome}/SJ.out.tab"
    params:
        whitelist = config["umi_whitelist"],
        prefix = "results/starsolo/{dataset}_{genome}/",
        fastq_str = lambda wc, input: ",".join(input["R2_fastqs"]) + " " + ",".join(input["R1_fastqs"]),
        CBstart = config["soloCBstart"],
        CBlen = config["soloCBlen"],
        UMIstart = config["soloUMIstart"],
        UMIlen = config["soloUMIlen"],
        BCReadlen = config["soloBarcodeReadLength"]
    conda: "../envs/star-scte.yaml"
    threads: 12
    resources:
        mem_mb = 48000,
        disk_mb = 164000,
        walltime = 720
    shell:
        """
        STAR --runThreadN 48 \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist {params.whitelist} \
            --soloCBstart {params.CBstart} \
            --soloCBlen {params.CBlen} \
            --soloUMIstart {params.UMIstart} \
            --soloUMIlen {params.UMIlen} \
            --soloBarcodeReadLength {params.BCReadlen} \
            --genomeDir {input.genome} \
            --readFilesIn {params.fastq_str} \
            --readFilesCommand zcat \
            --outSAMattributes NH HI nM AS CR CY UR UY CB UB GX GN sS sQ sM \
            --outSAMtype BAM SortedByCoordinate \
            --soloUMIfiltering MultiGeneUMI \
            --soloCBmatchWLtype 1MM_multi_pseudocounts \
            --limitBAMsortRAM 16111457846 \
            --limitSjdbInsertNsj=2000000 \
            --outFilterMultimapNmax 100 \
            --winAnchorMultimapNmax 100 \
            --outSAMmultNmax 1 \
            --twopassMode Basic \
            --runRNGseed 42 \
            --outFileNamePrefix {params.prefix}
        """
