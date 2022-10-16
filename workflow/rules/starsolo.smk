rule star_genome_generate:
    input:
        fasta = config["genome"]["fasta"],
        gtf = config["genome"]["gtf"]
    output:
        directory("resources/star_genome_generate/{genome}")
    threads: 4
    conda: "../envs/star-scte.yaml"
    resources:
        mem_mb = 32000,
        disk_mb = 50000
    shell:
        """
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --genomeSAsparseD 3
        """

rule starsolo:
    input:
        genome = "resources/star_genome_generate/{genome}",
        R1_fastqs = config["R1_fastqs"],
        R2_fastqs = config["R2_fastqs"]
    output:
        # see STAR manual for additional output files
        "results/starsolo/{dataset}_{genome}/Aligned.sortedByCoord.out.bam",
        "results/starsolo/{dataset}_{genome}/SJ.out.tab"
    params:
        whitelist = config["umi_whitelist"],
        prefix = "results/starsolo/{dataset}_{genome}/",
        fastq_str = lambda wc: ",".join(config["R2_fastqs"]) + " " + ",".join(config["R1_fastqs"]),
        CBstart = config["soloCBstart"],
        CBlen = config["soloCBlen"],
        UMIstart = config["soloUMIstart"],
        UMIlen = config["soloUMIlen"]
    conda: "../envs/star-scte.yaml"
    threads: 24
    resources:
        mem_mb = 2000,
        disk_mb = 50000
    shell:
        """
        STAR --runThreadN 48 \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist {params.whitelist} \
            --soloCBstart {params.CBstart} \
            --soloCBlen {params.CBlen} \
            --soloUMIstart {params.UMIstart} \
            --soloUMIlen {params.UMIlen} \
            --genomeDir {input.genome} \
            --readFilesIn {params.fastq_str} \
            --readFilesCommand zcat \
            --outSAMattributes NH HI nM AS CR CY UR UY CB UB GX GN sS sQ sM \
            --outSAMtype BAM SortedByCoordinate \
            --soloUMIfiltering MultiGeneUMI \
            --soloCBmatchWLtype 1MM_multi_pseudocounts \
            --limitBAMsortRAM 16111457846 \
            --outFilterMultimapNmax 100 \
            --winAnchorMultimapNmax 100 \
            --outSAMmultNmax 1 \
            --twopassMode Basic \
            --runRNGseed 42 \
            --outFileNamePrefix {params.prefix}
        """
        
