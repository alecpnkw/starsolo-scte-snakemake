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
        fastqs = config["fastqs"]
    output:
        # see STAR manual for additional output files
        "results/starsolo/{dataset}_{genome}/Aligned.sortedByCoord.out.bam",
        "results/starsolo/{dataset}_{genome}/Aligned.toTranscriptome.out.bam",
        "results/starsolo/{dataset}_{genome}/SJ.out.tab",
	    "results/starsolo/{dataset}_{genome}/ReadsPerGene.out.tab"
    params:
        whitelist = config["umi_whitelist"],
        fastq_str = lambda wc: ",".join(config["fastqs"])
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
            --soloUMIlen 12 \
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
            --runRNGseed 42
        """
        