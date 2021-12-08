''' 
Snakemake rules for mapping and filtering 10x genomic 3' scRNA-seq 
input: FASTQ, alternative: "cellranger_count/{sample}/outs/possorted_genome_bam.bam"
files used in other parts: 
"out/{sample}/filtered/cbtags.bam"
'''

rule cellranger_count:
    ''' perform alignment and barcode and UMI correction 
    Note: only works on slurm cluster
    '''
    input:

    output:
        os.path.join(config['cellranger_dir'], 
            "{sample}", "outs", "possorted_genome_bam.bam")
    params:
        fastqs=lambda wildcards: os.path.abspath(samples.at[wildcards.sample, "fastqs"]),
        transcriptome=os.path.abspath(config["transcriptome"]),
        id=lambda wildcards: samples.at[wildcards.sample, "name"],
        out_path = config['cellranger_dir']
    threads: 8
    resources:
        mem_mb=65536,
        mem_gb=30
    shell: """
            cd {params.out_path}; \
            rm -r {wildcards.sample}/; \
            cellranger count \
                --id={wildcards.sample} \
                --transcriptome={params.transcriptome} \
                --fastqs={params.fastqs} \
                --sample={params.id} \
                --localcores={threads} \
                --localmem={resources.mem_gb} \
                --nosecondary
            """


rule filter_high_quality:
    ''' filter aligned reads by mapq quality scores '''
    input:
        os.path.join(config['cellranger_dir'], 
            "{sample}", "outs", "possorted_genome_bam.bam")
    output:
        temp(os.path.join(config['mapping_dir'], "{sample}", "filtered", "mapq.bam"))
    resources:
        mem_mb=2048
    params:
        mapqual=config["filter_high_quality"]["mapq"]
    conda: 
        os.path.join(workflow.basedir, config['envs_dir'],"samtools.yaml")
    shell: """
        samtools view -bh -q {params.mapqual} -o {output} {input}
        """


rule filter_non_CB_tags:
    ''' filter reads that don't have the CB tag '''
    input:
        os.path.join(config['mapping_dir'], "{sample}", "filtered", "mapq.bam")
    output:
        temp(os.path.join(config['mapping_dir'], "{sample}", "filtered", "cbtags_duplicated.bam"))
    params:
        tag = "CB"
    resources:
        mem_mb=2048
    conda: 
        os.path.join(workflow.basedir, config['envs_dir'],"samtools.yaml")
    shell: """
        samtools view -bh -d {params.tag} -o {output} {input}
        """


rule index_reads:
    input: 
        os.path.join(config['mapping_dir'], "{sample}", "filtered", "cbtags_duplicated.bam")
    output:
        temp(os.path.join(config['mapping_dir'], "{sample}", "filtered", "cbtags_duplicated.bam.bai"))
    resources:
        mem_mb=2048
    conda: 
        os.path.join(workflow.basedir, config['envs_dir'],"samtools.yaml")
    shell: """
        samtools index -b {input} > {output}
        """


rule de_duplicate:
    ''' de-duplicate the aligned reads with umi_tools; 
    remove reads that have the same cellular (CB) and unique (UB) tag and map to the same gene
    '''
    input:
        bam=os.path.join(config['mapping_dir'], "{sample}", "filtered", "cbtags_duplicated.bam"),
        bai=os.path.join(config['mapping_dir'], "{sample}", "filtered", "cbtags_duplicated.bam.bai")
    output:
        os.path.join(config['mapping_dir'], "{sample}", "filtered", "cbtags.bam")
    resources:
        mem_mb=lambda wildcards, input: (int(os.path.getsize(input.bam) / pow(1024,2) * 1.8) + 2000)
    log:
        os.path.join(workflow.basedir, config['local_log'], "de_duplicate", "{sample}.log")
    conda: 
        os.path.join(workflow.basedir, config['envs_dir'],"umi_tools.yaml")
    shell: """
        umi_tools dedup \
            --stdin={input.bam} \
            --log={log} \
            --extract-umi-method=tag \
            --umi-tag=UB --cell-tag=CB \
            --method=unique \
            --per-gene --gene-tag=GX --per-cell \
            > {output}
        """

