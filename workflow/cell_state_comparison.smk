'''
Module containing rules to perform AUC analysis for given cell types between sample origins (i.e. cell states)
'''


rule label_randomisation:
    '''
    Randomise the 3p end coverage between two samples
    memory: based on previous runs, made regression for MB usage of input file sizes.
    '''
    input:
        bed = os.path.join(config['mapping_dir'], "terminal_exons.bed"),
        one_plus = expand(os.path.join(config['out_dir'], "coverage_3p", 
          "{{cell_type}}_{origin}_plus.corr.bedgraph"),
          origin = SAMPLE_ORIGIN[0]),
        one_minus = expand(os.path.join(config['out_dir'], "coverage_3p", 
          "{{cell_type}}_{origin}_minus.corr.bedgraph"),
          origin = SAMPLE_ORIGIN[0]),
        two_plus = expand(os.path.join(config['out_dir'], "coverage_3p", 
          "{{cell_type}}_{origin}_plus.corr.bedgraph"),
          origin = SAMPLE_ORIGIN[1]),
        two_minus = expand(os.path.join(config['out_dir'], "coverage_3p", 
          "{{cell_type}}_{origin}_minus.corr.bedgraph"),
          origin = SAMPLE_ORIGIN[1])
    output:
        one_plus = expand(os.path.join(config['out_dir'], "randomised_coverage", 
          "{{cell_type}}_{origin}_plus.corr.bedgraph"),
          origin = SAMPLE_ORIGIN[0]),
        one_minus = expand(os.path.join(config['out_dir'], "randomised_coverage", 
          "{{cell_type}}_{origin}_minus.corr.bedgraph"),
          origin = SAMPLE_ORIGIN[0]),
        two_plus = expand(os.path.join(config['out_dir'], "randomised_coverage", 
          "{{cell_type}}_{origin}_plus.corr.bedgraph"),
          origin = SAMPLE_ORIGIN[1]),
        two_minus = expand(os.path.join(config['out_dir'], "randomised_coverage", 
          "{{cell_type}}_{origin}_minus.corr.bedgraph"),
          origin = SAMPLE_ORIGIN[1])
    params:
        script = os.path.join(config['script_dir'], "label_randomisation.py")
    threads: 8
    resources:
        mem_mb = lambda wildcards, input: int(sum([os.path.getsize(x) / pow(1024,2) for x in input]) * 20 + 1000)
    log: os.path.join(config['local_log'], "label_randomisation", "labels_{cell_type}.log")
    conda: os.path.join(workflow.basedir, config['envs_dir'], "python_basics.yaml")
    shell:
        """
        python {params.script} -bed {input.bed} \
          -s1p {input.one_plus} -s1m {input.one_minus} \
          -s2p {input.two_plus} -s2m {input.two_minus} \
          -o1p {output.one_plus} -o1m {output.one_minus} \
          -o2p {output.two_plus} -o2m {output.two_minus} \
          --cores {threads} --log-file {log}
        """

rule cumulative_coverage_3p_end:
    '''
    Compute reverse cumulative distribution in BED intervals from 3p end coverage.
    threads: if file < 80MB, use 1 thread, else 2.
    memory: based on previous runs, estimate MB usage from input file size.
    '''
    input:
        coverage = os.path.join(config['out_dir'], 
            "{cov_dir}", "{cell_type}_{origin}_{strand}.corr.bedgraph"),
        bed = os.path.join(config['mapping_dir'], "terminal_exons.bed")
    output:
        cumsum = os.path.join(config['out_dir'], 
            "{cov_dir}", "{cell_type}_{origin}_cumsum_{strand}.bedgraph")
    params:
        strand = "{strand}",
        script = os.path.join(config['script_dir'], "cumulative_coverage.py")
    threads: 
        lambda wildcards, input: 1 if (os.path.getsize(input.coverage) / 
          pow(1024,2)) < 50 else 3
    resources:
        mem_mb = lambda wildcards, input: 4096 if (os.path.getsize(input.coverage) / 
          pow(1024,2)) < 50 else 8192
    log: os.path.join(config['local_log'], "cumulative_coverage_3p_end", "{cov_dir}", "{cell_type}_{origin}_{strand}.log")
    conda: os.path.join(workflow.basedir, config['envs_dir'], "intervaltree.yaml")
    shell:
        """
        python {params.script} --bed-file {input.bed} \
          --coverage {input.coverage} \
          --output {output.cumsum} \
          --strand {params.strand} \
          --cores {threads} \
          --log-file {log}
        """

rule calculation_auc:
    '''
    Calculate the difference in cumulative distribution in each TE region 
    as the Area Under the Curve (AUC) between the sample origin for the same cell type.
    Set CPM threshold = 0 for returning all TEs in BED file. CPM = 1 returns all non-zero TEs.
    '''
    input:
        bed = os.path.join(config['mapping_dir'], "terminal_exons.bed"),
        one_plus = expand(os.path.join(config['out_dir'], "{{cov_dir}}", 
          "{{cell_type}}_{origin}_cumsum_plus.bedgraph"),
          origin = SAMPLE_ORIGIN[0]),
        one_minus = expand(os.path.join(config['out_dir'], "{{cov_dir}}", 
          "{{cell_type}}_{origin}_cumsum_minus.bedgraph"),
          origin = SAMPLE_ORIGIN[0]),
        two_plus = expand(os.path.join(config['out_dir'], "{{cov_dir}}", 
          "{{cell_type}}_{origin}_cumsum_plus.bedgraph"),
          origin = SAMPLE_ORIGIN[1]),
        two_minus = expand(os.path.join(config['out_dir'], "{{cov_dir}}", 
          "{{cell_type}}_{origin}_cumsum_minus.bedgraph"),
          origin = SAMPLE_ORIGIN[1]),
        script = os.path.join(config['script_dir'], "calculation_auc.py")
    output:
        auc_bed = os.path.join(config['out_dir'], 
          "auc", "{cov_dir}", "auc_{cell_type}.tsv")
    params:
        cpm_threshold = config['calculation_auc']['cpm_threshold'],
        quants = config['calculation_auc']['quantiles'],
        out_dir = os.path.join(config['out_dir'], "auc")
    resources:
        mem_mb=lambda wildcards, input: 3072 if sum([os.path.getsize(x) / pow(1024,2) for x in input]) < 100 else 8192
    log: os.path.join(config['local_log'], "calculation_auc", "{cov_dir}", "auc_{cell_type}.log")
    conda: os.path.join(workflow.basedir, config['envs_dir'], "python_basics.yaml")
    shell:
        """
        python {input.script} -bed {input.bed} \
          --cpm-threshold {params.cpm_threshold} \
          --one-plus {input.one_plus} --one-minus {input.one_minus} \
          --two-plus {input.two_plus} --two-minus {input.two_minus} \
          --out-file {output.auc_bed} \
          --quantiles {params.quants} \
          --log-file {log}
        """

rule auc_analysis:
    '''
    Perform AUC analysis and generate result plots
    '''
    input:
        auc_files = expand(os.path.join(config['out_dir'], 
          "auc", "{cov_dir}", "auc_{{cell_type}}.tsv"),
          cov_dir = ['coverage_3p', 'randomised_coverage']),
        script = os.path.join(config['script_dir'], "auc_analysis.R")
    output:
        shorter = os.path.join(config['out_dir'],
          "auc", "analysis_out", "{prefix}_TEs_shorter_{cell_type}.tsv"),
        longer = os.path.join(config['out_dir'],
          "auc", "analysis_out", "{prefix}_TEs_longer_{cell_type}.tsv")
    params:
        in_dir = os.path.join(config['out_dir'], "auc"),
        out_dir = lambda wildcards, output: os.path.dirname(output.shorter),
        cmprs = "{cell_type}",
        out_prefix = "{prefix}_",
        alpha = config['auc_analysis']['alpha'],
        cpm_threshold = config['auc_analysis']['cpm_threshold'],
        length_threshold = config['auc_analysis']['length_threshold'],
        iqr_threshold = config['auc_analysis']['iqr_threshold']
    conda: 
        os.path.join(workflow.basedir, config['envs_dir'], "r_plot.yaml")
    resources:
        mem_mb=2048
    shell:
        """
        Rscript {input.script} \
          --analysis-dir {params.in_dir} \
          --out-dir {params.out_dir} \
          --comparison {params.cmprs} \
          --cpm-threshold {params.cpm_threshold} \
          --length-threshold {params.length_threshold} \
          --iqr-threshold {params.iqr_threshold} \
          --ALPHA {params.alpha} \
          --out-prefix {params.out_prefix}
        """
