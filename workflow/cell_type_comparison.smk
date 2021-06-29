'''
Module containing rules to perform AUC analysis for given cell type comparisons
'''

rule copy_coverage_3p_end:
    'Copy files to proper directory structure'
    input:
        os.path.join(config['out_dir'], 
            "coverage_3p", "{cell_type}_{origin}_{strand}.corr.bedgraph")
    output:
        os.path.join(config['out_dir'], 
            "{comparison}", "coverage_3p", "{cell_type}_{origin}_{strand}.corr.bedgraph")
    resources:
        mem_mb = 1024
    shell: "cp {input} {output}"

rule copy_coverage_3p_end_ct:
    'Copy files to proper directory structure'
    input:
        os.path.join(config['out_dir'], 
            "coverage_3p", "{cell_type}_{strand}.corr.bedgraph")
    output:
        os.path.join(config['out_dir'], 
            "{comparison}", "coverage_3p", "{cell_type}_{strand}.corr.bedgraph")
    resources:
        mem_mb = 1024
    shell: "cp {input} {output}"

rule label_randomisation:
    '''
    Randomise the 3p end coverage between two samples
    memory: based on previous runs, made regression for MB usage of input file sizes.
    output: bit hacky by outputting all cell types, but in fact only files in params will be generated, the other files are empty.
    '''
    input:
        bed = os.path.join(config['mapping_dir'], "terminal_exons.bed"),
        one_plus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "coverage_3p", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct1'].values[0] + "_plus.corr.bedgraph"),
        one_minus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "coverage_3p", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct1'].values[0] + "_minus.corr.bedgraph"),
        two_plus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "coverage_3p", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct2'].values[0] + "_plus.corr.bedgraph"),
        two_minus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "coverage_3p", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct2'].values[0] + "_minus.corr.bedgraph")
    output:
        expand(os.path.join(config['out_dir'], "{{comparison}}", 
            "randomised_coverage", "{cell_type}_{origin}_{strand}.corr.bedgraph"),
            cell_type = CELL_TYPES,
            origin = SAMPLE_ORIGIN,
            strand = ['plus', 'minus'])
    params:
        one_plus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "randomised_coverage", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct1'].values[0] + "_plus.corr.bedgraph"),
        one_minus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "randomised_coverage", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct1'].values[0] + "_minus.corr.bedgraph"),
        two_plus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "randomised_coverage", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct2'].values[0] + "_plus.corr.bedgraph"),
        two_minus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "randomised_coverage", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct2'].values[0] + "_minus.corr.bedgraph"),
        script = os.path.join(config['script_dir'], "label_randomisation.py")
    threads: 8
    resources:
        mem_mb = lambda wildcards, input: int(sum([os.path.getsize(x) / pow(1024,2) for x in input]) * 20 + 1000)
    log: os.path.join(config['local_log'], "label_randomisation", "log_{comparison}.log")
    conda: os.path.join(workflow.basedir, config['envs_dir'], "python_basics.yaml")
    shell:
        """
        python {params.script} -bed {input.bed} \
          -s1p {input.one_plus} -s1m {input.one_minus} \
          -s2p {input.two_plus} -s2m {input.two_minus} \
          -o1p {params.one_plus} -o1m {params.one_minus} \
          -o2p {params.two_plus} -o2m {params.two_minus} \
          --cores {threads} --log-file {log} && \
        touch -a {output}
        """

rule label_randomisation_ct:
    '''
    Randomise the 3p end coverage between two samples
    memory: based on previous runs, made regression for MB usage of input file sizes.
    output: bit hacky by outputting all cell types, but in fact only files in params will be generated, the other files are empty.
    '''
    input:
        bed = os.path.join(config['mapping_dir'], "terminal_exons.bed"),
        one_plus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "coverage_3p", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct1'].values[0] + "_plus.corr.bedgraph"),
        one_minus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "coverage_3p", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct1'].values[0] + "_minus.corr.bedgraph"),
        two_plus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "coverage_3p", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct2'].values[0] + "_plus.corr.bedgraph"),
        two_minus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "coverage_3p", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct2'].values[0] + "_minus.corr.bedgraph")
    output:
        expand(os.path.join(config['out_dir'], "{{comparison}}", 
            "randomised_coverage", "{cell_type}_{strand}.corr.bedgraph"),
            cell_type = CELL_TYPES,
            strand = ['plus', 'minus'])
    params:
        one_plus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "randomised_coverage", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct1'].values[0] + "_plus.corr.bedgraph"),
        one_minus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "randomised_coverage", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct1'].values[0] + "_minus.corr.bedgraph"),
        two_plus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "randomised_coverage", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct2'].values[0] + "_plus.corr.bedgraph"),
        two_minus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, "randomised_coverage", 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct2'].values[0] + "_minus.corr.bedgraph"),
        script = os.path.join(config['script_dir'], "label_randomisation.py")
    threads: 8
    resources:
        mem_mb = lambda wildcards, input: int(sum([os.path.getsize(x) / pow(1024,2) for x in input]) * 20 + 1000)
    log: os.path.join(config['local_log'], "label_randomisation", "log_{comparison}.log")
    conda: os.path.join(workflow.basedir, config['envs_dir'], "python_basics.yaml")
    shell:
        """
        python {params.script} -bed {input.bed} \
          -s1p {input.one_plus} -s1m {input.one_minus} \
          -s2p {input.two_plus} -s2m {input.two_minus} \
          -o1p {params.one_plus} -o1m {params.one_minus} \
          -o2p {params.two_plus} -o2m {params.two_minus} \
          --cores {threads} --log-file {log} && \
        touch -a {output}
        """

# TODO: reduce number of rule executions by reverting to cell type specific result for coverage_3p
# e.g. compute cumulative sum for each cell type, then copy to comparison directory
rule cumulative_coverage_3p_end:
    '''
    Compute reverse cumulative distribution in BED intervals from 3p end coverage.
    threads: if file < 80MB, use 1 threads, else 2.
    memory: based on previous runs, estimate MB usage from input file size.
    '''
    input:
        coverage = os.path.join(config['out_dir'], 
            "{comparison}", "{cov_dir}", "{cell_type}_{origin}_{strand}.corr.bedgraph"),
        bed = os.path.join(config['mapping_dir'], "terminal_exons.bed")
    output:
        cumsum = os.path.join(config['out_dir'], 
            "{comparison}", "{cov_dir}", "{cell_type}_{origin}_cumsum_{strand}.bedgraph")
    params:
        strand = "{strand}",
        script = os.path.join(config['script_dir'], "cumulative_coverage.py")
    threads: 
        lambda wildcards, input: 1 if (os.path.getsize(input.coverage) / 
            pow(1024,2)) < 50 else 3
    resources:
        mem_mb = lambda wildcards, input: 4096 if (os.path.getsize(input.coverage) / 
            pow(1024,2)) < 50 else 8192
    log: os.path.join(config['local_log'], "cumulative_coverage_3p_end", "{comparison}", "{cov_dir}", "{cell_type}_{origin}_{strand}.log")
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


rule cumulative_coverage_3p_end_ct:
    '''
    Compute reverse cumulative distribution in BED intervals from 3p end coverage.
    threads: if file < 80MB, use 1 threads, else 2.
    memory: based on previous runs, estimate MB usage from input file size.
    '''
    input:
        coverage = os.path.join(config['out_dir'], 
            "{comparison}", "{cov_dir}", "{cell_type}_{strand}.corr.bedgraph"),
        bed = os.path.join(config['mapping_dir'], "terminal_exons.bed")
    output:
        cumsum = os.path.join(config['out_dir'], 
            "{comparison}", "{cov_dir}", "{cell_type}_cumsum_{strand}.bedgraph")
    params:
        strand = "{strand}",
        script = os.path.join(config['script_dir'], "cumulative_coverage.py")
    threads: 
        lambda wildcards, input: 1 if (os.path.getsize(input.coverage) / 
            pow(1024,2)) < 50 else 3
    resources:
        mem_mb = lambda wildcards, input: 4096 if (os.path.getsize(input.coverage) / 
            pow(1024,2)) < 50 else 8192
    log: os.path.join(config['local_log'], "cumulative_coverage_3p_end", "{comparison}", "{cov_dir}", "{cell_type}_{strand}.log")
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
        one_plus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, wildcards.cov_dir, 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct1'].values[0] + "_cumsum_plus.bedgraph"),
        one_minus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, wildcards.cov_dir, 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct1'].values[0] + "_cumsum_minus.bedgraph"),
        two_plus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, wildcards.cov_dir, 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct2'].values[0] + "_cumsum_plus.bedgraph"),
        two_minus = lambda wildcards: os.path.join(config['out_dir'], 
            wildcards.comparison, wildcards.cov_dir, 
            cmprs.loc[cmprs['comparison'] == wildcards.comparison, 'ct2'].values[0] + "_cumsum_minus.bedgraph"),
        script = os.path.join(config['script_dir'], "calculation_auc.py")
    output:
        auc_bed = os.path.join(config['out_dir'], 
          "auc_comparisons", "{cov_dir}", "auc_{comparison}.tsv")
    params:
        cpm_threshold = config['calculation_auc']['cpm_threshold'],
        quants = config['calculation_auc']['quantiles'],
        out_dir = os.path.join(config['out_dir'], "auc_comparisons")
    resources:
        mem_mb=lambda wildcards, input: 3072 if sum([os.path.getsize(x) / pow(1024,2) for x in input]) < 100 else 8192
    log: os.path.join(config['local_log'], "calculation_auc", "{cov_dir}", "auc_{comparison}.log")
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
        auc_bed = expand(os.path.join(config['out_dir'], 
          "auc_comparisons", "{cov_dir}", "auc_{{comparison}}.tsv"),
          cov_dir = ['coverage_3p', 'randomised_coverage']),
        script = os.path.join(config['script_dir'], "auc_analysis.R")
    output:
        shorter = os.path.join(config['out_dir'],
          "auc_comparisons", "analysis_out", "{prefix}_TEs_shorter_{comparison}.tsv"),
        longer = os.path.join(config['out_dir'],
          "auc_comparisons", "analysis_out", "{prefix}_TEs_longer_{comparison}.tsv")
    params:
        in_dir = os.path.join(config['out_dir'], "auc_comparisons"),
        out_dir = lambda wildcards, output: os.path.dirname(output.shorter),
        cmprs = "{comparison}",
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
