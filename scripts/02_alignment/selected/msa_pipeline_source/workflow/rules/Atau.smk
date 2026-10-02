# convert maf to fasta and run conservation tools 

from snakemake import available_cpu_count
from psutil import virtual_memory

available_mem_gb = lambda: '%dG' % (virtual_memory().available >> 30)
containerized: "docker://apscheben/msa_pipeline:latest"

SPECIES = config['species']
TREE = config['speciesTree']

rule call_conservation:
    input:
       'results/tree/gerp.complete',
       'results/tree/phylop.complete'

rule cds_fold:
    input:
       refFasta='data/{refname}.fa'.format(refname=config['refName'])
    output:
       txt=temp('results/tree/sitefold.txt'),
    params:
       gff='data/{ingff}'.format(ingff=config['refGFF'])
    conda:
      'workflow_envs'
    log:
      'logs/cds_fold_log.txt'
    threads: 1
    script:
      '../scripts/cds_fold.py'

rule sitefold2bed:
    input:
       'results/tree/sitefold.txt'
    output:
       'results/tree/sitefold.bed'
    log:
      'logs/sitefold2bed_log.txt'
    threads: 1
    script:
        '../scripts/sitefold2bed.py'

rule maf2fold:
    input:
       maf='results/roast/roast.maf',
       bed='results/tree/sitefold.bed'
    output:
       maf=temp('results/roast/roast.sitefold.maf')
    conda:
      'workflow_envs'
    log:
      'logs/maf2fold_log.txt'
    threads: 1
    shell:
      """
      source ~/miniforge3/bin/activate workflow_envs  && \
      mafsInRegion {input.bed} {output.maf} {input.maf} &> {log}
      """

rule clean_maf:
    input:
       maf='results/roast/roast.sitefold.maf'
    output:
       maf=temp('results/roast/roast.sitefold.clean.maf')
    threads: 1
    script:
      '../scripts/clean_mafsInRegion_output.py'

rule maf2fa:
    input:
       'results/roast/roast.maf' if not config['refGFF'] else 'results/roast/roast.sitefold.clean.maf'
    output:
       temp('results/tree/roast.fa')
    params:
       expand(['{species}','{refname}'],species=SPECIES,refname=config['refName'])
    log:
      'logs/maf2fa_log.txt'
    benchmark:
      'benchmark/maf2fa_bm.txt'
    conda:
      'workflow_envs'
    threads: 1
    script:
      '../scripts/maf2mfa.py'

rule filtN:
    input:
       'results/tree/roast.fa'
    output:
       temp('results/tree/roast_noN.fa')
    log:
      'logs/filtn_log.txt'
    benchmark:
      'benchmark/filtn_bm.txt'
    threads: 1
    shell:
      """
      sed "/^>/! s/N/-/g"  {input} > {output} 2> {log}
      """

rule filtcols:
    input:
       'results/tree/roast_noN.fa'
    output:
       'results/tree/roast_mincol.fa'
    log:
      'logs/filtcols_log.txt'
    benchmark:
      'benchmark/filtcols_bm.txt'
    conda:
      'workflow_envs'
    threads: 1
    shell:
      """
      source ~/miniforge3/bin/activate workflow_envs
      trimal -in {input} -out {output} -fasta -gt 0.9 > {log} && [ -f {output} ] || printf "ERROR: Insufficient sites with >0.9 sample coverage in results/tree/roast_noN.fa.\nCheck coverage in multiple alignment results/roast/roast.maf.\nNote that some highly diverged genomes may require tuning/changing of the aligner to achieve adequate coverage.\n" && printf "ERROR: Insufficient sites with >0.9 sample coverage in results/tree/roast_noN.fa.\nCheck coverage in multiple alignment results/roast/roast.maf.\nNote that some highly diverged genomes may require tuning/changing of the aligner to achieve adequate coverage.\n" >> {log}
      """

rule trimcols:
    input:
       'results/tree/roast_mincol.fa'
    output:
       'results/tree/roast_mincol_maxlen.fa'
    params:
       maxSites=config["maxNeutralSites"]
    log:
      'logs/filtcols_log.txt'
    benchmark:
      'benchmark/filtcols_bm.txt'
    conda:
      'workflow_envs'
    script:
      '../scripts/subsample_sites_in_fasta_alingment.py'

rule raxml:
    input:
       'results/tree/roast_mincol_maxlen.fa'
    output:
       'results/tree/raxml/RAxML_bestTree.roast_mincol_maxlen'
    params:
       workDir='results/tree',
       prefix='roast_mincol_maxlen'
    log:
      '../../logs/raxml_log.txt'
    benchmark:
      'benchmark/raxml_bm.txt'
    conda:
      'workflow_envs'
    threads: 20
    shell:
      """
      source ~/miniforge3/bin/activate workflow_envs  && \
      cd {params.workDir} && \
      raxmlHPC-PTHREADS -T {threads} -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s {params.prefix}.fa -n {params.prefix} &> {log} && \
      mv RAxML_* raxml/
      """

rule phylofit:
    input:
       tree='results/tree/raxml/RAxML_bestTree.roast_mincol_maxlen',
       alignment='results/tree/roast_mincol_maxlen.fa'
    output:
       'results/tree/neutral.tre'
    params:
       prefix='results/tree/neutral'
    log:
      'logs/phylofit_log.txt'
    benchmark:
      'benchmark/phylofit_bm.txt'
    conda:
      'workflow_envs'
    threads: 1
    shell:
      """
      source ~/miniforge3/bin/activate workflow_envs  && \
      phyloFit --out-root {params.prefix} --tree {input.tree} --subst-mod REV --msa-format FASTA {input.alignment} &> {log} && \
      tree_doctor -t {params.prefix}.mod > {output} 2> {log}
      """ 

checkpoint rule maf_split:
    input:
       maf='results/roast/roast.maf',
       tree='results/tree/neutral.tre'
    output:
       dir=directory('results/conservation_raw')
    log:
      'logs/mafsplit_log.txt'
    benchmark:
      'benchmark/mafsplit_bm.txt'
    conda:
      'workflow_envs'
    threads: 1
    shell:
      """
      source ~/miniforge3/bin/activate workflow_envs  && \
      mkdir {output.dir} && \
      mafSplit -byTarget -useFullSequenceName splits.bed {output.dir}/ {input.maf} &> {log}
      """ 

rule gerp:
    input:
       maf='results/conservation_raw/{chr}.maf'
    output:
       temp=('results/conservation_raw/{chr}.maf.rates')
    params:
       refname=config['refName'],
       neutral='results/tree/neutral.tre'
    log:
      'logs/gerp_{chr}_log.txt'
    conda:
      'workflow_envs'
    threads: 1
    shell:
      """
      source ~/miniforge3/bin/activate workflow_envs  && \
      gerpcol -t {params.neutral} -f {input.maf} -e {params.refname} -j -z -x ".rates" &> {log}
      """

rule phylop:
    input:
       maf='results/conservation_raw/{chr}.maf'
    output:
       temp=('results/conservation_raw/{chr}.phylop.wig')
    params:
       neutral='results/tree/neutral.mod'
    log:
      'logs/phylop_{chr}_log.txt'
    conda:
      '../envs/phast.yaml'
    threads: 1
    shell:
      """
      source ~/miniforge3/bin/activate workflow_envs  && \
      phyloP --mode CONACC --method LRT --wig-scores {params.neutral} {input.maf} > {output} 2> {log}
      """

rule wig2bed:
    input:
       'results/conservation_raw/{chr}.phylop.wig'
    output:
       temp=('results/conservation_raw/{chr}.phylop.bed')
    log:
      'logs/wig2bed_{chr}_log.txt'
    threads: 1
    script:
      '../scripts/phylopwig2bed.py'

rule add_phylop_header:
    input:
       'results/conservation_raw/{chr}.phylop.bed'
    output:
       'results/phylop/{chr}.phylop.final.bed'
    threads: 1
    shell:
      """
      cat <(printf 'chr\tstart\tend\tphyloP_-log_pvalue\n') {input} > {output}
      """

rule gerp2bed:
    input:
       'results/conservation_raw/{chr}.maf.rates'
    output:
       temp('results/conservation_raw/{chr}.maf.rates.bed')
    threads: 1
    script:
       '../scripts/gerp_rates_add_coord.py'

rule maf2cov:
    input:
       'results/conservation_raw/{chr}.maf'
    output:
       temp('results/conservation_raw/{chr}.maf.cov.bed')
    conda:
      'workflow_envs'
    threads: 1
    script:
      '../scripts/maf2bedcov.py'

rule add_cov:
    input:
       cov='results/conservation_raw/{chr}.maf.cov.bed',
       rates='results/conservation_raw/{chr}.maf.rates.bed'
    output:
       temp('results/conservation_raw/{chr}.rates.cov.bed')
    conda:
      '../envs/bedtools.yaml'
    threads: 1
    shell:
      """
      source ~/miniforge3/bin/activate workflow_envs  && \
      bedtools intersect -loj -a {input.cov} -b {input.rates} | \
      awk -v OFS='\t' '{{print $5,$6,$7,$8,$9,$4}}' > {output} 
      """

rule add_gerp_header:
    input:
      'results/conservation_raw/{chr}.rates.cov.bed' 
    output:
       'results/gerp/{chr}.gerp.final.bed'
    conda:
      '../envs/bedtools.yaml'
    threads: 1
    shell:
      """
      source ~/miniforge3/bin/activate workflow_envs  && \
      cat <(printf "chr\tstart\tend\tGERP_ExpSubst\tGERP_RejSubstScore\tTaxaAligned\n") {input} > {output}
      """

def gerp_aggregate(wildcards):
    checkpoint_output = checkpoints.maf_split.get(**wildcards).output[0]
    return expand("results/gerp/{chr}.gerp.final.bed",
           chr=glob_wildcards(os.path.join(checkpoint_output, "{chr}.maf")).chr)

def phylop_aggregate(wildcards):
    checkpoint_output = checkpoints.maf_split.get(**wildcards).output[0]
    return expand("results/phylop/{chr}.phylop.final.bed",
           chr=glob_wildcards(os.path.join(checkpoint_output, "{chr}.maf")).chr)

rule gerp_clean:
    input:
       gerp_aggregate
    output:
       temp('results/tree/gerp.complete')
    conda:
      'workflow_envs'
    threads: 1
    shell:
      """
      source ~/miniforge3/bin/activate workflow_envs 
      touch {output}
      """

rule phylop_clean:
    input:
       phylop_aggregate
    output:
       temp('results/tree/phylop.complete')
    conda:
      'workflow_envs'
    threads: 1
    shell:
      """
      source ~/miniforge3/bin/activate workflow_envs 
      touch {output}
      """

