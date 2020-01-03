shell.executable('/bin/bash')

localrules: all, copy_assembly, compute_sig

import glob, os


r1_dict = {}
r1_list = glob.glob('reads/BCW-*_R1.fastq.gz')

for name in r1_list:
    i = os.path.basename(name).split('_')[0]
    r1_dict[i] = (name, name.replace('_R1', '_R2'))

def lookup_isolate(wildcards):
   return r1_dict[wildcards.isolate]

print(r1_dict)


outdir = config['outdir'] + 'data/'
logdir = config['outdir'] + 'log/'


rule all:
    input:
        expand(outdir + "{isolate}/sourmash/{isolate}-k31.sig", isolate=r1_dict.keys()),
        expand(outdir + "{isolate}/prokka", isolate=r1_dict.keys())

rule trim_isolate:
    input:
        lookup_isolate

    output:
        outdir + "{isolate}/{isolate}.trim_R1.fq.gz",
        outdir + "{isolate}/{isolate}.orphan_R1.fq.gz",
        outdir + "{isolate}/{isolate}.trim_R2.fq.gz",
        outdir + "{isolate}/{isolate}.orphan_R2.fq.gz"

    log:
        logdir + "trimmomatic/{isolate}.log"

    conda:
        "trimmomatic.yml"

    threads: 24

    shell:
        "trimmomatic PE \
        -threads {threads} \
        {input[0]} \
        {input[1]} \
        {output} \
        ILLUMINACLIP:TruSeq3-PE-2.fa:2:40:15:8:TRUE \
        LEADING:2 \
        TRAILING:2 \
        SLIDINGWINDOW:4:15 \
        MINLEN:50 \
        > /dev/null 2> {log}"

rule assemble_isolate:
    input: 
        outdir + "{isolate}/{isolate}.trim_R1.fq.gz",
        outdir + "{isolate}/{isolate}.trim_R2.fq.gz"

    output:
        outdir + "{isolate}/megahit"
    
    conda:
        "megahit.yml"
    
    threads: 24

    shell:
        "megahit \
        -1 {input[0]} \
        -2 {input[1]} \
        --out-dir {output} \
        -t {threads}"

rule copy_assembly:
    input:
        outdir + "{isolate}/megahit"

    output:
        outdir + "{isolate}/{isolate}.contigs.fa"

    shell:
        "cp {input}/final.contigs.fa {output}"

rule prokka_annotation:
    input:
        outdir + "{isolate}/{isolate}.contigs.fa"

    output:
        outdir + "{isolate}/prokka"

    conda:
        "prokka.yml"

    threads: 24

    shell:
        "prokka \
            --outdir {output} \
            --prefix {wildcards.isolate} \
            --locustag {wildcards.isolate} \
            --cpus 24 \
            {input}"

rule compute_sig:
    input:
        outdir + "{isolate}/{isolate}.contigs.fa"

    output:
        outdir + "{isolate}/sourmash/{isolate}-k31.sig"

    conda:
        "sourmash.yml"

    shell:
        "sourmash compute -k 31 --scaled 2000 {input} -o {output}"
