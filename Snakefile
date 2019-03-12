import glob, os

r1_dict = {}
r1_list = glob.glob('reads/*_*_*_R1_*.fastq.gz')

for name in r1_list:
    i = os.path.basename(name).split('_')[0]
    r1_dict[i] = (name, name.replace('_R1_', '_R2_'))

def lookup_isolate(wildcards):
   return r1_dict[wildcards.isolate]

print(r1_dict)

rule all:
    input:
        expand("{isolate}.fa.sig", isolate=r1_dict.keys()),
        "smash_k31-cmp.csv"

rule trim_isolate:
    input:
        lookup_isolate

    output:
        "{isolate}.trim_R1.fq.gz",
        "{isolate}.orphan_R1.fq.gz",
        "{isolate}.trim_R2.fq.gz",
        "{isolate}.orphan_R2.fq.gz"

    conda:
        "trimmomatic.yml"

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
        MINLEN:50"

rule assemble_isolate:
    input: 
        "{isolate}.trim_R1.fq.gz",
        "{isolate}.trim_R2.fq.gz"

    output:
        directory("{isolate}.assembly.d")
    
    conda:
        "megahit.yml"
    
    shell:
        "megahit -1 {input[0]} -2 {input[1]} -o {output}"

rule copy_assembly:
    input:
        "{isolate}.assembly.d"

    output:
        "{isolate}.fa"

    shell:
        "cp {wildcards.isolate}.assembly.d/final.contigs.fa {output}"

rule compute_sig:
    input:
        "{isolate}.fa"

    output:
        "{isolate}.fa.sig"

    conda:
        "sourmash.yml"

    shell:
        "sourmash compute -k 31 {input} -o {output}"

rule compare_sigs:
    input:
        expand("{isolate}.fa.sig", isolate=r1_dict.keys())
    output:
        "smash.k31.cmp",
        "smash.k31.cmp.labels.txt",
        "smash_k31-cmp.csv" 
    conda:
        "sourmash.yml"
    shell: 
        "sourmash compare -k 31 {input} -o {output[0]} --csv {output[2]}"
