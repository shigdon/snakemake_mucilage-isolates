import glob, os

r1_dict = {}
r1_list = glob.glob('reads/BCW-*_R1.fastq.gz')

for name in r1_list:
    i = os.path.basename(name).split('_')[0]
    r1_dict[i] = (name, name.replace('_R1', '_R2'))

def lookup_isolate(wildcards):
   return r1_dict[wildcards.isolate]

print(r1_dict)

rule all:
    input:
        expand("{isolate}/sourmash/{isolate}.fa.sig", isolate=r1_dict.keys())

rule make_isolate_directory:
    input:
        lookup_isolate
    
    output:
        directory("{isolate}")
    
    shell:
        "mkdir {wildcards.isolate}"

rule trim_isolate:
    input:
        lookup_isolate

    output:
        directory("{isolate}/trimmomatic"),
        "{isolate}/trimmomatic/{isolate}.trim_R1.fq.gz",
        "{isolate}/trimmomatic/{isolate}.orphan_R1.fq.gz",
        "{isolate}/trimmomatic/{isolate}.trim_R2.fq.gz",
        "{isolate}/trimmomatic/{isolate}.orphan_R2.fq.gz"

    conda:
        "trimmomatic.yml"

    threads: 24

    shell:
        "trimmomatic PE \
        -threads {threads} \
        {input[0]} \
        {input[1]} \
        {output[1]} \
        {output[2]} \
        {output[3]} \
        {output[4]} \
        ILLUMINACLIP:TruSeq3-PE-2.fa:2:40:15:8:TRUE \
        LEADING:2 \
        TRAILING:2 \
        SLIDINGWINDOW:4:15 \
        MINLEN:50"

rule assemble_isolate:
    input: 
        "{isolate}/trimmomatic/{isolate}.trim_R1.fq.gz",
        "{isolate}/trimmomatic/{isolate}.trim_R2.fq.gz"

    output:
        directory("{isolate}/megahit")
    
    conda:
        "megahit.yml"
    
    threads: 8

    shell:
        "megahit -1 {input[0]} -2 {input[1]} --out-dir {output} --out-prefix {isolate} -t {threads}"

rule compute_sig:
    input:
        "{isolate}/megahit/{isolate}.fa"

    output:
        directory("{isolate}/sourmash"),
        "{isolate}/sourmash/{isolate}.fa.sig"

    conda:
        "sourmash.yml"

    shell:
        "sourmash compute -k 31 --scaled 2000 {input} -o {output}"
