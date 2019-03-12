import glob

r1_dict = {}
r1_list = glob.glob('*_*_*_R1_*.fastq.gz')

for name in r1_list:
    isolate = name.split('_')[0]
    r1_dict[isolate] = (name, name.replace('_R1_', '_R2_'))

def lookup_isolate(wildcards):
   return r1_dict[wildcards.isolate]

print(r1_dict)

rule all:
    input:
        expand("{isolate}.fa", isolate=r1_dict.keys())

rule assemble_isolate:
    input: 
        lookup_isolate   

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
        "cp {isolate}.assembly.d/final.contigs.fa {output}"
