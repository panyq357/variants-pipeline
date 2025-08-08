import gzip

def get_chromosome_names(fasta_path):
    open_func = gzip.open if fasta_path.endswith('.gz') else open
    chrom_names = []

    with open_func(fasta_path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                chrom_name = line[1:].split()[0]  # take first word after '>'
                chrom_names.append(chrom_name)

    return chrom_names


def get_bam_list(w):
    return config["variant_calling_jobs"][w.job_name]


rule bcftools_call_by_chr:
    input:
        bam_list = get_bam_list,
        bai_list = lambda w: [f"{bam}.bai" for bam in get_bam_list(w)],
        genome = config["genome"]
    output:
        temp("resources/bcftools_call_by_chr/{job_name}.{chr}.vcf.gz")
    log:
        "resources/bcftools_call_by_chr/{job_name}.{chr}.vcf.gz.log"
    priority:
        10
    resources:
        io = lambda w: len(get_bam_list(w)) * 10 if len(get_bam_list(w)) * 10 < 100 else 100
    shell:
        '''
        bcftools mpileup -O u -a AD,DP -r {wildcards.chr} -f {input.genome} {input.bam_list} 2> {log} \
        | bcftools call -v -m -O z -a GQ -o {output} 2>> {log}
        '''


rule bcftools_concat:
    input:
        chr_vcf_list = expand("resources/bcftools_call_by_chr/{{job_name}}.{chrom}.vcf.gz", chrom=get_chromosome_names(config["genome"])),
    output:
        "results/bcftools_concat/{job_name}.vcf.gz",
    wildcard_constraints:
        grp = "[^.]+"
    log:
        "results/bcftools_concat/{job_name}.vcf.gz.log",
    priority:
        100
    shell:
        '''
        bcftools concat {input.chr_vcf_list} 2> {log} \
        | bgzip -c > {output} 2>> {log}
        '''

