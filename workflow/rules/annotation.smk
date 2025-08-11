rule filter_vcf:
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.filter-{filter}.vcf.gz"
    params:
        filter_expression = lambda w: config["bcftools"]["filters"][w.filter]
    log:
        "{prefix}.filter-{filter}.vcf.gz.log"
    wildcard_constraints:
        filter = "[^.]+"
    shell:
        '''
        bcftools view -e {params.filter_expression:q} {input} 2> {log} \
        | bgzip > {output} 2>> {log}
        '''


rule tabix:
    input:
        "{prefix}.gz"
    output:
        "{prefix}.gz.tbi"
    shell:
        "tabix {input}"


rule sort_gtf:
    input:
        gtf = config["gtf"]
    output:
        sorted_gtf = "resources/sorted.gtf.gz"
    shell:
        '''
        (zcat < {input.gtf} || cat < {input.gtf}) \
        | grep -v '^#' \
        | sort -k1,1 -k4,4n -k5,5n -t$'\\t' \
        | bgzip -c > {output.sorted_gtf}
        tabix -p gff {output.sorted_gtf}
        '''


rule vep_anno:
    input:
        vcf = "{prefix}.vcf.gz",
        genome = config["genome"],
        gtf = "resources/sorted.gtf.gz"
    output:
        vcf = "{prefix}.vep_anno.vcf.gz",
        summary = "{prefix}.vep_anno.html"
    log:
        "{prefix}.vep_anno.log"
    threads:
        4
    params:
        fields = config["vep"]["fields"]
    shell:
        '''
        bcftools view -G {input.vcf} \
        | bcftools annotate -x INFO \
        | vep \
            --fork {threads} \
            --gtf {input.gtf} \
            --fasta {input.genome} \
            --format vcf \
            --output_file STDOUT \
            --vcf \
            --force_overwrite \
            --stats_file {output.summary} \
            --per_gene \
            --fields "{params.fields}" \
            2> {log} \
        | bgzip -c > {output.vcf} 2>> {log}
        '''

