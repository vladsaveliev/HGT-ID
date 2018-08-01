
rule to_viral_sortname:
    input:   join(viral_mapping, 'partially_mapped_back_to_virus.bam')
    output:  join(viral_mapping, 'partially_mapped_back_to_virus.sort.bam')
    threads: config['threads']
    shell:   sortname_shell


rule to_viral_merge:
    input:   join(viral_mapping, 'partially_mapped_to_virus_mate_unmapped.bam'),
             join(viral_mapping, 'partially_mapped_to_virus_mate_mapped.bam')
    output:  join(viral_mapping, 'partially_mapped_back_to_virus.bam')
    shell:   merge_shell


rule to_viral_from_unmapped_sortname:
    input:   join(viral_mapping, 'unmapped_to_virus.bam')
    output:  join(viral_mapping, 'unmapped_to_virus.sort.bam')
    threads: config['threads']
    shell:   sortname_shell + ' && samtools index {output}'


rule to_viral_extract_mate_unmapped:
    input:  join(viral_mapping, 'partially_mapped_to_virus.bam')
    output: join(viral_mapping, 'partially_mapped_to_virus_mate_unmapped.bam')
    shell:  extract_mate_unmapped_shell


rule to_viral_extract_mate_mapped:
    input:  join(viral_mapping, 'partially_mapped_to_virus.bam')
    output: join(viral_mapping, 'partially_mapped_to_virus_mate_mapped.bam')
    shell:  extract_mate_mapped_shell


rule align_unmapped_to_viral:
    input:   fq1 = join(scoring_dir, 'read1.fq'),
             fq2 = join(scoring_dir, 'read2.fq'),
             sname = join(work_dir, 'sample_name')
    output:  join(scoring_dir, 'unmapped_to_virus.bam')
    threads: config['threads']
    params:  bwa_idx = config['virus_bwa']
    run:     get_bwa_run(input, output, threads, params)


rule align_human_to_viral:
    input:   fq1 = join(human_mapping_again, 'read1.fq'),
             fq2 = join(human_mapping_again, 'read2.fq'),
             sname = join(work_dir, 'sample_name')
    output:  join(viral_mapping, 'partially_mapped_to_virus.bam')
    threads: config['threads']
    params:  bwa_idx = config['virus_bwa']
    run:     get_bwa_run(input, output, threads, params)
