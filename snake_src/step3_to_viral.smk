rule _step3:
    input:   join(viral_mapping, 'virus.sortname.bam')


rule step3_to_viral_sortname:
    input:   join(viral_mapping, 'virus.bam')
    output:  join(viral_mapping, 'virus.sortname.bam')
    threads: config['threads']
    shell:   sortname_shell


rule step3_to_viral_merge:
    input:   join(viral_mapping, 'human_mapping_again_to_virus_mate_unmapped.bam'),
             join(viral_mapping, 'human_mapping_again_to_virus_mate_mapped.bam')
    output:  join(viral_mapping, 'virus.bam')
    shell:   merge_shell


rule step3_to_viral_from_unmapped_sortname:
    input:   join(viral_mapping, 'unmapped_to_virus.bam')
    output:  join(viral_mapping, 'unmapped_to_virus.sortname.bam')
    threads: config['threads']
    shell:   sortname_shell


rule step3_to_viral_extract_mate_unmapped:
    input:  join(viral_mapping, 'human_mapping_again_to_virus.bam')
    output: join(viral_mapping, 'human_mapping_again_to_virus_mate_unmapped.bam')
    shell:  extract_mate_unmapped_shell


rule step3_to_viral_extract_mate_mapped:
    input:  join(viral_mapping, 'human_mapping_again_to_virus.bam')
    output: join(viral_mapping, 'human_mapping_again_to_virus_mate_mapped.bam')
    shell:  extract_mate_mapped_shell


rule step3_align_unmapped_to_viral:
    input:   fq1 = join(scoring_dir, 'read1.fq'),
             fq2 = join(scoring_dir, 'read2.fq'),
             sname = join(work_dir, 'sample_name')
    output:  join(scoring_dir, 'unmapped_to_virus.bam')
    threads: config['threads']
    params:  bwa_idx = config['virus_bwa']
    run:     shell(get_bwa_run(input))


rule step3_align_human_to_viral:
    input:   fq1 = join(human_mapping_again, 'read1.fq'),
             fq2 = join(human_mapping_again, 'read2.fq'),
             sname = join(work_dir, 'sample_name')
    output:  join(viral_mapping, 'human_mapping_again_to_virus.bam')
    threads: config['threads']
    params:  bwa_idx = config['virus_bwa']
    run:     shell(get_bwa_run(input))
