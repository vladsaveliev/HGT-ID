rule _step2_back_to_human:
    input: read1 = join(human_mapping_again, 'read1.fq'),
           read2 = join(human_mapping_again, 'read2.fq')


rule step2_back_to_human_to_fastq:
    input:   join(human_mapping_again, 'human.sortname.bam')
    output:  read1 = join(human_mapping_again, 'read1.fq'),
             read2 = join(human_mapping_again, 'read2.fq')
    shell:   to_fastq_shell


rule step2_back_to_human_sortname:
    input:   join(human_mapping_again, 'human.bam')
    output:  join(human_mapping_again, 'human.sortname.bam')
    threads: config['threads']
    shell:   sortname_shell


rule step2_back_to_human_merge:
    input:   join(human_mapping_again, 'mate_unmapped.bam'),
             join(human_mapping_again, 'mate_mapped.bam')
    output:  join(human_mapping_again, 'human.bam')
    shell:   merge_shell


rule step2_back_to_human_extract_mate_unmapped:
    input:   join(human_mapping_again, 'human_back_from_first_pass.bam')
    output:  join(human_mapping_again, 'mate_unmapped.bam')
    shell:   extract_mate_unmapped_shell


rule step2_back_to_human_extract_mate_mapped:
    input:   join(human_mapping_again, 'human_back_from_first_pass.bam')
    output:  join(human_mapping_again, 'mate_mapped.bam')
    shell:   extract_mate_mapped_shell


rule step2_back_to_human_align:
    input:   fq1 = join(human_mapping, 'read1.fq'),
             fq2 = join(human_mapping, 'read2.fq'),
             sname = join(work_dir, 'sample_name')
    output:  join(human_mapping_again, 'human_back_from_first_pass.bam')
    threads: config['threads']
    params:  bwa_idx = config['human_bwa']
    run:     shell(get_bwa_run(input))
