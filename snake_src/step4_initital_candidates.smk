rule _step4_initialcandidates:
    input:  join(step4_dir, 'virus.human.sortpos.filt.bam'),
            join(step4_dir, 'chromosomes_to_keep.txt')


rule step4_initialcandidates_get_well_covered_regions:
    """ Filter chromosomes with the number of reads mapped or partially 
        mapped of at least `virus_read_support`
    """
    input:  join(step4_dir, 'virus.human.sortpos.filt.bam')
    output: join(step4_dir, 'chromosomes_to_keep.txt')
    params: virus_read_support = config['virus_read_support']
    shell:  'samtools idxstats {input} | '
            'awk -v depth={params.virus_read_support} \'$NF+$(NF-1)>virus_read_support\' | '
            'cut -f1 > {output}'


rule step4_initialcandidates_filter_with_seq_complexity_to_bam:
    input:  candidates      = join(step4_dir, 'HGT.candidates.filt.txt'),
            viral_human_bam = join(step4_dir, 'virus.human.sortpos.bam')
    output: join(step4_dir, 'virus.human.sortpos.filt.bam')
    log:    join(logs_dir, '6.FilterSamReads.log')
    # shell:  'picard FilterSamReads '
    #         'I={input.viral_human_bam} '
    #         'READ_LIST_FILE={input.candidates} '
    #         'FILTER=includeReadList O={output} '
    #         'WRITE_READS_FILES=false '
    #         'VALIDATION_STRINGENCY=SILENT '
    #         'CREATE_INDEX=TRUE '
    #         '> {log}'
    shell: 'samtools view -h {input.viral_human_bam} | ' \
           'grep -vf <(cut -f1 {input.candidates}) | ' \
           'samtools view -bS - -o {output}'


rule step4_initialcandidates_filter_with_seq_complexity:
    input:   candidates = join(step4_dir, 'HGT.candidates.txt'),
             script     = join(config['scripts_dir'], 'filter.pl')
    output:  join(step4_dir, 'HGT.candidates.filt.txt')
    params:  sequence_complexity = config['sequence_complexity']
    shell:   'cat {input.candidates} | '
             'perl {input.script} | '
             'awk \'NR>1\' | '
             'awk -v complex={params.sequence_complexity} \'$NF>=complex\' > {output}'
integration.pl

rule step4_initialcandidates_virushuman_sortpos:
    input:   join(step4_dir, 'virus.human.bam')
    output:  join(step4_dir, 'virus.human.sortpos.bam')
    threads: config['threads']
    shell:   sortpos_shell


rule step4_initialcandidates_virushuman_tobam:
    input:   sam = join(step4_dir, 'virus.human.sam'),
             fai =  config['virus_human_fa'] + '.fai'
    output:  join(step4_dir, 'virus.human.bam')
    shell:  "cat {input.sam} | sed '/^$/d' | samtools view -bt {input.fai} - > {output}"


rule step4_initialcandidates_map_virushuman:
    input:   bam             = join(human_mapping_again,   'human.sortname.bam'),
             viral_bam       = join(viral_mapping,         'virus.sortname.bam'),
             script          = join(config['scripts_dir'], 'map.virus.human.pl')
    output:  candidates      = join(step4_dir,             'HGT.candidates.txt'),
             virus_human_sam = join(step4_dir,             'virus.human.sam')
    log:     join(logs_dir, '5.map.virus.human.log')
    shell:  'perl {input.script} {input.bam} {input.viral_bam} {output.candidates} {output.virus_human_sam} > {log}'


