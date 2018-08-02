
rule initialcandidates_map_virushuman:
    input:   bam             = join(human_mapping_again, 'to_human.sort.bam'),
             script          = join(config['scripts_dir'], 'map.virus.human.pl')
    output:  candidates      = join(output_dir, 'HGT.candidates.txt'),
             virus_human_sam = join(output_dir, 'virus.human.sam')
    log:     join(logs_dir, '5.map.virus.human.log')
    shell:  'perl {input.script} {input.bam} {input.viral_bam} {output.candidates} {output.virus_human_sam} > {log} 2>&1'


rule initialcandidates_virushuman_tobam:
    input:   sam = join(output_dir, 'virus.human.sam')
             fa =  config['virus_human_fa']
    output:  join(output_dir, 'virus.human.bam')
    shell:  "cat {input.sam} | sed '/^$/d' | samtools view -bt {input.fa}.fai - > {output}"


rule initialcandidates_virushuman_sortname:
    input:   join(output_dir, 'virus.human.bam')
    output:  join(output_dir, 'virus.human.sort.bam')
    threads: config['threads']
    shell:   sortname_shell + ' && samtools index {output}'


rule initialcandidates_filter_with_seq_complexity:
    input:   candidates = join(output_dir, 'HGT.candidates.txt'),
             script     = join(config['scripts_dir'], 'filter.pl')
    output:  join(output_dir, 'HGT.candidates.filt.txt')
    params:  sequence_complexity = config['sequence_complexity']
    shell:   """
    cat {input.candidates} | 
    perl {input.script} | 
    awk 'NR>1' | 
    awk -v complex={params.sequence_complexity} '$NF>=complex' > {output}
    """


rule initialcandidates_filter_with_seq_complexity_to_bam:
    input:  candidates      = join(output_dir, 'HGT.candidates.filt.txt'),
            viral_human_bam = join(output_dir, 'virus.human.sort.bam')
    output: join(output_dir, 'virus.human.sort.filt.bam')
    log:    join(logs_dir, '6.FilterSamReads.log')
    shell:  """picard FilterSamReads \
               I={input.viral_human_bam} \
               RLF={input.candidates} \
               FILTER=includeReadList O={output} \
               WRITE_READS_FILES=false \
               VALIDATION_STRINGENCY=SILENT \
               CREATE_INDEX=TRUE \
               > {log} 2>&1
    """


rule initialcandidates_get_well_covered_regions:
    """ Filter chromosomes with the number of reads mapped or partially mapped of at least `virus_read_support`
    """
    input:  join(output_dir, 'virus.human.sort.filt.bam')
    output: join(output_dir, 'chromosomes2keep.txt')
    params: virus_read_support = config['virus_read_support']
    shell:  """
        samtools idxstats {input} 
        | awk -v depth={params.depth} '$NF+$(NF-1)>virus_read_support' 
        | cut -f1 > {output}
    """

