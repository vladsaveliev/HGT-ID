rule _step1_human_firstpass:
    input: read1 = join(human_mapping, 'read1.fq'),
           read2 = join(human_mapping, 'read2.fq')


rule step1_human_firstpass_unmapped_to_fastq:
    input:   join(scoring_dir, 'both_unmapped.bam')
    output:  read1 = join(scoring_dir, 'read1.fq'),
             read2 = join(scoring_dir, 'read2.fq')
    shell:   to_fastq_shell


rule step1_human_firstpass_to_fastq:
    input:   join(human_mapping, 'human.sortname.bam')
    output:  read1 = join(human_mapping, 'read1.fq'),
             read2 = join(human_mapping, 'read2.fq')
    shell:   to_fastq_shell


rule step1_human_firstpass_sortname:
    input:   join(human_mapping, 'human.bam')
    output:  join(human_mapping, 'human.sortname.bam')
    threads: config['threads']
    shell:   sortname_shell


rule step1_human_firstpass_merge:
    input:   join(human_mapping, 'mate_unmapped.bam'),
             join(human_mapping, 'mate_mapped.bam'),
             join(human_mapping, 'soft.bam')
    output:  join(human_mapping, 'human.bam')
    shell:   merge_shell


rule step1_human_firstpass_soft_bam:
    input:   bam = join(human_mapping, 'softclip.bam'),
             fai = config['human_fa'] + '.fai',
             split_reads_pl = join(config['scripts_dir'], 'splitReads.pl')
    output:  join(human_mapping, 'soft.bam')
    log:     join(logs_dir, '0.soft_bam.log')
    shell:   """
    samtools view {input.bam} | 
    perl {input.split_reads_pl} | 
    samtools view -bt {input.fai} - > {output} 2>{log}
    """


rule step1_human_firstpass_softclip_bam:
    input:     bam = config['bam'],
               ids = join(human_mapping, 'IDS')
    output:    join(human_mapping, 'softclip.bam')
    log:       join(logs_dir, '0.FilterSamReads.log')
    resources: mem_mb = 6000
    shell:     """picard FilterSamReads \
    I={input.bam} \
    WRITE_READS_FILES=false \
    RLF={input.ids} \
    FILTER=includeReadList \
    O={output} \
    SORT_ORDER=queryname &> {log}
    """


rule step1_human_firstpass_all_read_ids:
    input:  join(human_mapping, 'soft_clipped_read1.fq'),
            join(human_mapping, 'soft_clipped_read2.fq')
    output: join(human_mapping, 'IDS')
    shell: 'cat {input} | cut -f1 | sort -T $PWD | uniq > {output}'


rule step1_human_firstpass_softclip_reads:
    """ Take mate mapped & mapped & paired & first in pair,
        Take with soft clipped in CIGAR,
        Skip with insertions in CIGAR,
        Skip with deletions in CIGAR,
        Take reads with a S series in CIGAR longer than len/2
    """
    input:  bam = config['bam'],
            filt_pl = join(config['scripts_dir'], 'filter_soft_clipped.pl')
    output: join(human_mapping, 'soft_clipped_read{num_in_pair}.fq')
    params: flag = lambda wc: {'1': 64, '2': 128}[wc.num_in_pair]
    shell:  """
    samtools view -F8 -f2 -F4 -f{params.flag} {input.bam} |
    awk '$6 ~ /S/' |
    awk '$6 !~ /I/' |
    awk '$6 !~ /D/' |
    perl -an {input.filt_pl} > {output}
    """


rule step1_human_firstpass_mate_unmapped:
    input:  config['bam']
    output: join(human_mapping, 'mate_unmapped.bam')
    shell:  extract_mate_unmapped_shell


rule step1_human_firstpass_mate_mapped:
    input:  config['bam']
    output: join(human_mapping, 'mate_mapped.bam')
    shell:  extract_mate_mapped_shell


rule step1_human_firstpass_both_unmapped:
    # unmapped & mate unmapped
    input:  config['bam']
    output: join(scoring_dir, 'both_unmapped.bam')
    shell: 'samtools view -u -b -f 12 {input} > {output}'

