""" Creating regions to go back to original BAM file to extract reads for to find soft clipping...
"""

rule step5:
    input:   join(output_dir, 'forcalling.bam')



rule back_to_original_sortname:
    input:   join(candidates_dir, 'merged.cyto.bam')
    output:  join(output_dir, 'forcalling.bam')
    threads: config['threads']
    shell:   sortpos_shell



rule back_to_original_cyto:
    input:  bam =  join(candidates_dir, 'merged.bam'),
            cyto = config['human_cytoband']
    output:        join(candidates_dir, 'merged.cyto.bam')
    shell: 'zcat {input.cyto} | '
           'grep acen | '
           'awk \'{{ if($4 ~ /^p/) {{ print $1"\\t"$2-50000"\\t"$3 }} else {{ print $1"\\t"$2"\\t"$3+50000 }} }}\' | '
           'intersectBed -abam {input.bam} -b stdin -v > {output}'


rule back_to_original_merge:
    input:  join(candidates_dir, 'virus.human.sortpos.filt.bam'),
            join(candidates_dir, 'virus.org.bam'),
            join(candidates_dir, 'human.again.bam'),
            join(candidates_dir, 'virus.org.bam')
    output: join(candidates_dir, 'merged.bam')
    shell:  merge_shell
"merge the BAM file to find the integration point...\n";


# 256=primary and 1024=not_duplicate
subset_bam_shell = 'samtools view -b -F 256 -F 1024 -L {input.regions} {input.bam} > {output}'


rule back_to_original_extract_nearby_reads_viral:
    input:  bam = join(viral_mapping, 'human_mapping_again_to_virus.bam'),
            regions = join(candidates_dir, 'regions.to_keep.bed')
    output: join(candidates_dir, 'virus.org.bam')
    shell:  subset_bam_shell


rule back_to_original_extract_nearby_reads_humanagain:
    input:  bam = join(human_mapping_again, 'human_back_from_first_pass.bam'),
            regions = join(candidates_dir, 'regions.to_keep.bed')
    output: join(candidates_dir, 'human.again.bam')
    shell:  subset_bam_shell


rule back_to_original_extract_nearby_reads:
    input:  bam = config['bam'],
            regions = join(candidates_dir, 'regions.to_keep.bed')
    output: join(candidates_dir, 'human.org.bam')
    shell:  subset_bam_shell


# rule chr_to_look:
#     input: fai = config['human_fa'] + '.fai'
#            bed = join(candidates_dir, 'regions.to_keep.bed')
#     shell: """cat {input.bed} | cut -f1 | sort | uniq | '
#               grep -w -f <cut -f1 {input.fai}) |
#               tr '\\n' ' ' | sed \'s/\\s*\$\/\/g\'"""


rule back_to_original_regions_in_chroms_to_keep:
    input:  chroms_to_keep = join(candidates_dir, 'chromosomes_to_keep.txt'),
            regions = join(candidates_dir, 'regions.bed')
    output: join(candidates_dir, 'regions.to_keep.bed')
    shell: 'cat {input.regions} | grep -w -f {input.chroms_to_keep} > {output}'


rule back_to_original_create_regions:
    input:  bam = join(candidates_dir, 'virus.human.sortpos.filt.bam'),
            fai = config['virus_human_fa'] + '.fai',
            read_len = join(candidates_dir, 'readlen.txt'),
            lib_size = join(candidates_dir, 'libsize.txt')
    output: join(candidates_dir, 'regions.bed')
    run:
        read_len = int(open(input.read_len).read().strip())
        lib_size = int(open(input.lib_size).read().strip())
        slop_len = read_len + lib_size
        shell("""
    bamToBed -i {input.bam} | 
    sortBed -i stdin | 
    mergeBed -i stdin | 
    slopBed -i stdin -g {input.fai} -l {slop_len} -r {slop_len} | 
    mergeBed -i stdin > {output}
    """)


rule lib_size:
    """ Take mapped (-f2) reads with quality above 20 (-q20),
        Remove alignments with gaps,
        Remove non-paired reads,
        Take reads with insert size below 1000,
        Based on first 10k of them,
        Calculate the average insert absolute size.
    """
    input:   bam = config['bam'],
             script = join(config['snake_src_dir'], 'lib_size_from_bam.sh')
    output:  join(candidates_dir, 'libsize.txt')
    shell:  'bash {input.script} {input.bam} > {output}'


rule read_len:
    input:   config['bam']
    output:  join(candidates_dir, 'readlen.txt')
    shell:  "samtools view {input} | head -1 | awk '{{print length($10)}}' > {output}"

