""" Creating regions to go back to original BAM file to extract reads for to find soft clipping...
"""
import math

rule _step5_back_to_original:
    input:   join(step5_dir, 'forcalling.bam')



rule step5_forcalling_sortname:
    input:   join(step5_dir, 'merged.cyto.bam')
    output:  join(step5_dir, 'forcalling.bam')
    threads: config['threads']
    shell:   sortpos_shell



rule step5_forcalling_cyto:
    input:  bam =  join(step5_dir, 'merged.bam'),
            cyto = config['human_cytoband']
    output:        join(step5_dir, 'merged.cyto.bam')
    shell: 'zcat {input.cyto} | '
           'grep acen | '
           'awk \'{{ if($4 ~ /^p/) {{ print $1"\\t"$2-50000"\\t"$3 }} else {{ print $1"\\t"$2"\\t"$3+50000 }} }}\' | '
           'intersectBed -abam {input.bam} -b stdin -v > {output}'


rule step5_forcalling_merge:
    input:  join(step4_dir, 'virus.human.sortpos.filt.bam'),
            join(step5_dir, 'virus.org.bam'),
            join(step5_dir, 'human.again.bam'),
            join(step5_dir, 'human.org.bam')
    output: join(step5_dir, 'merged.bam')
    shell:  merge_shell


# 256=primary and 1024=not_duplicate
subset_bam_shell = 'samtools view -b -F 256 -F 1024 -L {input.regions} {input.bam} > {output}'


rule step5_extract_nearby_reads_viral:
    input:  bam     = join(viral_mapping, 'human_mapping_again_to_virus.bam'),
            regions = join(step5_dir, 'regions.to_keep.bed')
    output:           join(step5_dir, 'virus.org.bam')
    shell:  subset_bam_shell


rule step5_extract_nearby_reads_humanagain:
    input:  bam     = join(human_mapping_again, 'human_back_from_first_pass.bam'),
            regions = join(step5_dir, 'regions.to_keep.bed')
    output:           join(step5_dir, 'human.again.bam')
    shell:  subset_bam_shell


rule step5_extract_nearby_reads:
    input:  bam     = config['bam'],
            regions = join(step5_dir, 'regions.to_keep.bed')
    output:           join(step5_dir, 'human.org.bam')
    shell:  subset_bam_shell


# rule chr_to_look:
#     input: fai = config['human_fa'] + '.fai'
#            bed = join(step4_dir, 'regions.to_keep.bed')
#     shell: """cat {input.bed} | cut -f1 | sort | uniq | '
#               grep -w -f <cut -f1 {input.fai}) |
#               tr '\\n' ' ' | sed \'s/\\s*\$\/\/g\'"""


rule step5_regions_in_chroms_to_keep:
    input:  chroms_to_keep = join(step4_dir, 'chromosomes_to_keep.txt'),
            regions        = join(step5_dir, 'regions.bed')
    output:                  join(step5_dir, 'regions.to_keep.bed')
    shell: 'cat {input.regions} | grep -w -f {input.chroms_to_keep} > {output}'


rule step5_create_regions:
    input:  bam      = join(step4_dir, 'virus.human.sortpos.filt.bam'),
            fai      = config['virus_human_fa'] + '.fai',
            read_len = join(work_dir,  'readlen.txt'),
            lib_size = join(work_dir,  'libsize.txt')
    output:            join(step5_dir, 'regions.bed')
    run:
        read_len = int(math.floor(float(open(input.read_len).read().strip())))
        lib_size = int(math.floor(float(open(input.lib_size).read().strip())))
        slop_len = read_len + lib_size
        shell("""
    bamToBed -i {input.bam} | 
    sortBed -i stdin | 
    mergeBed -i stdin | 
    slopBed -i stdin -g {input.fai} -l {slop_len} -r {slop_len} | 
    mergeBed -i stdin > {output}
    """)

