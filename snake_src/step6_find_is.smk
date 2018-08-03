""" Creating regions to go back to original BAM file to extract reads for to find soft clipping...
"""

rule step6_find_is:
    input: join(find_is_dir, 'output.filt.gene.ann.txt')


rule annotate_with_genes:
    input:  file     = join(find_is_dir, 'output.filt.txt'),
            gene_bed = join(find_is_dir, 'gene.bed')
    output:            join(find_is_dir, 'output.filt.gene.txt')
    shell: 'cat {input.file} | awk \'NR>1 {{ print $1"\\t"$2"\\t"$2}\' | '
           'closestBed -t first -a stdin -d -b {input.gene_bed} | '
           'awk \'{{ print $1"\\t"$2"\\t"$(NF-1)"\\t"$NF }}\' > {output}'


rule annotate_with_genes_2:
    input:  file      = join(find_is_dir, 'output.filt.txt'),
            gene_file = join(find_is_dir, 'output.filt.gene.txt'),
            script    = join(config['scripts_dir'], 'map.hgt.annot.pl')
    output:             join(find_is_dir, 'output.filt.gene.ann.txt')
    shell: 'perl {input.script} {input.gene_file} {input.file} | grep -v NA > {output}'


rule create_gene_bed_to_annotate:
    input:   ref_flat = config['ref_flat'],
             script   = join(config['scripts_dir'], 'uniqgene.pl')
    output:  join(find_is_dir, 'gene.bed')
    shell:  'zcat {input.ref_flat} | '
            'grep -v random | '
            'grep -v chrUn | '
            'grep -v hap | '
            'awk \'{{print $3"\\t"$5"\\t"$6"\\t"$1 }}\' | '
            'sortBed -i stdin | '
            'perl {input.script} > {output}'


rule filter_centromere:
    # filter the file to remove the reads very close to 5KB centromere using the cytoband file
    input:  bed = join(find_is_dir, 'output.bed'),
            cyto = config['human_cytoband']
    output: join(find_is_dir, 'output.filt.txt')
    shell:  'zcat {cyto} | '
            'grep acen | '
            'awk \'{{ if($4 ~ /^p/) {{ print $1"\\t"$2-50000"\\t"$3 }} else {{ print $1"\\t"$2"\\t"$3+50000 }} }}\' | '
            'intersectBed -a {input.bed} -b stdin -v -header | '
            'cut -f4- > {output}'


rule filter_centromere_to_bed:
    input:  join(find_is_dir, 'output.txt')
    output: join(find_is_dir, 'output.bed')
    shell: 'cat {input} | awk \'{{ if (NR==1) {{print "#chr\\tstart\\tend\\t"$0 }} else {{ print $1"\\t"$2"\\t"$2"\\t"$0 }} }}\' > {output}'


rule find_is:
    input:  bam             = join(output_dir, 'forcalling.bam'),
            script          = join(config['scripts_dir'], 'integration.pl'),
            human_viral_bwa = config['virus_human_bwa'],
            human_bwa       = config['human_bwa'],
            lib_size        = join(candidates_dir, 'libsize.txt')
    output: join(find_is_dir, 'output.txt')
    params: min_pr   = config['min_pr'],
            min_soft = config['min_soft']
    run:
        lib_size = int(open(input.lib_size).read().strip())
        shell('perl {input.script} -b {input.bam} -f {input.human_viral_bwa} {output} -d {lib_size} '
              '-m {params.min_pr} -l {params.min_soft} -r {input.human_bwa} -t -v')


