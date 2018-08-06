""" Creating regions to go back to original BAM file to extract reads for to find soft clipping...
"""
from ngs_utils.file_utils import safe_mkdir


rule _step6_find_is:
    input: join(step6_dir, 'output.filt.gene.ann.viral_names.txt'),
           join(step6_dir, 'coverage.txt')


rule step6_viral_coverage:
    input:   file       = join(step6_dir,     'output.filt.gene.ann.viral_names.txt'),
             virus_fai  = config['virus_fa'] + '.fai',
             cov_script = join(config['scripts_dir'], 'coverage_split.pl'),
             bam        = join(viral_mapping, 'virus.sortname.bam')
    output:  coverage   = join(step6_dir,     'coverage.txt')
    run:
        lines_cnt = len(open(input.file).readlines())
        if lines_cnt > 1:
            shell('cat {input.virus_fai} | ' 
                  'cut -f1,2 | '
                  'perl {input.cov_script} | '
                  'coverageBed -abam {input.bam} -b stdin > {output}')
        else:
            shell('touch {output}')


# #### get the viral coverage
# my $coverage=$OUTPUT ."/.coverage";
# mkdir $coverage;
# my $count = `wc -l < $OUTNAME`;
# sleep 10;
#
# if ( $count > 1 )	{
# 	open FH, ">$autocode/coverage.sh" or die "can't open the script to write\n";
# 	print FH "for i in `cat $OUTNAME | awk 'NR>1' | cut -f8 | sort | uniq`\n";
# 	print FH "do\n";
# 	print FH "cat $VIRUS_database.fai | grep -w \$i | cut -f1,2 | perl $SCRIPT_DIR/coverage_split.pl | coverageBed -abam $viral_mapping/virus.fix.sort.bam -b stdin > $coverage/\$i.coverage.out \&\n";
# 	print FH "pid=\$!\n";
# 	print FH "cat $VIRUS_database.fai | grep -w \$i | cut -f1,2 | perl $SCRIPT_DIR/coverage_split.pl | coverageBed -abam $scoring/virus.fromUnmapped.sort.bam -b stdin > $coverage/\$i.proper.coverage.out \&\n";
# 	print FH "pid1=\$!\n";
# 	print FH "wait \$pid \$pid1\n";
# 	print FH "done\n";
# 	close FH;
# 	$command=join("","chmod 777 ",$autocode,"/coverage.sh");
# 	submit($command);
# 	$command=join ("","sh ",$autocode,"/coverage.sh");
# 	submit($command);
# 	sleep 10;
# 	submit($command);
# }
# sleep 10;


rule step6_filter_with_viral_fai:
    input:  virus_fai = config['virus_fa'] + '.fai',
            file      = join(step6_dir, 'output.filt.gene.ann.txt')
    output:             join(step6_dir, 'output.filt.gene.ann.viral_names.txt')
    run:
        virus_names = {l.strip.split()[0] for l in open(input.virus_fai)}
        with open(input.file) as inp, open(output, 'w') as out:
            for i, l in enumerate(inp):
               if i != 0:
                   fs = l.split('\t')
                   if fs[2] in virus_db or fs[7] in virus_db:
                       out.write(l)


rule step6_annotate_with_genes:
    input:  file     = join(step6_dir, 'output.filt.txt'),
            gene_bed = join(step6_dir, 'gene.bed')
    output:            join(step6_dir, 'output.filt.gene.txt')
    shell: 'cat {input.file} | awk \'NR>1 {{ print $1"\\t"$2"\\t"$2 }}\' | '
           'closestBed -t first -a stdin -d -b {input.gene_bed} | '
           'awk \'{{ print $1"\\t"$2"\\t"$(NF-1)"\\t"$NF }}\' > {output}'


rule step6_annotate_with_genes_2:
    input:  file      = join(step6_dir, 'output.filt.txt'),
            gene_file = join(step6_dir, 'output.filt.gene.txt'),
            script    = join(config['scripts_dir'], 'map.hgt.annot.pl')
    output:             join(step6_dir, 'output.filt.gene.ann.txt')
    shell: 'perl {input.script} {input.gene_file} {input.file} | grep -v NA > {output}'


rule step6_create_gene_bed_to_annotate:
    input:   ref_flat = config['ref_flat'],
             script   = join(config['scripts_dir'], 'uniqgene.pl')
    output:             join(step6_dir, 'gene.bed')
    shell:  'zcat {input.ref_flat} | '
            'grep -v random | '
            'grep -v chrUn | '
            'grep -v hap | '
            'awk \'{{print $3"\\t"$5"\\t"$6"\\t"$1 }}\' | '
            'sortBed -i stdin | '
            'perl {input.script} > {output}'


rule step6_filter_centromere:
    # filter the file to remove the reads very close to 5KB centromere using the cytoband file
    input:  bed  = join(step6_dir, 'output.bed'),
            cyto = config['human_cytoband']
    output:        join(step6_dir, 'output.filt.txt')
    shell:  'zcat {input.cyto} | '
            'grep acen | '
            'awk \'{{ if($4 ~ /^p/) {{ print $1"\\t"$2-50000"\\t"$3 }} else {{ print $1"\\t"$2"\\t"$3+50000 }} }}\' | '
            'intersectBed -a {input.bed} -b stdin -v -header | '
            'cut -f4- > {output}'


rule step6_filter_centromere_to_bed:
    input:  join(step6_dir, 'output.txt')
    output: join(step6_dir, 'output.bed')
    shell: 'cat {input} | awk \'{{ if (NR==1) {{print "#chr\\tstart\\tend\\t"$0 }} else {{ print $1"\\t"$2"\\t"$2"\\t"$0 }} }}\' > {output}'


rule step6_find_is_raw:
    input:  bam             = join(step5_dir, 'forcalling.bam'),
            script          = join(config['scripts_dir'], 'integration.pl'),
            lib_size        = join(work_dir,  'libsize.txt')
    output:                   join(step6_dir, 'output.txt')
    params: min_pr          = config['min_pr'],
            min_soft        = config['min_soft'],
            human_bwa       = config['human_bwa'],
            human_viral_fa  = config['virus_human_fa'],
            tmp_dir         = safe_mkdir(join(step6_dir, 'integration_tmp'))
    log:    join(logs_dir, 'integration.pl.log')
    run:
        lib_size = int(open(input.lib_size).read().strip())
        shell('perl {input.script} -b {input.bam} -f {params.human_viral_fa} -o {output} '
              '-d {lib_size} -m {params.min_pr} -l {params.min_soft} -r {params.human_bwa} '
              '-t {params.tmp_dir} -v 2>{log} && [[ $(wc -l <{output}) -gt 1 ]]')


