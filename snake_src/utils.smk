sortname_shell = 'samtools sort -n -@ {threads} {input} -Obam > {output}'
sortpos_shell  = 'samtools sort -@ {threads} {input} -Obam > {output} && samtools index {output}'
merge_shell    = 'samtools merge -u -n {output} {input}'
to_fastq_shell = 'samtools fastq {input} -1 {output.read1} -2 {output.read2}'

extract_mate_unmapped_shell = 'samtools view -u -b -f 8 -F 260 {input} > {output}' # mate unmapped & primary
extract_mate_mapped_shell   = 'samtools view -u -b -f 4 -F 264 {input} > {output}' # only mate mapped & primary


def get_bwa_run(input, output, threads, params):
    sname = open(input.sname).read().strip()
    rg = f'@RG\\tID:{sname}\\tSM:{sname}'
    shell('bwa mem -t {threads} -M -R "{rg}" '
          '{params.bwa_idx} {input.fq1} {input.fq2} | '
          'samtools view -u -bS - > {output}')


rule sample_name:
    """ Read the sample name from the RG tag """
    input:   bam = config['bam'],
             script = join(config['snake_src_dir'], 'sample_name_from_bam.sh')
    output:  join(work_dir, 'sample_name')
    shell:  'bash {input.script} {input.bam} > {output}'
