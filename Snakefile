from os.path import join


output_dir          = config.get('output', 'output')
output_file         = join(output_dir, 'output.txt')
work_dir            = join(output_dir, 'work')
logs_dir            = join(work_dir, 'logs')
autocode_dir        = join(work_dir, 'autocode')
human_mapping       = join(work_dir, 'human_mapping')
human_mapping_again = join(work_dir, 'human_mapping_again')
viral_mapping       = join(work_dir, 'viral_mapping')
scoring_dir         = join(work_dir, 'scoring')


rule all:
    input:  join(viral_mapping, 'partially_mapped_back_to_virus.sort.bam')


include: 'snake_src/utils.smk'
include: 'snake_src/step1_human_firstpass.smk'
include: 'snake_src/step2_back_to_human.smk'
include: 'snake_src/step3_to_viral.smk'

