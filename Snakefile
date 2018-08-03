from os.path import join


shell.executable(os.environ.get('SHELL', "bash"))
shell.prefix("")


output_dir          = config.get('output', 'output')
output_file         = join(output_dir, 'output.txt')
work_dir            = join(output_dir, 'work')
logs_dir            = join(work_dir, 'logs')
autocode_dir        = join(work_dir, 'autocode')
human_mapping       = join(work_dir, 'human_mapping')
human_mapping_again = join(work_dir, 'human_mapping_again')
viral_mapping       = join(work_dir, 'viral_mapping')
scoring_dir         = join(work_dir, 'scoring')
candidates_dir      = join(work_dir, 'candidates')
find_is_dir         = join(work_dir, 'find_is')


rule all:
    input: join(find_is_dir, 'output.filt.gene.ann.txt')


include: 'snake_src/utils.smk'
include: 'snake_src/step1_human_firstpass.smk'
include: 'snake_src/step2_back_to_human.smk'
include: 'snake_src/step3_to_viral.smk'
include: 'snake_src/step4_initital_candidates.smk'
include: 'snake_src/step5_back_to_original_bam.smk'
include: 'snake_src/step6_find_is.smk'

