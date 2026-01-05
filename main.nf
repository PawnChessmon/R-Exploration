nextflow.enable.dsl = 2

workflow {
  DESEQ2_ANALYSIS( file(params.input_dir) )
}

process DESEQ2_ANALYSIS {
  tag "deseq2"
  publishDir params.outdir, mode: 'copy'

  input:
  path input_dir

  output:
  path "output_step2a"
  path "output_step2b", optional: true
  path "output_step2c", optional: true
  path "output_step2d", optional: true

  script:
  """
  Rscript ${projectDir}/scripts/run_deseq2.R \
    --input_dir ${input_dir} \
    --output_base . \
    --input_type ${params.input_type} \
    --de_tables ${params.de_tables} \
    --pcol ${params.pcol} \
    --pthresh ${params.pthresh} \
    --lfc ${params.lfc}
  """
}
