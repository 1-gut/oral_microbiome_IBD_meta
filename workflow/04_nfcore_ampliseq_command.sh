# nf-core ampliseq pipeline command
nextflow run nf-core/ampliseq \
    -r 2.14.0 \
    -profile eddie \
    --input_folder "SRA" \
    --multiple_sequencing_runs \
    --input "./oral_samplesheet.tsv" \
    --metadata "./oral_metadata.tsv" \
    --skip_cutadapt \
    --trunc_qmin 30 \
    --trunc_rmin 0.8 \
    --ignore_empty_input_files \
    --skip_dada_addspecies \
    --skip_ancom \
    --skip_qiime \
    --metadata_category "IBD_group_name" \
    --picrust \
    --outdir "./ampliseq_out"