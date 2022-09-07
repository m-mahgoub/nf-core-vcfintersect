# nextflow run main.nf  -profile test,docker -resume --outdir test
# nextflow run main.nf  -profile test,conda -resume --outdir test
rm -rf null
rm -rf work
rm -rf .nextflow
rm -rf test
rm -rf .nextflow.log*