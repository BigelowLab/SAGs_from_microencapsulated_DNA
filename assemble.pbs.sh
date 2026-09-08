#PBS -N spc_github
#PBS -m bea -j oe -o status_of_assembly.pbs
#PBS -l ncpus=1,mem=50MB,walltime=120:00:00
module purge
module load singularity
module load nextflow
cd /mnt/scgc/scgc_temp/spc_github/
nextflow run /mnt/scgc/scgc_temp/spc_github/main.nf -profile charlie -resume \
  -name SAGs_from_microcapsules_$(date +%Y%m%d_%H%M%S) \
  --contam_ref_fasta /mnt/scgc_nfs/ref/assembly_protocol/GRCh38_AG665_mm10.fa --viral

