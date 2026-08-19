#PBS -N spc_github
#PBS -M ggavelis@bigelow.org
#PBS -m bea -j oe -o status_of_assembly.pbs
#PBS -l ncpus=1,mem=50MB,walltime=120:00:00
module purge
module load singularity
module load nextflow
cd /mnt/scgc/scgc_temp/spc_github/
NXF_VER=25.04.6 nextflow run /mnt/scgc/scgc_temp/spc_github/main.nf -profile charlie -resume \
  -name SAGs_from_microcapsules_$(date +%Y%m%d_%H%M%S) \
  --contam_ref_fasta /mnt/scgc_nfs/ref/assembly_protocol/GRCh38_AG665_mm10.fa --viral \
 --target_IDs "1_00230416,1_00130812,1_03162203,1_03022317,1_01171107,1_05160404,1_00222022,1_00180013,1_01002215,1_06000010,1_01160416,1_07072103,1_03052311,1_00150716,1_06001123,1_02151502,1_00181702,1_03102109,1_01040211,1_03100501,1_08081014,1_01161907,1_02220315,1_02131510,1_01181315,1_10220915,1_04161020,1_01122206,1_10210321,1_01040313,1_07191113,1_06010902,1_07060107,1_01051222,1_09090907,1_10100807,1_07201808,1_01120318"

