# Prerequisites

1. Download and install the RFMix2
2.  Download genetic maps `zip` file for build **GRCh38** from `https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/` into `Genetic_maps` directory, unzip and re-format:
```
cd Genetic_maps
wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
unzip plink.GRCh38.map.zip
for i in {1..22}; do awk '{print "chr"$1,$4,$3}' OFS='\t' plink.chr${i}.GRCh38.map > chr${i}.GRCh38.map; done
cd ..
```
3. Edit the `nextflow.config` file accordingly.
4. Load `bcftools` and `nextflow` modules:
```
module load bcftools
module load nextflow
```
5. Load Python3 virtual environment with Pandas
6. Run nextflow e.g.:
  ```
  nextflow run Pipeline.nf
  ```
