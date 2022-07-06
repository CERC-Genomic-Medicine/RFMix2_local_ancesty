## About
This is a Nextflow pipeline built around the RFMix2 software for the local genetic ancestry estimation using phased data and reference panel.

## Prerequisites

1. Download and install the RFMix2 from `https://github.com/slowkoni/rfmix`.
2. Clone this repository:
    ```
    git clone https://github.com/CERC-Genomic-Medicine/RFMix2_local_ancesty.git
    cd RFMix2_local_ancesty
    ```
2.  Download genetic maps `zip` file for build **GRCh38** from `https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/` into `Genetic_maps` directory, unzip and re-format:
    ```
    mkdir Genetic_maps
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
5. Load your Python 3 virtual environment and make sure that `Pandas` is installed.
6. Run nextflow e.g.:
    ```
    nextflow run Pipeline.nf
    ```
