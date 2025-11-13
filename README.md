# Prod-Seq
This repository hosts the source code for the data processing and analysis model of Prod-seq. Please see below for the installation instructions, basic usage, and the documentation of the importable functions.

## Installation
1. Create a conda environment using the prodseq_analysis_environment.yml file (note: this step may take some time for conda to resolve the environment); conda will create an environment named “prodseq”
```
conda env create -f prodseq_analysis_environment.yml
```
**Note:** for users with experience with conda, the packages listed in the .yml file can also be installed manually (the packages are all commonly-used packages)

**Note:** if a different environment name is desired, the user can change the first line of the .yml file to “name: desired_environment_name”

2. Download _ProdSeq_Data_Analysis.py_, _PQ_Data_Analysis.py_ and _ProdSeqAnalysis_utils.py_ to the desired directory.
3. Make the python scripts executable
   ```
   conda activate <environment_name>
   chmod +x <full_download_path>/ProdSeq_Data_Analysis.py
   chmod +x <full_download_path>/PQSeq_Data_Analysis.py
   chmod +x <full_download_path>/ProdSeqAnalysis_utils.py
   ```

4. To check that the scripts are ready to execute, run the following command. Correct setup should output the helper message specifying the usage of the script and arguments and options:
   ```
   <full_download_path>/ProdSeq_Data_Analysis.py [--barcodetsv BARCODETSV] sampletsv output_prefix
   ```


## Basic Usage

### ProdSeq_Data_Analysis.py

This is the main analysis script and the command-line executable.

Usage: 

```
<full_download_path>/ProdSeq_Data_Analysis.py [--barcodetsv BARCODETSV] sampletsv output_prefix
```


**Positional arguments:**

  * _sampletsv_  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Path to the sample specification tsv file

    This tsv file should contain three columns, where the first column dictates a sample’s name, and the second and third columns of the corresponding row specify the path to the sample’s read1 and read2 files (as fastq.gz files). 
      
  * _output_prefix_  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Output prefix

    This specifies the prefix of the output files (see “Output files” for more information).

**Optional arguments:**

  * _--barcodetsv BARCODETSV_  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Path to the barcoce list tsv file

    This optional argument allows the user to specify a custom set of protein barcodes that is different from the default Prod-seq set of 12 Ab-oligo and/or specify a different set of control Ab-oligos used as independent variables in the two-round GLM estimation. 
The first column should contain the target/protein name of the protein barcodes.
The second column of the corresponding lines should contain the reverse complement of the barcode sequences for each protein target. For example, in the default Prod-seq Ab-oligo set, the Ab-oligo sequence of SUZ12 is:

    CCTTGAACCACTTCTCTAAATCGACTCANNNNNNNNNNNNNNNgcttacaaccagactg

    so the second column of SUZ12 should enter the reverse complement of AATCGACTCA, which is TGAGTCGATT.

    The third column should be “control” for each control Ab-oligo, and empty for the other Ab-oligos. Note: it is required to have at least one control Ab-oligo.


