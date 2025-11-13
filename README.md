# Prod-Seq
This repository hosts the source code for the data processing and analysis model of Prod-seq. Please see below for the installation instructions, basic usage, and the documentation of the importable functions.

## Installation
1. Create a conda environment using the prodseq_analysis_environment.yml file (note: this step may take some time for conda to resolve the environment); conda will create an environment named “prodseq”
```
conda env create -f prodseq_analysis_environment.yml
```
**Note:** for users with experience with conda, the packages listed in the .yml file can also be installed manually (the packages are all commonly-used packages)

**Note:** if a different environment name is desired, the user can change the first line of the .yml file to “name: desired_environment_name”

2. Download _ProdSeq_Data_Analysis.py_, _PQSeq_Data_Analysis.py_ and _ProdSeqAnalysis_utils.py_ to the desired directory.
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

<br>

<br>

## Basic Usage

### ProdSeq_Data_Analysis.py

This is the main analysis script and the command-line executable.

Usage: 

```
<full_download_path>/ProdSeq_Data_Analysis.py [--barcodetsv BARCODETSV] sampletsv output_prefix
```

<br>

**<ins>Positional arguments:</ins>**

  * _sampletsv_  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Path to the sample specification tsv file

    This tsv file should contain three columns, where the first column dictates a sample’s name, and the second and third columns of the corresponding row specify the path to the sample’s read1 and read2 files (as fastq.gz files). 
      
  * _output_prefix_  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Output prefix

    This specifies the prefix of the output files (see “Output files” for more information).

<br>

**<ins>Optional arguments:</ins>**

  * _--barcodetsv BARCODETSV_  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Path to the barcoce list tsv file

    This optional argument allows the user to specify a custom set of protein barcodes that is different from the default Prod-seq set of 12 Ab-oligo and/or specify a different set of control Ab-oligos used as independent variables in the two-round GLM estimation.
    
    The first column should contain the target/protein name of the protein barcodes.
    
    The second column of the corresponding lines should contain the reverse complement of the barcode sequences for each protein target. For example, in the default Prod-seq Ab-oligo set, the Ab-oligo sequence of SUZ12 is $${\color{black}CCTTGAACCACTTCTCTA} {\color{green}AATCGACTCA} {\color{black}NNNNNNNNNNNNNNNgcttacaaccagactg}$$, so the second column of SUZ12 should enter the reverse complement of $${\color{green}AATCGACTCA}$$, which is $${\color{green}TGAGTCGATT}$$.

    The third column should be “control” for each control Ab-oligo, and empty for the other Ab-oligos.

    **Note:** it is required to have at least one control ab-oligo.

<br>

**<ins>Output files:</ins>**

All the output files are tab-delimited files where the first row is the column headers, and each of the remaining rows corresponds to one Prod-seq sample. The entry values of the output tsv files differ (detailed below).

**Note:** as the two-round GLM uses the non-present control ab-oligo(s) (default: HA-tag) in the prediction of the background levels, output files that store post-GLM outputs (_PredictedBackground_, _PPIEnrichmentLevels_) do not contain values related to the non-present control ab-oligos(s).

  * **_output_prefix.QCCnts.tsv_**
    
    For each sample, the read depth, number of detected Prod-seq product read pairs, UMI-deduplicated number of product read pairs, and the number of self-self byproduct read pairs for each protein barcode.

  * **_output_prefix.RawUMIPairs.tsv_**

    For each sample, the number of UMI pairs for each pairwise protein combination (excluding the self-self byproducts).

  * **_output_prefix.DepthNormedUMIPairs.tsv_**

    For each sample, the number of UMI pairs for each pairwise protein combination (excluding the self-self byproducts) after normalizing by the total number of UMI combinations in all the samples listed in this file.

  * **_output_prefix.PredictedBackground.tsv_**

    The level of background predicted by the two-round GLM for each PPI in each sample.

  * **_output_prefix.PPIEnrichment.tsv_**

    The enrichment value of each PPI in each sample calculated using the predicted background levels.


**Note:** the script runs a two-round GLM to estimate the PPI general background distribution, where the first round is only to exclude extreme positive outliers. In other words, it is expected that the first round of maximum likelihood estimation does not fully converge and it is normal behavior that the first round of GLM outputs the following warning messages:
```
RuntimeWarning: invalid value encountered in subtract
np.max(np.abs(fsim[0] - fsim[1:])) <= fatol):

message: Maximum number of iterations has been exceeded.
success: False
```

The second round of GLM is expected to terminate successfully; example terminal message:

```
message: Optimization terminated successfully.
success: True
```


<br>

<br>

### PQSeq_Data_Analysis.py

This is the command-line executable for read processing of protein quantification subroutine.

Usage: 

```
<full_download_path>/PQSeq_Data_Analysis.py [--barcodetsv BARCODETSV] sampletsv output_prefix
```

<br>

**<ins>Positional arguments:</ins>**

  * _sampletsv_  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Path to the sample specification tsv file

    This tsv file should contain two columns, where the first column dictates a sample’s name, and the second column specifies the path to the sample’s read1 fastq.gz files. 
      
  * _output_prefix_  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Output prefix

    This specifies the prefix of the output files (see “Output files” for more information).

<br>

**<ins>Optional arguments:</ins>**

The same as **ProdSeq_Data_Analysis.py**.

<br>

**<ins>Output files:</ins>**

All the output files are tab-delimited files where the first row is the column headers, and each of the remaining rows corresponds to one Prod-seq protein quantification sample. The entry values of the output tsv files differ (detailed below).

  * **_output_prefix.PQQCCnts.tsv_**
    
    For each sample, the read depth, number of detected product reads, UMI-deduplicated number of product reads.

  * **_output_prefix.PQRawUMIs.tsv_**

    For each sample, the number of UMIs for each protein target.


<br>

<br>

### ProdSeqAnalysis_utils.py

This file contains the data processing and analysis functions used in **ProdSeq_Data_Analysis.py** and **PQSeq_Data_Analysis.py**. In addition, this file provides useful file reading and plotting functions (documentation in the section below). To import functions from this script, include the following lines:

```
import sys
sys.path.insert(1, '<path_to_ProdSeqAnalysis_utils.py>')
from ProdSeqAnalysis_utils import *
```

**<ins>Documentation of additional functions:</ins>**

  * **_ReadProdTSVFile_**

    Read TSV output files into sample names (row names), column names, and values as lists.


  * **_GenerateProdQCPlots_**

    Generate QC plots for Prod-seq PPI sequencing result. 


  * **_GroupedPPIHeatmap_**

    Generate PPI enrichment heatmaps for Prod-seq data visualization.


  * **_CalcPQFromUMICnts_**

    Convert Prod-seq protein quantification (PQ) readout into (each sample separately) UMI proportions and UMI counts normalized by positive control protein barcodes.


  * **_PQHeatmap_**

    Generate PQ value heatmaps for PQ data visualization.





