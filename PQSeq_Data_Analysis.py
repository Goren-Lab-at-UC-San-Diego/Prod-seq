#!/usr/bin/env python
# coding: utf-8

# In[2]:


#!/usr/bin/python3
import argparse
from ProdSeqAnalysis_utils import *
import sys


# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument("sampletsv", help = "Path to the sample specification tsv file")
parser.add_argument("output_prefix", help = "Output prefix")
parser.add_argument("--barcodetsv", default = "", help = "Path to the barcoce list tsv file")


# In[ ]:


args = parser.parse_args()
sample_tsv = args.sampletsv
out_prefix = args.output_prefix
bc_tsv = args.barcodetsv


# In[ ]:


if (bc_tsv == ""):
    H4_bc = "TGTATCAGTT"
    H3K4me3_bc = "GTAGTGGCAT"
    H3K27ac_bc = "GTTATTAGGC"
    H3K27me3_bc = "TAACATGCGG"
    EZH2_bc = "TGGCTAATGT"
    EGFR_bc = "TGACCTTATG"
    HAtag_bc = "GATTGTCCGC"
    H3K27M_bc = "CACGATTGTT"
    AEBP2_bc = "TTGCATGGTA"
    EED_bc = "TTGAGTAACC"
    SUZ12_bc = "TGAGTCGATT"
    MED12_bc = "CTATGTTGGT"
    CycC_bc = "GAGGATAAGT"

    bc_pool = [
        H3K4me3_bc, 
        H3K27ac_bc, 
        H3K27me3_bc, EZH2_bc, EED_bc, SUZ12_bc, AEBP2_bc, 
        MED12_bc, CycC_bc, 
        H3K27M_bc,
        EGFR_bc, HAtag_bc]
    bc_nms = [
        "H3K4me3", 
        "H3K27ac", 
        "H3K27me3", "EZH2", "EED", "SUZ12", "AEBP2", 
        "MED12", "CycC", 
        "H3K27M",
        "EGFR", "HA-Tag"
    ]

    control_bc_idxs = [11]

else:
    bc_stream = (open(bc_tsv, "r")).readlines()
    bc_lines = [(line.strip()).split("\t") for line in bc_stream if line.strip()]
    bc_nms = [x[0] for x in bc_lines]
    bc_pool = [x[1] for x in bc_lines]

    control_bc_idxs = []
    for idx in range(len(bc_lines)):
        line_split = bc_lines[idx]
        if (len(line_split) >= 3 and line_split[2] == "control"):
            control_bc_idxs.append(idx)

    if (len(control_bc_idxs) == 0):
        print("Error: control barcode not specified. Exiting...")
        sys.exit()
        


# In[ ]:


sample_spec_stream = (open(sample_tsv, "r")).readlines()
sample_spec_lines = [(line.strip()).split("\t") for line in sample_spec_stream if line.strip()]
sample_nms = [x[0] for x in sample_spec_lines]
sample_R1s = [x[1] for x in sample_spec_lines]


# In[ ]:


# UMI counting
prod_midseq = "TAGAGAAG"
prod_arm_umi_len = 15
PQ_umi_cnts = []
PQ_QC_cnts = []
for idx in range(len(sample_spec_lines)):
    curr_umi_cnts, curr_QC_cnts = CountPQUMIs(sample_R1s[idx], bc_pool, prod_midseq, prod_arm_umi_len)
    PQ_umi_cnts.append(curr_umi_cnts)
    PQ_QC_cnts.append(curr_QC_cnts)


# In[ ]:


# Output results to files
QC_out = open(out_prefix + ".PQQCCnts.tsv", "w")
raw_umi_out = open(out_prefix + ".PQRawUMIs.tsv", "w")

QC_out.write("#Sample\tRead_Count\tPQ_Product_Read_Count\tUMI_Count_Total\n")
raw_umi_out.write("#Sample\t")

for i in range(len(bc_pool)):
    raw_umi_out.write(bc_nms[i])
    if (i == len(bc_pool) - 1):
        raw_umi_out.write("\n")
    else:
        raw_umi_out.write("\t")


for k in range(len(sample_spec_lines)):
    QC_out.write(sample_nms[k] + "\t")
    QC_out.write(str((PQ_QC_cnts[k])[0]) + "\t")
    QC_out.write(str((PQ_QC_cnts[k])[1]) + "\t")
    QC_out.write(str(sum(PQ_umi_cnts[k])) + "\n")
    

for k in range(len(sample_spec_lines)):
    raw_umi_out.write(sample_nms[k] + "\t")
    
    for i in range(len(bc_pool)):
        raw_umi_out.write(str(PQ_umi_cnts[k][i]))
        if (i == len(bc_pool) - 1):
            raw_umi_out.write("\n")
        else:
            raw_umi_out.write("\t")


QC_out.close()
raw_umi_out.close()

