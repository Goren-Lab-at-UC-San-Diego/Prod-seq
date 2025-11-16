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
sample_R2s = [x[2] for x in sample_spec_lines]


# In[ ]:


# UMI pair counting
prod_midseq = "TAGAGAAG"
prod_arm_umi_len = 15
prod_umipair_cnts = []
prod_QC_cnts = []
for idx in range(len(sample_spec_lines)):
    curr_umi_cnts, curr_QC_cnts = CountUMIPairs(sample_R1s[idx], sample_R2s[idx], bc_pool, prod_midseq, prod_arm_umi_len)
    prod_umipair_cnts.append(curr_umi_cnts)
    prod_QC_cnts.append(curr_QC_cnts)

prod_umipair_normed = ScaleProdData(prod_umipair_cnts)


# In[ ]:


# Background inference
pmats, noise_preds = CalcPPIBackground(prod_umipair_normed, control_bc_idxs = control_bc_idxs, bc_nms = bc_nms)


# In[ ]:


# Output results to files
QC_out = open(out_prefix + ".QCCnts.tsv", "w")
raw_umi_out = open(out_prefix + ".RawUMIPairs.tsv", "w")
norm_umi_out = open(out_prefix + ".DepthNormedUMIPairs.tsv", "w")
noise_out = open(out_prefix + ".PredictedBackground.tsv", "w")
signal_out = open(out_prefix + ".PPIEnrichment.tsv", "w")

QC_out.write("#Sample\tRead_Count\tProd_Product_Read_Count\tUMI_Combination_Count\t")
raw_umi_out.write("#Sample\t")
norm_umi_out.write("#Sample\t")
noise_out.write("#Sample\t")
signal_out.write("#Sample\t")

for i in range(len(bc_pool)):
    QC_out.write(bc_nms[i] + "_Selfself_Count")
    if (i == len(bc_pool) - 1):
        QC_out.write("\n")
    else:
        QC_out.write("\t")

for i in range(len(bc_pool)):
    for j in range(i + 1, len(bc_pool)):
        raw_umi_out.write(bc_nms[i] + "&" + bc_nms[j])
        norm_umi_out.write(bc_nms[i] + "&" + bc_nms[j])
        
        if (i == len(bc_pool) - 2 and j == len(bc_pool) - 1):
            raw_umi_out.write("\n")
            norm_umi_out.write("\n")
        else:
            raw_umi_out.write("\t")
            norm_umi_out.write("\t")


for k in range(len(sample_spec_lines)):
    QC_out.write(sample_nms[k] + "\t")
    QC_out.write(str((prod_QC_cnts[k])[0]) + "\t")
    QC_out.write(str((prod_QC_cnts[k])[1]) + "\t")
    QC_out.write(str((prod_QC_cnts[k])[2]) + "\t")

    for i in range(len(bc_pool)):
        QC_out.write(str((prod_umipair_cnts[k])[i][i]))
        if (i == len(bc_pool) - 1):
            QC_out.write("\n")
        else:
            QC_out.write("\t")

for k in range(len(sample_spec_lines)):
    raw_umi_out.write(sample_nms[k] + "\t")
    norm_umi_out.write(sample_nms[k] + "\t")
    
    for i in range(len(bc_pool)):
        for j in range(i + 1, len(bc_pool)):
            raw_umi_out.write(str(prod_umipair_cnts[k][i][j]))
            norm_umi_out.write(str(prod_umipair_normed[k][i][j]))

            
            if (i == len(bc_pool) - 2 and j == len(bc_pool) - 1):
                raw_umi_out.write("\n")
                norm_umi_out.write("\n")
            else:
                raw_umi_out.write("\t")
                norm_umi_out.write("\t")

target_bc_idxs = []
for idx in range(len(bc_pool)):
    if (idx not in control_bc_idxs):
        target_bc_idxs.append(idx)

for i in range(len(target_bc_idxs)):
    for j in range(i + 1, len(target_bc_idxs)):
        noise_out.write(bc_nms[target_bc_idxs[i]] + "&" + bc_nms[target_bc_idxs[j]])
        signal_out.write(bc_nms[target_bc_idxs[i]] + "&" + bc_nms[target_bc_idxs[j]])
        
        if (i == len(target_bc_idxs) - 2 and j == len(target_bc_idxs) - 1):
            noise_out.write("\n")
            signal_out.write("\n")
        else:
            noise_out.write("\t")
            signal_out.write("\t")

for k in range(len(sample_spec_lines)):
    noise_out.write(sample_nms[k] + "\t")
    signal_out.write(sample_nms[k] + "\t")
    
    for i in range(len(target_bc_idxs)):
        for j in range(i + 1, len(target_bc_idxs)):
            
            noise_out.write(str(noise_preds[k][i][j]))

            curr_sig_diff = prod_umipair_normed[k][target_bc_idxs[i]][target_bc_idxs[j]]
            curr_signal_norm = curr_sig_diff / noise_preds[k][i][j]
            signal_out.write(str(curr_signal_norm))

            if (i == len(target_bc_idxs) - 2 and j == len(target_bc_idxs) - 1):
                noise_out.write("\n")
                signal_out.write("\n")
            else:
                noise_out.write("\t")
                signal_out.write("\t")

QC_out.close()
raw_umi_out.close()
norm_umi_out.close()
noise_out.close()
signal_out.close()

