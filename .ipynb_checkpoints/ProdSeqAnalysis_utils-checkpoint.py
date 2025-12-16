#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/python3
import gzip
import math
import numpy as np
import os
from tqdm import tqdm
import glob

from scipy.optimize import minimize
import seaborn as sns
from scipy.stats import betabinom
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from matplotlib import rc
plt.rcParams["font.family"] = "Arial"

from scipy.stats import mannwhitneyu


# In[ ]:


bc_nms_default = ["H3K4me3", "H3K27ac", "H3K27me3", "EZH2", "EED", "SUZ12", "AEBP2", "MED12", "CycC", "H3K27M", "EGFR", "HA-Tag"]


# In[ ]:


def ReadFileSeqs(file_path):
    print("Reading file: " + file_path)
    fastq_file = gzip.open(file_path, "r")
    lines = []
    ln_idx = -1
    for curr_line in tqdm(fastq_file):
        ln_idx += 1
        if (ln_idx % 4 == 1):
            lines.append(curr_line.decode("utf-8").strip("\n"))
    
    fastq_file.close()
    return lines

def RevComp(seq):
    res = ""
    for i in range(len(seq)):
        if (seq[i] == "A"):
            res = "T" + res
        elif(seq[i] == "C"):
            res = "G" + res
        elif(seq[i] == "G"):
            res = "C" + res
        elif(seq[i] == "T"):
            res = "A" + res
    return res

def SingleArmCheckStruct(arm_seq, list_barcodes, arm_find_seq, UMI_len):
    if (arm_find_seq in arm_seq):
        find_pos = arm_seq.find(arm_find_seq)
        arm_barcode = arm_seq[find_pos - 10 : find_pos]
        arm_umi = arm_seq[find_pos - 10 - UMI_len : find_pos - 10]
        arm_pre_umi = arm_seq[find_pos - 15 - UMI_len : find_pos - 10 - UMI_len]
        if (arm_pre_umi == "TAAGC" and arm_barcode in list_barcodes):
            return list_barcodes.index(arm_barcode), arm_umi
        
    return -1, ""

def CountUMIPairs(R1_file, R2_file, barcodes, arm_find_seq, UMI_len):
    R1_lines = ReadFileSeqs(R1_file)
    R2_lines = ReadFileSeqs(R2_file)
    all_UMI_combs = set()
    UMIs_cnts = np.zeros((len(barcodes), len(barcodes))) # matrix format of interaction UMI pair counts
    QC_cnts = [0, 0, 0] # Number of reads, number of Prod-seq product reads, UMI dedup number of product reads
    
    for i in tqdm(range(len(R1_lines))):
        QC_cnts[0] += 1
        if (arm_find_seq in R1_lines[i] and arm_find_seq in R2_lines[i]):
            arm1_idx, arm1_umi = SingleArmCheckStruct(R1_lines[i], barcodes, arm_find_seq, UMI_len)
            arm2_idx, arm2_umi = SingleArmCheckStruct(R2_lines[i], barcodes, arm_find_seq, UMI_len)
            
            if ((arm1_idx != -1) and (arm2_idx != -1)):
                QC_cnts[1] += 1
                if ((arm1_umi + arm2_umi not in all_UMI_combs) and (arm2_umi + arm1_umi not in all_UMI_combs)
                   and (arm1_umi != arm2_umi)):
                    # New pair of unique UMIs
                    all_UMI_combs.add(arm1_umi + arm2_umi)
                    QC_cnts[2] += 1
                    
                    if (arm1_idx == arm2_idx):
                        UMIs_cnts[arm1_idx][arm1_idx] += 1
                    else:
                        UMIs_cnts[arm1_idx][arm2_idx] += 1
                        UMIs_cnts[arm2_idx][arm1_idx] += 1
    return UMIs_cnts, QC_cnts

def CountPQUMIs(R1_file, barcodes, arm_find_seq, UMI_len):
    R1_lines = ReadFileSeqs(R1_file)
    read_cnts = [0, 0] # [num all read pairs, num correct read pairs found product]
    
    single_UMIs = []
    for i in range(len(barcodes)):
        single_UMIs.append(set())
    
    for i in range(len(R1_lines)):
        read_cnts[0] += 1
        single_protein_idx, single_protein_umi = SingleArmCheckStruct(R1_lines[i], barcodes, arm_find_seq, UMI_len)
        if (single_protein_idx != -1):
            read_cnts[1] += 1
            (single_UMIs[single_protein_idx]).add(single_protein_umi)
    
    single_protein_counts = []
    for i in range(len(barcodes)):
        single_protein_counts.append(len(single_UMIs[i]))
    return single_protein_counts, read_cnts


# In[ ]:


def ScaleProdData(prod_xl_original_prenorm):
    prod_xl_noselfself = np.array(prod_xl_original_prenorm)
    for k in range(len(prod_xl_original_prenorm)):
        np.fill_diagonal(prod_xl_noselfself[k], 0)
    min_depth = np.min([np.sum(prod_xl_noselfself[x]) for x in range(len(prod_xl_original_prenorm))])
    res = [prod_xl_noselfself[x] * min_depth / np.sum(prod_xl_noselfself[x]) for x in range(len(prod_xl_original_prenorm))]
    res = np.round(res)
    return res

def GLMMean(param_vec, x_vecs):
    a = param_vec[0]
    b = param_vec[1]
    param_vec_use = param_vec[2 : ]
    res = np.matmul(param_vec_use, np.transpose(x_vecs))
    return res
    

def CalcPPIBackground(prod_xl_original, control_bc_idxs = [11], bc_nms = bc_nms_default):
    # Fit the background distribution
    prod_xl_fitdata_r1 = []
    prod_xl_fitX_r1 = []
    linear_idxes_r1 = []

    target_bc_idxs = []
    for i in range(len(prod_xl_original[0])):
        if (i not in control_bc_idxs):
            target_bc_idxs.append(i)
    
    for k in range(len(prod_xl_original)):
        for idx_i in range(len(target_bc_idxs)):
            i = target_bc_idxs[idx_i]
            for idx_j in range(idx_i + 1, len(target_bc_idxs)):
                j = target_bc_idxs[idx_j]
                prod_xl_fitdata_r1.append((prod_xl_original[k][i][j]))
                curr_fitX = [1]

                for ctrl_idx in control_bc_idxs:
                    curr_fitX.append(prod_xl_original[k][i][ctrl_idx] + prod_xl_original[k][j][ctrl_idx])
                prod_xl_fitX_r1.append(curr_fitX)
                linear_idxes_r1.append([k, i, j])

    prod_xl_fitX_r1 = np.array(prod_xl_fitX_r1)
    prod_xl_fitdata_r1 = np.array(prod_xl_fitdata_r1)

    res_r1 = GLMFit(prod_xl_fitX_r1, prod_xl_fitdata_r1, n_ite = 500)
    pred_res_r1 = GLMMean(res_r1.x, prod_xl_fitX_r1)
    
    r1_a = (res_r1.x)[0]
    r1_b = (res_r1.x)[1]

    prod_xl_fitdata_r2 = []
    prod_xl_fitX_r2 = []
    
    for pred_idx in range(len(pred_res_r1)):
        curr_k = (linear_idxes_r1[pred_idx])[0]
        curr_i = (linear_idxes_r1[pred_idx])[1]
        curr_j = (linear_idxes_r1[pred_idx])[2]
        curr_n = round((pred_res_r1[pred_idx] * (r1_a + r1_b)) / r1_a)
        
        curr_observed = prod_xl_original[curr_k][curr_i][curr_j]
        curr_pval = betabinom.sf(k = curr_observed, n = curr_n, a = r1_a, b = r1_a)

        if (curr_pval >= 0.025):
            prod_xl_fitdata_r2.append(prod_xl_fitdata_r1[pred_idx])
            prod_xl_fitX_r2.append(prod_xl_fitX_r1[pred_idx])

    prod_xl_fitX_r2 = np.array(prod_xl_fitX_r2)
    prod_xl_fitdata_r2 = np.array(prod_xl_fitdata_r2)

    res_r2 = GLMFit(prod_xl_fitX_r2, prod_xl_fitdata_r2, n_ite = 1000, param_init = res_r1.x)
    pred_res_r2 = GLMMean(res_r2.x, prod_xl_fitX_r1)
    
    p_mat = np.ones((len(prod_xl_original),
                    len(target_bc_idxs),
                    len(target_bc_idxs)))
    xl_pred_r2 = np.zeros((len(prod_xl_original),
                           len(target_bc_idxs),
                           len(target_bc_idxs)))
    r2_a = (res_r2.x)[0]
    r2_b = (res_r2.x)[1]

    for r2_idx in range(len(pred_res_r2)):
        curr_k = (linear_idxes_r1[r2_idx])[0]
        curr_i = (linear_idxes_r1[r2_idx])[1]
        curr_j = (linear_idxes_r1[r2_idx])[2]
        curr_observed = prod_xl_original[curr_k][curr_i][curr_j]

        curr_n = round((pred_res_r2[r2_idx] * (r2_a + r2_b)) / r2_a)
        curr_pval = betabinom.sf(curr_observed, n = curr_n, a = (res_r2.x)[0], b = (res_r2.x)[1])

        p_kidx = curr_k
        p_iidx = target_bc_idxs.index(curr_i)
        p_jidx = target_bc_idxs.index(curr_j)

        p_mat[p_kidx][p_iidx][p_jidx] = curr_pval
        xl_pred_r2[p_kidx][p_iidx][p_jidx] = pred_res_r2[r2_idx]

    return p_mat, xl_pred_r2

def GLMFit(x_vecs, y_vec, n_ite = 200, param_init = []):

    zero_tolerance = 1e-8
    bnds_l = []
    bnds_l.append((zero_tolerance, np.inf))
    bnds_l.append((zero_tolerance, np.inf))
    for i in range(int(len(x_vecs[0] - 2))):
        bnds_l.append((-np.inf, np.inf))
    bnds = tuple(bnds_l)

    def CalcNegLogLikelihood(param_vec, x_vecs, y_vec):
        # Beta-Binom: a, b, betas for GLM, const for GLM
        res = 0
        a = param_vec[0]
        b = param_vec[1]
        param_vec_use = param_vec[2 : ]
        mean_vec = np.matmul(param_vec_use, np.transpose(x_vecs))
        n_vec = np.round(mean_vec * (a + b) / a)
        res = np.sum([-betabinom.logpmf(k = round(y_vec[idx]), n = n_vec[idx], a = a, b = b) for idx in range(len(n_vec))])

        return res
    
    FunctionToFit = lambda x : CalcNegLogLikelihood(x, x_vecs, y_vec)
    if (len(param_init) == 0):
        init_guess = [0.5, 0.5]
        for i in range(len(bnds_l) - 2):
            init_guess.append(np.sum(y_vec) / np.sum([curr_x[i] for curr_x in x_vecs]))
    else:
        init_guess = param_init

        
    res = minimize(FunctionToFit, x0 = init_guess, method = "Nelder-Mead", bounds = bnds, options = {"maxiter" : n_ite})
    print(res)
    return res


# In[ ]:


def ReadProdTSVFile(file_path):
    f_lines = open(file_path).readlines()
    f_split = [line.strip().split("\t") for line in f_lines if line.strip()]
    col_nms = (f_split[0])[1:]
    sample_nms = [x[0] for x in f_split[1:]]
    values = [np.float64(x[1:]) for x in f_split[1:]]

    return values, sample_nms, col_nms


def GenerateProdQCPlots(QC_tsv_path):
    QC_lines = [x.strip().split("\t") for x in ((open(QC_tsv_path)).readlines())[1:] ]
    QC_split = []

    for i in range(len(QC_lines)):
        curr_line = [QC_lines[i][0]]
        for j in range(1, len(QC_lines[0])):
            curr_line.append(float(QC_lines[i][j]))
        QC_split.append(curr_line)

    correct_portion = [(x[2] / x[1]) for x in QC_split]
    unique_portion = [(x[3] / x[1]) for x in QC_split]
    expr_nms = [x[0] for x in QC_split]
    all_selfself_ratios = [sum(x[4: ]) / x[3] for x in QC_split]
    
    fig1 = plt.figure(figsize = (12, 12))
    ax = fig1.add_subplot(111)
    plt.ylim(top = 1)
    p1 = ax.bar([str(x) for x in range(len(expr_nms))], correct_portion)
    p2 = ax.bar([str(x) for x in range(len(expr_nms))], unique_portion)
    plt.legend((p2[0], p1[0]), ("Unique UMI pairs", "Duplicates"))
    plt.title("Proportion of reads where Prod-seq product is found")
    x_legend = '\n'.join([(str(n) + " - " + str(nm)) for n,nm in enumerate(expr_nms)])
    ax.text(.7, .2, x_legend, transform = ax.figure.transFigure)
    fig1.subplots_adjust(right = .65)

    for i in range(len(expr_nms)):
        plt.text(i, correct_portion[i], "{:.1f}".format((QC_split[i])[2] / 1000) + "K reads", ha = 'left', va = "bottom",
                 fontsize = 12, rotation_mode = "anchor", rotation = 60)
    
    fig2 = plt.figure(figsize = (8, 8))
    ax = fig2.add_subplot(111)
    plt.ylim(top = 1)
    
    p1 = ax.bar([str(x) for x in range(len(expr_nms))], all_selfself_ratios)
    plt.title("Proportion of UMI pairs that are self-self")
    x_legend = '\n'.join([(str(n) + " - " + str(nm)) for n,nm in enumerate(expr_nms)])
    ax.text(.7, .2, x_legend, transform = ax.figure.transFigure)
    fig2.subplots_adjust(right = .65)

    for i in range(len(expr_nms)):
        plt.text(i, all_selfself_ratios[i], "{:.3f}".format(all_selfself_ratios[i]), ha = 'left', va = "bottom",
                 fontsize = 12, rotation_mode = "anchor", rotation = 60)

    return fig1, fig2


def GroupedPPIHeatmap(PPI_signals, sample_nms, PPI_pairnms, bc_pair_groups, plot_vmax,
                      plot_size = (6, 6),
                      plot_vmin = [],
                      plot_cmaps = [],
                      sample_group_sizes = [],
                      cbar_ax_locs = []
                      ):

    if len(plot_vmin) == 0:
        plot_vmin = [0 for x in range(len(plot_vmax))]
    signal_dict = [{} for x in range(len(sample_nms))]
    for k in range(len(sample_nms)):
        for i in range(len(PPI_pairnms)):
            (signal_dict[k])[PPI_pairnms[i]] = (PPI_signals[k])[i]
            pairnm_split = (PPI_pairnms[i]).split("&")
            pairnm_rev = pairnm_split[1] + "&" + pairnm_split[0]
            (signal_dict[k])[pairnm_rev] = (PPI_signals[k])[i]
    
    if (len(plot_vmax) == 0):
        plot_vmax = [None for x in range(len(bc_pair_groups))]

    if (len(plot_vmin) == 0):
        plot_vmin = [None for x in range(len(bc_pair_groups))]

    if (len(plot_cmaps) == 0):
        plot_cmaps = [sns.color_palette("YlOrBr", as_cmap=True) for x in range(len(bc_pair_groups))]
        
    fig = plt.figure(figsize = plot_size)
    num_row = len(sample_nms)
    num_col = sum([len(x) for x in bc_pair_groups])
    
    for i in range(len(bc_pair_groups)):
        curr_pairs = bc_pair_groups[i]
        curr_data = []
        for pair in curr_pairs:
            curr_data.append([signal_dict[k][pair] for k in range(len(signal_dict))])
        curr_data = np.transpose(np.array(curr_data))

        curr_ax = plt.subplot2grid((num_row, num_col), (1, sum([len(x) for x in bc_pair_groups[0 : i]])),
                                       rowspan = num_row, colspan = len(bc_pair_groups[i]))
        if (len(cbar_ax_locs) > 0):
            cbar_ax = fig.add_axes(cbar_ax_locs[i])
            
            if (i == 0):
                
                ax = sns.heatmap(curr_data, ax=curr_ax, vmax = plot_vmax[i], vmin = plot_vmin[i], cmap = plot_cmaps[i],
                                 yticklabels = sample_nms, linewidths = .5,
                                 cbar_ax = cbar_ax,
                                 cbar_kws = {"orientation": "horizontal", "location": "top",
                                             "ticks": [plot_vmin[i], (plot_vmin[i] + plot_vmax[i]) / 2, plot_vmax[i]],
                                             'label': 'PPI enrichment'
                                            },
                                 xticklabels = bc_pair_groups[i], linecolor='lightgrey')
            else:
                
                
                ax = sns.heatmap(curr_data, ax=curr_ax, vmax = plot_vmax[i], vmin = plot_vmin[i], cmap = plot_cmaps[i],
                                 yticklabels = False, linewidths = .5,
                                 cbar_ax = cbar_ax,
                                 cbar_kws = {"orientation": "horizontal", "location": "top",
                                             "ticks": [plot_vmin[i], (plot_vmin[i] + plot_vmax[i]) / 2, plot_vmax[i]],
                                             'label': 'PPI enrichment'
                                            },
                                 xticklabels = bc_pair_groups[i], linecolor='lightgrey')

        else:
            if (i == 0):
                
                ax = sns.heatmap(curr_data, ax=curr_ax, vmax = plot_vmax[i], vmin = plot_vmin[i], cmap = plot_cmaps[i],
                                 yticklabels = sample_nms, linewidths = .5,
                                 cbar_kws = {"shrink": 0.8, "orientation": "horizontal", "location": "top",
                                             "ticks": [plot_vmin[i], (plot_vmin[i] + plot_vmax[i]) / 2, plot_vmax[i]]},
                                 xticklabels = bc_pair_groups[i], linecolor='lightgrey')
            else:
                
                
                ax = sns.heatmap(curr_data, ax=curr_ax, vmax = plot_vmax[i], vmin = plot_vmin[i], cmap = plot_cmaps[i],
                                 yticklabels = False, linewidths = .5,
                                 cbar_kws = {"shrink": 0.8, "orientation": "horizontal", "location": "top",
                                             "ticks": [plot_vmin[i], (plot_vmin[i] + plot_vmax[i]) / 2, plot_vmax[i]]},
                                 xticklabels = bc_pair_groups[i], linecolor='lightgrey')

        ax.set_xticklabels(ax.get_xticklabels(), va = "top", ha = 'center', rotation=90, fontsize = 15) 
        ax.set_yticklabels(ax.get_yticklabels(), fontsize = 15) 
        cbar = ax.collections[0].colorbar
        cbar.ax.tick_params(labelsize = 12)
        cbar.set_label('PPI enrichment', size = 15)

        outbox_color = "k"
        outbox_linewidth = 4
        inner_break_color = "k"
        inner_break_linewidth = 2
        
        if (len(sample_group_sizes) > 0):
            for k in range(1, len(sample_group_sizes)):
                ax.axhline(y = sum([x for x in sample_group_sizes[0 : k]]), color = inner_break_color,
                           linewidth = inner_break_linewidth)
        ax.axhline(y = 0, color = outbox_color, linewidth = outbox_linewidth)
        ax.axhline(y = len(sample_nms), color = outbox_color, linewidth = outbox_linewidth)
        ax.axvline(x = 0, color = outbox_color, linewidth = outbox_linewidth)
        ax.axvline(x = len(bc_pair_groups[i]), color = outbox_color, linewidth = outbox_linewidth)
        
    return fig


# In[ ]:


def CalcPQFromUMICnts(raw_cnts, control_idx = [7, 8]):
    PQ_prop = []
    PQ_Ctrl_normed = []
    for i in range(len(raw_cnts)):
        PQ_prop.append([x / sum(raw_cnts[i]) for x in raw_cnts[i]])
        curr_ctrl_sum = sum([raw_cnts[i][x] for x in control_idx])
        PQ_Ctrl_normed.append([x / curr_ctrl_sum for x in raw_cnts[i]])

    return PQ_prop, PQ_Ctrl_normed

def PQHeatmap(plot_data, sample_labels, plot_barcodes = bc_nms_default, cmap_use = "Greys"):
    fig, ax = plt.subplots(figsize = (20, len(plot_data)))
    
    im = ax.imshow(plot_data, cmap = cmap_use)
    cbar = fig.colorbar(im, ax = ax)
    cbar.ax.tick_params(labelsize = 15)
    
    ax.set_xticks(np.arange(len(plot_barcodes)))
    ax.set_yticks([x for x in range(len(plot_data))])
    ax.set_xticklabels(plot_barcodes, fontsize = 20)
    ax.set_yticklabels(sample_labels, fontsize = 15)
    
    ax.xaxis.tick_top()
    plt.setp(ax.get_xticklabels(), rotation = 45, ha = "left",
         rotation_mode = "anchor")
    
    for i in range(len(plot_barcodes) - 1):
        ax.plot([i + 0.5, i + 0.5], [-0.5, -0.5 + len(plot_data)], color = "black", lw = 1.5)
    
    for i in range(len(plot_data) - 1):
        ax.plot([-0.5, -0.5 + len(plot_barcodes)], [i + 0.5, i + 0.5], color = "black", lw = 1.5)
    
    return fig


# In[ ]:





# In[ ]:




