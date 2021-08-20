#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 13:19:09 2021

@author: anthony.cesnik
"""

from SingleCellProteogenomics import (FucciPseudotime,
                                      RNACellCycleDependence,
                                      RNADataPreparation,
                                      RNAVelocity,
                                      utils)
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import shutil
import os
import argparse

# Make PDF text readable
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.dpi"] = 300

bioccd = np.genfromtxt(
    "input/ProteinData/BiologicallyDefinedCCD.txt", dtype="str"
)  # from mitotic structures
wp_ensg = np.load("output/pickles/wp_ensg.npy", allow_pickle=True)
ccd_comp = np.load("output/pickles/ccd_comp.npy", allow_pickle=True)
nonccd_comp = np.load("output/pickles/nonccd_comp.npy", allow_pickle=True)
u_plates = ["355","356","357"]

adata, phases_filt = RNADataPreparation.read_counts_and_phases(
    "Counts", False, "protein_coding", u_plates
)
adata = RNADataPreparation.zero_center_fucci(adata)
FucciPseudotime.pseudotime_rna(adata)
valuetype, use_spikeins, biotype_to_use = "Tpms", False, ""
adata, phases = RNADataPreparation.read_counts_and_phases(
    valuetype, use_spikeins, biotype_to_use, u_plates
)
adata, phasesfilt = RNADataPreparation.qc_filtering(
    adata, do_log_normalize=True, do_remove_blob=True
)
adata = RNADataPreparation.zero_center_fucci(adata)

ccd_regev_filtered, ccd_filtered, nonccd_filtered = utils.ccd_gene_lists(adata)
adata_ccdprotein, adata_nonccdprotein, adata_regevccdgenes = RNADataPreparation.is_ccd(
    adata, wp_ensg, ccd_comp, nonccd_comp, bioccd, ccd_regev_filtered
)

# Generate plots with expression of genes overlayed
expression_data = adata.X
normalized_exp_data = (expression_data.T / np.max(expression_data, axis=0)[:, None]).T

# Cluster the expression into phases and analyze it that way
bulk_phase_tests = RNACellCycleDependence.analyze_ccd_variation_by_phase_rna(
    adata, normalized_exp_data, biotype_to_use
)

# Isoform analysis with all isoforms and make dataframe
doAllPlotsAndAnalyses=False
adata_isoform, ccdtranscript_isoform = RNACellCycleDependence.analyze_isoforms(
    adata, wp_ensg, ccd_comp, nonccd_comp, u_plates, make_mvavg_plots_isoforms=doAllPlotsAndAnalyses, biotype_to_use=""
)

