#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 18:14:21 2023

@author: huk
"""

import pickle
import pandas as pd
import argparse
import warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description='SUMA: A Lightweigt model of Number of Neighbour Selection')
parser.add_argument('-c','--cells', help='Number of Cells', required=True)
parser.add_argument('-p','--pc', help='Number of PCs', required=True)
parser.add_argument('-e','--exp_type', help='Experiment Type (1: Droplet, 0: Spike (ERCC)', required=True)
parser.add_argument('-g','--hvgs', help='Number of Highly Variant Genes used for PCA', required=True)
parser.add_argument('-v','--variance', help='Percentage of Variants Explained With PCA', required=True)
args = parser.parse_args()


rf_model = pickle.load(open('SUMA.sav', 'rb'))

number_of_cells = args.cells
number_of_pc = args.pc
extype = args.exp_type
number_of_hvgs = args.hvgs
var_perc = args.variance

df = pd.DataFrame()
df["Number_of_Cells"] = [int(number_of_cells)] * 150
df["Number_of_Neighbours"] = list(range(1, 51)) * 3
df['Number_of_PCs'] = [int(number_of_pc)] * 150
df["Number_of_HVGs"] = [int(number_of_hvgs)] * 150
df["Algorithm"] = ['Walktrap'] * 50 + ['Leiden'] * 50 + ['Louvain'] * 50
df['Algorithm'] = df['Algorithm'].replace({'Leiden': 1, 'Louvain': 2, 'Walktrap': 3})
df["var_perc"] = [float(var_perc)] * 150
df["Exp_Type"] = [int(extype)] * 150

rf_prediction = rf_model.predict_proba(df)

df["RF_Pred"] = rf_prediction[: , 1]

df = df.sort_values(by = ["RF_Pred"] , ascending = [False])

df["Above_Moderate"] = df["RF_Pred"] > 0.50
df['Algorithm'] = df['Algorithm'].replace({1:'Leiden', 2:'Louvain', 3:'Walktrap'})

selected_results = df[0:10]

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)
pd.set_option('display.colheader_justify', 'center')
pd.set_option('display.precision', 2)
print(selected_results.to_string(index=False))
