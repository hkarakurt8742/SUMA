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
args = parser.parse_args()


rf_model = pickle.load(open('rf_model_65.sav', 'rb'))

number_of_cells = args.cells
number_of_pc = args.pc
extype = args.exp_type
number_of_hvgs = args.hvgs

df = pd.DataFrame()
df["Number_of_Neighbours"] = list(range(1, 51)) * 3
df['Number_of_PCs'] = [int(number_of_pc)] * 150
df["Number_of_Cells"] = [int(number_of_cells)] * 150
df["Number_of_HVGs"] = [int(number_of_hvgs)] * 150
df["Algorithm"] = ['Walktrap'] * 50 + ['Leiden'] * 50 + ['Louvain'] * 50

df = pd.get_dummies(df, columns = ['Algorithm'])
df.rename(columns = {'Algorithm_Leiden' : 'Algorithm.Leiden'} , inplace=True)
df.rename(columns = {'Algorithm_Louvain' : 'Algorithm.Louvain'} , inplace=True)
df.rename(columns = {'Algorithm_Walktrap' : 'Algorithm.Walktrap'} , inplace=True)
df["Exp_TypeDroplet"] = [int(extype)] * 150
df.rename(columns = {'Exp_TypeDroplet' : 'Exp_Type.Droplet'} , inplace=True)

rf_prediction = rf_model.predict(df)

df["RF_Pred"] = rf_prediction

df = df.sort_values(by = ["RF_Pred"] , ascending = False)

df["Above_Moderate"] = df["RF_Pred"] > 0.2226
df["Above_Good"] = df["RF_Pred"] > 0.8605
algorithm_results = df[['Algorithm.Leiden', 'Algorithm.Louvain' , 'Algorithm.Walktrap']]

df = df.drop(["Exp_Type.Droplet",
                           "Number_of_PCs" , "Number_of_Cells" , 
                           'Algorithm.Leiden', 'Algorithm.Louvain' , 'Algorithm.Walktrap'] , axis = 1)

df["Algorithm"] = algorithm_results.idxmax(1)
df["Algorithm"] = df["Algorithm"].str.replace("Algorithm.","")

selected_results = df[0:10]

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)
pd.set_option('display.colheader_justify', 'center')
pd.set_option('display.precision', 2)
print(selected_results.to_string(index=False))
