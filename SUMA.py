#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 18:14:21 2023

@author: huk
"""

import pickle
import pandas as pd
import sys

catboost_model = pickle.load(open('tuned_catboost.sav', 'rb'))
lgbm_model = pickle.load(open('tuned_lgbm.sav', 'rb'))

number_of_cells = int(sys.argv[1])
number_of_pc = int(sys.argv[2])
extype = int(sys.argv[3])

df = pd.DataFrame()
df["Number_of_Neighbours"] = list(range(1, 51)) * 12
df['Number_of_PCs'] = [int(number_of_pc)] * 600
df["Number_of_Cells"] = [int(number_of_cells)] * 600
df["Number_of_HVGs"] = [500] * 150 + [1000] * 150 + [1500] * 150 + [2000] * 150
df["Algorithm"] = ['Walktrap'] * 200 + ['Leiden'] * 200 + ['Louvain'] * 200
df["Exp_TypeDroplet"] = [int(extype)] * 600

df = pd.get_dummies(df, columns = ['Algorithm'])
df.rename(columns = {'Algorithm_Leiden' : 'AlgorithmLeiden'} , inplace=True)
df.rename(columns = {'Algorithm_Louvain' : 'AlgorithmLouvain'} , inplace=True)
df.rename(columns = {'Algorithm_Walktrap' : 'AlgorithmWalktrap'} , inplace=True)

catboost_prediction = catboost_model.predict(df)
lightgbm_prediction = lgbm_model.predict(df)

df["CatBoost_Pred"] = catboost_prediction
df["LightGBM_Pred"] = lightgbm_prediction

selected_results = (df["CatBoost_Pred"] == 1) & (df["LightGBM_Pred"] == 1)
selected_results = df[selected_results]
algorithm_results = selected_results[['AlgorithmLeiden', 'AlgorithmLouvain' , 'AlgorithmWalktrap']]

selected_results = selected_results.drop(["CatBoost_Pred" , "LightGBM_Pred" , "Exp_TypeDroplet",
                           "Number_of_PCs" , "Number_of_Cells" , 
                           'AlgorithmLeiden', 'AlgorithmLouvain' , 'AlgorithmWalktrap'] , axis = 1)

selected_results["Algorithm"] = algorithm_results.idxmax(1)
selected_results["Algorithm"] = selected_results["Algorithm"].astype(str).str.replace("Algorithm","")

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)
pd.set_option('display.colheader_justify', 'center')
pd.set_option('display.precision', 2)
print(selected_results)
