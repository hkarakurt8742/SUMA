This repository stores the codes of SUMA: A Machine-Learning Based Lighweight Application for scRNA-Seq Clustering


Preprocess and training folders involves R and Python script. R scripts are used to download and pre-process publicly available scRNA-Seq datasets.
Preprocessed and encoded output dataframe (Preprocessing_and_Training/Model_Training/encoded_dataset.txt) imported to Python.
2 Python Notebooks show the parameter tuning steps of LightGBM and CatBoost models.
Models folder includes the 2 models.
SUMA.py file is the main application file.

To run the SUMA, just download main file and models from models folder.


```bash
python app_v2.py "number_of_cells" "number_of_pcs_to_be_used" "the_platform_type(1 for Droplet, 0 for Spike-based)"
```

App will directly print the output to your terminal. It will take about 2 seconds.
App requires Pickle and Pandas packages. 
Note: Be sure app and models are in the same folder.
