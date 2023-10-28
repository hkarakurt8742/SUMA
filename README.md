This repository stores the codes of SUMA: A Machine-Learning Based Lightweight Application for prediction of optimal number of neighbours for scRNA-Seq Clustering


Functions and scripts folders involves R scripts for data download and preprocess. 
Models and Training folder involves Python notebook for model training and the random forest model.
Models folder includes a random forest model.
SUMA.py file is the main application file.

To run the SUMA, just download main file and models from models folder.


```bash
python SUMA.py [-h] -c "number_of_cells" -p "number_of_pcs" -e "experiment_type" -g "number_of_highly_variant_genes" -v "explained_variance_percentage"
```

App will directly print the output to your terminal. It will take about 2 seconds.
App requires Pickle and Pandas packages. 
Note: Be sure app and model are in the same folder.
