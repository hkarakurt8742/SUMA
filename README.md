This repository stores the codes of SUMA: A Machine-Learning Based Lightweight Application for prediction of optimal number of neighbours for scRNA-Seq Clustering


Functions and scripts folders involves R scripts for data download and preprocess. 
Models and Training folder involves Python notebook for model training and the random forest model.
Models folder includes a random forest model.
SUMA.py file is the main application file.

To run the SUMA, just download main file and models from models folder.


```bash
python SUMA.py -c "number of cells" -p "Number of PCs" -e "Experiment Type" -g "Number of Highly Variant Genes"
```

App will directly print the output to your terminal. It will take about 2 seconds.
App requires Pickle and Pandas packages. 
Note: Be sure app and model are in the same folder.
