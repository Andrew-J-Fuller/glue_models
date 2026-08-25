# NANOBretSeqr
The goal of this model is to predict the NANOBret signal ratio of the Cereblon E3Ligase and a protein of interest (POI) given only the amino acid sequence of each (including tags).

Currently, the model takes csv files which are named after the protein of interest (eg. "FOXA1.csv"), queries UNIPROT for the AA sequence, concatenates the nanoluc AA sequence, and then using scikit-learn's MLPRegressor, trains and predicts NANOBret signal ratios.

These CSV files are experimentally generated data from NANOBret assays, and come in the form of:

(E3-tag combo+NanoLuc-POI combo, replicate1 signal ratio, replicate2 signal ratio,...(8 replicates)).

## Featurizing the data
To better model the assay, the following implementations were considered to improve the model

1. Include tag-protein combinations.
  -  Since the cells are transfected with the tag's location in alternating possitions, it was thought that including the distinction in the model might better generalize, instead of overfitting to the tag sequence (which is ubiquitously included on every protein).
2. One-hot encoding + AA properties.
  -  Originally, each AA would get assigned a value 1->20. However, I was concerned with the model treating AA's as integers, i.e. _cystine_ as greater than, or later in sequence to _alanine_. This is solved using one-hot encoding. However, I also wasn't confident one-hot alone would be sufficient for the model to separate the AA's very distinct properties. Included then, is also the hydrophobicity using the Bandyopadhyay-Mahler hydrophobicity scale, the literal charge, and polarity of each AA. This way, each AA is utterly unique to another, and the model can exploit this.
3. Aggregation.
  - With the preceding features, this would create a 1D array length (n)24 where n is the number of amino acids in a sequence. This results in a maximum of 14,400 features per sequence!! (given a maximum sequence lenght of 600). To prevent n << f we average and sum the arrays down to a 1D array of length (n)4 where n is the number of amino acids per sequence.

## Training the model/Predicting signal ratios
The model was trained using 72 protein/ligase/tag combos with 20% of the data reserved for testing the model.

When trained in this fashion, and prompted to predict the NANOBret signal ratio of PIMT (a well described target of CRBN) the model resulted in:
- Test MSE: 2.74824
- Training R^2: 0.99323
- Predicted Ratio: 1.772

These results can be interpreted as demonstrating the models poor prediction power, and overfitting. (A proper positive hit in NANOBret might result in a ratio of 30 or higher)

## Issues
One of the main issues in this case is very small number of training samples. 72 samples are simply not enough data to properly train a MLPRegressor model. Additionally, no pretrained protein embeddings are used. One could imagine that protBERT might be used to provide the model much better featurization. Finally, the model is trained using noisy, experimental data lacking a large sample of positive hits preventing good generalizations.

## Next steps
In the future if this experimental data driven approach is pursued, I'd like to prioritize the denoising of the input data (preferably by wet lab assay optimization), the implementation of a caching system to prevent quering UNIPROT n times for every experiment, a far larger data set, and proper protein embeddings. I believe these could potentially increase the model's generalization, helping to focus the design, test, analyze loop of the wet lab on likely hits.  
