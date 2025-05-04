This is our DATA 2020 project repository. Here, we reproduce the work from the paper "Bayesian Regression Tree Models for Causal Inference” by Hahn et. al. We generate covariates and from there,  calculate propensity scores and assign treatments. We then simulate outcomes for four different generated datasets across three types of models: BCF, BART, and Lasso Regression. Through these simulations, we test model performance with a variety of evaluation metrics, and find that BCF performs better than BART and Lasso, obtaining similar results to the paper.  

Our BCF model took a significant amount of time to run on OSCAR, and crashed on several occasions. When it successfully ran, it took around 12 hours total. BART and Lasso ran well locally with no issues. 

FILE DESCRIPTIONS:
BART_model, BCF_model, Lasso_model: Code to generate models and get models. Data generation is included in these files except for BART_model.

Visualizations: Code to create graphs from our data and models.

Data_and_models: Creates data and combines code from all of the _model files.

Link to Paper: https://doi.org/10.48550/arXiv.1706.09523
