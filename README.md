This is the repository for the course AI foundation models in biomedicine.

# Project Structure

We created embeddings with scGPT for the data from the dataset [GSE158055](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158055). Then we created a baseline with the unsupervised clustering algorithm k-Means. After that we trained a feedforward neural network for classification of the severity and the disease stage. 

## Reconstruction

To reconstruct our findings, download the dataset. Depending on your hardware, you need to split it into chunks, because it won't fit into memory at once. We executed the scgpt_create_embeddings.ipynb on Kaggle with 25 chunks. If you get dependency issues for the scGPT library, don't worry, just try it anyway.
Then you can take the pickle file and use it for the different classifier notebooks. You might need to adjust the labels depending on your needs. 