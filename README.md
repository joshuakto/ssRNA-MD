# ssRNA-MD
## CVAE
Directory containing the source code for training neural networks to learn the distance array for atoms in MD simulation of RNA molecules in solutions with varying ions type and concentration
- Variational Autoencoder (latent space regularized, t-SNE plots of encoded data from different MD runs cluster by concentrations of ions and types of ions)
- Conditional Variational Autoencoder (work in progress)
  - Goal: generate prediction of distance array for out-of-distribution ions concentration
## VAE
Code archive: this model is not able to capture a regularized latent space for the distance array
## MD_sim
Directory containing the source code for running molecular dynamics simulation in varying ions concentration and types.
