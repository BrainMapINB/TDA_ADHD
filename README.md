# TDA ADHD
Scripts to implement TDA in the ADHD200-NYU dataset.

[01-Preprocessing](https://github.com/BrainMapINB/TDA_ADHD/tree/master/01-Preprocessing) contains the scripts to preprocess fMRI datasets and average time series to regions of interest (ROI) from four different brain atlases.

[02-TDA](https://github.com/BrainMapINB/TDA_ADHD/tree/master/02-TDA) contains the scripts to compute Topological Data Analysis metrics to the functional connectivity matrices.

[03-Inference](https://github.com/BrainMapINB/TDA_ADHD/tree/master/03-Inference) contains the scripts to compute sample inferences.

In order to just replicate the same manuscript results, then, download and uncompress this repository and run the [TDA_ADHD200_NYU_inference.R](https://github.com/BrainMapINB/TDA_ADHD/blob/master/03-Inference/TDA_ADHD200_NYU_inference.R) file.  
Note: open '[TDA_ADHD200_NYU_inference.R](https://github.com/BrainMapINB/TDA_ADHD/blob/master/03-Inference/TDA_ADHD200_NYU_inference.R)' in RStudio and source it (preferably without Echo) and it will pop out every result. Without RStudio, just R, remove the lines of code 6 to 11 and set the working directory to the path of this directory, and then source it.
