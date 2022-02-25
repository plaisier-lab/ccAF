# ccAF:  cell cycle ASU-Fred Hutch neural network based scRNA-seq cell cycle classifier
The ability to accurately assign a cell cycle phase based on a transcriptome profile has many potential uses in single cell studies and beyond. We have developed a cell cycle classifier based on a scRNA-seq optimized Neural Network (NN) based machine learning algorithm [ACTINN](https://pubmed.ncbi.nlm.nih.gov/31359028/). The ACTINN code was adapted from:  [https://github.com/mafeiyang/ACTINN](https://github.com/mafeiyang/ACTINN)

## Dependencies
There are four dependencies that must be met for ccAF to classify cell cycle states:
1. [numpy](https://numpy.org/) - ([install](https://numpy.org/install/))
2. [scipy](https://www.scipy.org/index.html) - ([install](https://www.scipy.org/install.html))
3. [scanpy](https://scanpy.readthedocs.io/en/latest/) - ([install](https://scanpy.readthedocs.io/en/latest/installation.html))
3. [tensorflow](https://www.tensorflow.org/) - ([install](https://www.tensorflow.org/install))

*Python dependency installation commands:*
> **NOTE!**  pip may need to be replaced with pip3 depending upon your setup.

```shell
pip3 install numpy scipy scanpy tensorflow
```

## Installation of ccAF classifier
The ccAF classifier can be installed with the following command:

```shell
pip install ccAF
```

## Alternatively use the ccAF Docker container
We facilitate the use of ccAF by providing a Docker Hub container [cplaisier/ccAF](https://hub.docker.com/r/cplaisier/ccAF) which has all the dependencies and libraries required to run the ccAF classifier. To see how the Docker container is configured plaese refer to the [Dockerfile](https://github.com/plaisier-lab/ccAF/blob/master/Dockerfile). Please [install Docker](https://docs.docker.com/get-docker/) and then from the command line run:

```shell
docker pull cplaisier/ccAF
```

Then run the Docker container using the following command (replace <path to scRNA-seq profiles directory> with the directory where you have the scRNA-seq data to be classified):

```shell
docker run -it -v '<path to scRNA-seq profiles directory>:/files' cplaisier/ccAF
```

This will start the Docker container in interactive mode and will leave you at a command prompt. You will then want to change directory to where you have your scRNA-seq or trasncriptome profiling data.

## Gene labels must be in human Gene Ensembl IDs to run ccAF
The data input into ccAF must use human Ensembl gene IDs (ENSG<#>), whithout the version number. If your data is not currenly labeled with Ensemble gene IDs you may try [mygene](https://docs.mygene.info/projects/mygene-py/en/latest/) or go to the [BioMart](http://uswest.ensembl.org/biomart/martview).
  
## Running ccAF against your scRNA-seq data
The first step in using ccAF is to import your scRNA-seq profiling data into scanpy. A scanpy data object is the expected input into the ccAF classifier:

```python
import scanpy
import ccAF

# Load WT U5 hNSC data used to train classifier as a loom file
set1_scanpy = sc.read_loom('data/WT.loom')

# Predict cell cycle phase labels
predictedLabels = ccAF.predict_labels(set1_scanpy)
```

More complete example is available as [test.py](https://github.com/plaisier-lab/ccAF/blob/master/tests/test.py) on the GitHub page.

## Contact
For issues or comments please contact:  [Chris Plaisier](mailto:plaisier@asu.edu)

## Citation
[Neural G0: a quiescent-like state found in neuroepithelial-derived cells and glioma.](https://doi.org/10.1101/446344) Samantha A. O'Connor, Heather M. Feldman, Chad M. Toledo, Sonali Arora, Pia Hoellerbauer, Philip Corrin, Lucas Carter, Megan Kufeld, Hamid Bolouri, Ryan Basom, Jeffrey Delrow, Jose L. McFaline-Figueroa, Cole Trapnell, Steven M. Pollard, Anoop Patel, Patrick J. Paddison, Christopher L. Plaisier. bioRxiv 446344; doi: [https://doi.org/10.1101/446344](https://doi.org/10.1101/446344)
