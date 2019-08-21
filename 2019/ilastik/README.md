# Shallow Segmentation  with ilastik

- Ilastik is a user-friendly program that helps segment images using non-deep methods e.g. Random Forests, Autocontext, Adaboost, LDA. In this tutorial we will explore how Ilastik can help us segment the nucleolus from a multi-dimensional stack of raw microscope images.
---
## Setup & Pre-Processing 

1. Install Anaconda/[Miniconda](https://docs.conda.io/en/latest/miniconda.html). If you already have conda, skip to 2.
2. Create Conda environment: `conda create --name ilastik_fun python=3.6`
3. Activate environment: `conda activate ilastik_fun`
3. Install required packages (this may take a while): `pip install -r requirements.txt`
4. Launch Jupyternotebook in this environment: `jupyter notebook`
5. PreProcess the data in RawData.
---
## Ilastik 
- Assumes you are using MacOSX, if not follow instructions here: https://www.ilastik.org/documentation/basics/installation.html

1. Download here https://www.ilastik.org/download.html
2. Navigate to the directory where you downloaded ilastik
3. Extract its contents: Double click on download or run `tar xjf ilastik-1.3.2post1-OSX.tar.bz2`
4. Run ilastik: Right click and open.
5. Start a new Pixel Classification Project. 
---
## Post-Processing

1. Export probability map from Pixel Classification workflow. 
2. Load it into an Object classification Project
3. See PostProcessing_Object Quantification.ipynb
--
## Reading

- ilastik labelling best practices: https://blog.cellprofiler.org/2017/01/19/cellprofiler-ilastik-superpowered-segmentation/
- Great visualization of decision trees (the unit of Random Forests): http://www.r2d3.us/visual-intro-to-machine-learning-part-1/
- UC Berkeley CS 189 (Intro to Machine Learning) notes on Random Forests: https://people.eecs.berkeley.edu/~jrs/189/lec/16.pdf 
---
## Similar Software:

- [Segmentify (CZI)](https://github.com/sofroniewn/segmentify)
- [Allen Cell Structure Segmenter](https://www.allencell.org/segmenter.html)
- [CellProfiler](https://cellprofiler.org/)
