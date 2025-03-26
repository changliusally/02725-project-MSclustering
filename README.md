This is for the project of course Computational methods for Proteomics and Metagenomics (02725), where we are studying and comparing possible methods to accelerate the process of Mass Spectrum data (mzXML) clustering, which is a key step in interpreting the Mass Spectrum experiment. Here is the code needed to compare the effect of sparse index on improving the clustering efficiency. 

**read_file.py** is written to read in different type of mass spectrometry data like mzXML or mzML, which is later integrated into the normal.py.

**normal.py** can do a normal hierarchical clustering on the spectrum data and record the time of this runnning.

**sparse_index.py** is an improved hierachical clustering with better data strcuture for efficient storage and calling of huge mass spectrometry data, and also record the running timing.
