scRNA-seq data analysis
===========

This process was applied to each predefined immune subset independently.
This section comprises the following major steps:
- Preprocessing of reads by mapping them to the reference genome and collapsing them into UMI count matrices (10x's cellranger)
- Donor deconvolution based on hashtag oligos (HTOs) and based on genetic variants (Demuxlet; depends on the output from the section above).
- Clustering and dimensionality reduction.
- Population definition
- Sex-biased transcript calling

The first step is described in the methods section of our manuscript in detail. Here we provide the code for the rest of the steps.
