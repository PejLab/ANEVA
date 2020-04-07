# ANEVA [![DOI](https://zenodo.org/badge/183066057.svg)](https://zenodo.org/badge/latestdoi/183066057)


ANalysis of Expression VAriance

Package developed by: Pejman Mohammadi (firstname@Scripps.edu). See our [paper](https://science.sciencemag.org/content/366/6463/351.abstract) for method descritions.

This Matlab package allows you to estimate the amount of genetic variation in regulation in a population (Vg) using a population AE data. `Run_Example.m` allows you to do a demo run of ANEVA on a set of 20 sample allelic expression data files. The package includes all inputs and output for the example data. Please feel free to reach out if you have trouble running this on your data we're always interested in calculating reference Vg for any dataset.

To run the demo you have to first decompress all data files (```gunzip -r ANEVA/01_Data```). Also you should include [Pejtools](https://github.com/PejLab/Pejtools) in your Matlab path. 


## Further notes:
- If you don't know how to **generate ASE data**, you can start [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0762-6). Note: This model assumes variant level ASE data and is not compatible with haplotypic read counts such as Phaser data.

- You can download our **pre-calculated Vg estimates** for available datasets at [Datasets/Reference_Vg_Estimates](https://github.com/PejLab/Datasets/tree/master/Reference_Vg_Estimates).
