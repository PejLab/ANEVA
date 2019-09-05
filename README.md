# ANEVA
ANalysis of Expression VAriance

Package developed by: Pejman Mohammadi (firstname@Scripps.edu). See our [preprint](https://www.biorxiv.org/content/biorxiv/early/2019/05/09/632794.full.pdf) for method descritions.

This Matlab package allows you to estimate the amount of genetic variation in regulation in a population (Vg) using a population AE data. `Run_Example.m` allows you to do a demo run of ANEVA on a set of 20 sample allelic expression data files. The package includes all inputs and output for the example data. 

To run the demo you have to first decompress all data files (```gunzip -r ANEVA/01_Data```). Also you should include [Pejtools](https://github.com/PejLab/Pejtools) in your Matlab path.

## Further notes:
- If you don't know how to **generate ASE data**, you can start [here](https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/) and [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0762-6). Note: This model assumes variant level ASE data and is not compatible with haplotypic read counts such as Phaser data.

- You can download our **pre-calculated Vg estimates** for available datasets at [Datasets/Reference_Vg_Estimates](https://github.com/PejLab/Datasets/tree/master/Reference_Vg_Estimates).
