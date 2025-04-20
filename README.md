# FNDIV_Final_Project_CBPT

This code is meant for ERP data preprocessed using BrainVision Analyzer v 2.3. The MLX version is a more user-friendly version of the code. It allows you to edit parameters and run the entire cluster based permutation testing locally in one file. However, the component .m files are also accessible. This includes separated function scripts, preprocessing script, processing script, and data visualization scripts that can be used independently.

# Requirements
- MATLAB R2024b
- ERP data preprocessed and exported from BVA (this includes interpolation of electrodes*, IIRC filters, ocular correction, baseline correction, artifact rejection, segmentation, and averaging ERPs across trials). *where necessary
	- This script assumes a sampling rate of 1000 Hz
- all required functions are contained in this repository. No additional packages are required to download in Matlab.
- 'Cond Key' that assigns group membership to each participant (in 1s and 2s).
- Folder with one csv file for each participant. These CSV files export from BVA and do not need any editing prior to loading them into the script.

# To get started
1. Clone the repository
``` bash
git clone https://github.com/malloryej/FNDIV_Final_Project.git
```

2. Open MLX file. Edit parameters in first cell.
	- NOTE: this script assumes a -200 to 800 milisecond ERP where -200 to 0 are treated as baseline. If using a different time scale, the linspace functions for data visualization will need to be updated to reflect the appropriate time window.

# Running the code
1. Raw data visualization includes average ERP waveforms with interquartile range shaded. Histograms are generated for hypothesized key timepoints for components N2 and P3. Script is set for N2 timepoints 200, 220, 250 and P3 timepoints 300, 375, and 450. This is editable. The code applies Sturge's formula to determine appropriate bin size. Histograms should show a normal distribution if there are no outliers.
2. Create a null distribution. The null distribution randomizes the data and searches for clusters. Timepoints are significant at a p-val <= 0.01. This is also editable (null_alpha in the parameters cell). The 95th percentile value is extracted and acts as the minimum cluster size against which to test true data clusters.
3. Generate final data visualization. Similar to the raw data visualization, this will generate an average ERP for both groups plotted on top of one another. The standard error is shaded around the lines and the time range for significant clusters is shaded in grey.

``` text
Project Structure │ 
├── data/ # Contains subfolders for each dataset with CSVs per participant. Also includes condition keys. 
├── datFiles_[DATASET]/ # Contains individual CSV files for each participant in a given experiment. 
├── MFiles/ # Independent .m files used in the MLX script (functions, processing scripts, etc.). 
├── CBPT_MLX.mlx # Main Live Script that runs the full cluster-based permutation testing pipeline. 
└── README.md # This file.
```