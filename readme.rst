Notos
=====

Notos is a suite that calculates CpN o/e ratios (e.g., the commonly used CpG o/e ratios) for a set of nucleotide sequences and uses Kernel Density Estimation (KDE) to model the obtained distribution.

It consists of two programs, CpGoe.pl is used to calculate the CpN o/e ratios and KDEanalysis.r estimates the model. 
In the following, these two programs are described briefly.

CpGoe.pl
--------


This program will calculate CpN o/e ratios on nucleotide multifasta files.
For each sequence that is found in the file it will output the sequence name followed by the CpN o/e ratio, where N can be any of the nucletides A, C, G or T, into a TAB separated file.

An example call would be:

    perl CpGoe.pl -f input_species.fasta -a 1 -c CpG -o input_species_cpgoe.csv -m 200
	

The available contexts (-c) are CpG, CpA, CpC, CpT. Default is CpG. 

The available algorithms (-a) for calculating the CpNo/e ratio are the following (here shown for CpG o/e)::

    1 => (CpG / (C * G)) * (L^2 / L-1)
    2 => (CpG / (C * G)) * L
    3 => (CpG / L) / ((C + G) / L)^2
    4 => (CpG / (C + G)/2)^2
		
Here L denotes the length of the sequence, CpG represents the count of CG dinucleotide, C and G represent the count for the respective bases.

KDEanalysis.r
-------------

This program carries out two steps.
First, the data preparation step, mainly to remove data artifacts.
Secondly, the mode detection step, which is baesd on a KDE modelling approach.

Example basic usage on command line:

    Rscript ~/src/github/notos/KDEanalysis.r "Input species" input_species_cpgoe.csv

	
In the above case "Input species" will be used to name the graphs that are generated as well as an identifier for each sample.
It has to be surrounded by " if the name of the species contains spaces.
The input of KDEanalysis.r is of the same format as the output of CpGoe.pl.

Any of the following parameters can be used

+--------+---------------------+-----------------------------------------------------------------------------------------+
| Option | Long option         | Description                                                                             |
+--------+---------------------+-----------------------------------------------------------------------------------------+
| -o     | --frac-outl         | maximum fraction of CpGo/e ratios excluded as outliers [default 0.01]                   |
+--------+---------------------+-----------------------------------------------------------------------------------------+
| -d     | --min-dist          | minimum distance between modes, modes that are closer are joined [default 0.2]          |
+--------+---------------------+-----------------------------------------------------------------------------------------+
| -c     | --conf-level        | level of the confidence intervals of the mode positions [default 0.95]                  |
+--------+---------------------+-----------------------------------------------------------------------------------------+
| -m     | --mode-mass         | minimum probability mass of a mode [default 0.05]                                       |
+--------+---------------------+-----------------------------------------------------------------------------------------+
| -b     | --band-width        | bandwidth constant for kernels [default 1.06]                                           |
+--------+---------------------+-----------------------------------------------------------------------------------------+
| -B     | --bootstrap         | calculate confidence intervals of mode positions using bootstrap.                       |
+--------+---------------------+-----------------------------------------------------------------------------------------+
| -r     | --bootstrap-reps    | number of bootstrap repetitions [default 1500]                                          |
+--------+---------------------+-----------------------------------------------------------------------------------------+
| -p     | --peak-file         | name of the output file describing the peaks of the KDE [default modes_basic_stats.csv] |
+--------+---------------------+-----------------------------------------------------------------------------------------+
| -s     | --bootstrap-file    | Name of the output file with bootstrap values [default "modes_bootstrap.csv"]           |
+--------+---------------------+-----------------------------------------------------------------------------------------+
| -H     | --outlier-hist-file | Outliers histogram file [default outliers_hist.pdf]                                     |
+--------+---------------------+-----------------------------------------------------------------------------------------+
| -C     | --cutoff-file       | Outliers cutoff file [default outliers_cutoff.csv]                                      |
+--------+---------------------+-----------------------------------------------------------------------------------------+
| -k     | --kde-file          | Kernel density estimation graph [default KDE.pdf]                                       |
+--------+---------------------+-----------------------------------------------------------------------------------------+

Of special interest is the -B parameter that will trigger the bootstrap calculations.
Default settings have been thoroughly calibrated through extensive testing, so we would advice to modify them only if you know what you are doing.

Output: Both the data preparation and the mode detection step return results in form of CSV files and figures to the user.
The two figures illustrate the results of the data cleaning and mode detection step, respectively.
The contents of the CSV files is described in the following.

1. outliers_cutoff.csv. The columns of this file contain

+---------------+----------------------------------------------------------------------------------------------------------------+
| Column        | description                                                                                                    |
+---------------+----------------------------------------------------------------------------------------------------------------+
| Name          | name of the file analyzed                                                                                      |
+---------------+----------------------------------------------------------------------------------------------------------------+
| prop.zero     | proportion of observations equal to zero excluded (relative to original sample)                                |
+---------------+----------------------------------------------------------------------------------------------------------------+
| prop.out.2iqr | proportion of values equal excluded if 2 * IQR was used, relative to sample after exclusion of zeros (0 - 100) |
+---------------+----------------------------------------------------------------------------------------------------------------+
| prop.out.3iqr | proportion of values equal excluded if 3 * IQR was used, relative to sample after exclusion of zeros (0 - 100) |
+---------------+----------------------------------------------------------------------------------------------------------------+
| prop.out.4iqr | proportion of values equal excluded if 4 * IQR was used, relative to sample after exclusion of zeros (0 - 100) |
+---------------+----------------------------------------------------------------------------------------------------------------+
| prop.out.5iqr | proportion of values equal excluded if 5 * IQR was used, relative to sample after exclusion of zeros (0 - 100) |
+---------------+----------------------------------------------------------------------------------------------------------------+
| used          | IQR used for exclusion of outliers / extreme values                                                            |
+---------------+----------------------------------------------------------------------------------------------------------------+
| no.obs.raw    | number of observations in the original sample                                                                  |
+---------------+----------------------------------------------------------------------------------------------------------------+
| no.obs.nozero | number of observations in sample after excluding values equal to zero                                          |
+---------------+----------------------------------------------------------------------------------------------------------------+
| no.obs.clean  | number of observations in sample after excluding outliers / extreme values                                     |
+---------------+----------------------------------------------------------------------------------------------------------------+

2. modes_basic_stats.csv. We use the following notation: sigma - standard deviation, mu - mean, nu - median, Mo - mode, Q_i - the i-th quartile, q_s - the s % quantile. The columns of this file contain

+-----------------------------------+------------------------------------------------------------------------------------+
| Column                            | description                                                                        |
+-----------------------------------+------------------------------------------------------------------------------------+
| Name                              | name of the file analyzed                                                          |
+-----------------------------------+------------------------------------------------------------------------------------+
| Number of modes                   | number of modes without applying any exclusion criterion                           |
+-----------------------------------+------------------------------------------------------------------------------------+
| Number of modes (5% excluded)     | number of modes after exclusion of those with less then 5% probability mass        |
+-----------------------------------+------------------------------------------------------------------------------------+
| Number of modes (10% excluded)    | number of modes after exclusion of those with less then 10% probability mass       |
+-----------------------------------+------------------------------------------------------------------------------------+
| Skewness                          | Pearson's moment coefficient of skewness E(X-mu/sigma)^3                           |
+-----------------------------------+------------------------------------------------------------------------------------+
| Mode skewness                     | Pearson's first skewness coefficient (mu - Mo)/sigma                               |
+-----------------------------------+------------------------------------------------------------------------------------+
| Nonparametric skew                | (mu - nu)/sigma                                                                    |
+-----------------------------------+------------------------------------------------------------------------------------+
| Q50 skewness                      | Bowley's measure of skewness / Yule's coefficient (Q_3 + Q_1 - 2Q_2) / (Q_3 - Q_1) |
+-----------------------------------+------------------------------------------------------------------------------------+
| Absolute Q50 mode skewness        | (Q_3 + Q_1) / 2 - Mo                                                               |
+-----------------------------------+------------------------------------------------------------------------------------+
| Absolute Q80 mode skewness        | (q_90 + q_10) / 2 - Mo                                                             |
+-----------------------------------+------------------------------------------------------------------------------------+
| Peak i, i = 1,..., 10             | location of peak i                                                                 |
+-----------------------------------+------------------------------------------------------------------------------------+
| Probability Mass i, i = 1,..., 10 | probability mass assigned to peak i                                                |
+-----------------------------------+------------------------------------------------------------------------------------+
| Warning close modes               | flag indicating that modes lie too close. The default threshold is 0.2             |
+-----------------------------------+------------------------------------------------------------------------------------+
| Number close modes                | number of modes lying too close, given the threshold                               |
+-----------------------------------+------------------------------------------------------------------------------------+
| Modes (close modes excluded)      | number of modes after exclusion of modes that are too close                        |
+-----------------------------------+------------------------------------------------------------------------------------+
| SD                                | sample standard deviation sigma                                                    |
+-----------------------------------+------------------------------------------------------------------------------------+
| IQR 80                            | 80% distance between the 90 % and 10 % quantile                                    |
+-----------------------------------+------------------------------------------------------------------------------------+
| IQR 90                            | 90% distance between the 95 % and 5 % quantile                                     |
+-----------------------------------+------------------------------------------------------------------------------------+
| Total number of sequences         | total number of sequences / CpG o/e ratios used for this analysis step             |
+-----------------------------------+------------------------------------------------------------------------------------+

3. modes_bootstrap.csv. The columns of this optional file resulting from the bootstrap procedure contains:

+-------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------+
| Column                                    | description                                                                                                                                |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------+
| Name                                      | name of the file analyzed                                                                                                                  |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------+
| Number of modes (NM)                      | number of modes detected for the original sample                                                                                           |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------+
| % of samples with same NM                 | proportion of bootstrap samples with the same number of modes (0 - 100)                                                                    |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------+
| % of samples with more NM                 | proportion of bootstrap samples a higher number of modes (0 - 100)                                                                         |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------+
| % of samples with less NM                 | proportion of bootstrap samples a lower number of modes (0 - 100)                                                                          |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------+
| no. of samples with same NM               | number of bootstrap samples with the same number of modes                                                                                  |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------+
| % BS samples excluded by prob.~mass crit. | proportion of bootstrap samples excluded due to strong deviations from the probability masses determined for the original sample (0 - 100) |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------+
