Notos
=====

Notos is a suite that calculates CpGo/e ratios for nucleotide sequences and uses Kernel Density Estimation to model the obtained distribution.

It consists of two programs, CpGoe.pl is used to prepare the CpG ratios and KDEanalysis.r estimates the model. 

CpGoe.pl
--------


This program will calculate CpN o/e ratios on nucleotide multifasta files. For each sequence that is found in the file it will output the sequence name followed by the CpN o/e ratio, where N can be any of the nucletides A, C, G or T, into a TAB separated file.

An example invocation would be:

    perl CpGoe.pl -f input_species.fasta -a 1 -c CpG -o input_species_cpgoe.csv -m 200
	

The available contexts (-c) are CpG, CpA, CpC, CpT. Default CpG. 

The available algorithms (-a) are the following::

    1 => (CpG / (C * G)) * (L^2 / L-1)
    2 => (CpG / (C * G)) * L
    3 => (CpG / L) / ((C + G) / L)^2
    4 => (CpG / (C + G)/2)^2
		
Where L represents the length of the sequence, CpG represents the count of CG dinucleotide, C and G represent the count for the respective bases.


KDEanalysis.r
-------------

Model the distribution of CpG o/e ratios using Kernel Density Estimation.

Example basic usage on command line:

    Rscript ~/src/github/notos/KDEanalysis.r "Input species" input_species_cpgoe.csv


In the above case "Input species" will be used to name the graphs that are produced as well as an identifier for each sample. It has to be surrounded by " if the name of the species contains spaces.

Any of the following parameters can be used::

  "-o", "--frac-outl"			maximum fraction of CpGo/e ratios excluded as outliers [default 0.01]
  "-d", "--min-dist"			minimum distance between modes, modes that are closer are joined [default 0.2]
  "-c", "--conf-level"			level of the confidence intervals of the mode positions [default 0.95]
  "-m", "--mode-mass"			minimum probability mass of a mode [default 0.05]
  "-b", "--band-width"			bandwidth constant for kernels [default 1.06]
  "-B", "--bootstrap"			calculate confidence intervals of mode positions using bootstrap.
  "-r", "--bootstrap-reps"		number of bootstrap repetitions [default 1500]
  "-p", "--peak-file"			name of the output file describing the peaks of the KDE [default modes_basic_stats.csv]
  "-s", "--bootstrap-file"		Name of the output file with bootstrap values [default "modes_bootstrap.csv"]
  "-H", "--outlier-hist-file"	Outliers histogram file [default outliers_hist.pdf]
  "-C", "--cutoff-file"			Outliers cutoff file [default outliers_cutoff.csv]
  "-k", "--kde-file"			Kernel density estimation graph [default KDE.pdf]

Of special interest is the -B parameter that will trigger the bootstrap calculations.
