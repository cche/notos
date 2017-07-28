Notos
=====

Notos is a suite that calculates CpGo/e ratios for nucleotide sequences and uses Kernel Density Estimation to model the distribution obtained.

It consists of two programs, CpGoe.pl is used to prepare the CpG ratios and KDEanalysis.r estimates the model. 

CpGoe.pl
--------


CpGoe.pl -f fasta-file -a algo -c context -o output_file -m minlength
	
Reads multi fasta files, writes name, length, CpGs and GpCs, CpGo/e ratios and other quantities into TAB file for each sequence.

The available contexts (-c) are CpG, CpA, CpC, CpT. Default CpG. 

The algorithms (-a) available for -a are the following::

    1 => (CpG / (C * G)) * (L^2 / L-1)
    2 => (CpG / (C * G)) * L
    3 => (CpG / L) / ((C + G) / L)^2
    4 => (CpG / (C + G)/2)^2
		
	Where L represents the length of the sequence, CpG represents the count of CG dinucleotide, C and G represent the count for the respective bases and TG represents the number of TG dinucleotides.

Using the CpGo/e ratio to infer the type of methylation of a genome, whether it is global or mosaic.



KDEanalysis.r

-------------


Model the distribution of CpG o/e ratios using Kernel Density Estimation.

Parameters::

  "-o", "--frac-outl"			maximum fraction of CpGo/e ratios excluded as outliers [default 0.01]
  "-d", "--min-dist"			minimum distance between modes, modes that are closer are joined [default 0.2]
  "-c", "--conf-level"			level of the confidence intervals of the mode positions [default 0.95]
  "-m", "--mode-mass"			minimum probability mass of a mode [default 0.05]
  "-b", "--band-width"			bandwidth constant for kernels [default 1.06]
  "-B", "--bootstrap"			calculate confidence intervals of mode positions using bootstrap.
  "-r", "--bootstrap-reps"		number of bootstrap repetitions [default 1500]
  "-p", "--peak-file"			name of the output file describing the peaks of the KDE [default modes_basic_stats.cs]
  "-s", "--bootstrap-file"		Name of the output file with bootstrap values [default "modes_bootstrap.csv"]
  "-H", "--outlier-hist-file"	Outliers histogram file [default outliers_hist.pdf]
  "-C", "--cutoff-file"			Outliers cutoff file [default outliers_cutoff.csv]
  "-k", "--kde-file"			Kernel density estimation graph [default KDE.pdf]


