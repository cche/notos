#! /usr/bin/perl

####################################
#
# CpGoe.pl -f fasta-file -o output_file -m minlength
#
# reads concatanated fasta files, writes name, length, CpGs and GpCs, CpGo/e ratios and other quantities into TAB file for each sequence
#####################################

use diagnostics;
use strict;
use Carp;
use FileHandle;
use File::Path;
use File::Basename;
use Data::Dumper;
use Getopt::Std;

my $VERSION = "1.0";
$Getopt::Std::STANDARD_HELP_VERSION = 1;


# Called when --help set as flag
sub HELP_MESSAGE
{
  print "Description: CpGoe processes a FASTA file and outputs the CpGo/e ratios and - if specified - further quantities\n\n" .
        "Usage: CpGoe.pl [OPTION] [FASTA_FILE] \n\n" .
        "Options:\n" .
        "-o OUT_FILE\n" .
        "      name of output file containing the results\n" .
        "-m MIN_LEN\n" .
        "      minimum length of sequences, shorter sequences are discarded\n" .
        "-d\n" .
        "      detailed output, providing other quantities additional to the CpGo/e ratios\n" .
        "-h\n" .
        "      output header line\n";
        "-v\n" .
        "      verbose messages\n" .
  exit 0;
}



# Called when --version set as flag
sub VERSION_MESSAGE
{
  print "CpGoe $VERSION\n";
}



# Command line parsing
# ... read argument
my %opts;
getopts('f:o:m:dvh', \%opts);
#if ($#ARGV != 0) {
#  print STDERR "Exactly one argument has to be provided, the name of the input FASTA file.\n" .
#               "Moreover, the options must be listed first, then the name of the input FASTA file.\n";
#  exit 1;
#}
my $fasta_fname = $opts{'f'};

# ... read options
my $out_fname;
my $has_file_output;
if (exists($opts{'o'})) {
  $out_fname = $opts{'o'};
  $has_file_output = 1;
} 
else {
  $has_file_output = 0;
}

my $min_len; 
if (exists($opts{'m'})) {
  $min_len = $opts{'m'}
}
else {
  $min_len = 1;
}

my $is_verbose = exists($opts{'v'});
my $has_header = exists($opts{'h'});
my $is_detailed = exists($opts{'d'});



# read input file and split into fasta sequences on the fly
# ... check whether input FASTA file exists
my $FASTA;
if (-e $fasta_fname) {
  if (-f $fasta_fname) {
    my $res = open($FASTA, $fasta_fname);
    if (!$res) {
      print STDERR "could not open $fasta_fname\n";
      exit 1;
    }
  } 
  else {
    print STDERR "$fasta_fname is not a file\n";
    exit 1;
  }
}
else {
  print STDERR "cannot open file $fasta_fname\n";
  exit 1;
}


# ... determine which linebreak is used (linux / windows / mac)
my $found_n = 0;
my $found_r = 0;
my $found_rn = 0;
while ( defined( my $ch = getc($FASTA) ) ) {
  if ($ch eq "\n") {
    if ($found_r) {
      $found_rn = 1;
      $found_r = 0;
    } else {
      $found_n = 1;
    }
    last;
  } elsif ($ch eq "\r") {
    $found_r = 1;
  } else {
    if ($found_r) {
      last;
    }
  }
}
close($FASTA);
if ($found_r + $found_n + $found_rn != 1) {
  print STDERR "something went wrong determining the linebreaks used in $fasta_fname\n";
}


# ... read in sequences
my $old_linebreak = $/;
if ($found_n) {
  $/ = "\n";
} elsif ($found_r) {
  $/ = "\r";
} else {
  $/ = "\r\n";
}
my $res = open($FASTA, $fasta_fname);
if (!$res) {
  print STDERR "could not open $fasta_fname\n";
  exit 1;
}

my @names = ();   # names of the sequences
my @seqs = ();   # sequences
my $is_first = 1;
while ( <$FASTA> ) {
  chomp;
  if (/^[^#]/) {
    if ( /^>(\S+)/) {
      s/^>|\s+$//g;   # remove leading '>' and trailing whitespaces
      s/\s/_/g;   # replace spaces by underscores
      push(@names, $_);
      push(@seqs, "");
      $is_first = 0;
    } 
    else {
      if ($is_first) {
        print STDERR "first non-comment line of FASTA file " . $fasta_fname . " does not start with \'>\'\n";
        exit 1;
      }
      else {
        s/[\-\.]*//g;   # remove dashes and dots
        $seqs[-1] .= $_;
      }
    }
  }
}
$res = close($FASTA);
if (!$res) {
  print STDERR "could not close $fasta_fname\n";
  exit 1;
}
$/ = $old_linebreak;

# print Dumper(@names) . "\n";
# print Dumper(@seqs) . "\n";


# ... check sequences
# ... ... are there any sequence names?
if ($#names < 0) {
  print STDERR "FASTA file $fasta_fname is empty\n";
  exit 1;
}

# ... ... are there empty sequences?}
my $str = "";
my $err = 0;
my $num = 0;
my $MAX = 50;   # maximum number of notifications about an empty sequence
for (my $i = 0; $i <= $#names; ++$i) {
  if ($seqs[$i] eq "") {
    if ($num < $MAX) {
      $str .= "Sequence " . $names[$i] . " in FASTA file $fasta_fname is empty\n";
    }
    $err = 1;
    ++$num;     
  }
}
if ($err) {
  print STDERR "$str";
  if ($num > $MAX) {
    print STDERR "$num empty sequences in total in FASTA file $fasta_fname \n";
  }
  exit 1;
}


# ... ... check for illegal characters in sequences
for (my $i = 0; $i <= $#names; ++$i) {
  my $str = $seqs[$i];
  $str =~ s/[abcdghkmnrstuvwy]//gi;
  if (length($str) > 0) {
    $str = join '', sort { $a cmp $b } split(//, $str);
    print STDERR "Sequence " . $names[$i] . " of FASTA file " . $fasta_fname . " contains illegal characters: ";
    for (my $j = 0; $j <= length($str); ++$j) {
      my $out = ($j == 0);
      my $curr = substr($str, $j, 1);
      if (!$out) {
        my $prev = substr($str, $j - 1, 1) ;
        $out = ($curr ne $prev);
      }
      if ($out) {
        print STDERR $curr;
      }
    }
    print STDERR "\n";
    exit 1;
  }
}


# ... output
if ($is_verbose) {
  print $#names + 1 . " sequence(s) read.\n";
}



# output quantities
# ... open output file
my $OUT;
if ($has_file_output) {
  if ((-e $out_fname) && !(-f $out_fname)) {
    print STDERR "$out_fname exists and is not a file\n";
    exit 1;
  }
  if (!open($OUT, ">$out_fname")) {
    print STDERR "cannot create file $out_fname\n";
    exit 1;
  }
} else {
  $OUT = *STDOUT;
}

# ... print header
if ($has_header) {
  print $OUT "#name\tlength\tCpGs\tGpCs\tCs\tGs\tNs\tCpG o\/e\n"; 
}



# ... for each sequence calculate CpGo/e ratios and related quantities:
# - length of the sequence
# - CpGs present in the sequence
# - GpCs present in the sequence
# - Cs the number of C present in the sequence
# - Gs the number of G present in the sequence
# - CpG o/e ratio of the sequence
my $num_short = 0;   # number of sequences which are too short
for my $i (0 .. $#names) {
  my @ar = ();
  @ar = split('\|', $names[$i]);
  my $seqname = $ar[1];
  my $num_N = () = ( $seqs[$i] =~ m/N/gi );
  my $len = length($seqs[$i]);
  my $l = $len - $num_N;
  if ($l >= $min_len) { 
    my $num_G = () = ( $seqs[$i] =~ m/G/gi );
    my $num_C = () = ( $seqs[$i] =~ m/C/gi );
    my $num_CG = () = ( $seqs[$i] =~ m/CG/gi );
    my $num_GC = () = ( $seqs[$i] =~ m/GC/gi );
    my $CpGoe;
    if ( ($num_G == 0) || ($num_C == 0) || ($l == 1) ) { 
      $CpGoe = 0;
    }   
    else {
      my $x = $num_CG / ($num_C * $num_G);
      my $y = $l**2 / ($l - 1);
      $CpGoe = $x * $y;
    }  
    print $OUT $names[$i] . "\t";
    if ($is_detailed) {
      print $OUT $len . "\t" . $num_CG . "\t" . $num_GC . "\t" .$num_C. "\t" .$num_G. "\t" .$num_N. "\t";
    }
    print $OUT $CpGoe . "\n";
  } else {
    ++$num_short;
  }
}
if ($is_verbose) {
  print $num_short . " sequence(s) discarded due to short length.\n";
}

if ($has_file_output) {
  my $res = close($OUT);
  if (!$res) {
    print STDERR "could not close $out_fname\n";
    exit 1;
  }
}

