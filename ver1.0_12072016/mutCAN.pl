#!/usr/local/bin/perl -w
# Bookmarks: 0,0 0,0 0,0 16,117 0,0 0,0 0,0 0,0 0,0 16,117
#  Program Description
#    This program processes the alignment of query files of mutant & wild-type libraries to reference genome,
#  exports the coverage depth & nucleotide frequency information for each single genomic position, and calls
#  candidate mutations by comparing the difference between the two libraries, and also provides reads at the
#  candidate sites for further confirmation.

#  Platform & Requirement
#    Platform: Linux system
#    Requirement: Bowtie

#  Contact & Version
#    Author: Jinfeng Lu, jinfeng.lu@ucr.edu
#    Contact author: Shou-wei Ding, shou-wei.ding@ucr.edu
#    Version: 1.0, Date: 08-23-2016

#  Usage
#    $ perl mutCAN.pl --generation [options]             alignment & generate files of coverage depth & nucleotide frequency (.can) in /CAN/ folder
#    $ perl mutCAN.pl --comparison [options]             call mutation candidates by comparing the two .can files in /CAN/ folder
#    $ perl mutCAN.pl --confirmation [options]           retrieve reads at the candidate site
#    $ perl mutCAN.pl -h/--help                          print this usage message

#  Example
#    $ perl mutCAN.pl --generation -r TAIR.fa -v 3 -c Chr2 -m mutant.fastq -w wtpair1.fastq wtpair2.fastq --wX 600
#    $ perl mutCAN.pl --comparison -a TAIR.gff --cor correction.txt --dmc 0.1 --dwc 15 --smc 10 --swc 11 --swf 0.65 0.67 --smf 1 0.99 --sjud 0
#    $ perl mutCAN.pl --confirmation -r TAIR.fa -m mutant.fastq -w wtpair1.fastq wtpair2.fastq -c Chr2 --reg 10080 10080 -f 35

#  Command-line Option for --generation
#    -r/--reference <file>                    genome reference file (fasta format) for alignment
#    -m/--mutant <file>                       single-read file (fastq format) of the mutant library
#              or
#    -m/--mutant <file1> <file2>              space-separated paired-end files (fw/rev direction, fastq format) of the mutant library
#    --mX <int>                               maximum insert size for paired-end alignment for the mutant library (default: 250)
#    -w/--wildtype <file>                     single-read file (fastq format) of the wild-type library
#              or
#    -w/--wildtype <file1> <file2>            space-separated paired-end files (fw/rev direction, fastq format) of the wild-type library
#    --wX <int>                               maximum insert size for paired-end alignment for the wild-type library (default: 250)
#    -c/--chromosome <chrID>                  specify the chromosome <chrID> to analyze (default: export information for all chromosomes)
#                                             <chrID> must be identical to the name in reference file (found after '>' or '>gi|')
#    --reg <int1> <int2>                      specify the genomic position/region to analyze if <chrID> is specified, offset=1 (default: export information for full length)
#    -v <int>                                 set the maximum mismatch number (0-3) for alignment (default: 2)
#    --dr <int>                               divide the query file by <int> reads per subfile to reduce memory usage (default: 100000)
#    --dm <int>                               divide the alignment result by <int> genomic positions per block to reduce memory usage (default: 5000000)

#  Command-line Option for --comparison
#    -a/--annotation <file>                   genome annotation file (gff format) for candidate annotation (default: no annotation will be made)
#    --cor <file>                             a file containing lines to correct those chromosome names if they are different between the reference & the annotation file
#                                             format for each line: <chrname_annotation> <chrname_reference>
#    --smc <dec>                              minimum coverage depth in mutant library for candidate SNP calling (default: candidate SNPs will not be called)
#    --smf <dec1> <dec2>                      set the frequency range for top nucleotide in mutant library for candidate SNP calling (default: candidate SNPs will not be called)
#    --swc <dec>                              minimum coverage depth in wild-type library for candidate SNP calling (default: candidate SNPs will not be called)
#    --swf <dec1> <dec2>                      set the frequency range for top nucleotide in wild-type library for candidate SNP calling (default: candidate SNPs will not be called)
#    --sjud <int>                             the top nucleotides from the two libraries are different (<int> =0) or same (<int> !=0) (default: candidate SNPs will not be called)
#    --dmc <dec>                              maximum coverage depth in mutant library for candidate deletion calling (default: candidate deletions will not be called)
#    --dwc <dec>                              minimum coverage depth in wild-type library for candidate deletion calling (default: candidate deletions will not be called)
#    --imc <dec>                              maximum coverage depth in mutant library for candidate insertion calling (default: candidate insertions will not be called)
#    --iwc <dec>                              minimum coverage depth in wild-type library for candidate insertion calling (default: candidate insertions will not be called)

#  Command-line Option for --confirmation
#    -r/--reference <file>                    genome reference file (fasta format)
#    -m/--mutant <file>                       single-read file (fastq format) of the mutant library
#              or
#    -m/--mutant <file1> <file2>              space-separated paired-end files (fastq format) of the mutant library
#    -w/--wildtype <file>                     single-read file (fastq format) of the wild-type library
#              or
#    -w/--wildtype <file1> <file2>            space-separated paired-end files (fastq format) of the wild-type library
#    -c/--chromosome <chrID>                  the chromosome <chrID> containing the candidate site
#                                             <chrID> must be identical to the name in reference file (found after '>' or '>gi|')
#    --reg <int1> <int2>                      specify the candidate site as a single nucleotide position or a region (offset=1)
#                                             if the region size > 100bp, reads at either end will be shown separately
#    -d/--distance <int>                      set the distance (0-50) between the candidate site & the flanking sequences to retrieve reads (default: 30)
#    -f/--flank <int>                         set the length (10-50) of the flanking sequences to retrieve reads (default: 25)

use strict;
my $error="no";
my $help="no";
my $command="";
my $i=0;
my $setcount=0;
my %set=();
my $modcount=0;
my @setload=();
my @setlist=("-r","--reference","-m","--mutant","-w","--wildtype","-c","--chromosome","--reg","-d","--distance","-f","--flank","--mX","--wX","-v","--dr","--dm","-a","--annotation","--cor","--smc","--smf","--swc","--swf","--sjud","--dmc","--dwc","--imc","--iwc");
my $isset="";
my @tempdata=();
my $temp="";
my $missing;
if ($#ARGV==-1) {
      $help="yes";
} else {
      for ($i=0;$i<=$#ARGV;$i++) {
            if (($ARGV[$i] eq "-h") or ($ARGV[$i] eq "--help")) {
                  $help="yes";
            } elsif (($ARGV[$i] eq "--generation") or ($ARGV[$i] eq "--comparison") or ($ARGV[$i] eq "--confirmation")) {
                  $set{"module"}=$ARGV[$i];
                  $modcount++;
            } else {
                  $isset="no";
                  foreach (@setlist) {
                       if ($_ eq $ARGV[$i]) {
                            $isset="yes";
                       }
                  }
                  if ($isset eq "yes") {
                       $setcount++;
                       $setload[$setcount*2-1]=$ARGV[$i];
                  } else {
                       if (!(defined($setload[$setcount*2]))) {
                            $setload[$setcount*2]=$ARGV[$i];
                       } else {
                            $setload[$setcount*2]=join " ",($setload[$setcount*2],$ARGV[$i]);
                       }
                  }
            }
      }
}
if (($help eq "yes") or ($modcount==0)) {
     print "Usage\n";
     print "  \$ perl mutCAN.pl --generation [options]             alignment & generate files of coverage depth & nucleotide frequency (.can) in /CAN/ folder\n";
     print "  \$ perl mutCAN.pl --comparison [options]             call mutation candidates by comparing the two .can files in /CAN/ folder\n";
     print "  \$ perl mutCAN.pl --confirmation [options]           retrieve reads at the candidate site\n";
     print "  \$ perl mutCAN.pl -h/--help                          print this usage message\n";
     print "\n";
     print "Example\n";
     print "  \$ perl mutCAN.pl --generation -r TAIR.fa -v 3 -c Chr2 -m mutant.fastq -w wtpair1.fastq wtpair2.fastq --wX 600\n";
     print "  \$ perl mutCAN.pl --comparison -a TAIR.gff --cor correction.txt --dmc 0.1 --dwc 15 --smc 10 --swc 11 --swf 0.65 0.67 --smf 1 0.99 --sjud 0\n";
     print "  \$ perl mutCAN.pl --confirmation -r TAIR.fa -m mutant.fastq -w wtpair1.fastq wtpair2.fastq -c Chr2 --reg 10080 10080 -f 35\n";
     print "\n";
     print "Command-line Option for --generation\n";
     print "  -r/--reference <file>                    genome reference file (fasta format) for alignment\n";
     print "  -m/--mutant <file>                       single-read file (fastq format) of the mutant library\n";
     print "            or\n";
     print "  -m/--mutant <file1> <file2>              space-separated paired-end files (fw/rev direction, fastq format) of the mutant library\n";
     print "  --mX <int>                               maximum insert size for paired-end alignment for the mutant library (default: 250)\n";
     print "  -w/--wildtype <file>                     single-read file (fastq format) of the wild-type library\n";
     print "            or\n";
     print "  -w/--wildtype <file1> <file2>            space-separated paired-end files (fw/rev direction, fastq format) of the wild-type library\n";
     print "  --wX <int>                               maximum insert size for paired-end alignment for the wild-type library (default: 250)\n";
     print "  -c/--chromosome <chrID>                  specify the chromosome <chrID> to analyze (default: export information for all chromosomes)\n";
     print "                                           <chrID> must be identical to the name in reference file (found after '>' or '>gi|')\n";
     print "  --reg <int1> <int2>                      specify the genomic position/region to analyze if <chrID> is specified, offset=1 (default: export information for full length)\n";
     print "  -v <int>                                 set the maximum mismatch number (0-3) for alignment (default: 2)\n";
     print "  --dr <int>                               divide the query file by <int> reads per subfile to reduce memory usage (default: 100000)\n";
     print "  --dm <int>                               divide the alignment result by <int> genomic positions per block to reduce memory usage (default: 5000000)\n";
     print "\n";
     print "Command-line Option for --comparison\n";
     print "  -a/--annotation <file>                   genome annotation file (gff format) for candidate annotation (default: no annotation will be made)\n";
     print "  --cor <file>                             a file containing lines to correct those chromosome names if they are different between the reference & the annotation file\n";
     print "                                           format for each line: <chrname_annotation> <chrname_reference>\n";
     print "  --smc <dec>                              minimum coverage depth in mutant library for candidate SNP calling (default: candidate SNPs will not be called)\n";
     print "  --smf <dec1> <dec2>                      set the frequency range for top nucleotide in mutant library for candidate SNP calling (default: candidate SNPs will not be called)\n";
     print "  --swc <dec>                              minimum coverage depth in wild-type library for candidate SNP calling (default: candidate SNPs will not be called)\n";
     print "  --swf <dec1> <dec2>                      set the frequency range for top nucleotide in wild-type library for candidate SNP calling (default: candidate SNPs will not be called)\n";
     print "  --sjud <int>                             the top nucleotides from the two libraries are different (<int> =0) or same (<int> !=0) (default: candidate SNPs will not be called)\n";
     print "  --dmc <dec>                              maximum coverage depth in mutant library for candidate deletion calling (default: candidate deletions will not be called)\n";
     print "  --dwc <dec>                              minimum coverage depth in wild-type library for candidate deletion calling (default: candidate deletions will not be called)\n";
     print "  --imc <dec>                              maximum coverage depth in mutant library for candidate insertion calling (default: candidate insertions will not be called)\n";
     print "  --iwc <dec>                              minimum coverage depth in wild-type library for candidate insertion calling (default: candidate insertions will not be called)\n";
     print "\n";
     print "Command-line Option for --confirmation\n";
     print "  -r/--reference <file>                    genome reference file (fasta format)\n";
     print "  -m/--mutant <file>                       single-read file (fastq format) of the mutant library\n";
     print "            or\n";
     print "  -m/--mutant <file1> <file2>              space-separated paired-end files (fastq format) of the mutant library\n";
     print "  -w/--wildtype <file>                     single-read file (fastq format) of the wild-type library\n";
     print "            or\n";
     print "  -w/--wildtype <file1> <file2>            space-separated paired-end files (fastq format) of the wild-type library\n";
     print "  -c/--chromosome <chrID>                  the chromosome <chrID> containing the candidate site\n";
     print "                                           <chrID> must be identical to the name in reference file (found after '>' or '>gi|')\n";
     print "  --reg <int1> <int2>                      specify the candidate site as a single nucleotide position or a region (offset=1)\n";
     print "                                           if the region size > 100bp, reads at either end will be shown separately\n";
     print "  -d/--distance <int>                      set the distance (0-50) between the candidate site & the flanking sequences to retrieve reads (default: 30)\n";
     print "  -f/--flank <int>                         set the length (10-50) of the flanking sequences to retrieve reads (default: 25)\n";
} elsif ($modcount>1) {
     print "Error: multiple --generation/--comparison/--confirmation\n";
} else {
     for ($i=1;$i<=$setcount;$i++) {
           if ($setload[$i*2-1] eq "--reference") {
                $setload[$i*2-1]="-r";
           }
           if ($setload[$i*2-1] eq "--mutant") {
                $setload[$i*2-1]="-m";
           }
           if ($setload[$i*2-1] eq "--wildtype") {
                $setload[$i*2-1]="-w";
           }
           if ($setload[$i*2-1] eq "--chromosome") {
                $setload[$i*2-1]="-c";
           }
           if ($setload[$i*2-1] eq "--annotation") {
                $setload[$i*2-1]="-a";
           }
           if ($setload[$i*2-1] eq "--distance") {
                $setload[$i*2-1]="-d";
           }
           if ($setload[$i*2-1] eq "--flank") {
                $setload[$i*2-1]="-f";
           }
           if (defined($set{$setload[$i*2-1]})) {
                print "Error: duplicate $setload[$i*2-1]\n";
                $error="yes";
           } elsif (defined($setload[$i*2])) {
                $set{$setload[$i*2-1]}=$setload[$i*2];
           }
     }
     if ($error eq "no") {
           if ($set{"module"} eq "--generation") {
                $command="perl can_generation.pl";
                if (defined($set{"-r"})) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"-r"};
                    $command=join " ",($command,$tempdata[0]);
                } else {
                    $error="yes";
                    if (defined($missing)) {
                          $missing=join " ",($missing,"-r");
                    } else {
                          $missing="-r";
                    }
                }
                if (defined($set{"-v"})) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"-v"};
                    if (($tempdata[0]<0) or ($tempdata[0]>3)) {
                          print "-v out of range, using default value\n";
                          $command=join " ",($command,2);
                    } else {
                          $command=join " ",($command,int($tempdata[0]));
                    }
                } else {
                    $command=join " ",($command,2);
                }
                if (defined($set{"--dr"})) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"--dr"};
                    $command=join " ",($command,int($tempdata[0]));
                } else {
                    $command=join " ",($command,100000);
                }
                if (defined($set{"--dm"})) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"--dm"};
                    $command=join " ",($command,int($tempdata[0]));
                } else {
                    $command=join " ",($command,5000000);
                }
                if (defined($set{"-m"})) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"-m"};
                    if ($#tempdata==0) {
                          $command=join " ",($command,"s",$tempdata[0]);
                    } else {
                          $command=join " ",($command,"p",$tempdata[0],$tempdata[1]);
                          if (defined($set{"--mX"})) {
                              @tempdata=();
                              @tempdata=split /\s+/, $set{"--mX"};
                              $command=join " ",($command,int($tempdata[0]));
                          } else {
                              $command=join " ",($command,250);
                          }
                    }
                } else {
                    $error="yes";
                    if (defined($missing)) {
                          $missing=join " ",($missing,"-m");
                    } else {
                          $missing="-m";
                    }
                }
                if (defined($set{"-w"})) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"-w"};
                    if ($#tempdata==0) {
                          $command=join " ",($command,"s",$tempdata[0]);
                    } else {
                          $command=join " ",($command,"p",$tempdata[0],$tempdata[1]);
                          if (defined($set{"--wX"})) {
                              @tempdata=();
                              @tempdata=split /\s+/, $set{"--wX"};
                              $command=join " ",($command,int($tempdata[0]));
                          } else {
                              $command=join " ",($command,250);
                          }
                    }
                } else {
                    $error="yes";
                    if (defined($missing)) {
                          $missing=join " ",($missing,"-w");
                    } else {
                          $missing="-w";
                    }
                }
                if (defined($set{"-c"})) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"-c"};
                    $command=join " ",($command,$tempdata[0]);
                    if (defined($set{"--reg"})) {
                         @tempdata=();
                         @tempdata=split /\s+/, $set{"--reg"};
                         if ($#tempdata==0) {
                               $command=join " ",($command,int($tempdata[0]),int($tempdata[0]));
                         } elsif ($tempdata[0]<$tempdata[1]) {
                               $command=join " ",($command,int($tempdata[0]),int($tempdata[1]));
                         } else {
                               $command=join " ",($command,int($tempdata[1]),int($tempdata[0]));
                         }
                    }
                }
           } elsif ($set{"module"} eq "--comparison") {
                $command="perl can_comparison.pl";
                if (defined($set{"-a"})) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"-a"};
                    $command=join " ",($command,"a",$tempdata[0]);
                    if (defined($set{"--cor"})) {
                        @tempdata=();
                        @tempdata=split /\s+/, $set{"--cor"};
                        $command=join " ",($command,$tempdata[0]);
                    }
                }
                if ((defined($set{"--dmc"})) and (defined($set{"--dwc"}))) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"--dmc"};
                    $command=join " ",($command,"d",$tempdata[0]);
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"--dwc"};
                    $command=join " ",($command,$tempdata[0]);
                }
                if ((defined($set{"--imc"})) and (defined($set{"--iwc"}))) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"--imc"};
                    $command=join " ",($command,"i",$tempdata[0]);
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"--iwc"};
                    $command=join " ",($command,$tempdata[0]);
                }
                if ((defined($set{"--smc"}))and(defined($set{"--smf"}))and(defined($set{"--swc"}))and(defined($set{"--swf"}))and(defined($set{"--sjud"}))) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"--smc"};
                    $command=join " ",($command,"m",$tempdata[0]);
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"--smf"};
                    if ($#tempdata==0) {
                          $command=join " ",($command,$tempdata[0],$tempdata[0]);
                    } elsif ($tempdata[0]<$tempdata[1]) {
                          $command=join " ",($command,$tempdata[0],$tempdata[1]);
                    } else {
                          $command=join " ",($command,$tempdata[1],$tempdata[0]);
                    }
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"--swc"};
                    $command=join " ",($command,$tempdata[0]);
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"--swf"};
                    if ($#tempdata==0) {
                          $command=join " ",($command,$tempdata[0],$tempdata[0]);
                    } elsif ($tempdata[0]<$tempdata[1]) {
                          $command=join " ",($command,$tempdata[0],$tempdata[1]);
                    } else {
                          $command=join " ",($command,$tempdata[1],$tempdata[0]);
                    }
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"--sjud"};
                    if ($tempdata[0]==0) {
                          $command=join " ",($command,"different");
                    } else {
                          $command=join " ",($command,"same");
                    }
                }
           } elsif ($set{"module"} eq "--confirmation") {
                $command="perl can_confirmation.pl";
                if (defined($set{"-r"})) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"-r"};
                    $command=join " ",($command,$tempdata[0]);
                } else {
                    $error="yes";
                    if (defined($missing)) {
                          $missing=join " ",($missing,"-r");
                    } else {
                          $missing="-r";
                    }
                }
                if (defined($set{"-c"})) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"-c"};
                    $command=join " ",($command,$tempdata[0]);
                } else {
                    $error="yes";
                    if (defined($missing)) {
                          $missing=join " ",($missing,"-c");
                    } else {
                          $missing="-c";
                    }
                }
                if (defined($set{"--reg"})) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"--reg"};
                    if ($#tempdata==0) {
                          $command=join " ",($command,int($tempdata[0]),int($tempdata[0]));
                    } else {
                          $tempdata[2]=$tempdata[1]-$tempdata[0];
                          if (($tempdata[2]>100) or ($tempdata[2]<-100)) {
                               $command=join " ",($command,int($tempdata[0]),int($tempdata[0]));
                               $set{"module"}="repeat";
                               $set{"--reg"}=$tempdata[1];
                          } elsif ($tempdata[2]>0) {
                               $command=join " ",($command,int($tempdata[0]),int($tempdata[1]));
                          } else {
                               $command=join " ",($command,int($tempdata[1]),int($tempdata[0]));
                          }
                    }
                } else {
                    $error="yes";
                    if (defined($missing)) {
                          $missing=join " ",($missing,"--reg");
                    } else {
                          $missing="--reg";
                    }
                }
                if (defined($set{"-d"})) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"-d"};
                    if (($tempdata[0]<0) or ($tempdata[0]>50)) {
                          print "-d out of range, using default value\n";
                          $command=join " ",($command,30);
                    } else {
                          $command=join " ",($command,int($tempdata[0]));
                    }
                } else {
                    $command=join " ",($command,30);
                }
                if (defined($set{"-f"})) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"-f"};
                    if (($tempdata[0]<10) or ($tempdata[0]>50)) {
                          print "-f out of range, using default value\n";
                          $command=join " ",($command,25);
                    } else {
                          $command=join " ",($command,int($tempdata[0]));
                    }
                } else {
                    $command=join " ",($command,25);
                }
                if ((defined($set{"-m"})) and (defined($set{"-w"}))) {
                    $command=join " ",($command,$set{"-m"},"z",$set{"-w"});
                } elsif (defined($set{"-m"})) {
                    $command=join " ",($command,$set{"-m"});
                } elsif (defined($set{"-w"})) {
                    $command=join " ",($command,$set{"-w"});
                } else {
                    $error="yes";
                    if (defined($missing)) {
                          $missing=join " ",($missing,"-m/-w");
                    } else {
                          $missing="-m/-w";
                    }
                }
           }
           if ($error eq "no") {
                system ($command);
           } else {
                $temp=join " ",("Error: missing",$missing,"for",$set{"module"});
                print "$temp\n";
           }
           if (($error eq "no") and ($set{"module"} eq "repeat")) {
                $command="perl can_confirmation.pl";
                @tempdata=();
                @tempdata=split /\s+/, $set{"-r"};
                $command=join " ",($command,$tempdata[0]);
                @tempdata=();
                @tempdata=split /\s+/, $set{"-c"};
                $command=join " ",($command,$tempdata[0]);
                @tempdata=();
                @tempdata=split /\s+/, $set{"--reg"};
                $command=join " ",($command,int($tempdata[0]),int($tempdata[0]));
                if (defined($set{"-d"})) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"-d"};
                    if (($tempdata[0]<0) or ($tempdata[0]>50)) {
                          $command=join " ",($command,30);
                    } else {
                          $command=join " ",($command,int($tempdata[0]));
                    }
                } else {
                    $command=join " ",($command,30);
                }
                if (defined($set{"-f"})) {
                    @tempdata=();
                    @tempdata=split /\s+/, $set{"-f"};
                    if (($tempdata[0]<10) or ($tempdata[0]>50)) {
                          $command=join " ",($command,25);
                    } else {
                          $command=join " ",($command,int($tempdata[0]));
                    }
                } else {
                    $command=join " ",($command,25);
                }
                if ((defined($set{"-m"})) and (defined($set{"-w"}))) {
                    $command=join " ",($command,$set{"-m"},"z",$set{"-w"});
                } elsif (defined($set{"-m"})) {
                    $command=join " ",($command,$set{"-m"});
                } else {
                    $command=join " ",($command,$set{"-w"});
                }
                system ($command);
           }
     }

}


exit;







