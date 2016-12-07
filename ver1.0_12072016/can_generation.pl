#!/usr/local/bin/perl -w
# Bookmarks: 0,0 0,0 0,0 0,71 0,0 0,0 0,0 0,0 0,0 0,71; CollapsedSubs: genome_size
### This program maps raw .fastq files and generates coverage&nucleotide_variation array files (./CAN/mutant.can OR wt.can)
### Requires buildcan.pl, mapping_single.pl, mapping_pair.pl
### Also generates a log.txt under ./CAN/
### *************
### Usage: perl can_generation.pl genome v dr dm s/p mutant_pair1 _pair2 (or mutant_single) size s/p wt_pair1 _pair2 (or wt_single) size chromosome start end
### *************
### "genome" is reference using ">" to indicate segments
### "v" is for "-v" in mapping_single/pair.pl
### "dr" is for $maxline in mapping_single/pair.pl
### "dm" is for $divide in buildcan.pl
### first "s/p" indicates whether the reads for mutant pool are single/pair ended
### second "s/p" is for wt pool
### "size" (can't be omitted if related "s/p"=p) is the value of -X (bowtie setting) for pair-end mapping
### "chromosome start end" specifies the region to export .can
### "chromosome" should be identical to the reference (the * of >* or >gi|*)
### "start end" starts from position 1
### If "chromosome" is defined but "start end" is undefined, the program will export results for the whole "chromosome"
### If "chromosome start end" is undefined, the program will export results for the whole genome

sub genome_size {
    print ("Preparing reference $_[0]\n");
    open (input_data, "<$_[0]") or die "Couldn't open: $!";
    my @size=();
    my @segmentID=();
    my $segmentnumber=0;
    my $temp="";
    my $i=0;

    while (1) {
           chomp ($temp=<input_data>);
           if ($temp eq "") {
                next;
           } elsif ($temp=~/^>/) {
                $segmentnumber++;
                $size[$segmentnumber]=0;
                $temp=$';
                if ($temp=~/^gi\|/) {
                     $temp=$';
                }
                if ($temp=~/\s+/) {
                     $temp=$`;
                }
                $segmentID[$segmentnumber]=$temp;
                next;
           } else {
                $size[$segmentnumber]=$size[$segmentnumber] + length ($temp);
           }
           if (eof) {
              last;
           }
    }
    close (input_data);

    open (output_log, ">./CAN/log.txt") or die "Couldn't open: $!";
    open (output_g, ">./genome_size.txt") or die "Couldn't open: $!";
    print output_log ">$_[0]\n";
    for ($i=1;$i<=$segmentnumber;$i++) {
          $temp=join "",($segmentID[$i]," = ",$size[$i]," nt");
          print output_log "$temp\n";
          print output_g "$temp\n";
          print ("$temp\n");
    }
    close (output_log);
    close (output_g);

    my $bowtiebuildcommand=join " ", ("bowtie-build --quiet", $_[0], "./genome.binary");
    system($bowtiebuildcommand);

}

### Main ###
system("mkdir ./CAN");
system("mkdir ./mapping");
&genome_size($ARGV[0]);
my $point=4;
my $command="";
my $chr="";
my $start=0;
my $end=0;
my $temp="";
if ($ARGV[$point] eq "s") {
     $point=$point+2;
} elsif ($ARGV[$point] eq "p") {
     $point=$point+4;
}
if ($ARGV[$point] eq "s") {
     $point=$point+2;
} elsif ($ARGV[$point] eq "p") {
     $point=$point+4;
}
if ($point==$#ARGV) {
     $chr=$ARGV[$point];
} elsif ($point<$#ARGV) {
     $chr=$ARGV[$point];
     $start=$ARGV[$point+1];
     $end=$ARGV[$point+2];
}

$point=4;
open (output_log, ">>./CAN/log.txt") or die "Couldn't open: $!";
print output_log "\n";
print output_log ">Mutant Pool\n";
print "Processing mutant pool\n";
close (output_log);
if ($ARGV[$point] eq "s") {
     $point++;
     $command=join " ",("perl mapping_single.pl",$ARGV[$point],$ARGV[0],$ARGV[1],$ARGV[2]);
     system("$command");
     $temp=$ARGV[$point];
     if ($temp=~/\.fastq$/) {
          $temp=$`;
     }
     while ($temp=~/\//) {
          $temp=$';
     }
     $command=$temp;
     $temp=$ARGV[0];
     while ($temp=~/\//) {
          $temp=$';
     }
     $command=join "",($command,"_MappingTo_",$temp);
     if ($chr eq "") {
          $command=join " ",("perl buildcan.pl",$command,"mutant.can",$ARGV[3]);
     } elsif (($start==0) and ($end==0)) {
          $command=join " ",("perl buildcan.pl",$command,"mutant.can",$ARGV[3],$chr);
     } else {
          $command=join " ",("perl buildcan.pl",$command,"mutant.can",$ARGV[3],$chr,$start,$end);
     }
     system("$command");
     $point++;
} elsif ($ARGV[$point] eq "p") {
     $point++;
     if ($ARGV[$point+2] eq "default") {
          $command=join " ",("perl mapping_pair.pl",$ARGV[$point],$ARGV[$point+1],$ARGV[0],250,$ARGV[1],$ARGV[2]);
     } else {
          $command=join " ",("perl mapping_pair.pl",$ARGV[$point],$ARGV[$point+1],$ARGV[0],$ARGV[$point+2],$ARGV[1],$ARGV[2]);
     }
     system("$command");
     $temp=$ARGV[$point];
     if ($temp=~/\.fastq$/) {
          $temp=$`;
     }
     while ($temp=~/\//) {
          $temp=$';
     }
     $command=$temp;
     $temp=$ARGV[$point+1];
     if ($temp=~/\.fastq$/) {
          $temp=$`;
     }
     while ($temp=~/\//) {
          $temp=$';
     }
     $command=join "",($command,"_",$temp);
     $temp=$ARGV[0];
     while ($temp=~/\//) {
          $temp=$';
     }
     $command=join "",($command,"_MappingTo_",$temp);
     if ($chr eq "") {
          $command=join " ",("perl buildcan.pl",$command,"mutant.can",$ARGV[3]);
     } elsif (($start==0) and ($end==0)) {
          $command=join " ",("perl buildcan.pl",$command,"mutant.can",$ARGV[3],$chr);
     } else {
          $command=join " ",("perl buildcan.pl",$command,"mutant.can",$ARGV[3],$chr,$start,$end);
     }
     system("$command");
     $point=$point+3;
}

open (output_log, ">>./CAN/log.txt") or die "Couldn't open: $!";
print output_log "\n";
print output_log ">Wild-type Pool\n";
print "Processing wild-type pool\n";
close (output_log);
if ($ARGV[$point] eq "s") {
     $point++;
     $command=join " ",("perl mapping_single.pl",$ARGV[$point],$ARGV[0],$ARGV[1],$ARGV[2]);
     system("$command");
     $temp=$ARGV[$point];
     if ($temp=~/\.fastq$/) {
          $temp=$`;
     }
     while ($temp=~/\//) {
          $temp=$';
     }
     $command=$temp;
     $temp=$ARGV[0];
     while ($temp=~/\//) {
          $temp=$';
     }
     $command=join "",($command,"_MappingTo_",$temp);
     if ($chr eq "") {
          $command=join " ",("perl buildcan.pl",$command,"wt.can",$ARGV[3]);
     } elsif (($start==0) and ($end==0)) {
          $command=join " ",("perl buildcan.pl",$command,"wt.can",$ARGV[3],$chr);
     } else {
          $command=join " ",("perl buildcan.pl",$command,"wt.can",$ARGV[3],$chr,$start,$end);
     }
     system("$command");
} elsif ($ARGV[$point] eq "p") {
     $point++;
     if ($ARGV[$point+2] eq "default") {
          $command=join " ",("perl mapping_pair.pl",$ARGV[$point],$ARGV[$point+1],$ARGV[0],250,$ARGV[1],$ARGV[2]);
     } else {
          $command=join " ",("perl mapping_pair.pl",$ARGV[$point],$ARGV[$point+1],$ARGV[0],$ARGV[$point+2],$ARGV[1],$ARGV[2]);
     }
     system("$command");
     $temp=$ARGV[$point];
     if ($temp=~/\.fastq$/) {
          $temp=$`;
     }
     while ($temp=~/\//) {
          $temp=$';
     }
     $command=$temp;
     $temp=$ARGV[$point+1];
     if ($temp=~/\.fastq$/) {
          $temp=$`;
     }
     while ($temp=~/\//) {
          $temp=$';
     }
     $command=join "",($command,"_",$temp);
     $temp=$ARGV[0];
     while ($temp=~/\//) {
          $temp=$';
     }
     $command=join "",($command,"_MappingTo_",$temp);
     if ($chr eq "") {
          $command=join " ",("perl buildcan.pl",$command,"wt.can",$ARGV[3]);
     } elsif (($start==0) and ($end==0)) {
          $command=join " ",("perl buildcan.pl",$command,"wt.can",$ARGV[3],$chr);
     } else {
          $command=join " ",("perl buildcan.pl",$command,"wt.can",$ARGV[3],$chr,$start,$end);
     }
     system("$command");
}

system ("rm ./genome.binary.1.ebwt");
system ("rm ./genome.binary.2.ebwt");
system ("rm ./genome.binary.3.ebwt");
system ("rm ./genome.binary.4.ebwt");
system ("rm ./genome.binary.rev.1.ebwt");
system ("rm ./genome.binary.rev.2.ebwt");
system ("rm -rf ./mapping");
system ("rm ./genome_size.txt");

############



exit;







