#!/usr/local/bin/perl -w
# Bookmarks: 0,0 0,0 0,0 21,173 0,0 0,0 0,0 0,0 0,0 21,173
### This program maps $ARGV[0] $ARGV[1] to ./genome.binary, allowing $ARGV[4] mismatches
### $ARGV[0,1] is raw fastq file
### $ARGV[2] is reference, only for information printing
### $ARGV[3] is the value for -X (bowtie setting)
### First divide $ARGV[0,1] by $ARGV[5] reads per subfile to reduce memory usage

sub fastqtofasta {
    print ("Preparing read pairs $_[0] & $_[1]\n");
    open (input_data1, "<$_[0]") or die "Couldn't open: $!";
    open (input_data2, "<$_[1]") or die "Couldn't open: $!";
    $totalreads=0;
    $filenum=0;
    $readstitle1=$_[0];
    $readstitle2=$_[1];
    if ($readstitle1=~/\.fastq$/) {
         $readstitle1=$`;
    }
    while ($readstitle1=~/\//) {
         $readstitle1=$';
    }
    if ($readstitle2=~/\.fastq$/) {
         $readstitle2=$`;
    }
    while ($readstitle2=~/\//) {
         $readstitle2=$';
    }
    my $temp="";
    my $id1="";
    my $seq1="";
    my $str1="";
    my $qual1="";
    my $id2="";
    my $seq2="";
    my $str2="";
    my $qual2="";
    my $count=0;
    my $maxline=$ARGV[5];

    $filenum++;
    $temp=join "",($readstitle1,"_f",$filenum,".fa");
    open (output_fasta1, ">./$temp") or die "Couldn't open: $!";
    $temp=join "",($readstitle2,"_f",$filenum,".fa");
    open (output_fasta2, ">./$temp") or die "Couldn't open: $!";

    while (1) {
           chomp ($id1=<input_data1>);
           chomp ($seq1=<input_data1>);
           chomp ($str1=<input_data1>);
           chomp ($qual1=<input_data1>);
           chomp ($id2=<input_data2>);
           chomp ($seq2=<input_data2>);
           chomp ($str2=<input_data2>);
           chomp ($qual2=<input_data2>);
           if ((!($seq1=~/N/)) and (!($seq2=~/N/))) {
               $count++;
               $totalreads++;
               if ($totalreads%10000000==0) {
                    print ("$totalreads\n");
               }
               $temp=join "",(">p1seq",$totalreads);
               print output_fasta1 "$temp\n";
               print output_fasta1 "$seq1\n";
               $temp=join "",(">p2seq",$totalreads);
               print output_fasta2 "$temp\n";
               print output_fasta2 "$seq2\n";
           }
           if (eof) {
               last;
           } elsif ($count==$maxline) {
               close (output_fasta1);
               close (output_fasta2);
               $count=0;
               $filenum++;
               $temp=join "",($readstitle1,"_f",$filenum,".fa");
               open (output_fasta1, ">./$temp") or die "Couldn't open: $!";
               $temp=join "",($readstitle2,"_f",$filenum,".fa");
               open (output_fasta2, ">./$temp") or die "Couldn't open: $!";
       }
    }
    print ("$totalreads pairs prepared\n");
    close (input_data1);
    close (input_data2);
    close (output_fasta1);
    close (output_fasta2);

}

sub bowtie_program {
    my $output_file_trunc1="";
    my $output_file_trunc2="";
    my $output_file_map="";
    my $i=0;
    my $bowtiecommand="";
    $bowtiesetting="";

    for ($i=1;$i<=$filenum;$i++) {
          print ("Mapping $ARGV[0] & $ARGV[1] to $ARGV[2] ($i/$filenum)\n");
          $output_file_trunc1=join "",($readstitle1,"_f",$i,".fa");
          $output_file_trunc2=join "",($readstitle2,"_f",$i,".fa");
          $output_file_map=join "",($readstitle1,"_",$readstitle2,"_f",$i,"_mapping");
          $bowtiecommand=join "", ("bowtie -f -a ./genome.binary -1 ./", $output_file_trunc1, " -2 ./", $output_file_trunc2, " -X $ARGV[3] -v $ARGV[4] -p 4 -e 140 --quiet --suppress 6,7 ./", $output_file_map);
          $bowtiesetting="-X $ARGV[3] -v $ARGV[4] -e 140";
          system ($bowtiecommand);
          system ("rm ./$output_file_trunc1");
          system ("rm ./$output_file_trunc2");
    }
}

sub combinemapping {
    my $i=0;
    my $j=0;
    my $mapped=0;
    my %hit=();
    my $temp="";
    my @id=();
    my @info=();
    my $output_file_map="";
    $temp=$ARGV[2];
    while ($temp=~/\//) {
          $temp=$';
    }
    my $output_map=join "",($readstitle1,"_",$readstitle2,"_MappingTo_",$temp);

    for ($i=1;$i<=$filenum;$i++) {
          print ("Generatig mapping result ($i/$filenum)\n");
          %hit=();
          @id=();
          @info=();
          $output_file_map=join "",($readstitle1,"_",$readstitle2,"_f",$i,"_mapping");
          open (output_file, "<./$output_file_map") or die "Couldn't open: $!";
          while (1) {
                 chomp ($temp=<output_file>);
                 $temp=~/\s+/;
                 $id[$#id+1]=$`;
                 if (defined($hit{$`})) {
                     $hit{$`}++;
                 } else {
                     $hit{$`}=1;
                     $mapped++;
                 }
                 $info[$#info+1]=join "",($&,$');
                 if (eof) {
                     last;
                 }
          }
          close (output_file);
          system ("rm ./$output_file_map");

          if ($i==1) {
               open (output_result, ">./mapping/$output_map") or die "Couldn't open: $!";
          } else {
               open (output_result, ">>./mapping/$output_map") or die "Couldn't open: $!";
          }
          for ($j=0;$j<=$#id;$j++) {
               $temp=join "",($id[$j],"_",$hit{$id[$j]},$info[$j]);
               print output_result "$temp\n";
          }
          close (output_result);
    }

    $mapped=$mapped/2;
    open (output_log, ">>./CAN/log.txt") or die "Couldn't open: $!";
    print output_log "Reads file: $ARGV[0] & $ARGV[1]\n";
    print output_log "Total pairs = $totalreads\n";
    print ("Total pairs = $totalreads\n");
    $temp=int($mapped/$totalreads*1000+0.5)/10;
    print output_log "Mapped pairs = $mapped ($temp%)\n";
    print ("Mapped pairs = $mapped ($temp%)\n");
    print output_log "Bowtie mapping settings: $bowtiesetting\n";
    close (output_log);
}

### Main ###
&fastqtofasta($ARGV[0],$ARGV[1]);
&bowtie_program;
&combinemapping;
############

exit;







