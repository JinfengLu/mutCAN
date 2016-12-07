#!/usr/local/bin/perl -w
# Bookmarks: 0,0 0,0 0,0 0,142 0,0 0,0 0,0 0,0 0,0 0,142
### This program maps ./$ARGV[0] to ./genome.binary, allowing $ARGV[2] mismatches
### $ARGV[0] is raw fastq file
### $ARGV[1] is reference, only for information printing
### First divide $ARGV[0] by $ARGV[3] reads per subfile to reduce memory usage

sub fastqtofasta {
    print ("Preparing reads $_[0]\n");
    open (input_data, "<$_[0]") or die "Couldn't open: $!";
    $totalreads=0;
    $filenum=0;
    $readstitle=$_[0];
    if ($readstitle=~/\.fastq$/) {
         $readstitle=$`;
    }
    while ($readstitle=~/\//) {
         $readstitle=$';
    }
    my $temp="";
    my $id="";
    my $seq="";
    my $str="";
    my $qual="";
    my $count=0;
    my $maxline=$ARGV[3];

    $filenum++;
    $temp=join "",($readstitle,"_f",$filenum,".fa");
    open (output_fasta, ">./$temp") or die "Couldn't open: $!";

    while (1) {
           chomp ($id=<input_data>);
           chomp ($seq=<input_data>);
           chomp ($str=<input_data>);
           chomp ($qual=<input_data>);
           if (!($seq=~/N/)) {
               $count++;
               $totalreads++;
               if ($totalreads%10000000==0) {
                    print ("$totalreads\n");
               }
               $temp=join "",(">seq",$totalreads);
               print output_fasta "$temp\n";
               print output_fasta "$seq\n";
           }
           if (eof) {
               last;
           } elsif ($count==$maxline) {
               close (output_fasta);
               $count=0;
               $filenum++;
               $temp=join "",($readstitle,"_f",$filenum,".fa");
               open (output_fasta, ">./$temp") or die "Couldn't open: $!";
           }
    }
    print ("$totalreads reads prepared\n");
    close (input_data);
    close (output_fasta);

}

sub bowtie_program {
    my $output_file_trunc="";
    my $output_file_map="";
    my $i=0;
    my $bowtiecommand="";
    $bowtiesetting="";

    for ($i=1;$i<=$filenum;$i++) {
          print ("Mapping $ARGV[0] to $ARGV[1] ($i/$filenum)\n");
          $output_file_trunc=join "",($readstitle,"_f",$i,".fa");
          $output_file_map=join "",($readstitle,"_f",$i,"_mapping");
          $bowtiecommand=join "", ("bowtie -f -a ./genome.binary ./", $output_file_trunc, " -v $ARGV[2] -p 4 -e 140 --quiet --suppress 6,7 ./", $output_file_map);
          $bowtiesetting="-v $ARGV[2] -e 140";
          system ($bowtiecommand);
          system ("rm ./$output_file_trunc");
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
    $temp=$ARGV[1];
    while ($temp=~/\//) {
          $temp=$';
    }
    my $output_map=join "",($readstitle,"_MappingTo_",$temp);

    for ($i=1;$i<=$filenum;$i++) {
          print ("Generatig mapping result ($i/$filenum)\n");
          %hit=();
          @id=();
          @info=();
          $output_file_map=join "",($readstitle,"_f",$i,"_mapping");
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

    open (output_log, ">>./CAN/log.txt") or die "Couldn't open: $!";
    print output_log "Reads file: $ARGV[0]\n";
    print output_log "Total reads = $totalreads\n";
    print ("Total reads = $totalreads\n");
    $temp=int($mapped/$totalreads*1000+0.5)/10;
    print output_log "Mapped reads = $mapped ($temp%)\n";
    print ("Mapped reads = $mapped ($temp%)\n");
    print output_log "Bowtie mapping settings: $bowtiesetting\n";
    close (output_log);
}

### Main ###
&fastqtofasta($ARGV[0]);
&bowtie_program;
&combinemapping;
############

exit;







