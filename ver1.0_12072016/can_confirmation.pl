#!/usr/local/bin/perl -w
# Bookmarks: 0,0 0,0 0,0 0,259 0,0 0,0 0,0 0,0 0,0 0,259; CollapsedSubs: complementary  reverse  processgenome  processreads  drawpic
### Usage: perl can_confirmation.pl genome chromosome start end distance flanklength file1 file2 z file3 file4 file5 z file6...

sub complementary {
    my $i=0;
    my $temp="";
    my $letter="";
    for ($i=0;$i<=length($_[0])-1;$i++) {
          $letter=substr ($_[0], $i, 1);
          if ($letter eq "A") {
               $temp=join "", ($temp, "T");
          } elsif ($letter eq "T") {
               $temp=join "", ($temp, "A");
          } elsif ($letter eq "G") {
               $temp=join "", ($temp, "C");
          } elsif ($letter eq "C") {
               $temp=join "", ($temp, "G");
          } else {
               $temp=join "", ($temp, $letter);
          }
    }
    return $temp;
}

sub reverse {
    my $i=0;
    my $temp="";
    my $letter="";
    for ($i=length($_[0])-1;$i>=0;$i--) {
          $letter=substr ($_[0], $i, 1);
          $temp=join "", ($temp, $letter);
    }
    return $temp;
}

sub processgenome {
    print "Processing genome $_[0]\n";
    open (input_data, "<$_[0]") or die "Couldn't open: $!";
    my $find="no";
    my $temp="";
    my @temp=();
    my $linestart=0;
    my $seq="";
    my $start=$_[2]-$_[4]-$flanklength;
    my $end=$_[3]+$_[4]+$flanklength;
    while (1) {
           chomp ($temp=<input_data>);
           if ($temp ne "") {
                if ($temp=~/^>/) {
                     $temp=$';
                     if ($temp=~/^gi\|/) {
                          $temp=$';
                     }
                     @temp=split /\s+/, $temp;
                     if ($temp[0] eq $_[1]) {
                          $find="yes";
                          $linestart=1;
                     }
                } elsif ($find eq "yes") {
                     if ($linestart>$end) {
                          last;
                     } elsif (($linestart+length($temp))<$start) {
                          $linestart=$linestart+length($temp);
                     } else {
                          $seq="";
                          if (($start<=$linestart) and ($end>=$linestart+length($temp))) {
                               $seq=$temp;
                          } elsif (($start>$linestart) and ($end<$linestart+length($temp))) {
                               $seq=substr($temp,$start-$linestart,$end-$start+1);
                          } elsif (($start<=$linestart) and ($end<$linestart+length($temp))) {
                               $seq=substr($temp,0,$end-$linestart+1);
                          } elsif (($start>$linestart) and ($end>=$linestart+length($temp))) {
                               $seq=substr($temp,$start-$linestart,length($temp)-$start+$linestart);
                          }
                          $refplus=join "",($refplus,$seq);
                          $linestart=$linestart+length($temp);
                     }

                }
           }
           if (eof) {
              last;
           }
    }
    close (input_data);
    $flank1=substr($refplus,0,$flanklength);
    $flank2=&complementary($flank1);
    $flank2=&reverse($flank2);
    $flank3=substr($refplus,length($refplus)-$flanklength,$flanklength);
    $flank4=&complementary($flank3);
    $flank4=&reverse($flank4);
    $refminus=&complementary($refplus);
}

sub processreads {
    open (input_data, "<$_[0]") or die "Couldn't open: $!";
    my $id1="";
    my $seq="";
    my $id2="";
    my $qual="";
    if (-s "./$_[0]") {
         while (1) {
                chomp ($id1=<input_data>);
                chomp ($seq=<input_data>);
                chomp ($id2=<input_data>);
                chomp ($qual=<input_data>);
                if ($seq=~/\Q$flank1\E/) {
                      $seq=join "",($flank1,$');
                      if (length($seq)>length($refplus)) {
                          $seq=substr($seq,0,length($refplus));
                      }
                      $reads1[$#reads1+1]=$seq;
                } elsif ($seq=~/\Q$flank2\E/) {
                      $seq=join "",($`,$flank2);
                      $seq=&reverse($seq);
                      if (length($seq)>length($refminus)) {
                          $seq=substr($seq,0,length($refminus));
                      }
                      $reads2[$#reads2+1]=$seq;
                } elsif ($seq=~/\Q$flank3\E/) {
                      $seq=join "",($`,$flank3);
                      if (length($seq)>length($refplus)) {
                          $seq=substr($seq,length($seq)-length($refplus),length($refplus));
                      }
                      $reads3[$#reads3+1]=$seq;
                } elsif ($seq=~/\Q$flank4\E/) {
                      $seq=join "",($flank4,$');
                      $seq=&reverse($seq);
                      if (length($seq)>length($refminus)) {
                          $seq=substr($seq,length($seq)-length($refminus),length($refminus));
                      }
                      $reads4[$#reads4+1]=$seq;
                }
                if (eof) {
                      last;
                }
         }
    }
    close input_data;
    print ("$_[0] reads checked\n");
}

sub drawpic {

    my $i=0;
    my $j=0;
    my $ref="";
    my $start=$_[2]-$_[4]-$flanklength;
    my $nt="";
    my $maxbarnumberplus=$#reads1+$#reads3+2;
    my $maxbarnumberminus=$#reads2+$#reads4+2;

    use GD;
    my $output_file=join "", ("$_[1] $_[2] $_[3] $_[0]",".png");

    if (($maxbarnumberplus+$maxbarnumberminus)>=300) {
          $output_file=~s/png$/txt/;
          open (output_result, ">./CAN/$output_file") or die "Couldn't open: $!";;
          print output_result "Plus reads number = $maxbarnumberplus\n";
          print output_result "Minus reads number = $maxbarnumberminus\n";
          close output_result;
          print ("Failed: Too many reads\n");
          return;
    } else {
          print ("Success\n");
    }

    open (output_result, ">./CAN/$output_file") or die "Couldn't open: $!";;
    my $letterwidth=7;
    my $barwidth=10;
    my $first="yes";
    my $img=GD::Image->new(100+length($refplus)*$letterwidth,20+20+($maxbarnumberplus+$maxbarnumberminus+4)*$barwidth);

    local $white=$img->colorAllocate(255,255,255);
    local $red=$img->colorAllocate(255,0,0);
    local $black=$img->colorAllocate(0,0,0);
    local $blue=$img->colorAllocate(0,0,255);
    $img->fill(1,1,$white);

    $img->filledRectangle(50,10+5+($maxbarnumberplus+1)*$barwidth,50+length($refplus)*$letterwidth,10+5+3+($maxbarnumberplus+1)*$barwidth,$black);
    $img->filledRectangle(50,10+10+2+($maxbarnumberplus+3)*$barwidth,50+length($refplus)*$letterwidth,10+10+2+3+($maxbarnumberplus+3)*$barwidth,$black);
    for ($i=1;$i<=length($refplus);$i++) {
          $ref=substr($refplus,$i-1,1);
          if (($start+$i-1>=$_[2]) and ($start+$i-1<=$_[3])) {
                 $img->string(gdSmallFont,50+($i-1)*$letterwidth,10+10+($maxbarnumberplus+1)*$barwidth,$ref,$blue);
          } else {
                 $img->string(gdSmallFont,50+($i-1)*$letterwidth,10+10+($maxbarnumberplus+1)*$barwidth,$ref,$black);
          }
          $ref=substr($refminus,$i-1,1);
          if (($start+$i-1>=$_[2]) and ($start+$i-1<=$_[3])) {
                 $img->string(gdSmallFont,50+($i-1)*$letterwidth,10+10+($maxbarnumberplus+2)*$barwidth,$ref,$blue);
          } else {
                 $img->string(gdSmallFont,50+($i-1)*$letterwidth,10+10+($maxbarnumberplus+2)*$barwidth,$ref,$black);
          }
          if (($start+$i-1)%10==0) {
               $img->line(50+($i-1+0.4)*$letterwidth,10+2+($maxbarnumberplus+1)*$barwidth,50+($i-1+0.4)*$letterwidth,10+5+($maxbarnumberplus+1)*$barwidth,$black);
               $img->line(50+($i-1+0.4)*$letterwidth,10+10+5+($maxbarnumberplus+3)*$barwidth,50+($i-1+0.4)*$letterwidth,10+10+5+3+($maxbarnumberplus+3)*$barwidth,$black);
               if ($first eq "yes") {
                    $img->string(gdSmallFont,50+($i-1-0.4)*$letterwidth,10+$maxbarnumberplus*$barwidth,$start+$i-1,$black);
                    $img->string(gdSmallFont,50+($i-1-0.4)*$letterwidth,10+20+($maxbarnumberplus+3)*$barwidth,$start+$i-1,$black);
                    $first="no";
               } else {
                    $img->string(gdSmallFont,50+($i-1-0.4)*$letterwidth,10+$maxbarnumberplus*$barwidth,($start+$i-1)%1000,$black);
                    $img->string(gdSmallFont,50+($i-1-0.4)*$letterwidth,10+20+($maxbarnumberplus+3)*$barwidth,($start+$i-1)%1000,$black);
               }
          }
    }

    for ($i=0;$i<=$#reads1;$i++) {
          for ($j=1;$j<=length($reads1[$i]);$j++) {
               $ref=substr($refplus,$j-1,1);
               $nt=substr($reads1[$i],$j-1,1);
               if ($ref eq $nt) {
                    $img->string(gdSmallFont,50+($j-1)*$letterwidth,10+($maxbarnumberplus-$i-1)*$barwidth,$nt,$black);
               } else {
                    $img->string(gdSmallFont,50+($j-1)*$letterwidth,10+($maxbarnumberplus-$i-1)*$barwidth,$nt,$red);
               }
          }
    }
    for ($i=0;$i<=$#reads2;$i++) {
          for ($j=1;$j<=length($reads2[$i]);$j++) {
               $ref=substr($refminus,$j-1,1);
               $nt=substr($reads2[$i],$j-1,1);
               if ($ref eq $nt) {
                    $img->string(gdSmallFont,50+($j-1)*$letterwidth,10+20+($maxbarnumberplus+4+$i)*$barwidth,$nt,$black);
               } else {
                    $img->string(gdSmallFont,50+($j-1)*$letterwidth,10+20+($maxbarnumberplus+4+$i)*$barwidth,$nt,$red);
               }
          }
    }
    for ($i=0;$i<=$#reads3;$i++) {
          for ($j=1;$j<=length($reads3[$i]);$j++) {
               $ref=substr($refplus,length($refplus)-length($reads3[$i])+$j-1,1);
               $nt=substr($reads3[$i],$j-1,1);
               if ($ref eq $nt) {
                    $img->string(gdSmallFont,50+(length($refplus)-length($reads3[$i])+$j-1)*$letterwidth,10+($maxbarnumberplus-$#reads1-$i-2)*$barwidth,$nt,$black);
               } else {
                    $img->string(gdSmallFont,50+(length($refplus)-length($reads3[$i])+$j-1)*$letterwidth,10+($maxbarnumberplus-$#reads1-$i-2)*$barwidth,$nt,$red);
               }
          }
    }
    for ($i=0;$i<=$#reads4;$i++) {
          for ($j=1;$j<=length($reads4[$i]);$j++) {
               $ref=substr($refminus,length($refminus)-length($reads4[$i])+$j-1,1);
               $nt=substr($reads4[$i],$j-1,1);
               if ($ref eq $nt) {
                    $img->string(gdSmallFont,50+(length($refminus)-length($reads4[$i])+$j-1)*$letterwidth,10+20+($maxbarnumberplus+4+$#reads2+$i+1)*$barwidth,$nt,$black);
               } else {
                    $img->string(gdSmallFont,50+(length($refminus)-length($reads4[$i])+$j-1)*$letterwidth,10+20+($maxbarnumberplus+4+$#reads2+$i+1)*$barwidth,$nt,$red);
               }
          }
    }

    binmode output_result;
    print output_result GD::Image::png($img);
    close (output_result);

}

$refplus=""; ### $flank1+"distance"+potential mutation region+"distance"+$flank3, 5' to 3' direction
$refminus=""; ### $refplus complementary, 3' to 5' direction
$flank1=""; ### 5' + upstream to the potential mutation region
$flank2=""; ### 5' - upstream to the potential mutation region
$flank3=""; ### 3' + downstream
$flank4=""; ### 3' - downstream

$flanklength=$ARGV[5];
&processgenome($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4]);
my $i=0;
my $filename="";
my $temp="";
@reads1=(); ### reads containing $flank1
@reads2=(); ### reads containing $flank2, sequence reversed
@reads3=(); ### reads containing $flank3
@reads4=(); ### reads containing $flank4, sequence reversed
for ($i=6;$i<=$#ARGV;$i++) {
      if ($ARGV[$i] eq "z") {
           &drawpic($filename,$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4]);
           @reads1=();
           @reads2=();
           @reads3=();
           @reads4=();
           $filename="";
      } else {
           $temp=$ARGV[$i];
           if ($temp=~/\.fastq$/) {
                $temp=$`;
           }
           while ($temp=~/\//) {
                $temp=$';
           }
           if ($filename eq "") {
                $filename=$temp;
           } else {
                $filename="$filename $temp";
           }
           &processreads($ARGV[$i]);
      }
}
&drawpic($filename,$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4]);

exit;







