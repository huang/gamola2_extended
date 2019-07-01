#Subroutines used by various CRITICA scripts. This file needs to be somewhere
#on your PERLLIB path and be named Critica.pm

package Critica;
require Exporter;

@ISA = qw(Exporter);
@EXPORT = qw(getScriptName getPrefixName flagArguments realArguments loadNextContig printFASTA createBlastFragments bestInits cdsSort byloc computeSDScores computePromScores computeInitScores);

# returns the file name of the running script.
sub getScriptName {
 my @temp=split("/",$0);
 my $scriptname=pop(@temp);
 return $scriptname;
}

# returns the part of a filename stripped of all directories and suffixes
sub getPrefixName {
  my ($file)=@_;
  my $name="";
  my @strings;

  if (index($file,"/")) {
    @strings=split('/',$file);
    $name=pop(@strings); # name is last part of filename
  }
  else {$name=$file;}
  if (index($name,".")) {
    @strings=split('[.]',$name);
    $name=shift(@strings); # name is first part of filename
  }
  return $name;
}

# return array of flag arguments (those that start with "-")
sub flagArguments {
  my $arg;
  my @flags=();

  foreach $arg (@_) {
    if (substr($arg,0,1) eq "-") {
      push(@flags,$arg);
    }
  }
  return @flags;
}

# return array of real arguments (those that are not flags)
sub realArguments {
  my $arg;
  my @flags=();

  foreach $arg (@_) {
    if (substr($arg,0,1) ne "-") {
      push(@flags,$arg);
    }
  }
  return @flags;
}

# returns array containing next sequence and the sequence header
sub loadNextContig {
  my $fname=shift;
  my ($seq, $header, $curpos);

  $seq="";
  $header="";
  $curpos=0;

  while(<$fname>) {
    chop($_);
    if (/^>/) {
      if ($header ne "") {last;}
      $header=substr($_,1);
      $header=~tr/ /\_/;
      $header=~tr/\t/\_/;
    }
    else{$seq.=uc($_);}
    $curpos = tell($fname);
  }
  seek($fname, $curpos, 0);
  $seq=~tr/ //d;
  $seq=~s/\s//gs;
  $header =~ s/\W/_/g;
  return ($seq,$header);
}

#prints (in FASTA format) a sequence to open(file) handle $file
sub printFASTA {
  my($file, $seq, $header) = @_;
  my ($i, $n, $ln);

  print $file ">$header\n";
  $n = length($seq);
  for ($i=0; ($i < $n); $i += 60)
  {
    if (($i + 60) <= $n) {
      $ln = substr($seq,$i,60);
    }
    else {
      $ln = substr($seq,$i,($n-$i));
    }
    print $file "$ln\n";
  }
}


#subroutine to create fragments of genome for blast-contigs
sub createBlastFragments {
  my ($dir, $fasta, $subseqlen, $minimum) = @_;
  my ($seq, $subseq, $name, $subname, $start, $end, $seqlen);

  open(SEQ,$fasta) || die ("Can't open $fasta!\n");
  if (!-e($dir)) {
    mkdir($dir,0777);
  }

  $seq=" ";

  while($seq ne "") {
    ($seq,$name)=loadNextContig(\*SEQ);
    $seqlen=length($seq);

    for($i=0;$i<$seqlen;$i+=($subseqlen)) {
      $subseq=substr($seq,$i,$subseqlen);
      $start=$i+1;
      $end=$start+length($subseq)-1;
      if ($seqlen-$end<$minimum) { # fragment of contig remaining is too short
        $subseq=substr($seq,$i);
        $end=$seqlen;
        $i=$seqlen+1;
      }
      $subname=$name."=".$start."=".$end;
      open(FILE,">$dir/$subname") || die ("Can't create $dir/$subname!\n");
      printFASTA(FILE, $subseq,$subname);
      close(FILE);
      if ($end>$seqlen) {
          last;
      }
    }
  }
  close(SEQ);
}

#subroutine to select CRITICA hits with the best inits possible for each ORF
sub bestInits {
  my ($input, $output)=@_;
  my ($name, $start, $p, $matrix, $orfkey, $len);
  my (%bestpm, %bestline, %bestlen);

  open(INPUT,$input) || die("Can't open $input in bestInits\n");
  open(OUTPUT,">$output") || die("Can't open $output in bestInits\n");

  while(<INPUT>) {
    ($name,$start,$end,$p, $matrix)=split(" ",$_);
    $orfkey=$name." ".$end;
    $len=abs($start-$end);
    if  (!defined($bestline{$orfkey})) {
      $bestp{$orfkey}=1.0;
      $bestlen{$orfkey}=0;
    }
    if (($p<$bestp{$orfkey}) ||
            (($p==$bestp{$orfkey}) && ($len>$bestlen{$orfkey}))) {
      $bestp{$orfkey}=$p;
      $bestlen{$orfkey}=$len;
      $bestline{$orfkey}=$_;
    }
  }
  foreach $orfkey (keys %bestline) {
    print OUTPUT $bestline{$orfkey};
  }
  close(INPUT);
  close(OUTPUT);
}

# subroutine to sort ORFs by end position
sub cdsSort {
  my ($input)=@_;
  my ($name, $start, $end, $orfkey, %line);
  my $output=$input.".tmp";

  open(INPUT,$input) || die("Can't open $input in cdsSort\n");
  open(OUTPUT,">$output") || die("Can't open $output in cdsSort\n");

  while(<INPUT>) {
    ($name,$start,$end)=split(" ",$_);
    $orfkey=$name." ".$end;
    $line{$orfkey}=$_;
  }
  foreach $orfkey (sort byloc (keys %line)) {
    print OUTPUT $line{$orfkey};
  }
  close(INPUT);
  close(OUTPUT);
  rename($output, $input);
}

# subroutine called by cdsSort to do the sorting
sub byloc {
  my ($contiga, $enda)=split(" ",$a);
  my ($contigb, $endb)=split(" ",$b);

  if ($contiga gt $contigb) {return 1;}
  elsif ($contiga lt $contigb) {return -1;}
  else {
    if ($enda>$endb) {return 1;}
    elsif ($enda<$endb) {return -1;}
    else {return 0;}
  }
}

# subroutine to compute SD bonus scores
sub computeSDScores {
  my ($input, $output)=@_;
  my ($name, $start, $end, $p, $matrix, $comp, $di, $ibonus, $init, $slen);
  my ($sdbonus, $sstart, $sdseq, $orfkey, %bestslen, %bestscore, %bestline);
  my (%random, %real, $totrand, $totreal, $freqrand, $freqreal, $rl, $rand);
  my ($score);

  $totrand=0;
  $totreal=0;

  open(INPUT,$input) || die("Can't open $input in computeSDScores\n");
  open(OUTPUT,">$output") || die("Can't open $output in computeSDScores\n");

  while(<INPUT>) {
    ($name,$start, $end, $p, $matrix, $comp,$di,$ibonus,$init,
     $sdbonus,$sstart,$sdseq)=split(" ",$_);
    $orfkey=$name." ".$end;
    if (!defined($bestslen{$orfkey})) {
      $bestslen{$orfkey}=0;
      $bestscore{$orfkey}=0;
    }
    $slen=length($sdseq);
    if ($start<4) {
      $slen=5;
    }
    if ($slen==0) {$slen=1;}
    if ($slen>=$bestslen{$orfkey}) {
      if (($slen>$bestslen{$orfkey}) || ($comp+$di>$bestscore{$orfkey})) {
        $bestslen{$orfkey}=$slen;
        $bestline{$orfkey}=$_;
        $bestscore{$orfkey}=$comp+$di;
      }
    }
    $random{$sdseq}++;
    $totrand++;
  }

  foreach $orfkey (sort(keys %bestline)) {
    ($name,$start, $end, $p, $matrix, $comp,$di,$ibonus,$init,
     $sdbonus,$sstart,$sdseq)=split(" ",$bestline{$orfkey});
    $real{$sdseq}++;
    $totreal++;
  }

  while(($sdseq,$rl)=each %real) {
    $random{$sdseq}-=$real{$sdseq};
    $totrand-=$real{$sdseq};
  }

  while(($sdseq,$rand)=each %random) {
    if (defined($real{$sdseq})) {
      $rl=$real{$sdseq};
    }
    else {$rl=0;}
    $freqreal=$rl/$totreal;
    $freqrand=$rand/$totrand;
    if ($freqrand==0) {
      $freqrand=0.5/$totrand;
    }
    if($freqreal>0) {
      $score=log($freqreal/$freqrand);
    }
    else {
      $score=0;
    }
    printf(OUTPUT "%11s %4d %4d %8.3f %8.3f %8.3f\n",$sdseq,$rl,$rand,
         $freqreal,$freqrand,$score);
  }
  close(INPUT);
  close(OUTPUT);
}

# subroutine to compute promoter bonus scores
sub computePromScores {
  my ($input, $output)=@_;
  my ($name, $start, $end, $p, $matrix, $comp, $di, $ibonus, $init);
  my ($sdbonus, $sstart, $sdseq, $orfkey, %bestplen, %bestscore, %bestline);
  my (%random, %real, $totrand, $totreal, $freqrand, $freqreal, $rl, $rand);
  my ($prombonus,$promstart,$promseq,$promerr);
  my ($score);

  $totrand=0;
  $totreal=0;

  open(INPUT,$input) || die("Can't open $input in computeSDScores\n");
  open(OUTPUT,">$output") || die("Can't open $output in computeSDScores\n");

  while(<INPUT>) {
    ($name,$start, $end, $p, $matrix, $comp,$di,$ibonus,$init,
     $sdbonus,$sstart,$sdseq,$prombonus,$promstart,$promseq,$promerr)
       =split(" ",$_);
    $orfkey=$name." ".$end;
    if (!defined($bestperr{$orfkey})) {
      $bestperr{$orfkey}=9;
      $bestscore{$orfkey}=0;
    }
    if ($promerr<=$bestperr{$orfkey}) {
      if (($promerr<$bestperr{$orfkey}) || ($comp+$di>$bestscore{$orfkey})) {
        $bestperr{$orfkey}=$promerr;
        $bestline{$orfkey}=$_;
        $bestscore{$orfkey}=$comp+$di;
      }
    }
    $random{$promerr}++;
    $totrand++;
  }

  foreach $orfkey (sort(keys %bestline)) {
    ($name,$start, $end, $p, $matrix, $comp,$di,$ibonus,$init,
     $sdbonus,$sstart,$sdseq,$prombonus,$promstart,$promseq,$promerr)
       =split(" ",$bestline{$orfkey});
    $real{$promerr}++;
    $totreal++;
  }

  while(($promerr,$rl)=each %real) {
    $random{$promerr}-=$real{$promerr};
    $totrand-=$real{$promerr};
  }

  while(($promerr,$rand)=each %random) {
    if (defined($real{$promerr})) {
      $rl=$real{$promerr};
    }
    else {$rl=0;}
    $freqreal=$rl/$totreal;
    $freqrand=$rand/$totrand;
    if ($freqrand==0) {
      $freqrand=0.5/$totrand;
    }
    if($freqreal>0) {
      $score=log($freqreal/$freqrand);
    }
    else {
      $score=0;
    }
    printf(OUTPUT "%11s %4d %4d %8.3f %8.3f %8.3f\n",$promerr,$rl,$rand,
         $freqreal,$freqrand,$score);
  }
  close(INPUT);
  close(OUTPUT);
}

# subroutine to compute initiator bonus scores
sub computeInitScores {
  my($cds,$seq,$output)=@_;
  my (@fields, @random, $init, %real, $totreal, $cmd, %random, $totrand);
  my ($freqreal, $freqrand, $score);

  $totreal=0;
  $totrand=0;

  open(CDS,$cds) || die("Can't open $cds in computeInitScores\n");
  open(OUTPUT,">$output") || die("Can't open $output in computeInitScores\n");

  while(<CDS>) {
    @fields=split(" ",$_);
    $init=$fields[8];
    if ($init eq "off") {next;}
    $real{$init}++;
    $totreal++;
  }

  foreach $init (keys %real) {
    $cmd="motiffind $init $seq";
    @random=`$cmd`;
    $random{$init}=@random;
    $totrand+=$random{$init};
  }

  foreach $init (sort (keys %real)) {
    $freqreal=$real{$init}/$totreal;
    $freqrand=$random{$init}/$totrand;
    if($freqreal>0) {
      $score=log($freqreal/$freqrand);
    }
    else {
      $score=0;
    }
    printf(OUTPUT "%3s %5d %7d %8.3f %8.3f %8.3f\n",
           $init,$real{$init},$random{$init},$freqreal,$freqrand,$score);
  }
  close(CDS);
  close(OUTPUT);
}

return 1;

