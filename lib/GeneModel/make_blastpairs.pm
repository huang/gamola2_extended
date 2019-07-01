#!/opt/ActivePerl-5.8/bin/perl

# Module to take output from blast-contigs and turn it into the "blast-pair"
# format readable by critica.
#input arguments: progress_bar, directory, filename, ini_ref, auto_ini_ref
#output: one blast-pair result file in input_directory "xxx.blast-pair"

package GeneModel::make_blastpairs;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&blastpair_critica);
use vars qw();

use initialise::read_me qw(:DEFAULT);
use GeneModel::blast_contigs;
use GeneModel::setup;

#local variables
my (%args, @data, $seq1, $seq2, $start1, $start2, $gstart, $gend, $exclude,
    $acc, $id, $ident, $pvalue, $q, $noprint, $aln, $e, $s, $end1, $end2, $v);

sub blastpair_critica {
   my %args = @_;
   $seq1    = ''; $seq2   = '';
   $start1  = 0 ; $start2 = 0;
   $gstart  = 0 ; $gend   = 0;
   $exclude = "nothing at all";

   ${$args{progress_bar}}->configure(-label=>"Generating Critica Blast pairs for $args{filename}");
   ${$args{main_window}}->update;

   #initialise input and output files
   open READ,   '<'.$args{directory}.'/'.$args{filename}.'.blast';
   open WRITE, '+>'.$args{directory}.'/'.$args{filename}.'.blast.pairs';

   #open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
   #print ERRORLOG "Error_make_blastpairs: .$args{directory}.'/'.$args{filename}.'_blastn and .blast.pairs\n\n";
   #close ERRORLOG;
   while(<READ>) {
     if (/^Query=/) {
       print_pairs($acc,$id,$ident,$pvalue,$start1,$end1,$start2,$end2,
                   $seq1,$seq2);
       ($q,$acc)=split(" ",$_);
       $_=<READ>;
       if (substr($_,0,1) ne " ") {
         $acc.=$_; # BLAST stupidly inserts carriage returns when names are long
       }
       ($acc,$gstart,$gend)=split(/\=/,$acc);
       $gstart--;
     }
     elsif (/^>/) {
       ($id)=split(" ",$_);
       ($q,$id)=split(">",$id);
       $noprint=0;
       if (/$exclude/) {$noprint=1;}
       $_=<READ>; # maybe species on next line
       if (/$exclude/) {$noprint=1;}
       $_=<READ>; # or even the next
       if (/$exclude/) {$noprint=1;}
     }
     elsif (/^Query:/) {
       ($q,$s,$aln,$e)=split(" ",$_);
       if ($start1==0) {
         $start1=$gstart+$s;
       }
       $end1=$gstart+$e;
       $seq1=$seq1.$aln."\n";
     }
     elsif (/^Sbjct:/) {
       ($q,$s,$aln,$e)=split(" ",$_);
       if ($start2==0) {$start2=$s;}
       $end2=$e;
       $seq2=$seq2.$aln."\n";
     }
     elsif ((/Expect =/) || (/Expect\([0-9]*\) =/))  {
       print_pairs($acc,$id,$ident,$pvalue,$start1,$end1,$start2,$end2,
                   $seq1,$seq2);
       @data=split(" ",$_);
       $pvalue=$data[$#data];
     }
     elsif (/Identities/) {
       ($q,$e,$v)=split(" ",$_);
       $ident=int(eval($v)*100);
     }

   }
   print_pairs($acc,$id,$ident,$pvalue,$start1,$end1,$start2,$end2,
               $seq1,$seq2);

   close READ;
   close WRITE;
   return (1);
}

sub print_pairs {
  ($acc,$id,$ident,$pvalue,$start1,$end1,$start2,$end2,$seq1,$seq2)=@_;
  if (($start1==0) || ($ident>97) || ($noprint)) {
    $seq1="";$seq2="";
    $start1=0;$start2=0;
    return;
  }
  printf WRITE (">$acc $start1-$end1 $ident $pvalue\n");
  print  WRITE uc($seq1);
  printf WRITE (">$id $start2-$end2 $ident $pvalue\n");
  print  WRITE uc($seq2);

  $seq1="";$seq2="";
  $start1=0;$start2=0;
}

1;