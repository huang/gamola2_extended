#!/opt/ActivePerl-5.8/bin/perl

#reformat sequences to fasta and translate nt to aa
#input arguments: main_window, progress_bar, sequence, left_bd, right_bd, orientation, ini_ref, auto_ini_ref

#nt2seq returns 1
#nt2aa  returns aa_sequence or 0

package ProgrammeModules::sequence;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&seq2fasta &nt2aa &nt2aa_unchecked &orient_sequence);
use vars qw();

use Basics::MesgBox;
use initialise::read_me qw(:DEFAULT);

#local variables
my (%args, %codon2aa, $fasta_seq, @temp_array, $codon, $header, $i);

#translation hash
#codons and corresponding amino acids
BEGIN {
   $codon2aa{"ttt"}= "F"; $codon2aa{"ttc"}= "F"; $codon2aa{"tta"}= "L";
   $codon2aa{"ttg"}= "L"; $codon2aa{"tct"}= "S"; $codon2aa{"tcc"}= "S";
   $codon2aa{"tca"}= "S"; $codon2aa{"tcg"}= "S"; $codon2aa{"tat"}= "Y";
   $codon2aa{"tac"}= "Y"; $codon2aa{"tgt"}= "C"; $codon2aa{"tgc"}= "C";
   $codon2aa{"ctt"}= "L"; $codon2aa{"ctc"}= "L"; $codon2aa{"cta"}= "L";
   $codon2aa{"ctg"}= "L"; $codon2aa{"cct"}= "P"; $codon2aa{"ccc"}= "P";
   $codon2aa{"cca"}= "P"; $codon2aa{"ccg"}= "P"; $codon2aa{"cat"}= "H";
   $codon2aa{"cac"}= "H"; $codon2aa{"caa"}= "Q"; $codon2aa{"cag"}= "Q";
   $codon2aa{"cgt"}= "R"; $codon2aa{"cgc"}= "R"; $codon2aa{"cga"}= "R";
   $codon2aa{"cgg"}= "R"; $codon2aa{"att"}= "I"; $codon2aa{"atc"}= "I";
   $codon2aa{"ata"}= "I"; $codon2aa{"atg"}= "M"; $codon2aa{"act"}= "T";
   $codon2aa{"acc"}= "T"; $codon2aa{"aca"}= "T"; $codon2aa{"acg"}= "T";
   $codon2aa{"aat"}= "N"; $codon2aa{"aac"}= "N"; $codon2aa{"aaa"}= "K";
   $codon2aa{"aag"}= "K"; $codon2aa{"agt"}= "S"; $codon2aa{"agc"}= "S";
   $codon2aa{"aga"}= "R"; $codon2aa{"agg"}= "R"; $codon2aa{"gtt"}= "V";
   $codon2aa{"gtc"}= "V"; $codon2aa{"gta"}= "V"; $codon2aa{"gtg"}= "V";
   $codon2aa{"gct"}= "A"; $codon2aa{"gcc"}= "A"; $codon2aa{"gca"}= "A";
   $codon2aa{"gcg"}= "A"; $codon2aa{"gat"}= "D"; $codon2aa{"gac"}= "D";
   $codon2aa{"gaa"}= "E"; $codon2aa{"gag"}= "E"; $codon2aa{"ggt"}= "G";
   $codon2aa{"ggc"}= "G"; $codon2aa{"gga"}= "G"; $codon2aa{"ggg"}= "G";
   $codon2aa{"tgg"}= "W";
   $codon2aa{"tag"}= "*";
   $codon2aa{"tga"}= "*";
   $codon2aa{"taa"}= "*";
}

sub seq2fasta {
   my %args = (header      => 'anonymous',
               orientation => 'sense',
               @_);

   #fasta format?
   if ($args{sequence} =~ /^\>/) {
      $args{sequence}  =~ s/\>([^\n]*?)\n//;
      $args{header}    = $1;
      $args{sequence}  =~ s/\s//gs;

   } else {
      $args{sequence}  =~ s/\s//gs;
   }

   #format it nicely
   for($i=0;$i<length($args{sequence});$i+=50){
      $fasta_seq.=(substr($args{sequence},$i,50))."\n";
   }

   $args{sequence} = '>'.$args{header}."\n".$fasta_seq;
   $fasta_seq = "";

   return (\$args{sequence});
}

sub nt2aa {
   my %args = (filename => '',
               module   => '',
               @_
              );
   my ($start_codon, $stop_codon, $aa_seq, $count_stop, $count_X);
   my (@temp_array, @codons);
   $count_stop = 0;
   $count_X = 0;

   if ($args{sequence} !~ /\S+/ || length ($args{sequence}) < 3 || !defined $args{sequence}) {
      ${$args{main_window}}->messageBox(
                                         -title   => 'Error',
                                         -message => "Empty or insufficient sequence string delivered for $args{left_bd}\nfrom module $args{module}.",
                                         -icon    => 'question',
                                         -type    => 'OK'
                                        );

      return (0);
   }

   #convert to lowercase
   $args{sequence} = lc ($args{sequence});

   #orient sequence
   $args{sequence} = &orient_sequence (seq         => $args{sequence},
                                       orientation => $args{orientation});

   #generate codon array
   @codons = ($args{sequence} =~ m/(\w\w\w)/g);

   #check for correct start codon  MODIFIED: initially allowed start_codons atg|gtg|ttg --> tolerated added start_codons ata|atc|att|ctg
   unless ($codons[0] =~ /(atg|gtg|ttg|ata|atc|att|ctg)/) {
      ${$args{progress_bar}}->configure(-label =>"Error in file $args{filename} from module $args{module}.".
                                                 "\nDNA sequence does not start with a recognized start codon: $codons[0]".
                                                 "\nIt is advisable to check the gene model at positions $args{left_bd} and $args{right_bd}.".
                                                 "\nAnalysis continues\n"
                                        );
      ${$args{main_window}}->update; #sleep(1);
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in file $args{filename} from module $args{module}.".
                  "\nDNA sequence does not start with a recognized start codon: $codons[0]".
                  "\nIt is advisable to check the gene model at positions $args{left_bd} and $args{right_bd}.\n\n";
      close WRITE;
   }
   #check for correct stop codon -> non-critical error
   unless ($codons[-1] =~ /(tag|tga|taa)/) {
      ${$args{progress_bar}}->configure(-label =>"Error in file $args{filename} from module $args{module}.".
                                                 "\nDNA sequence does not end with a recognized stop codon: $codons[-1]".
                                                 "\nIt is advisable to check the gene model at positions $args{left_bd} and $args{right_bd}.".
                                                 "\nAnalysis continues\n"
                                         );
      ${$args{main_window}}->update;#sleep(1);
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in file $args{filename} from module $args{module}.".
                  "\nDNA sequence does not end with a recognized stop codon: $codons[-1]".
                  "\nIt is advisable to check the gene model at positions $args{left_bd} and $args{right_bd}.\n\n";
      close WRITE;
   }

   foreach my $codon (@codons){
      if (defined $codon2aa{$codon}){
         $aa_seq.=$codon2aa{$codon};
      } else {
         $aa_seq.="X"; #undecipherable nucleic acid
      }
   }

   #check for correct frame: if there are too many 'x' and more than on '*', its the wrong frame
   $count_stop = ($aa_seq =~ tr/\*//);
   $count_X = ($aa_seq =~ tr/[Xx]//);

   if ($count_stop > 1) {
      ${$args{progress_bar}}->configure(-label =>"Error in file $args{filename} from module $args{module}.".
                                                 "\nThe deduced amino acid sequence contains more than one stop codon".
                                                 "\nThis could be a critical error".
                                                 "\nPlease verify the gene model at positions $args{left_bd} and $args{right_bd}"

                                         );
      ${$args{main_window}}->update;#sleep(1);
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in file $args{filename} from module $args{module}.".
                  "\nThe deduced amino acid sequence contains more than one stop codon".
                  "\nThis could be a critical error".
                  "\nPlease verify the gene model at positions $args{left_bd} and $args{right_bd}".
                  "\n$args{orientation}  Sequence: $aa_seq\n\n";
      close WRITE;
      #${$args{main_window}}->messageBox(
      #                                   -title   => 'Error',
      #                                   -message => "Error in file $args{filename}".
      #                                               "\nThe deduced amino acid sequence contains more than one stop codon".
      #                                               "\nThis could be a critical error".
      #                                               "\nPlease verify the gene model at positions $args{left_bd} and $args{right_bd}".
      #                                               "\n$args{orientation}  Sequence: $aa_seq",
      #                                   -icon    => 'question',
      #                                   -type    => 'OK'
      #                                  );
   }

   if ($count_X / (length($aa_seq)) > 0.10 ) { #more than 10% undefined aa?
      ${$args{progress_bar}}->configure(-label =>"Error in file $args{filename} from module $args{module}.".
                                                 "\nThe deduced amino acid sequence contains ".
                                                 (($count_X * 100) / (length($aa_seq))).' percent undefined codon(s)'.
                                                 "\nThis could be a critical error".
                                                 "\nPlease verify the gene model at positions $args{left_bd} and $args{right_bd}"
                                         );
      ${$args{main_window}}->update;#sleep(1);
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in file $args{filename} from module $args{module}.".
                  "\nThe deduced amino acid sequence contains ".
                  (($count_X * 100) / (length($aa_seq))).' percent undefined codon(s)'.
                  "\nThis could be a critical error".
                  "\nPlease verify the gene model at positions $args{left_bd} and $args{right_bd}\n\n";
      close WRITE;

      #${$args{main_window}}->messageBox(
      #                                   -title   => 'Error',
      #                                   -message => "Error in file $args{filename}".
      #                                               "\nThe deduced amino acid sequence contains ".
      #                                               (($count_X * 100) / (length($aa_seq))).' percent undefined codon(s)'.
      #                                               "\nThis could be a critical error".
      #                                               "\nPlease verify the gene model at positions $args{left_bd} and $args{right_bd}",
      #                                   -icon    => 'question',
      #                                   -type    => 'OK'
      #                                  );
   } elsif ($count_X / (length($aa_seq)) <= 0.1 && $count_X > 0) {
      ${$args{progress_bar}}->configure(-label =>"Error in file $args{filename} from module $args{module}.".
                                                 "\nThe deduced amino acid sequence contains ".
                                                 $count_X.' undefined codon(s)'.
                                                 "\nThis is a non-critical error at positions $args{left_bd} and $args{right_bd}".
                                                 "\nAnalysis will continue"
                                         );
      ${$args{main_window}}->update;#sleep(1);
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in file $args{filename} from module $args{module}.".
                  "\nThe deduced amino acid sequence contains ".
                  $count_X.' undefined codon(s)'.
                  "\nThis is a non-critical error at positions $args{left_bd} and $args{right_bd}\n\n";
      close WRITE;
   }


   #replace the first start codon with 'M' in case of 'L' or 'V'
   if ($aa_seq =~ /^L/) {
      $aa_seq =~  s /L/M/;
   } elsif ($aa_seq =~ /^V/) {
      $aa_seq =~ s /V/M/;
   }
   return ($aa_seq);
}

sub nt2aa_unchecked {
   my %args = (filename    => '',
               orientation => '',
               @_);
   my (@codons, $aa_seq);

   if ($args{sequence} eq "" || $args{sequence} =~ /^\s+$/) {
      return (0);
   }

   #convert to lowercase
   $args{sequence} = lc ($args{sequence});

   #invert sequence if selected
   if ($args{orientation} =~ /(sense|antisense)/) {
      $args{sequence} = &orient_sequence (seq         => $args{sequence},
                                          orientation => $args{orientation});
   }

   #generate codon array
   @codons = ($args{sequence} =~ m/(\w\w\w)/g);

   #translate
   foreach my $codon (@codons){
      if (defined $codon2aa{$codon}){
         $aa_seq.=$codon2aa{$codon};
      } else {
         $aa_seq.="X"; #undecipherable nucleic acid
      }
   }

   return ($aa_seq);
}

sub orient_sequence {
   my %args = (spacer => 'false',
               @_
              );

   #remove whitespaces
   if ($args{spacer} eq 'false') {
      $args{seq}=~ s/\s//gs;
   }

   if ($args{orientation} eq 'antisense') {
      $args{seq} = reverse ($args{seq});

      # IUPAC complemetary codes
      # A=T, C=G, G=C, T=A, M=K, K=M, R=Y, W=W, S=S, Y=R, B=V, V=B, D=H, H=D, N=N

      $args{seq} =~ tr/acgtmkrybvdhACGTMKRYBVDH/tgcakmyrvbhdTGCAKMYRVBHD/;
      $args{seq} = lc $args{seq};
   } else {
      $args{seq} =~ s/[^acgtmkrybvdhwsn ]/n/igs;
      $args{seq} = lc $args{seq};
   }
   return ($args{seq});
}






1;
