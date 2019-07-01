#!/opt/ActivePerl-5.8/bin/perl

#Verify final gene model before analyses start
#input arguments: progress_bar, directory, filename, ini_ref, auto_ini_ref


package GeneModel::verify_genemodel;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&verify_gene_model);
use vars qw();

use Basics::MesgBox;
use initialise::read_me qw(:DEFAULT);

#local variables
my (%args, );

sub verify_gene_model {
   #read sequence input directory
   $log_entry = "\nVerifying gene model".
                "\nReading sequence input files";
   print $log_entry;
   print WRITE_LOG $log_entry; $log_entry="";
   my (@seqfile) = &read_dir($seqinput);

   #make sure each input file has a GAMOLA gene model
   $log_entry = "\nReading gene models for each input file";
   print $log_entry;
   print WRITE_LOG $log_entry; $log_entry="";
   chdir $gloutput or die;
   foreach my $input (@seqfile) {
      unless (-e $input.'.glcoord.mod') {
         $log_entry = "\nGene model not found for $input";
         print $log_entry;
         print WRITE_LOG $log_entry; $log_entry="";
         return (0);
      }
   }

   #processing each sequence file
   foreach my $sequence (@seqfile) {
      my $fasta = "";
      #read input sequence file
      $fasta = &slurp($seqinput, $sequence);

      #clean up FASTA sequence
      $fasta =~ s/^\>[^\n]*?\n//s;
      $fasta =~ s/\s//gs;
      $fasta =~ s/[^acgtnACGTN]/n/gs;

      #read gene model
      my $genemodel = "";
      my @gene_model = ();
      $genemodel = &slurp($gloutput, $sequence.'.glcoord.mod');

      @gene_model = split/\n/,$genemodel;

      #check each ORF
      foreach my $orf (@gene_model) {
         if ($orf !~ /\w+/) {next}; #skip empty lines
         my $number = "";
         my $start = "";
         my $stop = "";
         my $aa_seq = "";
         $orf =~ m/\s+(\d+)\s+(\d+)\s+(\d+)\s+/;
         $number = $1;
         $start = $2;
         $stop = $3;

         if ($start < $stop) {
            my $seq_start = $start - 1;
            my $seq_stop = $stop + 3;
            $aa_seq = &seq2protein(substr($fasta,$seq_start,$seq_stop-$seq_start), $start, $stop);
         } elsif ($start > $stop) {
            my $seq_start = $start;
            my $seq_stop = $stop-4;
            $aa_seq = &seq2protein(substr($fasta,$seq_stop,$seq_start-$seq_stop), $start, $stop);
         } else {
            $log_entry = "\n\nCould not grab start and stop positions from\n$orf\n";
            print $log_entry;
            print WRITE_LOG $log_entry; $log_entry="";
            return (0);
         }
         if ($aa_seq eq '0') {
            $log_entry = "\nProblem with ORF $orf\n\n";
            print $log_entry;
            print WRITE_LOG $log_entry; $log_entry="";
            return (0);
         }
      }
   }
   return (1);
}
