#!/opt/ActivePerl-5.8/bin/perl

#Prodigal program calls
#input arguments: progress_bar, directory, filename, ini_ref, auto_ini_ref


package GeneModel::prodigal;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&prodigal);
use vars qw();

use Basics::MesgBox;
use initialise::read_me qw(:DEFAULT);
use Cwd;
use File::Find;

#local variables
my (%args, $tl);




sub prodigal {
   my %args = @_;
   my $dna_seq;
   my $prodigal_coord      = ${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.prodigal.coord';
   my $prodigal_parameters = '-q ';
   my $str                 = '';
   my $header              = '';
   my $res                 = '';
   my $executable;

format FASTA =
^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<~~
$dna_seq
.

   #generating status message
   ${$args{progress_bar}}->configure(-label=>"Generating Prodigal gene model for $args{input_file}.");
   ${$args{main_window}}->update;

   #setup optional parameters
   if (${$args{auto_ini_ref}}{prodigal_closed_ends} == 1) {
      $prodigal_parameters .= ' -c';
   }
   #
   'reset' =~ m/reset/;
   ${$args{auto_ini_ref}}{prodigal_translation_table} =~ m/^(\d+)/;
   $prodigal_parameters .= ' -g '.$1;
   #
   if (${$args{auto_ini_ref}}{prodigal_masked_seq} == 1) {
      $prodigal_parameters .= ' -m';
   }
   #
   if (${$args{auto_ini_ref}}{prodigal_bypass} == 1) {
      $prodigal_parameters .= ' -n';
   }
   #
   if (${$args{auto_ini_ref}}{prodigal_procedure} eq 'Single') {
      $prodigal_parameters .= ' -p single';
   } elsif (${$args{auto_ini_ref}}{prodigal_procedure} eq 'Meta') {
      $prodigal_parameters .= ' -p meta';
   }


   #reading input file
   my ($file_ref) = &slurp(main_window => \$args{progress_bar},
                           directory   => $args{directory},
                           filename    => $args{input_file}
                          );
   #parsing input file
   if (${$file_ref} =~ /^\>/) {
      ${$file_ref} =~ m/^(\>[^\n]+)\n(.+)/s;
      $header = $1;
      $dna_seq = $2;
   } else {
      $header = '>'.$args{input_file};
      $dna_seq = ${$file_ref};
   }

   #cleanup and format DNA sequence
   $dna_seq =~ s/\s//gs;
   $dna_seq =~ s/[^acgtmkrybvdhwsn]/n/igs;
   $dna_seq = lc($dna_seq);

   open  FASTA, '>', \$str or die;
   write FASTA;
   close FASTA;

   #write tempfile for glimmer
   open  WRITE, "+>".${$args{ini_ref}}{tempfile} or return (0);
   print WRITE $header."\n".$str."\n";
   close WRITE;

   #setup override for Prodigal if selected
   if (length($dna_seq) < 100000 && ${$args{auto_ini_ref}}{prodigal_override} == 1) {
      $prodigal_parameters =~ s/ -p single/ -p meta/;
   }

   #prepare Prodigal run
   find (sub {$executable = $File::Find::name if (-f && $_ =~ m/^prodigal/i) }, ${$args{ini_ref}}{prodigal});
   `$executable $prodigal_parameters -o $prodigal_coord -f sco -i ${$args{ini_ref}}{tempfile}`;

   if (-s $prodigal_coord < 1) {
      my $error_msg = ${$args{main_window}}->Dialog(-title  => 'Error',
                                                    -text   => "Error while creating Prodigal output file for $args{input_file}".
                                                               "\nin directory ${$args{ini_ref}}{genemodel_output}",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg->Show();
      open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error while creating Prodigal output for file $args{input_file}".
                     "\nin directory ${$args{ini_ref}}{genemodel_output}\n\n";
      close ERRORLOG;
      return (1);
   }

   #reformat Prodigal output
   &prodigal_format(main_window   => $args{main_window},
                    auto_ini_ref  => $args{auto_ini_ref},
                    ini_ref       => $args{ini_ref},
                    directory     => ${$args{ini_ref}}{genemodel_output},
                    filename      => $args{input_file}.'.prodigal.coord'
                   );



   unlink ${$args{ini_ref}}{tempfile};
   return (1);
}

sub prodigal_format {
   my %args                = @_;
   my $prodigal_mod        = $args{directory}.'/'.$args{filename}.'.mod';
   my $increase_orf_number = 1;
   my @sq_array            = ();

   my ($gene_model_ref) = &slurp(main_window   => $args{main_window},
                                 directory     => $args{directory},
                                 filename      => $args{filename}
                                );


   #test if putative genes were predicted, if not, create empty file and return
   unless (${$gene_model_ref} =~/\>\d/s) {
      open WRITE, "+>$prodigal_mod";
      close WRITE;
      return;
   }

   #make the ORF numbers sequentially
   @sq_array = split/\n/,${$gene_model_ref};
   open WRITE, "+>$prodigal_mod";
   foreach (@sq_array) {
      next if ($_ !~ /^\>/);
      #get all data
      'reset' =~ m/reset/;
      m/^\>\d+_(\d+)_(\d+)_([\+\-])/;
      my ($start, $stop, $frame) = ($1, $2, $3);
      if ($frame eq '+') {
         print WRITE '  '.$increase_orf_number.'  '.$start.'  '.($stop - 3).'  ['.$frame.'1  '."\n";
      } else {
         print WRITE '  '.$increase_orf_number.'  '.$stop.'  '.($start + 3).'  ['.$frame.'1  '."\n"
      }
      $increase_orf_number++;
   }
   close WRITE;
}


1;