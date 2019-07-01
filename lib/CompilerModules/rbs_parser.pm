#!/opt/ActivePerl-5.8/bin/perl

#RBSfinder result parser
#input arguments: main_window, progress_bar, gb_file, gene_model
#takes the RBSfinder output file and adds RBSs as features



package CompilerModules::rbs_parser;
use strict;
use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);

use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&rbs_parser);
use vars qw();


#local variables
my (%args, %seen, $value);

format ENTRY =
~~                   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$value
.

sub rbs_parser {
   my %args             = @_;
   my $counter          = 1;
   my ($feature_boundary, $key_rbs, $colour, $feature_product, $feature_rbs);

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling RBSfinder results",
                   label        => ''
                  );
   &show_pbar_2;

   #slurp RBSfinder results file
   my $seq_ref = &slurp(main_window  => $args{main_window},
                        progress_bar => $args{progress_bar},
                        auto_ini_ref => $args{auto_ini_ref},
                        ini_ref      => $args{ini_ref},
                        directory    => ${$args{ini_ref}}{genemodel_output},
                        filename     => $args{filename}.'.rbsfinder'
                      );
   if ($seq_ref eq '0') {
      &hide_pbar_2;
      return (0);
   }

   #parse each line
   my @RBS_line = split/\n/, $$seq_ref;
   foreach my $line (@RBS_line) {
      my ($ORF_start, $RBS_left_bd, $RBS_pattern);
      $line =~ m/^\s+\S+\s+(\d+)\s+\d+\s+(\S+)\s+(\S+)/;
      ($ORF_start, $RBS_pattern, $RBS_left_bd) = ($1, $2, $3);

      next unless (defined $RBS_left_bd);
      next if ($RBS_pattern eq '---');

      #feature boundary
      if ($ORF_start > $RBS_left_bd) { #sense
         $feature_boundary = $RBS_left_bd.'..'.($RBS_left_bd + length($RBS_pattern) - 1)."\n";
      } else { #antisense
         $feature_boundary = 'complement('.($RBS_left_bd - length($RBS_pattern) + 1).'..'.$RBS_left_bd.')'."\n";
      }
      #create note
      $value = '/note="'.$RBS_pattern.'"'."\n";
      open  ENTRY, '>', \$feature_product;
      write ENTRY;
      close ENTRY;
      #create colour
      $colour = "                     \/colour\=13\n";
      #combine feature
      $feature_rbs =  '     RBS             '.$feature_boundary.
                      $feature_product.
                      $colour;
      $key_rbs     = $RBS_left_bd.'_'.($RBS_left_bd + length($RBS_pattern)).'_rbs_match_'.$counter;

      #increase internal counter, will later be replaed with global feature counter $args{counter}
      $counter++;

      ${$args{genbank_ref}}{$key_rbs} = $feature_rbs;
      push (@{$args{feature_list_ref}}, $key_rbs);
   }

   &hide_pbar_2;
   return (1);

}

1;

