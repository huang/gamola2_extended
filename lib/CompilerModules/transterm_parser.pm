#!/opt/ActivePerl-5.8/bin/perl

#Transterm parser calls and reformatting of results
#input arguments: main_window, progress_bar, auto_ini_ref, ini_ref,
#                 directory, input_array_ref, output_folder, input_type

package CompilerModules::transterm_parser;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&transterm_parser);
use vars qw();

use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);
use Cwd;

#local variables
my (%args, $value);

format ENTRY =
~~                   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$value
.

sub transterm_parser {
   my %args = @_;
   my (@terminators);

   #skip if no results
   unless (-e ${$args{ini_ref}}{transterm_results}.'/'.$args{filename}.'_transterm') {return (1)};

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling Transterm results",
                   label        => 'Generating terminator model'
                  );
   &show_pbar_2;

   #slurp up result file
   my ($seq_ref) = &slurp(main_window  => $args{main_window},
                          progress_bar => $args{progress_bar},
                          auto_ini_ref => $args{auto_ini_ref},
                          ini_ref      => $args{ini_ref},
                          directory    => ${$args{ini_ref}}{transterm_results},
                          filename     => $args{filename}.'_transterm'
                         );
   if ($seq_ref eq '0') {
      &hide_pbar_2;
      return (0);
   }

   #cleanup results
   ${$seq_ref} =~ s/^TransTerm v2\.0.*?SEQUENCE.*?\n\n//s;
   #build terminator array
   @terminators = (${$seq_ref} =~ m/  TERM([^\n]+\n[^\n]+)\n/gs);
   undef $seq_ref;

   #generate terminator entries and push into start array
   my $counter = 1;
   foreach my $term (@terminators) {
      #update status bar
      if (($counter % 20) == 0) {
         &update_pbar_2(title        => "Compiling Transterm results",
                        label        => "Generating terminator model",
                        progress     => ($counter / $#terminators) * 100,
                       );
      }
      $counter++;
      my ($term_feature, $format, $term_number, $start, $stop, $orientation, $class,
          $score, $feature, $term_seq, $term_lb, $term_rb,
          $boundary_terminator, $label, $term_score, $note, $seq, $colour,
          $feature_transterm, $key_transterm);
      $term =~ m/\s*(\d+)\s+(\d+)\s+-\s+(\d+)\s+([\+-])\s+(\w+)\s+(\d+)[^\|]+\|([^\n]*?)\n(.*)/s;
      $term_number = $1;
      $start       = $2;
      $stop        = $3;
      $orientation = $4;
      $class       = $5;
      $score       = $6;
      $feature     = $7;
      $term_seq    = $8;

      $term_seq =~ s/ {2,}/ /g;
      $term_seq =~ s/^\s+//;
      $term_seq =~ s/\s+$//;
      $feature  =~ s/gap 1/gaps in hairpin-stem/;
      $feature  =~ s/^\s+//;

      #skip terminators residing in ORFs
      if ($class =~ /[gG]/) {next};

      #increase ID
      $args{counter}++;

      #create boundary entry
      if ($orientation eq '-') {
         $term_lb = $stop;
         $term_rb = $start;
         $boundary_terminator = '     terminator      complement('.$stop.'..'.$start."\)\n";
         #also reverse terminator sequence
         $term_seq = &orient_sequence(seq         => $term_seq,
                                      orientation => 'antisense',
                                      spacer      => 'true'
                                     );
      } elsif ($orientation eq '+') {
         $term_lb = $start;
         $term_rb = $stop;
         $boundary_terminator = '     terminator      '.$start.'..'.$stop."\n";
      }

      #create label
      $value = '/label=TERM'.$term_number."\n";
      open  ENTRY, '>', \$label;
      write ENTRY;
      close ENTRY;

      #create score
      $value = '/score='.$score."\n";
      open  ENTRY, '>', \$term_score;
      write ENTRY;
      close ENTRY;

      #create note
      if ($feature =~ /\w+/) {
         $value = '/note="'.$feature.'"'."\n";
         open  ENTRY, '>', \$note;
         write ENTRY;
         close ENTRY;
      } else {
         $note = '';
      }

      #create term_seq
      if ($term_seq =~ /\w+/) {
         $value = '/note="'.$term_seq.'"'."\n";
         open  ENTRY, '>', \$seq;
         write ENTRY;
         close ENTRY;
      } else {
         $seq = '';
      }

      #create colour
      $colour = "                     \/colour=6\n";

      #combine feature
      $feature_transterm =  $boundary_terminator.
                            $label.
                            $term_score.
                            $note.
                            $seq.
                            $colour;
      $key_transterm     =  $term_lb.'_'.$term_rb.'_terminator_'.$args{counter}; #create  uniqe ID instead of ORFnumber

      ${$args{genbank_ref}}{$key_transterm} = $feature_transterm;
      push (@{$args{feature_list_ref}}, $key_transterm);

   }
   &hide_pbar_2;
   return (1);
}

1;