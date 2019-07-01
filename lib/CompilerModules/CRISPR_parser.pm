#!/opt/ActivePerl-5.8/bin/perl

#CRISPR result parser
#input arguments: main_window, progress_bar, gb_file, gene_model



package CompilerModules::CRISPR_parser;
use strict;
use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);

use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&CRISPR_parser);
use vars qw();


#local variables
my (%args, %seen, $value);

format ENTRY =
~~                   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$value
.

sub CRISPR_parser {
   my %args = @_;
   my (@CRISPR, @structure, $CRISPR_structure, $CRISPR_numerical);

   #skip if no results
   unless (-e ${$args{ini_ref}}{CRISPR_results}.'/'.$args{filename}.'.CRISPR') {return (1)};

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling CRISPR results",
                   label        => 'Generating CRISPR model'
                  );
   &show_pbar_2;

   #slurp up result file
   my ($file_ref) = &slurp(main_window  => $args{main_window},
                           progress_bar => $args{progress_bar},
                           auto_ini_ref => $args{auto_ini_ref},
                           ini_ref      => $args{ini_ref},
                           directory    => ${$args{ini_ref}}{CRISPR_results},
                           filename     => $args{filename}.'.CRISPR'
                          );
   if ($file_ref eq '0') {
      &hide_pbar_2;
      return (0);
   }

   #cleanup results
   ${$file_ref} =~ s/^ORGANISM:.*?Bases:\s+\d+.*?\n\n/\n/s;
   ${$file_ref} =~ s/Time to find repeats.*//s;

   #build CRISPR array
   @CRISPR = (${$file_ref} =~ m/\n\n(.+?)\n\n/gs);

   #generate CRISPR entries and push into start array
   my $counter = 1;
   foreach my $entry (@CRISPR) {
      next if ($entry !~ /\w+/);
      my ($boundary_CRISPR);
      my ($overview, $note1, $note2, $feature_CRISPR, $key_CRISPR, $colour);
      #update status bar
      if (($counter % 20) == 0) {
         &update_pbar_2(title        => "Compiling CRISPR results",
                        label        => "Generating CRISPR model",
                        progress     => ($counter / $#CRISPR) * 100,
                       );
      }
      $counter++;

      #capture global boundaries
      'reset' =~ m/reset/;
      $entry  =~ m/Range\:\s+(\d+)\s+\-\s+(\d+).+Repeats\:\s+(\d+)\s+Average\s+Length\:\s+(\d+)\s+Average\s+Length\:\s+(\d+)\s*\n?\r?/s;
      my ($left_bd, $right_bd, $repeats, $av_repeat_length, $av_spacer_length) = ($1, $2, $3, $4, $5);

      unless (defined $av_spacer_length) {
         ${$args{main_window}}->messageBox(-title   => 'Error Parsing',
                                           -message => "Could not parse CRISPR entry $entry",
                                           -type    => 'OK',
                                           -icon    => 'info');
         &hide_pbar_2;
         return (0);
      }

      #capture local boundaries
      @structure = ();
      @structure = ($entry =~ /\n(\d+[^\n]+)/gs);

      $CRISPR_structure = '';
      $CRISPR_numerical = '';
      foreach my $structure (@structure) {
         next if ($structure !~ /\w+/);
         my ($location, $length_repeat, $length_spacer, $spacer, $repeat);
         'reset' =~ /reset/;
         $structure =~ m/\s*(\d+)\s+(\S+)\s+(\S+)/;
         ($location, $repeat, $spacer) = ($1, $2, $3);
         $length_repeat = length($repeat);
         $length_spacer = length($spacer);

         #last repeat?
         unless (defined $spacer && $structure =~ m/\[/) {
            'reset' =~ /reset/;
            $structure =~ m/\s*(\d+)\s+(\S+)/;
            ($location, $repeat) = ($1, $2);
            if (defined ($repeat)) {
               $spacer        = 0;
               $length_spacer = 0;
               $length_repeat = length($repeat);
            }
         }

         #if not, assume error
         unless (defined $length_spacer) {
            ${$args{main_window}}->messageBox(-title   => 'Error Parsing',
                                              -message => "Could not parse CRISPR structure $structure in $args{filename}",
                                              -type    => 'OK',
                                              -icon    => 'info');
            next;
         }

         $CRISPR_structure .= $location.'-'.$length_repeat.'-'.$length_spacer.'-'.$repeat.'-'.$spacer.'__';
         $CRISPR_numerical .= $location.'-Rep:'.$length_repeat.'-Spac:'.$length_spacer.'__';
      }
      $CRISPR_structure =~ s/__$//;
      undef @structure;

      ($boundary_CRISPR) = &create_domains(main_window  => $args{main_window},
                                           progress_bar => $args{progress_bar},
                                           auto_ini_ref => $args{auto_ini_ref},
                                           ini_ref      => $args{ini_ref},
                                           structure    => $CRISPR_structure,
                                          );
      if ($boundary_CRISPR eq '0') {
         ${$args{main_window}}->messageBox(-title   => 'Error Parsing',
                                           -message => "Could not parse CRISPR structure $CRISPR_structure in $entry",
                                           -type    => 'OK',
                                           -icon    => 'info');
         next;
      }

      #increase ID
      $args{counter}++;

      #create overview structure
      $CRISPR_numerical =~ s/\_\_$//;
      $value = '/CRISPR='.$CRISPR_numerical."\n";
      open  ENTRY, '>', \$overview;
      write ENTRY;
      close ENTRY;

      #create note1: statistics
      $value = '/note="Overall CRISPR Range: '.$left_bd.' to '.$right_bd.'; '.
               'Number of Repeats: '.$repeats.'; '.
               'Average Repeat Length: '.$av_repeat_length.'; '.
               'Average Spacer Length: '.$av_spacer_length.
               '"'."\n";
      open  ENTRY, '>', \$note1;
      write ENTRY;
      close ENTRY;

      #create note2: complete Model
      if ($CRISPR_structure =~ /\w+/) {
         #cleanup entry
         $CRISPR_structure =~ s/-0//g;
         $CRISPR_structure =~ s/\d+\-\d+\-(\d+\-)?//g;
         $value = '/note="'.$CRISPR_structure.'"'."\n";
         open  ENTRY, '>', \$note2;
         write ENTRY;
         close ENTRY;
      } else {
         $note2 = '';
      }

      #create colour
      $colour = "                     \/colour\=9\n";

      #combine feature
      $feature_CRISPR = $boundary_CRISPR.
                        $overview.
                        $note1.
                        $note2.
                        $colour;
      $key_CRISPR     = $left_bd.'_'.$right_bd.'_CRISPR_'.$args{counter}; #create  uniqe ID instead of ORFnumber

      ${$args{genbank_ref}}{$key_CRISPR} = $feature_CRISPR;
      push (@{$args{feature_list_ref}}, $key_CRISPR);
      undef $CRISPR_structure;
      undef $CRISPR_numerical;

   }
   undef $file_ref;
   &hide_pbar_2;
   return(1);
}

sub create_domains {
   my %args = @_;
   my (@individual, $boundary, $CRISPR_domain_bd);

   #split individual CRISPR's
   $args{structure} = '__'.$args{structure}.'__';
   @individual = ();
   @individual = ($args{structure} =~ m/__([^\_]+)/gs);

   #define initial boundary setup
   $value    = 'join(';

   foreach my $entry (@individual) {
      my ($lb, $repeat_length, $CRISPR_rb, $CRISPR_lb);
      'reset'        =~ m/reset/;
      $entry         =~ m/(\d+)-(\d+)/;
      $lb            = $1;
      $repeat_length = $2;

      if ($repeat_length == 0) {
      }

      unless (defined $repeat_length) {
            ${$args{main_window}}->messageBox(-title   => 'Error Parsing',
                                              -message => "Could not parse TMH $entry in $args{structure}",
                                              -type    => 'OK',
                                              -icon    => 'info');
            return (0);
      }

      $CRISPR_rb  = $lb;
      $CRISPR_lb  = $lb + $repeat_length;
      $value .= $CRISPR_rb.'..'.$CRISPR_lb.',';
   }

   #close boundary setup
   $value =~ s/\,$//;
   $value .= "\)\n";

   #format nicely
   open  ENTRY, '>', \$CRISPR_domain_bd;
   write ENTRY;
   close ENTRY;
   $CRISPR_domain_bd =~ s/^\s+/     CRISPR          /;

   return ($CRISPR_domain_bd);
}

1;

