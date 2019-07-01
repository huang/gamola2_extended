#!/opt/ActivePerl-5.8/bin/perl

#TMHMM result parser
#input arguments: main_window, progress_bar, gb_file, gene_model



package CompilerModules::tmhmm_parser;
use strict;
use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);

use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&tmhmm_parser);
use vars qw();


#local variables
my (%args, %seen, $value);

format ENTRY =
~~                   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$value
.

sub tmhmm_parser {
   my %args = @_;
   my (@inputfile, $max_number, $counter);
   #get result files
   opendir SEQINPUT, ${$args{ini_ref}}{tmhmm_results}.'/TMHMM_'.$args{filename} or do {
      ${$args{main_window}}->messageBox(-title   => 'Error Opening',
                                        -message => "Could not open directory ",
                                        -type    => 'OK',
                                        -icon    => 'info');
      return;
   };
   @inputfile = grep /_tmhmm$/, readdir(SEQINPUT);
   closedir SEQINPUT;

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling TMHMM results",
                   label        => ''
                  );
   &show_pbar_2;

   #max number of result files
   $max_number = $#inputfile;
   $counter = 0;

   #parse each result file
   foreach my $entry (@inputfile) {
      my ($id, $left_bd, $right_bd, $orientation, @structure,
          $tmhmm_lb, $tmhmm_rb, $tmhmm_length, $tmhmm_number, $tmhmm_prob,
          $tmhmm_signal, $tmhmm_structure, $tmhmm_start,
          $boundary_tmhmm, $helix_tmhmm, $note1, $note2, $note3, $colour,
          $feature_tmhmm, $key_tmhmm);
      #update progress bar
      $counter++;
      if (($counter % 200) == 0) {
         &update_pbar_2(title        => "Compiling TMHMM results",
                        label        => "Parsing TMHMM results for $args{filename}",
                        progress     => ($counter / $max_number) * 100,
                       );
      }

      #get boundaries
      'reset' =~ /reset/;
      $entry =~ m/(\d+)_(\d+)_(\d+)_(sense|antisense)_/;
      $id          = $1;
      $left_bd     = $2;
      $right_bd    = $3;
      $orientation = $4;

      unless (defined $orientation) {
         ${$args{main_window}}->messageBox(-title   => 'Error Parsing',
                                           -message => "Could not parse TMHMM entry $entry",
                                           -type    => 'OK',
                                           -icon    => 'info');
         &hide_pbar_2;
         return;
      }

      #get entry
      my ($file_ref) = &slurp(main_window  => $args{main_window},
                              progress_bar => $args{progress_bar},
                              auto_ini_ref => $args{auto_ini_ref},
                              ini_ref      => $args{ini_ref},
                              directory    => ${$args{ini_ref}}{tmhmm_results}.'/TMHMM_'.$args{filename},
                              filename     => $entry
                             );
      if ($file_ref eq '0') {
         ${$args{progress_bar}}->configure(-label =>"Could not grab TMHMM result $entry");
         ${$args{main_window}}->update;sleep(1);
         open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
         print WRITE "Error in TMHMM analysis->parsing".
                     "\nCould not grab TMHMM result $entry\n\n";
         close WRITE;
         next;
      }

      #get length of tmhmm region
      'reset' =~ /reset/;
      ${$file_ref} =~ m/Length:\s+(\d+)/s;
      $tmhmm_length = $1;

      #get number of tmhmm's
      'reset' =~ /reset/;
      ${$file_ref} =~ m/Number of predicted TMHs:\s+(\d+)/s;
      $tmhmm_number = $1;

      ######
      #ignore if number == 0
      next if ($tmhmm_number == 0);

      #get N-terminal probability
      'reset' =~ /reset/;
      ${$file_ref} =~ m/Total prob of N-in:\s+(\S+)\s*\n\r?/s;
      $tmhmm_prob = $1;

      #possible signal sequence?
      if (${$file_ref} =~ /POSSIBLE N-term signal sequence/) {
         $tmhmm_signal = 'POSSIBLE N-term signal sequence';
      } else {
         $tmhmm_signal = '';
      }

      #get TMHMM structural prediction
      ${$file_ref} =~ s/^.*?------\n\r?//s;
      ${$file_ref} =~ s/\#.*?\n\r?//gs;
      ${$file_ref} =~ s/$id\_$left_bd\_$right_bd\_$orientation\s+\S+\s+//gs;

      @structure = (${$file_ref} =~ /(.*?)\n\r?/gs);

      foreach my $loop (@structure) {
         next if ($loop !~ /\w+/);
         my ($location, $start, $stop);
         'reset' =~ /reset/;
         $loop =~ m/\s*(\S+)\s+(\d+)\s+(\d+)/;
         $location = $1;
         $start    = $2;
         $stop     = $3;

         unless (defined $stop) {
            ${$args{main_window}}->messageBox(-title   => 'Error Parsing',
                                              -message => "Could not parse TMHMM structure $loop in $entry",
                                              -type    => 'OK',
                                              -icon    => 'info');
            next;
         }

         unless (defined $tmhmm_start) {$tmhmm_start = $start};

         $tmhmm_structure .= $start.'-'.$location.'-'.$stop.'__';
      }
      $tmhmm_structure =~ s/__$//;
      undef @structure;

      #create boundary entry
      if (${$args{auto_ini_ref}}{tmhmm_domains} == 0) {
         if ($orientation eq 'sense') {
            $tmhmm_rb = $left_bd + ($tmhmm_length * 3) - 1;
            $tmhmm_lb = $left_bd + (($tmhmm_start) * 3) - 3;
            $boundary_tmhmm = '     TMM             '.$tmhmm_lb.'..'.$tmhmm_rb."\n";
         } elsif ($orientation eq 'antisense') {
            $tmhmm_rb = $right_bd - ($tmhmm_start * 3) + 3;
            $tmhmm_lb = $tmhmm_rb - ($tmhmm_length * 3) - 2;
            $boundary_tmhmm = '     TMM             complement('.$tmhmm_lb.'..'.$tmhmm_rb."\)\n";
         }
      } elsif (${$args{auto_ini_ref}}{tmhmm_domains} == 1) {
         ($boundary_tmhmm, $tmhmm_lb, $tmhmm_rb) = &create_domains(main_window  => $args{main_window},
                                                                   progress_bar => $args{progress_bar},
                                                                   auto_ini_ref => $args{auto_ini_ref},
                                                                   ini_ref      => $args{ini_ref},
                                                                   structure    => $tmhmm_structure,
                                                                   left_bd      => $left_bd,
                                                                   right_bd     => $right_bd,
                                                                   orientation  => $orientation
                                                                  );
         if ($boundary_tmhmm eq '0') {
            ${$args{main_window}}->messageBox(-title   => 'Error Parsing',
                                              -message => "Could not parse TMHMM structure $tmhmm_structure in $entry",
                                              -type    => 'OK',
                                              -icon    => 'info');
            next;
         }

      }

      #create helix structure
      $value = '/tmhelix='.$tmhmm_structure."\n";
      open  ENTRY, '>', \$helix_tmhmm;
      write ENTRY;
      close ENTRY;

      #create note1
      if ($tmhmm_signal =~ /\w+/) {
         $value = '/note="'.$tmhmm_signal.'"'."\n";
         open  ENTRY, '>', \$note1;
         write ENTRY;
         close ENTRY;
      } else {
         $note1 = '';
      }

      #create note2
      $value = '/note="Probability for N-terminus inside: '.$tmhmm_prob.'"'."\n";
      open  ENTRY, '>', \$note2;
      write ENTRY;
      close ENTRY;

      #create note3
      $value = '/note="Number of predicted TMHs: '.$tmhmm_number.'"'."\n";
      open  ENTRY, '>', \$note3;
      write ENTRY;
      close ENTRY;


      #create colour
      $colour = "                     \/colour\=13\n";

      #combine feature
      $feature_tmhmm =  $boundary_tmhmm.
                        $helix_tmhmm.
                        $note1.
                        $note2.
                        $note3.
                        $colour;
      $key_tmhmm     =  $tmhmm_lb.'_'.$tmhmm_rb.'_TMM_'.$id; #maintain original uniqe ID instead of ORFnumber

      ${$args{genbank_ref}}{$key_tmhmm} = $feature_tmhmm;
      push (@{$args{feature_list_ref}}, $key_tmhmm);

   }
   &hide_pbar_2;
   return(1);
}

sub create_domains {
   my %args = @_;
   my (@individual, $boundary, $tmhmm_rb, $tmhmm_lb, $counter, $left_bd, $right_bd,
       $tmhmm_domain_bd);

   #split individual tmh's
   #$args{structure} .= '__';
   @individual = ();
   @individual = ($args{structure} =~ m/(\d+-TMhelix-\d+)__/gs);

   #define initial boundary setup
   if ($#individual == 0) {
      if ($args{orientation} eq 'sense') {
         #$boundary = '     TMM             ';
         $value    = '';
      } elsif ($args{orientation} eq 'antisense') {
         #$boundary = '     TMM             ';#complement(';
         $value    = 'complement(';
      }
   } elsif ($#individual > 0) {
      if ($args{orientation} eq 'sense') {
         #$boundary = '     TMM             ';#join(';
         $value    = 'join(';
      } elsif ($args{orientation} eq 'antisense') {
         #$boundary = '     TMM             ';#complement(join(';
         $value    = 'complement(join(';
      }
   } else {
      return (0);
   }

   $counter = 0;
   foreach my $tmh (@individual) {
      my ($tmm_lb, $tmm_rb);
      'reset' =~ m/reset/;
      $tmh =~ m/(\d+)-TMhelix-(\d+)/;
      $tmm_lb = $1;
      $tmm_rb = $2;
      unless (defined $tmm_rb) {
         ${$args{main_window}}->messageBox(-title   => 'Error Parsing',
                                           -message => "Could not parse TMH $tmh in $args{structure}",
                                           -type    => 'OK',
                                           -icon    => 'info');
         return (0);
      }


      if ($args{orientation} eq 'sense') {
         $tmhmm_rb  = $args{left_bd} + ($tmm_rb * 3) - 1;
         $tmhmm_lb  = $args{left_bd} + ($tmm_lb * 3) - 3;
         $value .= $tmhmm_lb.'..'.$tmhmm_rb.',';
         #get left_bd
         if ($counter == 0) {$left_bd = $tmhmm_lb};
      } elsif ($args{orientation} eq 'antisense') {
         $tmhmm_rb  = $args{right_bd} - ($tmm_lb * 3) + 3;
         $tmhmm_lb  = $args{right_bd} - ($tmm_rb * 3) + 1;
         $value .= $tmhmm_lb.'..'.$tmhmm_rb.',';
         #get right_bd
         if ($counter == 0) {$right_bd = $tmhmm_rb};
      }
   }

   #get other boundary
   if ($args{orientation} eq 'sense') {
      $right_bd = $tmhmm_rb;
   } elsif ($args{orientation} eq 'antisense') {
      $left_bd = $tmhmm_lb;
   }

   #close boundary setup
   $value =~ s/\,$//;
   if ($#individual == 0) {
      if ($args{orientation} eq 'sense') {
         $value .= "\n";
      } elsif ($args{orientation} eq 'antisense') {
         $value .= "\)\n";
      }
   } elsif ($#individual > 0) {
      if ($args{orientation} eq 'sense') {
         $value .= "\)\n";
      } elsif ($args{orientation} eq 'antisense') {
         $value .= "\)\)\n";
      }
   }

   #format nicely
   open  ENTRY, '>', \$tmhmm_domain_bd;
   write ENTRY;
   close ENTRY;
   $tmhmm_domain_bd =~ s/^\s+/     TMM             /;

   return ($tmhmm_domain_bd, $tmhmm_lb, $tmhmm_rb);
}

1;

