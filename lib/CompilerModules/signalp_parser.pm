#!/opt/ActivePerl-5.8/bin/perl

#SignalP result parser
#input arguments: main_window, progress_bar, ini_ref, auto_ini_ref, p_bar, genbank_ref, feature_lsit_ref, filename
#SignalP summary file contains results for ALL entry files.



package CompilerModules::signalp_parser;
use strict;
use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);

use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&signalp_parser);
use vars qw();


#local variables
my (%args, %seen, $value);

format ENTRY =
~~                   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$value
.

sub signalp_parser {
   my %args = @_;
   my ($total_newlines, $counter);
   #test if summary file is present
   unless (-e ${$args{ini_ref}}{signalp_results}.'/SignalP_summary.txt') {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'No file found',
                                                    -text    => "SignalP summary file is not present in folder ${$args{ini_ref}}{signalp_results}",
                                                    -bitmap  => 'info',
                                                    -buttons => ['OK'],
                                                   );
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "SignalP:  SignalP summary file is not present in folder ${$args{ini_ref}}{signalp_results}".
                     "\nReturn to main\n\n";
      close ERRORLOG;
      return;
   }

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling SignalP results",
                   label        => ''
                  );
   &show_pbar_2;

   #determine length of Summary file
   $total_newlines = `wc -l ${$args{ini_ref}}{signalp_results}\/SignalP_summary.txt`;

   #open summary file
   $counter = 0;
   open READ, ${$args{ini_ref}}{signalp_results}.'/SignalP_summary.txt';
   my $signalP_hash; #hash to check for gram+ and gram- entries for the same left_bd and orientation and choose the better one.
   while (<READ>) {
      my ($sig_prob, $cleav_prob, $cleav_position, $left_bd, $right_bd, $orientation, $source,
          $ID, $filename, $D_score, $D_cutoff, $networks);

      $counter++;
      #update status bar
      if (($counter % 20) == 0) {
         &update_pbar_2(title        => "Compiling SignalP results",
                        label        => "Parsing SignalP results for $args{filename}",
                        progress     => ($counter / $total_newlines) * 100,
                       );
      }

      #skip no hits
      next if ($_ =~ /---\t---/ || $_ =~ m/---\tNO/);
      #skip header line
      next if ($_ =~ /Entry name\s+Signal peptide probability/ || $_ =~ /Entry name\s+D-score/);
      #skip unrelated entries
      next unless ($_ =~ /^$args{filename}/);

      #start parsing
      'reset' =~ /reset/;
      if ($_ =~ m/SignalP-/) {
         next unless ($_ =~ m/YES\s*$/);
         m/(.*?)___(\d+)___(\d+)___(\d+)___(sense|antisense)\t([^\t]*?)\t([^\t]*?)\t([^\t]*?)\t([^\t]*?)\t/;
         my ($filename, $ID, $left_bd, $right_bd, $orientation, $D_score, $D_threshold, $networks, $cleav_position) = ($1, $2, $3, $4, $5, $6, $7, $8, $9);
         unless (defined $cleav_position) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error Parsing',
                                                          -text    => "Could not parse SignalP4 entry $_",
                                                          -bitmap  => 'info',
                                                          -buttons => ['OK'],
                                                         );
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "SignalP4:  Could not parse SignalP4 entry $_".
                           "\nReturn to main\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return;
         }
         $signalP_hash->{$left_bd.'_'.$orientation.'_'.$networks} = {('filename'       => $filename,
                                                                      'ID'             => $ID,
                                                                      'left_bd'        => $left_bd,
                                                                      'right_bd'       => $right_bd,
                                                                      'orientation'    => $orientation,
                                                                      'D_score'        => $D_score,
                                                                      'D_threshold'    => $D_threshold,
                                                                      'cleav_position' => $cleav_position,
                                                                      'networks'       => $networks,
                                                                      'version'        => '4'
                                                                    )};
      } else {
         m/(.*?)___(\d+)___(\d+)___(\d+)___(sense|antisense)\t([^\t]*?)\t([^\t]*?)\t([^\t]+?)\t(\S+)/;
         $filename       = $1;
         $ID             = $2;
         $left_bd        = $3;
         $right_bd       = $4;
         $orientation    = $5;
         $sig_prob       = $6;
         $cleav_prob     = $7;
         $cleav_position = $8;
         $source         = $9;

         unless (defined $source) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error Parsing',
                                                          -text    => "Could not parse SignalP3 entry $_",
                                                          -bitmap  => 'info',
                                                          -buttons => ['OK'],
                                                         );
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "SignalP3:  Could not parse SignalP3 entry $_".
                           "\nReturn to main\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return;
         }
         $source =~ s/\s+//gs;

         $signalP_hash->{$left_bd.'_'.$orientation.'_'.$source} = {('filename'       => $filename,
                                                                    'ID'             => $ID,
                                                                    'left_bd'        => $left_bd,
                                                                    'right_bd'       => $right_bd,
                                                                    'orientation'    => $orientation,
                                                                    'sig_prob'       => $sig_prob,
                                                                    'cleav_prob'     => $cleav_prob,
                                                                    'cleav_position' => $cleav_position,
                                                                    'source'         => $source,
                                                                    'version'        => '3'
                                                                   )};
      }
   }
   close READ;

   #iterate though hash and find best hit if 'both' was selected
   my %seen = ();
   foreach my $key (keys %{ $signalP_hash }) {
      next if (exists $seen{$key}); #skip processed entries

      my ($sigp_lb, $sigp_rb, $query_seq_nt, $query_seq_aa,
             $boundary_signalp, $cleavage, $note, $signal, $colour, $probability,
             $feature_signalP, $key_signalP, $alternative_source,
             $D_score);

      'reset' =~ m/reset/;
      $key =~ m/(\d+)_(\w+?)_(.+)/;
      my ($left_bd, $orientation, $source) = ($1, $2, $3); #source and networks are the same

      #process Signalp3
      if ($signalP_hash->{$key}->{'version'} == 3) {
         if ($source eq 'Gram-positive') {
            $alternative_source = 'Gram-negative';
         } else {
            $alternative_source = 'Gram-positive';
         }

         #is there an alternative source? which is is the better hit?
         if (exists $signalP_hash->{ $left_bd.'_'.$orientation.'_'.$alternative_source }) {
            if ($signalP_hash->{ $left_bd.'_'.$orientation.'_'.$alternative_source }->{ 'sig_prob' } == $signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'sig_prob' }) {
               if ($signalP_hash->{ $left_bd.'_'.$orientation.'_'.$alternative_source }->{ 'cleav_prob' } > $signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'cleav_prob' }) {
                  $source = $alternative_source;
               }
            } elsif ($signalP_hash->{ $left_bd.'_'.$orientation.'_'.$alternative_source }->{ 'sig_prob' } > $signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'sig_prob' }) {
               $source = $alternative_source;
            }
         }
         #add to seen hash
         $seen{$left_bd.'_'.$orientation.'_'.$alternative_source}++;
      }

      #add to seen hash
      $seen{$left_bd.'_'.$orientation.'_'.$source}++;

      #get nt sequence
      my ($seq_ref) = &slurp(main_window  => $args{main_window},
                             progress_bar => $args{progress_bar},
                             auto_ini_ref => $args{auto_ini_ref},
                             ini_ref      => $args{ini_ref},
                             directory    => ${$args{ini_ref}}{input_files},
                             filename     => $signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'filename' }
                            );
      if ($seq_ref eq '0') {
         &hide_pbar_2;
         return (0);
      }
      #clean up sequence
      ${$seq_ref} =~ s/^\>[^\n]*?\n//;
      ${$seq_ref} =~ s/\s//gs;

      #extract signalP nt sequence from ORF
      if ($orientation eq 'sense') {
         $query_seq_nt = substr(${$seq_ref}, ($signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'left_bd' } - 1), (($signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'cleav_position' } * 3) + 1));
      } elsif ($orientation eq 'antisense') {
         $query_seq_nt = substr(${$seq_ref}, ($signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'right_bd' } - ($signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'cleav_position' } * 3) - 1), (($signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'cleav_position' } * 3) + 1));
      }

      #translate nt sequence to aa if selected
      $query_seq_aa = nt2aa_unchecked(orientation  => $orientation,
                                      sequence     => $query_seq_nt,
                                     );

      #create boundary entry
      if ($orientation eq 'sense') {
         $sigp_lb = $signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'left_bd' };
         $sigp_rb = $signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'left_bd' } + ($signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'cleav_position' } * 3) - 1;
         $boundary_signalp = '     signalP         '.$sigp_lb.'..'.$sigp_rb."\n";
      } elsif ($orientation eq 'antisense') {
         $sigp_lb = $signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'right_bd' } - ($signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'cleav_position' } * 3) + 3;
         $sigp_rb = $signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'right_bd' };
         $boundary_signalp = '     signalP         complement('.$sigp_lb.'..'.$sigp_rb."\)\n";
      }

      #create cleavage
      $value = '/cleavage='.$signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'cleav_position' }."\n";
      open  ENTRY, '>', \$cleavage;
      write ENTRY;
      close ENTRY;

      #create note
      #process Signalp3
      if ($signalP_hash->{$key}->{'version'} == 3) {
         $value = '/note="Source type: '. $source .'; SignalP Probability: '.$signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'sig_prob' }.'; Cleavage probability: '.$signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'cleav_prob' }.'"'."\n";
      } else {
        $value = '/note="Networks: '. $source .'; D-score: '.$signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'D_score' }.'; D-score threshold: '.$signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'D_threshold' }.'"'."\n";
      }
      open  ENTRY, '>', \$note;
      write ENTRY;
      close ENTRY;

      #create signal
      $value = '/signal="'.$query_seq_aa.'"'."\n";
      open  ENTRY, '>', \$signal;
      write ENTRY;
      close ENTRY;

      #create colour
      $colour = "                     \/colour\=1\n";

      #combine feature
      $feature_signalP =  $boundary_signalp.
                          $cleavage.
                          $note.
                          $signal.
                          $colour;
      $key_signalP     =  $sigp_lb.'_'.$sigp_rb.'_signalP_'.$signalP_hash->{ $left_bd.'_'.$orientation.'_'.$source }->{ 'ID' }; #maintain original uniqe ID instead of ORFnumber

      ${$args{genbank_ref}}{$key_signalP} = $feature_signalP;
      push (@{$args{feature_list_ref}}, $key_signalP);
   }
   &hide_pbar_2;
   return;
}

1;

