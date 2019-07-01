#!/opt/ActivePerl-5.8/bin/perl

#rRNA result parser
#input arguments: main_window, progress_bar, gb_file, gene_model

package CompilerModules::rrna_parser;
use strict;
use initialise::read_me  qw(:DEFAULT);
use Basics::progress_bar qw(:DEFAULT);
use File::Copy;

use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&rrna_parser);
use vars qw();


#local variables
my (%args, $value);

format ENTRY =
~~                   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$value
.

sub rrna_parser {
   my %args = @_;
   my (@res_files, $max_count, $max_seq_length, $count, $file, $file_ref, $local_threshold, $domains);
   my ($rRNA);

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling ribosomal RNA results",
                   label        => ''
                  );
   &show_pbar_2;

   #define max seq length
   $file_ref = slurp(main_window => $args{main_window},
                     directory   => ${$args{ini_ref}}{input_files},
                     filename    => $args{filename}
                    );
   if ($file_ref eq '0') {return (0)};
   ${$file_ref} =~ s/^\>[^\n]*?\n//s;
   $max_seq_length = length(${$file_ref});
   undef $file_ref;

   #get all result files for input file
   opendir SEQINPUT, ${$args{ini_ref}}{rrna_results};
   @res_files = grep /^$args{filename}\-.*?\.bl2seq$/, readdir(SEQINPUT);
   closedir SEQINPUT;
   $max_count = $#res_files + 2;
   $count = 1;

   #parse through rRNA results
   foreach my $res_file (@res_files) {
      my ($rRNA_db, $query_position, $qu_left_bd, $qu_right_bd, $file_ref, $length,
          $truncated);
      my (@candidates, @iterate);

      &update_pbar_2(title        => "Compiling ribosomal RNA results",
                     label        => "Parsing result files $count of $max_count",
                     progress     => ($count / $max_count) * 100,
                    );
      $count++;

      #capture rRNA db and start position
      'reset' =~ m/reset/;
      $res_file =~ m/^$args{filename}\-(.+?)\-(\d+)\.bl2seq/;
      if ($2) {($rRNA_db, $query_position) = ($1, $2)} else {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse $res_file for rRNA.",
                                                       -bitmap  => 'info',
                                                       -buttons => ['OK']
                                                       );
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing rRNA results for file $res_file".
                        "\nCould not parse database and query position.\n\n";
         close ERRORLOG;

         next;
      };

      #skip if selected
      next if (${$args{auto_ini_ref}}{rrna_5S}  == 0 && $rRNA_db eq '5S_rRNA' );
      next if (${$args{auto_ini_ref}}{rrna_16S} == 0 && $rRNA_db eq '16S_rRNA');
      next if (${$args{auto_ini_ref}}{rrna_23S} == 0 && $rRNA_db eq '23S_rRNA');

      #read input files and create main hash
      $file_ref = slurp(main_window => $args{main_window},
                        directory   => ${$args{ini_ref}}{rrna_results},
                        filename    => $res_file
                       );
      if ($file_ref eq '0') {return (0)};

      #next if no hits
      next if (${$file_ref} =~ /\*\*\*\*\* No hits found \*\*\*\*\*\*/s);

      #get query left and right boundary
      'reset'      =~ m/reset/;
      ${$file_ref} =~ m/^\s*Query=(.+?)\(.+letters\)/s;
      my $header   = $1;
      $header      =~ s/[\n\r]//gs;
      'reset'      =~ m/reset/;
      $header      =~ m/$args{filename}.+?(\d+)\-(\d+)\s*$/s;

      if ($2) {($qu_left_bd, $qu_right_bd) = ($1, $2)} else {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not grab right boundary from $res_file for rRNA.",
                                                       -bitmap  => 'info',
                                                       -buttons => ['OK']
                                                       );
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing rRNA results for file $res_file".
                        "\nCould not grab right boundary from \n...$header ...\n${$file_ref}.\n\n";
         close ERRORLOG;

         next;
      };
      undef $header;

      #get subject length
      'reset' =~ m/reset/;
      ${$file_ref} =~ m/>.*?Length = (\d+)/s;
      if ($1) {$length = $1} else {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not grab subject length from $res_file for rRNA.",
                                                       -bitmap  => 'info',
                                                       -buttons => ['OK']
                                                       );
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing rRNA results for file $res_file".
                        "\nCould not grab subject length from \n${$file_ref}.\n\n";
         close ERRORLOG;

         next;
      };

      #get all candidate hits (there could be adjacent tandem rRNAs)
      @candidates = split/Score = /, ${$file_ref};
      shift @candidates;

      #iterate through all candidate hits
      foreach my $entry (@candidates) {
         my ($score, $evalue, $alignment_length, $identities, $orientation, $qu_local, $sb_local);
         #get Length, score, evalue, alignment_length, identities [%], orientation and begin of local alignment
         'reset' =~ m/reset/;
         $entry =~ m/\s*(\d+).*?Expect = (\S+).*?Identities = \d+\/(\d+)\s+\((\d+)\%\).*?Strand = (\S+? \/ \S+?)\s.*?Query:\s+(\d+)\s.*?Sbjct:\s+(\d+)\s/s;
         if ($7) {
            ($score, $evalue, $alignment_length, $identities, $orientation, $qu_local, $sb_local) = ($1, $2, $3, $4, $5, $6, $7, $8);
         } else {
             my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                           -text    => "Could not grab parameters from $res_file for rRNA.",
                                                           -bitmap  => 'info',
                                                           -buttons => ['OK']
                                                           );
             $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
             $error_msg-> Show();
             open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
             print ERRORLOG "Error parsing rRNA results for file $res_file".
                            "\nCould not grab parameters from $entry.\n\n";
             close ERRORLOG;
            next;
         };

         #test for sufficient alignment length/score for at least 5S and 16S rRNA

         next if ($rRNA_db eq '23S_rRNA' && $alignment_length < 250);
         next if ($rRNA_db eq '16S_rRNA' && $alignment_length < 100);
         next if ($rRNA_db eq '5S_rRNA'  && $score            < 40);

         #define a few paramters
         if ($orientation =~ /minus/i) {
            $orientation = 'antisense';
            $qu_left_bd  = $qu_left_bd + $qu_local - ($length - $sb_local) - 6;
            $qu_right_bd = $qu_left_bd + $length - 6;
         } else {
            $orientation = 'sense';
            $qu_left_bd  = $qu_left_bd + $qu_local - $sb_local + 1;
            $qu_right_bd = $qu_left_bd + $length + 1;
         }

         #test if crossing sequence boundaries
         $truncated = 0;
         if ($qu_left_bd < 1) {
            $qu_left_bd = 1;
            $truncated  = 1;
         }

         if ($qu_left_bd > $max_seq_length) {
            $qu_left_bd = $max_seq_length - 3;
            $truncated = 1;
         }

         if ($qu_right_bd > $max_seq_length) {
            $qu_right_bd = $max_seq_length - 3;
            $truncated = 1;
         }

         #skip if feature is too small
         next if ($qu_right_bd <= $qu_left_bd);

         #skip hit if better hit for same rRNA db with same start position already exists
         next if (exists $rRNA->{$rRNA_db}->{$query_position} && $rRNA->{$rRNA_db}->{$query_position}->{'score'} >= $score);

         #test if overlapping tandem hit exists
         my $skip = 0;
         foreach my $test (@iterate) {
            my ($test_left, $test_right, $test_score);
            'reset' =~ m/reset/;
            $test =~ m/(\d+)_(\d+)_(.+)/;
            ($test_left, $test_right, $test_score) = ($1, $2, $3);
            if (($test_left < $qu_left_bd  && $test_right > $qu_left_bd ) ||
                ($test_left < $qu_right_bd && $test_right > $qu_right_bd) ||
                ($test_left < $qu_left_bd  && $test_right > $qu_right_bd)) {
               #keep better hit
               if ($test_score < $score) {
                  #delete
                  delete $rRNA->{ $rRNA_db }->{ $test_left }
               } elsif ($test_score >= $score) {
                  $skip = 1;
                  last;
               }
            }
         }
         next if ($skip == 1);
         #no skip, add to array
         push (@iterate, $qu_left_bd.'_'.$qu_right_bd.'_'.$score);

         #set up main hash
         $rRNA->{ $rRNA_db }->{ $qu_left_bd } = { 'db'               => $rRNA_db,
                                                  'score'            => $score,
                                                  'evalue'           => $evalue,
                                                  'alignment_length' => $alignment_length,
                                                  'subject_length'   => $length,
                                                  'identities'       => $identities,
                                                  'orientation'      => $orientation,
                                                  'left_bd'          => $qu_left_bd,
                                                  'right_bd'         => $qu_right_bd,
                                                  'truncated'        => $truncated
                                               };
      }
      undef @iterate;
      undef @candidates;
   }

   #remove overlapping rRNAs
   &update_pbar_2(title        => "ribosomal RNA analysis",
                  label        => "Removing overlapping rRNAs",
                  progress     => 1,
                 );
   $max_count = &remove_global_overlaps(rRNA_ref => \$rRNA,
                                        counter  => $1);
   if ($max_count < 1) {$max_count = 1};

   &update_pbar_2(title        => "ribosomal RNA analysis",
                  label        => "Compiling ribosomal RNA results",
                  progress     => 1,
                 );

   #generate new entries
   $count = 1;
   foreach my $rRNA_db (keys %{ $rRNA } ) {
      foreach my $qu_left_bd (keys %{ $rRNA->{$rRNA_db} } ) {
         #skip hash references
         next if ($qu_left_bd =~ /\D+/);

         #update status bar
         if (($count % 10) == 0) {
            &update_pbar_2(title        => "Compiling ribosomal RNA results",
                           label        => "Compiling ribosomal RNA results for $args{filename}",
                           progress     => ($count / $max_count) * 100,
                          );
         }
         $count++;

         my ($boundary_rrna, $label, $score, $colour, $note,
             $feature_rRNA, $key_rRNA, $rrna_product);

         #make sure entry is defined
         if (defined $rRNA->{$rRNA_db}->{$qu_left_bd}->{'orientation'}) {
            #increase ID
            $args{counter}++;

            #create boundary entry
            if ($rRNA->{$rRNA_db}->{$qu_left_bd}->{'orientation'} eq 'sense') {
               $boundary_rrna  = '     rRNA            '.$rRNA->{$rRNA_db}->{$qu_left_bd}->{'left_bd'}.'..'.$rRNA->{$rRNA_db}->{$qu_left_bd}->{'right_bd'}."\n";
            } elsif ($rRNA->{$rRNA_db}->{$qu_left_bd}->{'orientation'} eq 'antisense') {
               $boundary_rrna  = '     rRNA            complement('.$rRNA->{$rRNA_db}->{$qu_left_bd}->{'left_bd'}.'..'.$rRNA->{$rRNA_db}->{$qu_left_bd}->{'right_bd'}."\)\n";
            } else {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Could not set boundaries for rRNA $rRNA_db at position $qu_left_bd in file $args{filename}.",
                                                             -bitmap  => 'info',
                                                             -buttons => ['OK']
                                                             );
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error parsing rRNA results for file $args{filename}".
                              "\nCould not set boundaries for rRNA $rRNA_db at position $qu_left_bd in file $args{filename}.\n\n";
               close ERRORLOG;

               &hide_pbar_2;
               return (0);
            };

            #create rRNA label
            my $db = $rRNA->{$rRNA_db}->{$qu_left_bd}->{'db'};
            $db =~ s/_rRNA//;
            $value         = '/label='.$db."\n";
            open  ENTRY, '>', \$label;
            write ENTRY;
            close ENTRY;

            #add product
            if ($rRNA->{$rRNA_db}->{$qu_left_bd}->{'truncated'} == 1) {
               $value = '/product="'.$db.'_truncated'."\"\n";
            } else {
               $value = '/product="'.$db."\"\n";
            }
            open  ENTRY, '>', \$rrna_product;
            write ENTRY;
            close ENTRY;

            #create score
            $value = '/score='.$rRNA->{$rRNA_db}->{$qu_left_bd}->{'score'}."\n";
            open  ENTRY, '>', \$score;
            write ENTRY;
            close ENTRY;

            #set notes
            $value = '/note="Alignment length: '.$rRNA->{$rRNA_db}->{$qu_left_bd}->{'alignment_length'}.'nt out of '.$rRNA->{$rRNA_db}->{$qu_left_bd}->{'subject_length'}.'nt with '.
                     $rRNA->{$rRNA_db}->{$qu_left_bd}->{'identities'}.'% identity"'."\n";
            open  ENTRY, '>', \$note;
            write ENTRY;
            close ENTRY;

            #set colour
            $value = '/colour=7'."\n";
            open  ENTRY, '>', \$colour;
            write ENTRY;
            close ENTRY;

            #combine feature
            $feature_rRNA = $boundary_rrna.
                            $label.
                            $rrna_product.
                            $note.
                            $colour;
            $key_rRNA     = $rRNA->{$rRNA_db}->{$qu_left_bd}->{'left_bd'}.'_'.$rRNA->{$rRNA_db}->{$qu_left_bd}->{'right_bd'}.'_rRNA_'.$args{counter}; #create  uniqe ID instead of ORFnumber

            ${$args{genbank_ref}}{$key_rRNA} = $feature_rRNA;
            push (@{$args{feature_list_ref}}, $key_rRNA);

            #add to summeray file
            open WRITE, ">>${$args{ini_ref}}{rrna_results}\/rRNA_summary.txt" or do {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Could not append to rRNA summary file.\nAborting at file $args{filename}.",
                                                             -bitmap  => 'info',
                                                             -buttons => ['OK']
                                                             );
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error parsing rRNA results for file $args{filename}".
                              "\nCould not append to rRNA summary file.\nAborting at file $args{filename}.\n\n";
               close ERRORLOG;

               &hide_pbar_2;
               return (0);
            };

            print WRITE "Input file: $args{filename}\tType: $db\tLeft boundary: $rRNA->{$rRNA_db}->{$qu_left_bd}->{'left_bd'}\tRight boundary: $rRNA->{$rRNA_db}->{$qu_left_bd}->{'right_bd'}\tOrientation: $rRNA->{$rRNA_db}->{$qu_left_bd}->{'orientation'}\t".
                        "Score: $rRNA->{$rRNA_db}->{$qu_left_bd}->{'score'}\n";
            close WRITE;
         }
      }
   }
   &hide_pbar_2;
   return;
}

sub remove_global_overlaps {
   my %args = @_;
   my $iterate;
   foreach my $qu_rRNA_db (keys %{ ${$args{rRNA_ref}} }) {
      next if ($qu_rRNA_db =~ /\W+/);
      foreach my $qu_left_bd (keys %{ ${$args{rRNA_ref}} -> { $qu_rRNA_db } }) {
         next if ($qu_left_bd =~ /\D+/);
         next if (${$args{rRNA_ref}} -> { $qu_rRNA_db } -> { $qu_left_bd } -> { 'score' } !~ /\w+/); #ignore if deleted
         $args{counter}++;
         my $qu_score    = ${$args{rRNA_ref}} -> { $qu_rRNA_db } -> { $qu_left_bd } -> { 'score' };
         my $qu_right_bd = ${$args{rRNA_ref}} -> { $qu_rRNA_db } -> { $qu_left_bd } -> { 'right_bd' };
         foreach my $sb_rRNA_db (keys %{ ${$args{rRNA_ref}} }) {
            next if ($sb_rRNA_db =~ /\W+/);
            foreach my $sb_left_bd (keys %{ ${$args{rRNA_ref}} -> { $sb_rRNA_db } }) {
               next if ($sb_left_bd =~ /\D+/);
               next if ($qu_rRNA_db eq $sb_rRNA_db && $qu_left_bd == $sb_left_bd); #ignore self
               next if (${$args{rRNA_ref}} -> { $sb_rRNA_db } -> { $sb_left_bd } -> { 'score' } !~ /\w+/); #ignore if deleted
               my $sb_score    = ${$args{rRNA_ref}} -> { $sb_rRNA_db } -> { $sb_left_bd } -> { 'score' };
               my $sb_right_bd = ${$args{rRNA_ref}} -> { $sb_rRNA_db } -> { $sb_left_bd } -> { 'right_bd' };

               if (($sb_left_bd  >= $qu_left_bd  && $sb_left_bd  <= $qu_right_bd) ||
                   ($sb_right_bd >= $qu_left_bd  && $sb_right_bd <= $qu_right_bd) ||
                   ($sb_left_bd  <= $qu_left_bd  && $sb_right_bd >= $qu_right_bd)) {
                  $iterate = 1;
                  #first, remove weaker link
                  if ($sb_rRNA_db eq '16S_rRNA' && ${$args{rRNA_ref}} -> { $sb_rRNA_db } -> { $sb_left_bd } -> {'alignment_length'} < 500) {
                     delete ${$args{rRNA_ref}} -> { $sb_rRNA_db } -> { $sb_left_bd };
                     next;
                  }
                  if ($qu_rRNA_db eq '16S_rRNA' && ${$args{rRNA_ref}} -> { $qu_rRNA_db } -> { $qu_left_bd } -> {'alignment_length'} < 500) {
                     delete ${$args{rRNA_ref}} -> { $qu_rRNA_db } -> { $qu_left_bd };
                     next;
                  }

                  if ($qu_score >= $sb_score) {
                     #print "\nDeleted 5S $sb_left_bd ...".${$args{rRNA_ref}}->{$sb_rRNA_db}->{$sb_left_bd}->{'orientation'}." ...\n" if ($sb_rRNA_db eq '5S_rRNA');
                     delete ${$args{rRNA_ref}} -> { $sb_rRNA_db } -> { $sb_left_bd };
                  } else {
                     #print "\nDeleted 5S $qu_left_bd ...".${$args{rRNA_ref}}->{$qu_rRNA_db}->{$qu_left_bd}->{'orientation'}." ...\n" if ($qu_rRNA_db eq '5S_rRNA');
                     delete ${$args{rRNA_ref}} -> { $qu_rRNA_db } -> { $qu_left_bd };
                  }
               }
            }
         }
         if ($iterate == 1) {
            #restart iteration
            &remove_global_overlaps(rRNA_ref => $args{rRNA_ref},
                                    counter  => 1
                                   );
         }
      }
   }
   return($args{counter});
}

1;

