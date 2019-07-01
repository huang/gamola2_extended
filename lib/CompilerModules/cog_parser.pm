#!/opt/ActivePerl-5.8/bin/perl

#Blast result parser
#input arguments: main_window, progress_bar, gb_file, gene_model

package CompilerModules::cog_parser;
use strict;
use initialise::read_me  qw(:DEFAULT);
use Basics::progress_bar qw(:DEFAULT);

use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&cog_parser &translate_COG2003 &translate_POG2013 &translate_COG2014 &translate_arCOG &translate_arCOG2014 &translate_COG2008 );
use vars qw();


#local variables
my (%args, %seen, $value);

format ENTRY =
~~                   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$value
.

sub cog_parser {
   my %args = @_;

   #which COG database to reformat?
   if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
      &cog2003_parser(main_window         => $args{main_window},
                      progress_bar        => $args{progress_bar},
                      auto_ini_ref        => $args{auto_ini_ref},
                      ini_ref             => $args{ini_ref},
                      genbank_ref         => $args{genbank_ref},
                      feature_list_ref    => $args{feature_list_ref},
                      orf_to_id           => $args{orf_to_id},
                      counter             => $args{counter},
                      COGcode_to_number   => $args{COGcode_to_number},
                      COGcode_to_header   => $args{COGcode_to_header},
                      COGcode_to_letter   => $args{COGcode_to_letter},
                      COGcode_to_phylo    => $args{COGcode_to_phylo},
                      COGletter_to_family => $args{COGletter_to_family}
                     );
   } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2008') {
      &cog2008_parser(     main_window         => $args{main_window},
                           progress_bar        => $args{progress_bar},
                           auto_ini_ref        => $args{auto_ini_ref},
                           ini_ref             => $args{ini_ref},
                           genbank_ref         => $args{genbank_ref},
                           feature_list_ref    => $args{feature_list_ref},
                           orf_to_id           => $args{orf_to_id},
                           counter             => $args{counter},
                           COGcode_to_number   => $args{COGcode_to_number},
                           COGcode_to_header   => $args{COGcode_to_header},
                           COGcode_to_letter   => $args{COGcode_to_letter},
                           COGletter_to_family => $args{COGletter_to_family},
                           COG2008genomes      => $args{COG2008genomes},
                           COG2008_to_def      => $args{COG2008_to_def},
                           COG2008_to_acc      => $args{COG2008_to_acc},
                          );
   } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2014') {
      &cog2014_parser(     main_window         => $args{main_window},
                           progress_bar        => $args{progress_bar},
                           auto_ini_ref        => $args{auto_ini_ref},
                           ini_ref             => $args{ini_ref},
                           genbank_ref         => $args{genbank_ref},
                           feature_list_ref    => $args{feature_list_ref},
                           orf_to_id           => $args{orf_to_id},
                           counter             => $args{counter},
                           COGcode_to_number   => $args{COGcode_to_number},
                           COGcode_to_header   => $args{COGcode_to_header},
                           COGcode_to_letter   => $args{COGcode_to_letter},
                           COGletter_to_family => $args{COGletter_to_family},
                           COG2014genomes      => $args{COG2014genomes},
                           COG2014_to_def      => $args{COG2014_to_def},
                           COG2014_to_acc      => $args{COG2014_to_acc},
                           COG2014_to_refseq   => $args{COG2014_to_refseq},
                           COGletter_to_family => $args{COGletter_to_family},
                          );
   } elsif (${$args{auto_ini_ref}}{COG_type} eq 'arCOG') {
      &arcog_parser(       main_window         => $args{main_window},
                           progress_bar        => $args{progress_bar},
                           auto_ini_ref        => $args{auto_ini_ref},
                           ini_ref             => $args{ini_ref},
                           genbank_ref         => $args{genbank_ref},
                           feature_list_ref    => $args{feature_list_ref},
                           orf_to_id           => $args{orf_to_id},
                           counter             => $args{counter},
                           COGcode_to_number   => $args{COGcode_to_number},
                           COGcode_to_header   => $args{COGcode_to_header},
                           COGcode_to_letter   => $args{COGcode_to_letter},
                           COGletter_to_family => $args{COGletter_to_family},
                           arCOGcode_to_COG    => $args{arCOGcode_to_COG},
                           arCOGacc_to_org     => $args{arCOGacc_to_org},
                          );
   } elsif (${$args{auto_ini_ref}}{COG_type} eq 'arCOG2014') {
      &arcog2014_parser(  main_window          => $args{main_window},
                          progress_bar         => $args{progress_bar},
                          auto_ini_ref         => $args{auto_ini_ref},
                          ini_ref              => $args{ini_ref},
                          genbank_ref          => $args{genbank_ref},
                          feature_list_ref     => $args{feature_list_ref},
                          orf_to_id            => $args{orf_to_id},
                          counter              => $args{counter},
                          arCOGcode2014_to_COG => $args{arCOGcode2014_to_COG},
                          arCOGacc2014_to_org  => $args{arCOGacc2014_to_org},
                          COGletter_to_family  => $args{COGletter_to_family},
                        );
   } elsif (${$args{auto_ini_ref}}{COG_type} eq 'POG2013') {
      &pog2013_parser(    main_window          => $args{main_window},
                          progress_bar         => $args{progress_bar},
                          auto_ini_ref         => $args{auto_ini_ref},
                          ini_ref              => $args{ini_ref},
                          genbank_ref          => $args{genbank_ref},
                          feature_list_ref     => $args{feature_list_ref},
                          orf_to_id            => $args{orf_to_id},
                          counter              => $args{counter},
                          POG2013_to_gi        => $args{POG2013_to_gi},
                        );
   } else {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Could not determine which COG database was selected\.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error parsing global COG result hash".
                     "\nCould not determine which COG database was selected\.\n\n";
      close ERRORLOG;
      return (0);
   }
   return;
}

sub cog2008_parser {
   my %args = @_;

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling COG2008 results",
                   label        => ''
                  );
   &show_pbar_2;

   #parse through COG2008s
   foreach (my $i = 1; $i < $args{counter}; $i++) {
      my ($ORF_key, $filename, $id, $left_bd, $right_bd, $orientation, $best_hit, $evalue, $length, $score,
          $COG2008_code, $match, $COG2008_left_bd, $COG2008_right_bd, $boundary_COG2008, $colour, $label, $product,
          $feature_COG2008, $key_COG2008, $first_COG2008_code, $first_evalue, $first_score, $first_length, $COG2008_header,
          $delete_COG2008, $COG2008_acc, $first_COG2008_acc);

      #get filename, left_bd, right_bd, orientation
      'reset' =~ m/reset/;
      $args{orf_to_id}->{$i} =~ m/^(.*?)___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
      $filename    = $1;
      $id          = $2;
      $left_bd     = $3;
      $right_bd    = $4;
      $orientation = $5;
      unless (defined $orientation) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse from entry ${$args{orf_to_id}}->{$i}\.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing global COG2008 result hash".
                        "\nCould not parse from entry ${$args{orf_to_id}}->{$i}\.\n\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }

      #update status bar
      if (($i % 20) == 0) {
         &update_pbar_2(title        => "Compiling COG2008 results",
                        label        => "Compiling COG2008 results for $filename",
                        progress     => ($i / $args{counter}) * 100,
                       );
      }

      #read input file
      (my $file_ref) = slurp(main_window => $args{main_window},
                             auto_ini_ref => $args{auto_ini_ref},
                             directory   => ${$args{ini_ref}}{COG_results},
                             filename    => $filename.'_COG_'.$id
                            );
      if ($file_ref eq '0') {
         &hide_pbar_2;
         return (0);
      }

      #skip if no COG results
      next if (${$file_ref} =~ /No COG Hits found/s);

      #get COG header
      ${$file_ref} =~ m/Blast overview(.+)Complete list of Blast results/si;
      $COG2008_header = $1."\n\>No further hits\n";

      #catch best defined COG2008 hit
      $COG2008_code = '';
      while ($COG2008_code !~ m/\w+/) {
         'reset' =~ m/reset/;
         $COG2008_header =~ m/\>(.*?)\n\r?/si;
         $best_hit = $1;

         unless (defined $best_hit) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse best hit from COG2008 entry $filename\_COG_$id.\n$args{orf_to_id}->{$i}",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing COG2008 results for file $filename".
                           "\nCould not parse best hit from COG2008 entry $filename\_COG_$id.\n$args{orf_to_id}->{$i}\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }

         #no more hits?
         if ($best_hit eq 'No further hits') {
            #reset COG_code and evalue to first hit
            $COG2008_code = 'unclassified COG2008';
            $COG2008_acc  = $first_COG2008_acc;
            $evalue       = $first_evalue;
            $score        = $first_score;
            $length       = $first_length;
            last;
         }

         #catch evalue, score and length
         'reset' =~ m/reset/;
         $best_hit =~ m/Length\=(\d+)\s+Score\=(\d+)\s+Expect\=(.+)$/;
         $length = $1;
         $score  = $2;
         $evalue = $3;
         unless (defined $evalue) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse the evalue from entry $best_hit in file $filename, ID $id\.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing COG2008 results for file $filename".
                           "\nCould not parse the evalue from entry $best_hit in file $filename, ID $id\.\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }
         if ($evalue =~ /^e/) {$evalue = '1'.$evalue};

         #catch COG2008 code and COG2008 acc
         'reset' =~ m/reset/;
         $best_hit =~ m/^(\S+)\s+COG2008: (.+);\s+Class:/;
         ($COG2008_acc, $COG2008_code) = ($1, $2);
         unless (defined $COG2008_code) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse the COG2008 code from entry $best_hit.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing COG2008 results for file $filename".
                           "\nCould not parse the COG2008 code from entry $best_hit.\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }

         #if undefined COG, change COG_acc value to 'NULL' to bring in line with result file
         if ($COG2008_acc eq 'Not def') {
            $COG2008_acc = 'NULL';
         }

         #define best COG2008code and evalue in case no classified hit can be found
         unless (defined $first_COG2008_code) {
            $first_COG2008_code = $COG2008_code;
            $first_COG2008_acc  = $COG2008_acc;
            $first_evalue       = $evalue;
            $first_score        = $score;
            $first_length       = $length;
         }

         #skip if above threshold, define COG2008code to escape while loop, entry will be skipped anyway
         if ($evalue > ${$args{ini_ref}}{COG_cut}) {
            #reset COG2008_code and evalue to first hit
            $COG2008_code = $first_COG2008_code;
            $COG2008_acc  = $first_COG2008_acc;
            $evalue     = $first_evalue;
            $score      = $first_score;
            $length     = $first_length;
            unless (defined $COG2008_code) {
               $COG2008_code = 'unclassified COG2008';
               $delete_COG2008 = 1;
            }
         }
      }
      undef $COG2008_header;

      #skip if above threshold
      next if ($evalue > ${$args{ini_ref}}{COG_cut});

      #catch COG2008 fragment boundaries
      ${$file_ref} =~ s/^.+?Complete list of Blast results//s;
      'reset' =~ /reset/;
      ${$file_ref} =~ m/\>\w+\|\Q$COG2008_acc\E.*?(Query.*?)[\n\r]+\>/s;
      $match = $1;

      #catch if only one COG2008 hit present
      unless (defined $match) {
         'reset' =~ /reset/;
         ${$file_ref} =~ m/\>\w+\|\Q$COG2008_acc\E.*?(Query.*?)\s+Database:/s;
         $match = $1;
      }
      unless (defined $match) {
         'reset' =~ /reset/;
         ${$file_ref} =~ m/\>\w+\|\Q$COG2008_acc\E.*?(Query.*?)\s+Lambda/s;
         $match = $1;
      }
      unless (defined $match) {
         'reset' =~ /reset/;
         ${$file_ref} =~ m/\>\Q$COG2008_acc\E.*?(Query.*?)\s+Lambda/s;
         $match = $1;
      }
      unless (defined $match) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse COG2008 alignment from entry $filename.\n\n$args{orf_to_id}->{$i}",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing COG2008 results Code:\'$COG2008_code\' Acc:\'$COG2008_acc\' for file $filename".
                        "\nCould not parse COG2008 alignment from entry $filename.\n$args{orf_to_id}->{$i}\n${$file_ref}\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }

      #catch COG2008 left_bd
      'reset' =~ /reset/;
      $match =~ m/Query\:\s+(\d+)\s+/s;
      $COG2008_left_bd = $1;
      unless (defined $COG2008_left_bd) {
         $match =~ m/Query\s+(\d+)\s+/s;
         $COG2008_left_bd = $1;
      }

      #catch COG2008 right_bd
      'reset' =~ /reset/;
      $match =~ m/.*Query\:\s+\d+\s+\D+\s+(\d+)/s;
      $COG2008_right_bd = $1;
      unless (defined $COG2008_right_bd) {
         $match =~ m/.*Query\s+\d+\s+\D+\s+(\d+)/s;
         $COG2008_right_bd = $1;
      }

      unless (defined $COG2008_left_bd && defined $COG2008_right_bd) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse COG2008 alignment boundaries from entry $filename",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing COG2008 results for file $filename".
                        "\nCould not parse COG2008 alignment boundaries from entry $match in entry ${$file_ref}.\n\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }
      #create boundary entry
      if ($orientation eq 'sense') {
         $COG2008_right_bd = $left_bd + ($COG2008_right_bd * 3) - 1;
         $COG2008_left_bd  = $left_bd + ($COG2008_left_bd * 3) - 3;

         #modify boundaries is necessary
         if ($COG2008_left_bd == $left_bd) {
            $COG2008_left_bd += 3;  #add one codon to avoid equal to start position
         }
         if ($COG2008_right_bd == $right_bd) {
            $COG2008_right_bd -= 3;  #remove one codon to avoid equal to stop position
         }

         $boundary_COG2008  = '     COG_match       '.$COG2008_left_bd.'..'.$COG2008_right_bd."\n";
      } elsif ($orientation eq 'antisense') {
         my $length = $COG2008_right_bd - $COG2008_left_bd;
         $COG2008_right_bd = $right_bd - ($COG2008_left_bd * 3) + 3;
         $COG2008_left_bd  = $COG2008_right_bd - ($length * 3) - 2;

         #modify boundaries is necessary
         if ($COG2008_left_bd == $left_bd) {
            $COG2008_left_bd += 3;  #add one codon to avoid equal to stop position
         }
         if ($COG2008_right_bd == $right_bd) {
            $COG2008_right_bd -= 3;  #remove one codon to avoid equal to start position
         }

         $boundary_COG2008  = '     COG_match       complement('.$COG2008_left_bd.'..'.$COG2008_right_bd."\)\n";
      }

      #create product entry
      if ($COG2008_code eq 'unclassified COG2008') {
         $value = '/product="['.$COG2008_code.'] '.
                  ' Length='.$length.' Score='.$score.' Expect='.$evalue.'"'."\n";
      } else {
         #clean up COG annotation
         ${$args{COG2008_to_def}}->{$COG2008_code}->{'annotation'} =~ s/\"/\'/gs;

         $value = '/product="['.${$args{COG2008_to_def}}->{$COG2008_code}->{'class'}.'] '.
                  'COG2008: '.$COG2008_code.' '.
                  ${$args{COG2008_to_def}}->{$COG2008_code}->{'annotation'}.
                  ' Length='.$length.' Score='.$score.' Expect='.$evalue.'"'."\n";
      }
      #no COG2008 id?
      if ($value =~ /\[\]/) {$value =~ s/\[\]/\[not assigned\]/};
      open  ENTRY, '>', \$product;
      write ENTRY;
      close ENTRY;

      #create COG2008 label
      $value = 'label= ['.${$args{COG2008_to_def}}->{$COG2008_code}->{'class'}.'] '.
               ${$args{COGletter_to_family}}->{${$args{COG2008_to_def}}->{$COG2008_code}->{'class'}}.
               "\n";
      #clean up COG2008 label
      $value =~ s/\//-/gs;
      $value = '/'.$value;

      #no COG2008 id?
      if ($value =~ /\[\]/) {$value =~ s/\[\]/\[unclassified COG2008\]/};
      open  ENTRY, '>', \$label;
      write ENTRY;
      close ENTRY;

      #create codon start and translation table tags
      $colour = '                     /colour=2'."\n";

      #combine feature
      $feature_COG2008 = $boundary_COG2008.
                         $label.
                         $colour.
                         $product;
      $key_COG2008     = $COG2008_left_bd.'_'.$COG2008_right_bd.'_COG_match_'.$id; #maintain original uniqe ID instead of ORFnumber

      ${$args{genbank_ref}}{$key_COG2008} = $feature_COG2008;
      push (@{$args{feature_list_ref}}, $key_COG2008);

      #delete unassigned COG code again, otherwise it's stored for next ORF
      if ($delete_COG2008) {
         $COG2008_code = '';
         undef $delete_COG2008;
      }
   }
   &hide_pbar_2;
   return;
}

sub cog2014_parser {
   my %args = @_;

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling COG2014 results",
                   label        => ''
                  );
   &show_pbar_2;

   #parse through COG2014s
   foreach (my $i = 1; $i < $args{counter}; $i++) {
      my ($ORF_key, $filename, $id, $left_bd, $right_bd, $orientation, $best_hit, $evalue, $length, $score,
          $COG2014_code, $match, $COG2014_left_bd, $COG2014_right_bd, $boundary_COG2014, $colour, $label, $product,
          $feature_COG2014, $key_COG2014, $COG2014_cl_verbose, $first_COG2014_cl_verbose, $COG2014_class, $first_COG2014_class, $first_COG2014_code, $COG2014_anno, $first_COG2014_anno, $first_evalue, $first_score, $first_length, $COG2014_header,
          $delete_COG2014, $COG2014_acc, $first_COG2014_acc);

      #get filename, left_bd, right_bd, orientation
      'reset' =~ m/reset/;
      $args{orf_to_id}->{$i} =~ m/^(.*?)___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
      $filename    = $1;
      $id          = $2;
      $left_bd     = $3;
      $right_bd    = $4;
      $orientation = $5;
      unless (defined $orientation) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse from entry ${$args{orf_to_id}}->{$i}\.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing global COG2014 result hash".
                        "\nCould not parse from entry ${$args{orf_to_id}}->{$i}\.\n\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }

      #update status bar
      if (($i % 20) == 0) {
         &update_pbar_2(title        => "Compiling COG2014 results",
                        label        => "Compiling COG2014 results for $filename",
                        progress     => ($i / $args{counter}) * 100,
                       );
      }

      #read input file
      (my $file_ref) = slurp(main_window => $args{main_window},
                             auto_ini_ref => $args{auto_ini_ref},
                             directory   => ${$args{ini_ref}}{COG_results},
                             filename    => $filename.'_COG_'.$id
                            );
      if ($file_ref eq '0') {
         &hide_pbar_2;
         return (0);
      }

      #skip if no COG results
      next if (${$file_ref} =~ /No COG Hits found/s);

      #get COG header
      ${$file_ref} =~ m/Blast overview(.+)Complete list of Blast results/si;
      $COG2014_header = $1."\n\>No further hits\n";

      #catch best defined COG2014 hit
      $COG2014_code = '';
      while ($COG2014_code !~ m/\w+/) {
         'reset' =~ m/reset/;
         $COG2014_header =~ m/\>(.*?)\n\r?/si;
         $best_hit = $1;

         unless (defined $best_hit) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse best hit from COG2014 entry $filename\_COG_$id.\n$args{orf_to_id}->{$i}",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing COG2014 results for file $filename".
                           "\nCould not parse best hit from COG2014 entry $filename\_COG_$id.\n$args{orf_to_id}->{$i}\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }

         #no more hits?
         if ($best_hit eq 'No further hits') {
            #reset COG_code and evalue to first hit
            $COG2014_code       = 'unclassified COG2014';
            $COG2014_acc        = $first_COG2014_acc;
            $evalue             = $first_evalue;
            $score              = $first_score;
            $length             = $first_length;
            $COG2014_class      = $first_COG2014_class;
            $COG2014_anno       = $first_COG2014_anno;
            $COG2014_cl_verbose = $first_COG2014_cl_verbose;
            last;
         }

         #catch evalue, score and length
         'reset' =~ m/reset/;
         $best_hit =~ m/Length\=(\d+)\s+Score\=(\d+)\s+Expect\=(.+)$/;
         ($length, $score, $evalue) = ($1, $2, $3);

         unless (defined $evalue) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse the evalue from COG2014 entry $best_hit in file $filename, ID $id\.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing COG2014 results for file $filename".
                           "\nCould not parse the evalue from entry $best_hit in file $filename, ID $id\.\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }
         if ($evalue =~ /^e/) {$evalue = '1'.$evalue};

         #catch COG2014 code and COG2014 acc
         'reset' =~ m/reset/;
         $best_hit =~ m/^(\S+)\s+COG2014: (COG\d+);\s+Class:\s+(\w+)\,\s+(.+)\;.+?Annotation:\s+(.+) Length\=\d+/;
         ($COG2014_acc, $COG2014_code, $COG2014_class, $COG2014_cl_verbose, $COG2014_anno) = ($1, $2, $3, $4, $5);
         unless (defined $COG2014_code) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse the COG2014 code from entry $best_hit.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing COG2014 results for file $filename".
                           "\nCould not parse the COG2014 code from entry $best_hit.\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }

         #if undefined COG, change COG_acc value to 'NULL' to bring in line with result file
         if ($COG2014_acc eq 'Not def') {
            $COG2014_acc = 'NULL';
         }

         #clean up annotation
         $COG2014_anno =~ s/~+//g;

         #define best COG2014code and evalue in case no classified hit can be found
         unless (defined $first_COG2014_code) {
            $first_COG2014_code       = $COG2014_code;
            $first_COG2014_acc        = $COG2014_acc;
            $first_evalue             = $evalue;
            $first_score              = $score;
            $first_length             = $length;
            $first_COG2014_class      = $COG2014_class;
            $first_COG2014_anno       = $COG2014_anno;
            $first_COG2014_cl_verbose = $COG2014_cl_verbose;
         }

         #skip if above threshold, define COG2014code to escape while loop, entry will be skipped anyway
         if ($evalue > ${$args{ini_ref}}{COG_cut}) {
            #reset COG2014_code and evalue to first hit
            $COG2014_code       = $first_COG2014_code;
            $COG2014_acc        = $first_COG2014_acc;
            $evalue             = $first_evalue;
            $score              = $first_score;
            $length             = $first_length;
            $COG2014_class      = $first_COG2014_class;
            $COG2014_anno       = $first_COG2014_anno;
            $COG2014_cl_verbose = $first_COG2014_cl_verbose;
            unless (defined $COG2014_code) {
               $COG2014_code = 'unclassified COG2014';
               $delete_COG2014 = 1;
            }
         }
      }
      undef $COG2014_header;

      #skip if above threshold
      next if ($evalue > ${$args{ini_ref}}{COG_cut});

      #clean up COG annotation
      $COG2014_anno =~ s/\"/\'/gs;

      #catch COG2014 fragment boundaries
      ${$file_ref} =~ s/^.+?Complete list of Blast results//s;
      'reset' =~ /reset/;
      ${$file_ref} =~ m/\>\Q$COG2014_acc\E.*?(Query.*?)[\n\r]+\>/s;
      $match = $1;

      #catch if only one COG2014 hit present
      unless (defined $match) {
         'reset' =~ /reset/;
         ${$file_ref} =~ m/\>\w+\|\Q$COG2014_acc\E.*?(Query.*?)\s+Database:/s;
         $match = $1;
      }
      unless (defined $match) {
         'reset' =~ /reset/;
         ${$file_ref} =~ m/\>\w+\|\Q$COG2014_acc\E.*?(Query.*?)\s+Lambda/s;
         $match = $1;
      }
      unless (defined $match) {
         'reset' =~ /reset/;
         ${$file_ref} =~ m/\>\Q$COG2014_acc\E.*?(Query.*?)\s+Lambda/s;
         $match = $1;
      }
      unless (defined $match) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse COG2014 alignment from entry $filename.\n\n$args{orf_to_id}->{$i}",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing COG2014 results Code:\'$COG2014_code\' Acc:\'$COG2014_acc\' for file $filename".
                        "\nCould not parse COG2014 alignment from entry $filename.\n$args{orf_to_id}->{$i}\n${$file_ref}\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }

      #catch COG2014 left_bd
      'reset' =~ /reset/;
      $match =~ m/Query\:\s+(\d+)\s+/s;
      $COG2014_left_bd = $1;
      unless (defined $COG2014_left_bd) {
         $match =~ m/Query\s+(\d+)\s+/s;
         $COG2014_left_bd = $1;
      }

      #catch COG2014 right_bd
      'reset' =~ /reset/;
      $match =~ m/.*Query\:\s+\d+\s+\D+\s+(\d+)/s;
      $COG2014_right_bd = $1;
      unless (defined $COG2014_right_bd) {
         $match =~ m/.*Query\s+\d+\s+\D+\s+(\d+)/s;
         $COG2014_right_bd = $1;
      }

      unless (defined $COG2014_left_bd && defined $COG2014_right_bd) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse COG2014 alignment boundaries from entry $filename",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing COG2014 results for file $filename".
                        "\nCould not parse COG2014 alignment boundaries from entry $match in entry ${$file_ref}.\n\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }
      #create boundary entry
      if ($orientation eq 'sense') {
         $COG2014_right_bd = $left_bd + ($COG2014_right_bd * 3) - 1;
         $COG2014_left_bd  = $left_bd + ($COG2014_left_bd  * 3) - 3;

         #modify boundaries is necessary
         if ($COG2014_left_bd == $left_bd) {
            $COG2014_left_bd += 3;  #add one codon to avoid equal to start position
         }
         if ($COG2014_right_bd == $right_bd) {
            $COG2014_right_bd -= 3;  #remove one codon to avoid equal to stop position
         }

         $boundary_COG2014  = '     COG_match       '.$COG2014_left_bd.'..'.$COG2014_right_bd."\n";
      } elsif ($orientation eq 'antisense') {
         my $length = $COG2014_right_bd - $COG2014_left_bd;
         $COG2014_right_bd = $right_bd - ($COG2014_left_bd * 3) + 3;
         $COG2014_left_bd  = $COG2014_right_bd - ($length * 3) - 2;

         #modify boundaries is necessary
         if ($COG2014_left_bd == $left_bd) {
            $COG2014_left_bd += 3;  #add one codon to avoid equal to stop position
         }
         if ($COG2014_right_bd == $right_bd) {
            $COG2014_right_bd -= 3;  #remove one codon to avoid equal to start position
         }

         $boundary_COG2014  = '     COG_match       complement('.$COG2014_left_bd.'..'.$COG2014_right_bd."\)\n";
      }

      #create product entry
      if ($COG2014_code eq 'unclassified COG2014') {
         $value = '/product="['.$COG2014_code.'] '.
                  ' Length='.$length.' Score='.$score.' Expect='.$evalue.'"'."\n";
      } else {
         my $COG_protID = $_; #ACC is different here
         $COG_protID    =~ s/^\S+?\|//;
         $COG_protID    =~ s/\.\d+[\|\s]//g;
         $value = "/product=\"\[$COG2014_class\]".
                  'COG2014: '.$COG2014_code.' '. $COG2014_anno.
                  ' Length='.$length.' Score='.$score.' Expect='.$evalue.'"'."\n";
      }
      #no COG2014 id?
      if ($value =~ /\[\]/) {$value =~ s/\[\]/\[not assigned\]/};
      open  ENTRY, '>', \$product;
      write ENTRY;
      close ENTRY;

      #create COG2014 label
      $value = 'label= ['.$COG2014_class.'] '.$COG2014_cl_verbose."\n";
      #clean up COG2014 label
      $value =~ s/\//-/gs;
      $value = '/'.$value;

      #no COG2014 id?
      if ($value =~ /\[\]/) {$value =~ s/\[\]/\[unclassified COG2014\]/};
      open  ENTRY, '>', \$label;
      write ENTRY;
      close ENTRY;

      #create codon start and translation table tags
      $colour = '                     /colour=2'."\n";

      #combine feature
      $feature_COG2014 = $boundary_COG2014.
                         $label.
                         $colour.
                         $product;
      $key_COG2014     = $COG2014_left_bd.'_'.$COG2014_right_bd.'_COG_match_'.$id; #maintain original uniqe ID instead of ORFnumber

      ${$args{genbank_ref}}{$key_COG2014} = $feature_COG2014;
      push (@{$args{feature_list_ref}}, $key_COG2014);

      #delete unassigned COG code again, otherwise it's stored for next ORF
      if ($delete_COG2014) {
         $COG2014_code = '';
         undef $delete_COG2014;
      }
   }
   &hide_pbar_2;
   return;
}

sub arcog_parser {
   my %args = @_;
   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling archaeal COG results",
                   label        => ''
                  );
   &show_pbar_2;

   #parse through arCOGs
   foreach (my $i = 1; $i < $args{counter}; $i++) {
      my ($ORF_key, $filename, $id, $left_bd, $right_bd, $orientation, $best_hit, $evalue, $length, $score,
          $arCOG_code, $match, $arCOG_left_bd, $arCOG_right_bd, $boundary_arCOG, $colour, $label, $product,
          $feature_arCOG, $key_arCOG, $first_arCOG_code, $first_evalue, $first_score, $first_length, $arCOG_header,
          $delete_arCOG, $arCOG_acc, $first_arCOG_acc, $arCOG_class, $first_arCOG_class, $arCOG_annotation,
          $first_arCOG_annotation);

      #get filename, left_bd, right_bd, orientation
      'reset' =~ m/reset/;
      $args{orf_to_id}->{$i} =~ m/^(.*?)___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
      $filename    = $1;
      $id          = $2;
      $left_bd     = $3;
      $right_bd    = $4;
      $orientation = $5;
      unless (defined $orientation) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse from entry ${$args{orf_to_id}}->{$i}\.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing global arCOG result hash".
                        "\nCould not parse from entry ${$args{orf_to_id}}->{$i}\.\n\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }

      #update status bar
      if (($i % 20) == 0) {
         &update_pbar_2(title        => "Compiling archaeal COG results",
                        label        => "Compiling arCOG results for $filename",
                        progress     => ($i / $args{counter}) * 100,
                       );
      }

      #read input file
      (my $file_ref) = slurp(main_window => $args{main_window},
                             auto_ini_ref => $args{auto_ini_ref},
                             directory   => ${$args{ini_ref}}{COG_results},
                             filename    => $filename.'_COG_'.$id
                            );
      if ($file_ref eq '0') {
         &hide_pbar_2;
         return (0);
      }

      #skip if no COG results
      next if (${$file_ref} =~ /No COG Hits found/s);

      #get COG header
      ${$file_ref} =~ m/Blast overview(.+)Complete list of Blast results/si;
      $arCOG_header = $1."\n\>No further hits\n";

      #catch best defined arCOG hit
      $arCOG_code = '';
      while ($arCOG_code !~ m/\w+/) {
         'reset'       =~ m/reset/;
         $arCOG_header =~ m/\>(.*?)\n\r?/si;
         $best_hit     = $1;

         unless (defined $best_hit) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse best hit from arCOG entry $filename\_COG_$id.\n$args{orf_to_id}->{$i}",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing arCOG results for file $filename".
                           "\nCould not parse best hit from arCOG entry $filename\_COG_$id.\n$args{orf_to_id}->{$i}\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }

         #no more hits?
         if ($best_hit eq 'No further hits') {
            #reset COG_code and evalue to first hit
            $arCOG_code       = 'unclassified archaeal COG';
            $arCOG_acc        = $first_arCOG_acc;
            $arCOG_class      = $first_arCOG_class;
            $arCOG_annotation = $first_arCOG_annotation;
            $evalue           = $first_evalue;
            $score            = $first_score;
            $length           = $first_length;
            last;
         }

         #catch evalue, score and length
         'reset'   =~ m/reset/;
         $best_hit =~ m/Length\=(\d+)\s+Score\=(\d+)\s+Expect\=(.+)$/;
         $length   = $1;
         $score    = $2;
         $evalue   = $3;
         unless (defined $evalue) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse the evalue from entry $best_hit in file $filename, ID $id\.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing arCOG results for file $filename".
                           "\nCould not parse the evalue from entry $best_hit in file $filename, ID $id\.\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }
         if ($evalue =~ /^e/) {$evalue = '1'.$evalue};

         #catch arCOG code and arCOG acc
         'reset' =~ m/reset/;
         if ($best_hit =~ m/^\w+\|([^\|]+)\|?/) {
            $best_hit =~ m/^\w+\|([^\|]+)\|?\s+arCOG: (.+) Class: (\w);.+?Annotation:\s([^\[]+)/;
            ($arCOG_acc, $arCOG_code, $arCOG_class, $arCOG_annotation) = ($1, $2, $3, $4);
         } else {
            $best_hit =~ m/^(\S+)\s+arCOG\:\s+([^\;]+)\;\s+Class\:\s+([^\;]+)\;.+?Annotation\:\s+([^\~]+)/;
            ($arCOG_acc, $arCOG_code, $arCOG_class, $arCOG_annotation) = ($1, $2, $3, $4);
         }
         unless (defined $arCOG_annotation) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse the arCOG code from entry $best_hit.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing arCOG results for file $filename".
                           "\nCould not parse the COG code from entry $best_hit.\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }
         $arCOG_code       =~ s/;//;
         $arCOG_annotation =~ s/\~+.+$//;

         #define best arCOGcode and evalue in case no classified hit can be found
         unless (defined $first_arCOG_code) {
            $first_arCOG_code       = $arCOG_code;
            $first_arCOG_acc        = $arCOG_acc;
            $first_arCOG_class      = $arCOG_class;
            $first_arCOG_annotation = $arCOG_annotation;
            $first_evalue           = $evalue;
            $first_score            = $score;
            $first_length           = $length;
         }

         #skip if above threshold, define arCOGcode to escape while loop, entry will be skipped anyway
         if ($evalue > ${$args{ini_ref}}{COG_cut}) {
            #reset arCOG_code and evalue to first hit
            $arCOG_code       = $first_arCOG_code;
            $arCOG_acc        = $first_arCOG_acc;
            $arCOG_class      = $first_arCOG_class;
            $arCOG_annotation = $first_arCOG_annotation;
            $evalue           = $first_evalue;
            $score            = $first_score;
            $length           = $first_length;
            unless (defined $arCOG_code) {
               $arCOG_code = 'unclassified archaeal COG';
               $delete_arCOG = 1;
            }
         }
      }
      undef $arCOG_header;

      #skip if above threshold
      next if ($evalue > ${$args{ini_ref}}{COG_cut});

      #cleanup annotation
      $arCOG_annotation =~ s/\"/\'/gs;

      #catch arCOG fragment boundaries
      ${$file_ref} =~ s/^.+?Complete list of Blast results//s;
      'reset' =~ /reset/;
      if (${$file_ref} =~ m/\>\w+\|/) {
         ${$file_ref} =~ m/\>\w+\|\Q$arCOG_acc\E.*?(Query.*?)[\n\r]+\>/s;
         $match = $1;
      } else {
         ${$file_ref} =~ m/\>\Q$arCOG_acc\E.*?(Query.*?)[\n\r]+\>/s;
         $match = $1;
      }

      #catch if only one arCOG hit present
      unless (defined $match) {
         'reset' =~ /reset/;
         if (${$file_ref} =~ m/\>\w+\|/) {
            ${$file_ref} =~ m/\>\w+\|\Q$arCOG_acc\E.*?(Query.*?)\s+Database:/s;
            $match = $1;
         } else {
            ${$file_ref} =~ m/\>\Q$arCOG_acc\E.*?(Query.*?)\s+Database:/s;
            $match = $1;
         }
      }
      unless (defined $match) {
         'reset' =~ /reset/;
         ${$file_ref} =~ m/\>\Q$arCOG_acc\E.*?(Query.*?)\s+Database:/s;
         $match = $1;
      }
      unless (defined $match) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse arCOG alignment from entry $filename.\n\n$args{orf_to_id}->{$i}",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing arCOG results $arCOG_code for file $filename".
                        "\nCould not parse arCOG alignment from entry $filename.\n$args{orf_to_id}->{$i}\n${$file_ref}\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }

      #catch arCOG left_bd
      'reset' =~ /reset/;
      $match =~ m/Query\:\s+(\d+)\s+/s;
      $arCOG_left_bd = $1;
      unless (defined $arCOG_left_bd) {
         $match =~ m/Query\s+(\d+)\s+/s;
         $arCOG_left_bd = $1;
      }

      #catch arCOG right_bd
      'reset' =~ /reset/;
      $match =~ m/.*Query\:\s+\d+\s+\D+\s+(\d+)/s;
      $arCOG_right_bd = $1;
      unless (defined $arCOG_right_bd) {
         $match =~ m/.*Query\s+\d+\s+\D+\s+(\d+)/s;
         $arCOG_right_bd = $1;
      }

      unless (defined $arCOG_left_bd && defined $arCOG_right_bd) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse arCOG alignment boundaries from entry $filename",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing arCOG results for file $filename".
                        "\nCould not parse arCOG alignment boundaries from entry $match in entry ${$file_ref}.\n\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }
      #create boundary entry
      if ($orientation eq 'sense') {
         $arCOG_right_bd = $left_bd + ($arCOG_right_bd * 3) - 1;
         $arCOG_left_bd  = $left_bd + ($arCOG_left_bd * 3) - 3;

         #modify boundaries is necessary
         if ($arCOG_left_bd == $left_bd) {
            $arCOG_left_bd += 3;  #add one codon to avoid equal to start position
         }
         if ($arCOG_right_bd == $right_bd) {
            $arCOG_right_bd -= 3;  #remove one codon to avoid equal to stop position
         }

         $boundary_arCOG  = '     COG_match       '.$arCOG_left_bd.'..'.$arCOG_right_bd."\n";
      } elsif ($orientation eq 'antisense') {
         my $length = $arCOG_right_bd - $arCOG_left_bd;
         $arCOG_right_bd = $right_bd - ($arCOG_left_bd * 3) + 3;
         $arCOG_left_bd  = $arCOG_right_bd - ($length * 3) - 2;

         #modify boundaries is necessary
         if ($arCOG_left_bd == $left_bd) {
            $arCOG_left_bd += 3;  #add one codon to avoid equal to stop position
         }
         if ($arCOG_right_bd == $right_bd) {
            $arCOG_right_bd -= 3;  #remove one codon to avoid equal to start position
         }

         $boundary_arCOG  = '     COG_match       complement('.$arCOG_left_bd.'..'.$arCOG_right_bd."\)\n";
      }

      #create product entry
      if ($arCOG_code eq 'unclassified archaeal COG') {
         $value = '/product="['.$arCOG_class.'] '.
                  ' Length='.$length.' Score='.$score.' Expect='.$evalue.'"'."\n";
      } else {
         $value = '/product="['.$arCOG_class.'] '.
                  'arCOG: '.$arCOG_code.', '.
                  'COG: '.${$args{arCOGcode_to_COG}}->{$arCOG_code}->{'COGcluster'}.' '.
                  $arCOG_annotation.
                  'Length='.$length.' Score='.$score.' Expect='.$evalue.'"'."\n";
      }
      #no arCOG id?
      if ($value =~ /\[\]/) {$value =~ s/\[\]/\[not assigned\]/};
      open  ENTRY, '>', \$product;
      write ENTRY;
      close ENTRY;

      #create arCOG label
      $value = 'label= ['.$arCOG_class.'] '.
               ${$args{COGletter_to_family}}->{${$args{arCOGcode_to_COG}}->{$arCOG_code}->{'COGclass'}}.
               "\n";

      #clean up arCOG label
      $value =~ s/\//-/gs;
      $value = '/'.$value;

      #no arCOG id?
      if ($value =~ /\[\]/) {$value =~ s/\[\]/\[unclassified arCOG\]/};
      open  ENTRY, '>', \$label;
      write ENTRY;
      close ENTRY;

      #create codon start and translation table tags
      $colour = '                     /colour=2'."\n";

      #combine feature
      $feature_arCOG = $boundary_arCOG.
                       $label.
                       $colour.
                       $product;
      $key_arCOG     = $arCOG_left_bd.'_'.$arCOG_right_bd.'_COG_match_'.$id; #maintain original uniqe ID instead of ORFnumber

      ${$args{genbank_ref}}{$key_arCOG} = $feature_arCOG;
      push (@{$args{feature_list_ref}}, $key_arCOG);

      #delete unassigned COG code again, otherwise it's stored for next ORF
      if ($delete_arCOG) {
         $arCOG_code = '';
         undef $delete_arCOG;
      }
   }
   &hide_pbar_2;
   return;
}

sub arcog2014_parser {
   my %args = @_;

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling 2014 archaeal COG results",
                   label        => ''
                  );
   &show_pbar_2;

   #parse through arCOG2014s
   foreach (my $i = 1; $i < $args{counter}; $i++) {
      my ($ORF_key, $filename, $id, $left_bd, $right_bd, $orientation, $best_hit, $evalue, $length, $score,
          $arCOG2014_code, $match, $arCOG2014_left_bd, $arCOG2014_right_bd, $boundary_arCOG2014, $colour, $label, $product,
          $feature_arCOG2014, $key_arCOG2014, $arCOG2014_cl_verbose, $first_arCOG2014_cl_verbose, $first_arCOG2014_code, $first_evalue, $first_score, $first_length, $arCOG2014_header,
          $delete_arCOG2014, $arCOG2014_acc, $first_arCOG2014_acc, $arCOG2014_class, $first_arCOG2014_class, $arCOG2014_annotation,
          $first_arCOG2014_annotation);

      #get filename, left_bd, right_bd, orientation
      'reset' =~ m/reset/;
      $args{orf_to_id}->{$i} =~ m/^(.*?)___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
      $filename    = $1;
      $id          = $2;
      $left_bd     = $3;
      $right_bd    = $4;
      $orientation = $5;
      unless (defined $orientation) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse from entry ${$args{orf_to_id}}->{$i}\.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing global arCOG2014 result hash".
                        "\nCould not parse from entry ${$args{orf_to_id}}->{$i}\.\n\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }

      #update status bar
      if (($i % 20) == 0) {
         &update_pbar_2(title        => "Compiling 2014 archaeal COG results",
                        label        => "Compiling arCOG2014 results for $filename",
                        progress     => ($i / $args{counter}) * 100,
                       );
      }

      #read input file
      (my $file_ref) = slurp(main_window  => $args{main_window},
                             auto_ini_ref => $args{auto_ini_ref},
                             directory    => ${$args{ini_ref}}{COG_results},
                             filename     => $filename.'_COG_'.$id
                            );
      if ($file_ref eq '0') {
         &hide_pbar_2;
         return (0);
      }

      #skip if no COG results
      next if (${$file_ref} =~ /No COG Hits found/s);

      #get COG header
      ${$file_ref}      =~ m/Blast overview(.+)Complete list of Blast results/si;
      $arCOG2014_header = $1."\n\>No further hits\n";

      #catch best defined arCOG2014 hit
      $arCOG2014_code = '';
      while ($arCOG2014_code !~ m/\w+/) {
         'reset'           =~ m/reset/;
         $arCOG2014_header =~ m/\>(.*?)\n\r?/si;
         $best_hit         = $1;

         unless (defined $best_hit) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse best hit from arCOG2014 entry $filename\_COG_$id.\n$args{orf_to_id}->{$i}",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing arCOG2014 results for file $filename".
                           "\nCould not parse best hit from arCOG2014 entry $filename\_COG_$id.\n$args{orf_to_id}->{$i}\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }

         #no more hits?
         if ($best_hit eq 'No further hits') {
            #reset COG_code and evalue to first hit
            $arCOG2014_code       = 'unclassified archaeal COG';
            $arCOG2014_acc        = $first_arCOG2014_acc;
            $arCOG2014_class      = $first_arCOG2014_class;
            $arCOG2014_annotation = $first_arCOG2014_annotation;
            $arCOG2014_cl_verbose = $first_arCOG2014_cl_verbose;
            $evalue               = $first_evalue;
            $score                = $first_score;
            $length               = $first_length;
            last;
         }

         #catch evalue, score and length
         'reset'   =~ m/reset/;
         $best_hit =~ m/Length\=(\d+)\s+Score\=(\d+)\s+Expect\=(.+)$/;
         $length   = $1;
         $score    = $2;
         $evalue   = $3;
         unless (defined $evalue) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse the evalue from entry $best_hit in file $filename, ID $id\.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing arCOG2014 results for file $filename".
                           "\nCould not parse the evalue from entry $best_hit in file $filename, ID $id\.\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }
         if ($evalue =~ /^e/) {$evalue = '1'.$evalue};

         #catch arCOG2014 code and arCOG2014 acc
         'reset'   =~ m/reset/;
         $best_hit =~ m/^(\S+)\s+arCOG2014: (.+?);\s+Class:\s+(\w+)\,?\s*(.*)\;.+?Annotation:\s+(.+) Length\=\d+/;
         ($arCOG2014_acc, $arCOG2014_code, $arCOG2014_class, $arCOG2014_cl_verbose, $arCOG2014_annotation) = ($1, $2, $3, $4, $5);
         unless (defined $arCOG2014_annotation) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse the arCOG2014 code from entry $best_hit.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing arCOG2014 results for file $filename".
                           "\nCould not parse the arCOG2014 code from entry $best_hit.\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }
         $arCOG2014_code       =~ s/;//;
         $arCOG2014_annotation =~ s/\~//g;
         #no arCOG2014 verbose description?
         unless (defined $arCOG2014_cl_verbose && $arCOG2014_cl_verbose =~ m/\w+/) {
            $arCOG2014_cl_verbose = '';
         }
         $arCOG2014_cl_verbose =~ s/\s*\;?\s*Genome\:.+//;

         #define best arCOG2014code and evalue in case no classified hit can be found
         unless (defined $first_arCOG2014_code) {
            $first_arCOG2014_code       = $arCOG2014_code;
            $first_arCOG2014_acc        = $arCOG2014_acc;
            $first_arCOG2014_class      = $arCOG2014_class;
            $first_arCOG2014_annotation = $arCOG2014_annotation;
            $first_arCOG2014_cl_verbose = $arCOG2014_cl_verbose;
            $first_evalue               = $evalue;
            $first_score                = $score;
            $first_length               = $length;
         }

         #skip if above threshold, define arCOG2014code to escape while loop, entry will be skipped anyway
         if ($evalue > ${$args{ini_ref}}{COG_cut}) {
            #reset arCOG2014_code and evalue to first hit
            $arCOG2014_code       = $first_arCOG2014_code;
            $arCOG2014_acc        = $first_arCOG2014_acc;
            $arCOG2014_class      = $first_arCOG2014_class;
            $arCOG2014_annotation = $first_arCOG2014_annotation;
            $arCOG2014_cl_verbose = $first_arCOG2014_cl_verbose;
            $evalue               = $first_evalue;
            $score                = $first_score;
            $length               = $first_length;
            unless (defined $arCOG2014_code) {
               $arCOG2014_code    = 'unclassified archaeal COG';
               $delete_arCOG2014  = 1;
            }
         }
      }
      undef $arCOG2014_header;

      #skip if above threshold
      next if ($evalue > ${$args{ini_ref}}{COG_cut});

      #cleanup annotation
      $arCOG2014_annotation =~ s/\"/\'/gs;

      #catch arCOG2014 fragment boundaries
      ${$file_ref} =~ s/^.+?Complete list of Blast results//s;
      'reset'      =~ /reset/;
      ${$file_ref} =~ m/\>\Q$arCOG2014_acc\E.*?(Query.*?)[\n\r]+\>/s;
      $match       = $1;

      #catch if only one arCOG2014 hit present
      unless (defined $match) {
         'reset'      =~ /reset/;
         ${$file_ref} =~ m/\>\Q$arCOG2014_acc\E.*?(Query.*?)\s+(Database:|Lambda)/s;
         $match       = $1;
      }
      unless (defined $match) {
         'reset'      =~ /reset/;
         ${$file_ref} =~ m/\>\w+|\Q$arCOG2014_acc\E.*?(Query.*?)\s+(Database:|Lambda)/s;
         $match       = $1;
      }
      unless (defined $match) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse arCOG2014 alignment from entry $filename.\n\n$args{orf_to_id}->{$i}",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing arCOG2014 results $arCOG2014_code for file $filename".
                        "\nCould not parse arCOG2014 alignment for arCOG2014 acc $arCOG2014_acc from entry $filename.\n$args{orf_to_id}->{$i}\n${$file_ref}\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }

      #catch arCOG2014 left_bd
      'reset' =~ /reset/;
      $match  =~ m/Query\:\s+(\d+)\s+/s;
      $arCOG2014_left_bd = $1;
      unless (defined $arCOG2014_left_bd) {
         $match =~ m/Query\s+(\d+)\s+/s;
         $arCOG2014_left_bd = $1;
      }

      #catch arCOG2014 right_bd
      'reset' =~ /reset/;
      $match  =~ m/.*Query\:\s+\d+\s+\D+\s+(\d+)/s;
      $arCOG2014_right_bd = $1;
      unless (defined $arCOG2014_right_bd) {
         $match =~ m/.*Query\s+\d+\s+\D+\s+(\d+)/s;
         $arCOG2014_right_bd = $1;
      }

      unless (defined $arCOG2014_left_bd && defined $arCOG2014_right_bd) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse arCOG2014 alignment boundaries from entry $filename",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing arCOG2014 results for file $filename".
                        "\nCould not parse arCOG2014 alignment boundaries from entry $match in entry ${$file_ref}.\n\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }
      #create boundary entry
      if ($orientation eq 'sense') {
         $arCOG2014_right_bd = $left_bd + ($arCOG2014_right_bd * 3) - 1;
         $arCOG2014_left_bd  = $left_bd + ($arCOG2014_left_bd * 3)  - 3;

         #modify boundaries is necessary
         if ($arCOG2014_left_bd == $left_bd) {
            $arCOG2014_left_bd += 3;  #add one codon to avoid equal to start position
         }
         if ($arCOG2014_right_bd == $right_bd) {
            $arCOG2014_right_bd -= 3;  #remove one codon to avoid equal to stop position
         }

         $boundary_arCOG2014  = '     COG_match       '.$arCOG2014_left_bd.'..'.$arCOG2014_right_bd."\n";
      } elsif ($orientation eq 'antisense') {
         my $length = $arCOG2014_right_bd - $arCOG2014_left_bd;
         $arCOG2014_right_bd = $right_bd - ($arCOG2014_left_bd * 3) + 3;
         $arCOG2014_left_bd  = $arCOG2014_right_bd - ($length * 3) - 2;

         #modify boundaries is necessary
         if ($arCOG2014_left_bd == $left_bd) {
            $arCOG2014_left_bd += 3;  #add one codon to avoid equal to stop position
         }
         if ($arCOG2014_right_bd == $right_bd) {
            $arCOG2014_right_bd -= 3;  #remove one codon to avoid equal to start position
         }

         $boundary_arCOG2014  = '     COG_match       complement('.$arCOG2014_left_bd.'..'.$arCOG2014_right_bd."\)\n";
      }

      #create product entry
      if ($arCOG2014_code eq 'unclassified archaeal COG') {
         $value = '/product="['.$arCOG2014_class.'] '.
                  ' Length='.$length.' Score='.$score.' Expect='.$evalue.'"'."\n";
      } else {
         $value = '/product="['.$arCOG2014_class.'] '.$arCOG2014_cl_verbose.', '.
                  'arCOG2014: '.$arCOG2014_code.', ';
         if (${$args{arCOG2014code2014_to_COG}}->{$arCOG2014_code}->{'COGcode'} =~ m/\w+/) {
            $value .= 'COG: '.${$args{arCOG2014code2014_to_COG}}->{$arCOG2014_code}->{'COGcode'}.', ';
         }
         $value .= $arCOG2014_annotation.
                   ' Length='.$length.' Score='.$score.' Expect='.$evalue.'"'."\n";
      }
      #no arCOG2014 id?
      if ($value =~ /\[\]/) {$value =~ s/\[\]/\[not assigned\]/};
      open  ENTRY, '>', \$product;
      write ENTRY;
      close ENTRY;

      #create arCOG2014 label
      $value = 'label= ['.$arCOG2014_class.'] '.$arCOG2014_cl_verbose."\n";

      #clean up arCOG2014 label
      $value =~ s/\//-/gs;
      $value = '/'.$value;

      #no arCOG2014 id?
      if ($value =~ /\[\]/) {$value =~ s/\[\]/\[unclassified arCOG2014\]/};
      open  ENTRY, '>', \$label;
      write ENTRY;
      close ENTRY;

      #create codon start and translation table tags
      $colour = '                     /colour=2'."\n";

      #combine feature
      $feature_arCOG2014 = $boundary_arCOG2014.
                           $label.
                           $colour.
                           $product;
      $key_arCOG2014     = $arCOG2014_left_bd.'_'.$arCOG2014_right_bd.'_COG_match_'.$id; #maintain original uniqe ID instead of ORFnumber

      ${$args{genbank_ref}}{$key_arCOG2014} = $feature_arCOG2014;
      push (@{$args{feature_list_ref}}, $key_arCOG2014);

      #delete unassigned COG code again, otherwise it's stored for next ORF
      if ($delete_arCOG2014) {
         $arCOG2014_code = '';
         undef $delete_arCOG2014;
      }
   }
   &hide_pbar_2;
   return;
}

sub pog2013_parser {
   my %args = @_;

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling 2013 phage COG results",
                   label        => ''
                  );
   &show_pbar_2;

   #parse through POG2013s
   foreach (my $i = 1; $i < $args{counter}; $i++) {
      my ($ORF_key, $filename, $id, $left_bd, $right_bd, $orientation, $best_hit, $evalue, $length, $score,
          $POG2013_code, $match, $POG2013_left_bd, $POG2013_right_bd, $boundary_POG2013, $colour, $label, $product,
          $feature_POG2013, $key_POG2013, $first_POG2013_code, $first_evalue, $first_score, $first_length, $POG2013_header,
          $delete_POG2013, $POG2013_acc, $first_POG2013_acc, $POG2013_phage, $first_POG2013_phage, $POG2013_annotation,
          $first_POG2013_annotation);

      #get filename, left_bd, right_bd, orientation
      'reset' =~ m/reset/;
      $args{orf_to_id}->{$i} =~ m/^(.*?)___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
      $filename    = $1;
      $id          = $2;
      $left_bd     = $3;
      $right_bd    = $4;
      $orientation = $5;
      unless (defined $orientation) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse from entry ${$args{orf_to_id}}->{$i}\.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing global POG2013 result hash".
                        "\nCould not parse from entry ${$args{orf_to_id}}->{$i}\.\n\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }

      #update status bar
      if (($i % 20) == 0) {
         &update_pbar_2(title        => "Compiling 2013 phage COG results",
                        label        => "Compiling POG2013 results for $filename",
                        progress     => ($i / $args{counter}) * 100,
                       );
      }

      #read input file
      (my $file_ref) = slurp(main_window  => $args{main_window},
                             auto_ini_ref => $args{auto_ini_ref},
                             directory    => ${$args{ini_ref}}{COG_results},
                             filename     => $filename.'_COG_'.$id
                            );
      if ($file_ref eq '0') {
         &hide_pbar_2;
         return (0);
      }

      #skip if no COG results
      next if (${$file_ref} =~ /No COG Hits found/s);

      #get COG header
      ${$file_ref}      =~ m/Blast overview(.+)Complete list of Blast results/si;
      $POG2013_header   = $1."\n\>No further hits\n";

      #catch best defined POG2013 hit
      $POG2013_code = '';
      while ($POG2013_code !~ m/\w+/) {
         'reset'           =~ m/reset/;
         $POG2013_header   =~ m/\>(.*?)\n\r?/si;
         $best_hit         = $1;

         unless (defined $best_hit) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse best hit from POG2013 entry $filename\_COG_$id.\n$args{orf_to_id}->{$i}",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing POG2013 results for file $filename".
                           "\nCould not parse best hit from POG2013 entry $filename\_COG_$id.\n$args{orf_to_id}->{$i}\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }

         #no more hits?
         if ($best_hit eq 'No further hits') {
            #reset COG_code and evalue to first hit
            $POG2013_code       = 'unclassified POG';
            $POG2013_acc        = $first_POG2013_acc;
            $POG2013_phage      = $first_POG2013_phage;
            $POG2013_annotation = $first_POG2013_annotation;
            $evalue             = $first_evalue;
            $score              = $first_score;
            $length             = $first_length;
            last;
         }

         #catch evalue, score and length
         'reset'   =~ m/reset/;
         $best_hit =~ m/Length\=(\d+)\s+Score\=(\d+)\s+Expect\=(.+)$/;
         $length   = $1;
         $score    = $2;
         $evalue   = $3;
         unless (defined $evalue) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse the evalue from entry $best_hit in file $filename, ID $id\.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing POG2013 results for file $filename".
                           "\nCould not parse the evalue from entry $best_hit in file $filename, ID $id\.\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }
         if ($evalue =~ /^e/) {$evalue = '1'.$evalue};

         #catch POG2013 code and POG2013 acc
         'reset'   =~ m/reset/;
         $best_hit =~ m/(\S+)\s+POG2013\:\s+(POG\d+)\;\s+Phage\:\s+(.+)\;\s+Annotation\:\s+(.+) Length\=\d+/;
         ($POG2013_acc, $POG2013_code, $POG2013_phage, $POG2013_annotation) = ($1, $2, $3, $4);
         #no defined?
         unless (defined $POG2013_annotation) {
            'reset'   =~ m/reset/;
            $best_hit =~ m/(\S+)\s+POG2013\:\s+([^\;]+?)\;\s+Phage\:\s+(.+)\;\s+Annotation\:\s+(.+) Length\=\d+/;
            ($POG2013_acc, $POG2013_code, $POG2013_phage, $POG2013_annotation) = ($1, $2, $3, $4);
         }

         unless (defined $POG2013_annotation) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse the POG2013 code from entry $best_hit.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing POG2013 results for file $filename".
                           "\nCould not parse the POG2013 code from entry $best_hit.\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }
         $POG2013_annotation =~ s/\~//g;
         $POG2013_annotation =~ s/\s+$//;
         $POG2013_phage      =~ s/\s+$//;

         #define best POG2013code and evalue in case no classified hit can be found
         unless (defined $first_POG2013_code) {
            $first_POG2013_code       = $POG2013_code;
            $first_POG2013_acc        = $POG2013_acc;
            $first_POG2013_phage      = $POG2013_phage;
            $first_POG2013_annotation = $POG2013_annotation;
            $first_evalue             = $evalue;
            $first_score              = $score;
            $first_length             = $length;
         }

         #skip if above threshold, define POG2013code to escape while loop, entry will be skipped anyway
         if ($evalue > ${$args{ini_ref}}{COG_cut}) {
            #reset POG2013_code and evalue to first hit
            $POG2013_code       = $first_POG2013_code;
            $POG2013_acc        = $first_POG2013_acc;
            $POG2013_phage      = $first_POG2013_phage;
            $POG2013_annotation = $first_POG2013_annotation;
            $evalue             = $first_evalue;
            $score              = $first_score;
            $length             = $first_length;
            unless (defined $POG2013_code) {
               $POG2013_code    = 'unclassified POG';
               $delete_POG2013  = 1;
            }
         }
      }
      undef $POG2013_header;

      #skip if above threshold
      next if ($evalue > ${$args{ini_ref}}{COG_cut});

      #cleanup annotation
      $POG2013_annotation =~ s/\"/\'/gs;

      #catch POG2013 fragment boundaries
      ${$file_ref} =~ s/^.+?Complete list of Blast results//s;
      'reset'      =~ /reset/;
      if (${$file_ref} =~ m/\>.+?\|/) {
         ${$file_ref} =~ m/\>.+?\|\Q$POG2013_acc\E.*?(Query.*?)[\n\r]+\>/s;
         $match       = $1;
      } else {
         ${$file_ref} =~ m/\>\Q$POG2013_acc\E.*?(Query.*?)[\n\r]+\>/s;
         $match       = $1;
      }

      #catch if only one POG2013 hit present
      unless (defined $match) {
         'reset'      =~ /reset/;
         if (${$file_ref} =~ m/\>.+?|/) {
            ${$file_ref} =~ m/\>.+?|\Q$POG2013_acc\E.*?(Query.*?)\s+Database:/s;
            $match       = $1;
         } else {
            ${$file_ref} =~ m/\>\Q$POG2013_acc\E.*?(Query.*?)\s+Database:/s;
            $match       = $1;
         }
      }
      #catch if only one POG2013 hit present and no terminal header present
      unless (defined $match) {
         'reset'      =~ /reset/;
         if (${$file_ref} =~ m/\>.+?|/) {
            ${$file_ref} =~ m/\>.+?|\Q$POG2013_acc\E.*?(Query.+)/s;
            $match       = $1;
         } else {
            ${$file_ref} =~ m/\>\Q$POG2013_acc\E.*?(Query.+)/s;
            $match       = $1;
         }
      }
      #try to catch anything
      unless (defined $match) {
         'reset'      =~ /reset/;
         ${$file_ref} =~ m/\>\Q$POG2013_acc\E.*?(Query.+)\n\n\n/s;
         $match       = $1;
      }
      unless (defined $match) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse POG2013 alignment from entry $filename.\n\n$args{orf_to_id}->{$i}",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing POG2013 results $POG2013_code, $POG2013_acc for file $filename".
                        "\nCould not parse POG2013 alignment from entry $filename.\n$args{orf_to_id}->{$i}\n${$file_ref}\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }

      #catch POG2013 left_bd
      'reset' =~ /reset/;
      $match  =~ m/Query\:\s+(\d+)\s+/s;
      $POG2013_left_bd = $1;
      unless (defined $POG2013_left_bd) {
         $match =~ m/Query\s+(\d+)\s+/s;
         $POG2013_left_bd = $1;
      }

      #catch POG2013 right_bd
      'reset' =~ /reset/;
      $match  =~ m/.*Query\:\s+\d+\s+\D+\s+(\d+)/s;
      $POG2013_right_bd = $1;
      unless (defined $POG2013_right_bd) {
         $match =~ m/.*Query\s+\d+\s+\D+\s+(\d+)/s;
         $POG2013_right_bd = $1;
      }

      unless (defined $POG2013_left_bd && defined $POG2013_right_bd) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse POG2013 alignment boundaries from entry $filename",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing POG2013 results for file $filename".
                        "\nCould not parse POG2013 alignment boundaries from entry $match in entry ${$file_ref}.\n\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }
      #create boundary entry
      if ($orientation eq 'sense') {
         $POG2013_right_bd = $left_bd + ($POG2013_right_bd * 3) - 1;
         $POG2013_left_bd  = $left_bd + ($POG2013_left_bd * 3)  - 3;

         #modify boundaries is necessary
         if ($POG2013_left_bd == $left_bd) {
            $POG2013_left_bd += 3;  #add one codon to avoid equal to start position
         }
         if ($POG2013_right_bd == $right_bd) {
            $POG2013_right_bd -= 3;  #remove one codon to avoid equal to stop position
         }

         $boundary_POG2013  = '     COG_match       '.$POG2013_left_bd.'..'.$POG2013_right_bd."\n";
      } elsif ($orientation eq 'antisense') {
         my $length = $POG2013_right_bd - $POG2013_left_bd;
         $POG2013_right_bd = $right_bd - ($POG2013_left_bd * 3) + 3;
         $POG2013_left_bd  = $POG2013_right_bd - ($length * 3) - 2;

         #modify boundaries is necessary
         if ($POG2013_left_bd == $left_bd) {
            $POG2013_left_bd += 3;  #add one codon to avoid equal to stop position
         }
         if ($POG2013_right_bd == $right_bd) {
            $POG2013_right_bd -= 3;  #remove one codon to avoid equal to start position
         }

         $boundary_POG2013  = '     COG_match       complement('.$POG2013_left_bd.'..'.$POG2013_right_bd."\)\n";
      }

      #create product entry
      if ($POG2013_code eq 'unclassified POG') {
         $value = '/product="'.$POG2013_phage.', '.
                  ' Length='.$length.' Score='.$score.' Expect='.$evalue.'"'."\n";
      } else {
         $value = '/product="'.$POG2013_phage.', '.
                  'POG2013: '.$POG2013_code.', '.
                  $POG2013_annotation.
                   ' Length='.$length.' Score='.$score.' Expect='.$evalue.'"'."\n";
      }
      #no POG2013 id?
      if ($value =~ /\[\]/) {$value =~ s/\[\]/\[not assigned\]/};
      open  ENTRY, '>', \$product;
      write ENTRY;
      close ENTRY;

      #create POG2013 label
      $value = 'label= ['.$POG2013_code.'] '.$POG2013_annotation."\n";

      #clean up POG2013 label
      $value =~ s/\//-/gs;
      $value = '/'.$value;

      #no POG2013 id?
      if ($value =~ /\[\]/) {$value =~ s/\[\]/\[unclassified POG2013\]/};
      open  ENTRY, '>', \$label;
      write ENTRY;
      close ENTRY;

      #create codon start and translation table tags
      $colour = '                     /colour=2'."\n";

      #combine feature
      $feature_POG2013 = $boundary_POG2013.
                         $label.
                         $colour.
                         $product;
      $key_POG2013     = $POG2013_left_bd.'_'.$POG2013_right_bd.'_COG_match_'.$id; #maintain original uniqe ID instead of ORFnumber

      ${$args{genbank_ref}}{$key_POG2013} = $feature_POG2013;
      push (@{$args{feature_list_ref}}, $key_POG2013);

      #delete unassigned COG code again, otherwise it's stored for next ORF
      if ($delete_POG2013) {
         $POG2013_code = '';
         undef $delete_POG2013;
      }
   }
   &hide_pbar_2;
   return;
}

sub cog2003_parser {
   my %args = @_;

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling COG results",
                   label        => ''
                  );
   &show_pbar_2;

   #parse through COGs
   foreach (my $i = 1; $i < $args{counter}; $i++) {
      my ($ORF_key, $filename, $id, $left_bd, $right_bd, $orientation, $best_hit, $evalue, $length, $score,
          $COG_code, $match, $COG_left_bd, $COG_right_bd, $boundary_COG, $colour, $label, $product,
          $feature_COG, $key_COG, $first_COG_code, $first_evalue, $first_score, $first_length, $COG_header,
          $delete_COG);

      #get filename, left_bd, right_bd, orientation
      'reset' =~ m/reset/;
      $args{orf_to_id}->{$i} =~ m/^(.*?)___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
      $filename    = $1;
      $id          = $2;
      $left_bd     = $3;
      $right_bd    = $4;
      $orientation = $5;
      unless (defined $orientation) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse from entry ${$args{orf_to_id}}->{$i}\.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing global COG result hash".
                        "\nCould not parse from entry ${$args{orf_to_id}}->{$i}\.\n\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }

      #update status bar
      if (($i % 20) == 0) {
         &update_pbar_2(title        => "Compiling COG results",
                        label        => "Compiling COG results for $filename",
                        progress     => ($i / $args{counter}) * 100,
                       );
      }

      #read input file
      (my $file_ref) = slurp(main_window => $args{main_window},
                             auto_ini_ref => $args{auto_ini_ref},
                             directory   => ${$args{ini_ref}}{COG_results},
                             filename    => $filename.'_COG_'.$id
                            );
      if ($file_ref eq '0') {
         &hide_pbar_2;
         return (0);
      }

      #skip if no COG results
      next if (${$file_ref} =~ /No COG Hits found/s);

      #get COG header
      ${$file_ref} =~ m/Blast overview(.+)Complete list of Blast results/si;
      $COG_header = $1."\n\>No further hits\n";

      #catch best defined COG hit
      $COG_code = '';
      delete ${$args{COGcode_to_header}}->{$COG_code};
      while (!defined ${$args{COGcode_to_header}}->{$COG_code}) {
         'reset' =~ m/reset/;
         $COG_header =~ m/\>(.*?)\n\r?/si;
         $best_hit = $1;

         unless (defined $best_hit) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse best hit from COG entry $filename\_COG_$id.\n$args{orf_to_id}->{$i}",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing COG results for file $filename".
                           "\nCould not parse best hit from COG entry $filename\_COG_$id.\n$args{orf_to_id}->{$i}\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }

         #no more hits?
         if ($best_hit eq 'No further hits') {
            #reset COG_code and evalue to first hit
            $COG_code = $first_COG_code;
            $evalue   = $first_evalue;
            $score    = $first_score;
            $length   = $first_length;
            unless (defined ${$args{COGcode_to_header}}->{$COG_code}) {
               ${$args{COGcode_to_header}}->{$COG_code} = 'unclassified COG';
               $delete_COG = 1;
            }
            last;
         }

         #catch evalue, score and length
         'reset' =~ m/reset/;
         $best_hit =~ m/Length\=(\d+)\s+Score\=(\d+)\s+Expect\=(.+)$/;
         $length = $1;
         $score  = $2;
         $evalue = $3;
         unless (defined $evalue) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse the evalue from entry $best_hit in file $filename, ID $id\.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing COG results for file $filename".
                           "\nCould not parse the evalue from entry $best_hit in file $filename, ID $id\.\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }
         if ($evalue =~ /^e/) {$evalue = '1'.$evalue};

         #catch COG code
         'reset' =~ m/reset/;
         $best_hit =~ m/^(\S+)\s/;
         $COG_code = $1;
         unless (defined $COG_code) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse the COG code from entry $best_hit.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing COG results for file $filename".
                           "\nCould not parse the COG code from entry $best_hit.\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }
         $COG_code =~ s/~$//g;

         #define best COGcode and evalue in case no classified hit can be found
         unless (defined $first_COG_code) {
            $first_COG_code = $COG_code;
            $first_evalue   = $evalue;
            $first_score    = $score;
            $first_length   = $length;
         }

         unless (defined ${$args{COGcode_to_header}}->{$COG_code}) {
            $COG_header =~ s/\>(lcl\|)?\Q$best_hit\E.*?\n\r?//s;
         }

         #skip if above threshold, define COGcode_to_header to escape while loop, entry will be skipped anyway
         if ($evalue > ${$args{ini_ref}}{COG_cut}) { # && !defined ${$args{COGcode_to_header}}->{$COG_code}) {
            #reset COG_code and evalue to first hit
            $COG_code = $first_COG_code;
            $evalue   = $first_evalue;
            $score    = $first_score;
            $length   = $first_length;
            unless (defined ${$args{COGcode_to_header}}->{$COG_code}) {
               ${$args{COGcode_to_header}}->{$COG_code} = 'unclassified COG';
               $delete_COG = 1;
            }
         }
      }
      undef $COG_header;

      #skip if above threshold
      next if ($evalue > ${$args{ini_ref}}{COG_cut});

      #catch COG fragment boundaries
      'reset' =~ /reset/;
      ${$file_ref} =~ m/\>\Q$COG_code\E.*?(Query.*?)[\n\r]+\>/s;
      $match = $1;
      unless (defined $match) {
         ${$file_ref} =~ m/\>lcl\|\Q$COG_code\E.*?(Query.*?)[\n\r]+\>/s;
         $match = $1;
      }

      #catch if only one COG hit present
      unless (defined $match) {
         'reset' =~ /reset/;
         ${$file_ref} =~ m/\>\Q$COG_code\E.*?(Query.*?)\s+Database:/s;
         $match = $1;
         unless (defined $match) {
            ${$file_ref} =~ m/\>lcl\|\Q$COG_code\E.*?(Query.*?)\s+Database:/s;
            $match = $1;
         }
         #if nothing captured yet, catch all
         unless (defined $match) {
            ${$file_ref} =~ m/\>(lcl\|)?\Q$COG_code\E.*?(Query.+)/s;
            $match = $2;
         }
      }
      unless (defined $match) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse COG alignment from entry $filename.\n\n$args{orf_to_id}->{$i}",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing COG results $COG_code for file $filename".
                        "\nCould not parse COG alignment from entry $filename.\n$args{orf_to_id}->{$i}\n${$file_ref}\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }

      #catch COG left_bd
      'reset' =~ /reset/;
      $match =~ m/Query\:\s+(\d+)\s+/s;
      $COG_left_bd = $1;
      unless (defined $COG_left_bd) {
         $match =~ m/Query\s+(\d+)\s+/s;
         $COG_left_bd = $1;
      }

      #catch COG right_bd
      'reset' =~ /reset/;
      $match =~ m/.*Query\:\s+\d+\s+\D+\s+(\d+)/s;
      $COG_right_bd = $1;
      unless (defined $COG_right_bd) {
         $match =~ m/.*Query\s+\d+\s+\D+\s+(\d+)/s;
         $COG_right_bd = $1;
      }

      unless (defined $COG_left_bd && defined $COG_right_bd) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse COG alignment boundaries from entry $filename",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing COG results for file $filename".
                        "\nCould not parse COG alignment boundaries from entry $match in entry ${$file_ref}.\n\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }
      #create boundary entry
      if ($orientation eq 'sense') {
         $COG_right_bd = $left_bd + ($COG_right_bd * 3) - 1;
         $COG_left_bd  = $left_bd + ($COG_left_bd * 3) - 3;

         #modify boundaries is necessary
         if ($COG_left_bd == $left_bd) {
            $COG_left_bd += 3;  #add one codon to avoid equal to start position
         }
         if ($COG_right_bd == $right_bd) {
            $COG_right_bd -= 3;  #remove one codon to avoid equal to stop position
         }

         $boundary_COG  = '     COG_match       '.$COG_left_bd.'..'.$COG_right_bd."\n";
      } elsif ($orientation eq 'antisense') {
         my $length = $COG_right_bd - $COG_left_bd;
         $COG_right_bd = $right_bd - ($COG_left_bd * 3) + 3;
         $COG_left_bd  = $COG_right_bd - ($length * 3) - 2;

         #modify boundaries is necessary
         if ($COG_left_bd == $left_bd) {
            $COG_left_bd += 3;  #add one codon to avoid equal to stop position
         }
         if ($COG_right_bd == $right_bd) {
            $COG_right_bd -= 3;  #remove one codon to avoid equal to start position
         }

         $boundary_COG  = '     COG_match       complement('.$COG_left_bd.'..'.$COG_right_bd."\)\n";
      }

      #create product entry
      ${$args{COGcode_to_header}}->{$COG_code} =~ s/\"/\'/gs;

      $value = '/product="['.${$args{COGcode_to_letter}}->{$COG_code}.'] '.
               ${$args{COGcode_to_number}}->{$COG_code}.' '.
               ${$args{COGcode_to_header}}->{$COG_code}.
               ' Length='.$length.' Score='.$score.' Expect='.$evalue.'"'."\n";
      #no COG id?
      if ($value =~ /\[\]/) {$value =~ s/\[\]/\[not assigned\]/};
      open  ENTRY, '>', \$product;
      write ENTRY;
      close ENTRY;

      #create COG label
      $value = 'label= ['.${$args{COGcode_to_letter}}->{$COG_code}.'] '.
               ${$args{COGletter_to_family}}->{${$args{COGcode_to_letter}}->{$COG_code}}.
               "\n";
      #clean up COG label
      $value =~ s/\//-/gs;
      $value = '/'.$value;

      #no COG id?
      if ($value =~ /\[\]/) {$value =~ s/\[\]/\[unclassified COG\]/};
      open  ENTRY, '>', \$label;
      write ENTRY;
      close ENTRY;

      #create codon start and translation table tags
      $colour = '                     /colour=2'."\n";

      #combine feature
      $feature_COG = $boundary_COG.
                     $label.
                     $colour.
                     $product;
      $key_COG     = $COG_left_bd.'_'.$COG_right_bd.'_COG_match_'.$id; #maintain original unique ID instead of ORFnumber

      ${$args{genbank_ref}}{$key_COG} = $feature_COG;
      push (@{$args{feature_list_ref}}, $key_COG);

      #delete unassigned COG code again, otherwise it's stored for next ORF
      if ($delete_COG) {
         delete ${$args{COGcode_to_header}}->{$COG_code};
         undef $delete_COG;
      }
   }
   &hide_pbar_2;
   return;
}

sub translate_arCOG {
   my %args = @_;

   my ($file_ref, $arCOG_dir, $arCOGdef, @arCOGs);

   #arCOGdef.tab contains a link between COG, arCOG, functional class and arCOGannotation
   #expected default location is in the arCOG directory
   ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
   $arCOG_dir = $1;
   $arCOGdef  = 'arCOGdef.tab';

   unless (-e $arCOG_dir.'/'.$arCOGdef) {
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $arCOG_dir,
                                                    -title      => 'Select arCOG def file',
                                                    -defaultextension => '',
                                                   );
      unless (defined $file) {return (0)};
      'reset'        =~ /reset/;
      $file          =~ m/(.+)\/(.+)$/;
      $arCOG_dir     = $1;
      $arCOGdef      = $2;
   }

   unless (-e $arCOG_dir.'/'.$arCOGdef) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No COG def file selected, aborting.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error selecting COG def file:".
                     "\ndir...$arCOG_dir ...file:...$arCOGdef ...\n\n";
      close ERRORLOG;
      return (0);
   }

   #read arCOG codes and defs
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $arCOG_dir,
                      filename     => $arCOGdef
                     );
   ${$file_ref} = "\n".${$file_ref};
   @arCOGs = split /\n/, ${$file_ref};
   undef ${$file_ref};
   foreach my $entry (@arCOGs) {
      next unless ($entry =~ m/\w+/);
      my ($arCOG, $COG, $class, $annotation);
      'reset' =~ m/reset/;
      $entry =~ m/^([^\t]+)\t([^\t]*)\t([^\t]*)\t(.+)/;
      ($arCOG, $COG, $class, $annotation) = ($1, $2, $3, $4);

      unless (defined $annotation) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing arCOG definition for entry $entry.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing arCOG definition for file $arCOGdef".
                        "\nError parsing arCOG definition for entry $entry.\n\n";
         close ERRORLOG;
      }
      if ($COG !~ /\w+/) {$COG = 'not def '};
      $annotation =~ s/\n//gs;

      ${$args{arCOGcode_to_COG}}->{$arCOG} = {(COGcluster      => $COG,
                                               COGclass        => $class,
                                               arCOGannotation => $annotation,
                                               organism        => ''
                                             )};
   }
   undef @arCOGs;

   #add acc to organism qualifier
   #arCOGcsv contains a link between domain ID and genome name
   #arCOG has been MODIFIED and the first ID has been replaced by ref acc
   #this CHANGE allows to properly parse COG Blast results
   #expected default location is in the arCOG directory
   ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
   $arCOG_dir = $1;
   $arCOGdef  = 'arCOG2.txt';

   unless (-e $arCOG_dir.'/'.$arCOGdef) {
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $arCOG_dir,
                                                    -title      => 'Select arCOG2 ID to genome file',
                                                    -defaultextension => '',
                                                   );
      unless (defined $file) {return (0)};
      'reset'        =~ /reset/;
      $file          =~ m/(.+)\/(.+)$/;
      $arCOG_dir     = $1;
      $arCOGdef      = $2;
   }

   unless (-e $arCOG_dir.'/'.$arCOGdef) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No arCOG2 ID file selected, aborting.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error selecting COG ID file:".
                     "\ndir...$arCOG_dir ...file:...$arCOGdef ...\n\n";
      close ERRORLOG;
      return (0);
   }

   #read arCOG IDs and genome names
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $arCOG_dir,
                      filename     => $arCOGdef
                     );
   ${$file_ref} = "\n".${$file_ref};
   @arCOGs = split /\n/, ${$file_ref};
   undef ${$file_ref};
   foreach my $entry (@arCOGs) {
      next unless ($entry =~ m/\w+/);
      my ($ref_id, $arCOG, $genome, $arCOGid);
      'reset' =~ m/reset/;
      $entry =~ m/^([^\t]+)\t([^\t]+)\t([^\t]+)\t[^\t]*\t[^\t]*\t[^\t]*\t([^\t]*)\t/;
      ($ref_id, $genome, $arCOGid, $arCOG) = ($1, $2, $3, $4);

      unless (defined $arCOG) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing arCOG2 genome to ID for entry $entry.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing arCOG2 genome to ID for file $arCOGdef".
                        "\nError parsing arCOG definition for entry $entry.\n\n";
         close ERRORLOG;
      }
      $ref_id =~ s/\.\d+\s*$//g;
      ${$args{arCOGacc_to_org}}->{$ref_id} = {(arCOGid => $arCOGid,
                                               genome  => $genome,
                                               arCOG   => $arCOG
                                             )};
   }
   undef @arCOGs;
   return (1);
}

sub translate_arCOG2014 {
   my %args = @_;

   my ($file_ref, $arCOG_dir, $arCOGdef, @arCOGs);

   #arCOGdef.tab contains a link between COG, arCOG, functional class and arCOGannotation
   #expected default location is in the arCOG directory
   ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
   $arCOG_dir = $1;
   $arCOGdef  = 'ar14.arCOGdef.tab';

   unless (-e $arCOG_dir.'/'.$arCOGdef) {
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $arCOG_dir,
                                                    -title      => 'Select arCOG2014 def file',
                                                    -defaultextension => '',
                                                   );
      unless (defined $file) {return (0)};
      'reset'        =~ /reset/;
      $file          =~ m/(.+)\/(.+)$/;
      $arCOG_dir     = $1;
      $arCOGdef      = $2;
   }

   unless (-e $arCOG_dir.'/'.$arCOGdef) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No arCOG2014 def file selected, aborting.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error selecting arCOG2014 def file:".
                     "\ndir...$arCOG_dir ...file:...$arCOGdef ...\n\n";
      close ERRORLOG;
      return (0);
   }

   #read arCOG codes and defs
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $arCOG_dir,
                      filename     => $arCOGdef
                     );
   ${$file_ref} = "\n".${$file_ref};
   @arCOGs = split /\n/, ${$file_ref};
   undef ${$file_ref};
   foreach my $entry (@arCOGs) {
      next unless ($entry =~ m/\w+/);
      my ($arCOG, $COG, $class, $annotation, $gene_name);
      'reset' =~ m/reset/;
      $entry =~ m/^([^\t]+)\t(\w+)\t([^\t]*)\t([^\t]+)\t([^\t]*)\t/;
      ($arCOG, $class, $gene_name, $annotation, $COG) = ($1, $2, $3, $4, $5);

      unless (defined $COG) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing arCOG2014 definition for entry $entry.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing arCOG2014 definition for file $arCOGdef".
                        "\nError parsing arCOG2014 definition for entry $entry.\n\n";
         close ERRORLOG;
      }
      if ($COG !~ /\w+/) {$COG = 'not def '};
      if ($gene_name eq '-') {$gene_name = ''};
      $annotation =~ s/\n//gs;

      ${$args{arCOGcode2014_to_COG}}->{$arCOG} = {(COGcode        => $COG,
                                                   arCOGclass      => $class,
                                                   arCOGannotation => $annotation,
                                                   gene_name       => $gene_name,
                                                   organism        => ''
                                                 )};
   }
   undef @arCOGs;

   #read prot ID, genome name, arCOGs and COG
   #expected default location is in the arCOG directory
   ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
   $arCOG_dir = $1;
   $arCOGdef  = 'arCOG14.txt';

   unless (-e $arCOG_dir.'/'.$arCOGdef) {
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $arCOG_dir,
                                                    -title      => 'Select arCOG2014 file (ar14.arCOG.csv)',
                                                    -defaultextension => '',
                                                   );
      unless (defined $file) {return (0)};
      'reset'        =~ /reset/;
      $file          =~ m/(.+)\/(.+)$/;
      $arCOG_dir     = $1;
      $arCOGdef      = $2;
   }

   unless (-e $arCOG_dir.'/'.$arCOGdef) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No arCOG2014 file (ar14.arCOG.csv) selected, aborting.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error selecting COG ID file (ar14.arCOG.csv):".
                     "\ndir...$arCOG_dir ...file:...$arCOGdef ...\n\n";
      close ERRORLOG;
      return (0);
   }

   #read arCOG IDs and genome names
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $arCOG_dir,
                      filename     => $arCOGdef
                     );
   ${$file_ref} = "\n".${$file_ref};
   @arCOGs = split /\n/, ${$file_ref};
   undef ${$file_ref};
   foreach my $entry (@arCOGs) {
      next unless ($entry =~ m/\w+/);
      my ($ref_id, $arCOG, $genome, $COG);
      my @tmp = split /\,/, $entry;

      #ignore if no arCOG defined
      next unless (defined $tmp[6]);

      $ref_id = $tmp[0];
      $genome = $tmp[1];
      $arCOG  = $tmp[6];
      if (defined $tmp[8]) {
         $COG = $tmp[8];
      } else {
         $COG = '';
      }

      unless (defined $arCOG) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing arCOG2014 genome to ID for entry $entry.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing arCOG2014 genome to ID for file $arCOGdef".
                        "\nError parsing arCOG2014 definition for entry $entry.\n\n";
         close ERRORLOG;
      }
      $genome =~ s/_uid\d+//;
      $ref_id =~ s/\.\d+$//;

      ${$args{arCOGacc2014_to_org}}->{$ref_id} = {(arCOG2014 => $arCOG,
                                                   genome    => $genome,
                                                   COG       => $arCOG
                                                 )};
   }
   undef @arCOGs;

   #read arCOG2014 assignment (funclass.tab)
   #expected default location is in the arCOG2014 directory
   ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
   $arCOG_dir  = $1;
   $arCOGdef = 'funclass.tab';

   unless (-e $arCOG_dir.'/'.$arCOGdef) {
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $arCOG_dir,
                                                    -title      => 'Select arCOG2014 fun file (funclass.tab)',
                                                    -defaultextension => '',
                                                   );
      unless (defined $file) {return (0)};
      'reset'        =~ /reset/;
      $file          =~ m/(.+)\/(.+)$/;
      $arCOG_dir     = $1;
      $arCOGdef      = $2;
   }

   unless (-e $arCOG_dir.'/'.$arCOGdef) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No arCOG2014 fun file selected, aborting.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error selecting arCOG2014 fun file:".
                     "\ndir...$arCOG_dir ...file:...$arCOGdef ...\n\n";
      close ERRORLOG;
      return (0);
   }

   #read COG function file
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $arCOG_dir,
                      filename     => $arCOGdef
                     );
   ${$file_ref} = "\n".${$file_ref};
   @arCOGs = split /\n/, ${$file_ref};

   undef ${$file_ref};

   foreach my $entry (@arCOGs) {
      next if ($entry =~ m/^\d/);
      my ($letter, $family);
      'reset' =~ m/reset/;
      $entry  =~ m/(\w)\t\d\t(.+)\n?\r?/;
      $letter = $1;
      $family = $2;
      ${$args{COGletter_to_family}}->{$letter} = $family;
   }
   undef @arCOGs;

   return (1);
}

sub translate_COG2008 {
   my %args = @_;

   my ($file_ref, $COG2008_dir, $COG2008_file, @COG2008);

   #genomes.txt contains a link between genome abbreviation and full genome name
   #expected default location is in the COG2008 directory
   ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
   $COG2008_dir  = $1;
   $COG2008_file = 'genomes.txt';

   unless (-e $COG2008_dir.'/'.$COG2008_file) {
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $COG2008_dir,
                                                    -title      => 'Select COG2008 genomes file',
                                                    -defaultextension => '',
                                                   );
      unless (defined $file) {return (0)};
      'reset'        =~ /reset/;
      $file          =~ m/(.+)\/(.+)$/;
      $COG2008_dir   = $1;
      $COG2008_file  = $2;
   }

   unless (-e $COG2008_dir.'/'.$COG2008_file) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No COG2008 genomes file selected, aborting.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error selecting COG208 genomes file:".
                     "\ndir...$COG2008_dir ...file:...$COG2008_file ...\n\n";
      close ERRORLOG;
      return (0);
   }

   #read COG2008 genomes
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $COG2008_dir,
                      filename     => $COG2008_file
                     );
   ${$file_ref} = "\n".${$file_ref};
   @COG2008     = split /\n/, ${$file_ref};
   undef ${$file_ref};
   foreach my $entry (@COG2008) {
      next unless ($entry =~ m/\w+/);
      my ($COG_abb, $COG_genome);
      'reset' =~ m/reset/;
      $entry  =~ m/^([^\t]+)\t[^\t]+\t([^\t]+)\t/;
      ($COG_abb, $COG_genome) = ($1, $2);

      unless (defined $COG_genome) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing COG2008 genomes for entry $entry.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing COG2008 genomes for file $COG2008_file".
                        "\nError parsing COG2008 genomes for entry $entry.\n\n";
         close ERRORLOG;
      }

      ${$args{COG2008genomes}}->{$COG_abb} = $COG_genome;
   }
   undef @COG2008;

   #link COG2008 ID to classes and annotation
   #expected default location is in the COG2008 directory
   ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
   $COG2008_dir  = $1;
   $COG2008_file = 'cogdef.txt';

   unless (-e $COG2008_dir.'/'.$COG2008_file) {
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $COG2008_dir,
                                                    -title      => 'Select COGdef file',
                                                    -defaultextension => '',
                                                   );
      unless (defined $file) {return (0)};
      'reset'          =~ /reset/;
      $file            =~ m/(.+)\/(.+)$/;
      $COG2008_dir     = $1;
      $COG2008_file    = $2;
   }

   unless (-e $COG2008_dir.'/'.$COG2008_file) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No COG2008def file selected, aborting.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error selecting COG2008 def file:".
                     "\ndir...$COG2008_dir ...file:...$COG2008_file ...\n\n";
      close ERRORLOG;
      return (0);
   }

   #read COG2008 IDs and classes
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $COG2008_dir,
                      filename     => $COG2008_file
                     );
   ${$file_ref} = "\n".${$file_ref};
   @COG2008     = split /\n/, ${$file_ref};
   undef ${$file_ref};
   foreach my $entry (@COG2008) {
      next unless ($entry =~ m/\w+/);
      my ($COG2008, $class, $name);
      'reset' =~ m/reset/;
      $entry  =~ m/^([^\t]+)\t([^\t]+)\t(.+)/;
      ($COG2008, $class, $name) = ($1, $2, $3);

      unless (defined $name) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing COG2008 COG to ID for entry $entry.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing COG2008 COG to ID for file $COG2008_file".
                        "\nError parsing COG2008 definition for entry $entry.\n\n";
         close ERRORLOG;
      }
      $name =~ s/\s+$//g;
      $name =~ s/\n//g;
      ${$args{COG2008_to_def}}->{$COG2008} = {(class      => $class,
                                               annotation => $name
                                             )};
   }
   undef @COG2008;

   #link accs to genomes, COG IDs and annotation
   #expected default location is in the COG2008 directory
   ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
   $COG2008_dir  = $1;
   $COG2008_file = 'domdat.txt';

   unless (-e $COG2008_dir.'/'.$COG2008_file) {
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $COG2008_dir,
                                                    -title      => 'Select COG2008domdat file',
                                                    -defaultextension => '',
                                                   );
      unless (defined $file) {return (0)};
      'reset'          =~ /reset/;
      $file            =~ m/(.+)\/(.+)$/;
      $COG2008_dir     = $1;
      $COG2008_file    = $2;
   }

   unless (-e $COG2008_dir.'/'.$COG2008_file) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No COG2008COGdomdat file selected, aborting.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error selecting COG2008 domdat file:".
                     "\ndir...$COG2008_dir ...file:...$COG2008_file ...\n\n";
      close ERRORLOG;
      return (0);
   }

   #read COGdomdat IDs and genome names
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $COG2008_dir,
                      filename     => $COG2008_file
                     );
   ${$file_ref} = "\n".${$file_ref};
   @COG2008     = split /\n/, ${$file_ref};
   undef ${$file_ref};
   foreach my $entry (@COG2008) {
      next unless ($entry =~ m/\w+/);
      my ($acc, $genome_abbr, $gb_id, $COG2008, $class);
      'reset' =~ m/reset/;
      $entry  =~ m/^([^\t]+)\t([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t(.+)/;
      ($acc, $genome_abbr, $gb_id, $COG2008, $class) = ($1, $2, $3, $4, $5);

      unless (defined $class) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing COG2008 acc to ID for entry $entry.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing COG2008 acc to ID for file $COG2008_file".
                        "\nError parsing COG2008 definition for entry $entry.\n\n";
         close ERRORLOG;
      }

      #if 'NULL', change to ''
      if ($COG2008 eq 'NULL') {$COG2008 = ''};
      if ($class   eq 'NULL') {$class   = ''};

      ${$args{COG2008_to_acc}}->{$acc} = {(genome_abbr     => $genome_abbr,
                                           genbank_id      => $gb_id,
                                           COG2008code     => $COG2008,
                                           COG2008class    => $class
                                          )};
   }
   undef @COG2008;
   return (1);
}

sub translate_COG2014 {
   my %args = @_;

   my ($file_ref, $COG2014_dir, $COG2014_file, @COG2014);

   #genomes2003-2014.tab contains a link between genome abbreviation and full genome name
   #expected default location is in the COG2014 directory
   ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
   $COG2014_dir  = $1;
   $COG2014_file = 'genomes2003-2014.tab';

   unless (-e $COG2014_dir.'/'.$COG2014_file) {
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $COG2014_dir,
                                                    -title      => 'Select COG2014 genomes file (genomes2003-2014.tab)',
                                                    -defaultextension => '',
                                                   );
      unless (defined $file) {return (0)};
      'reset'        =~ /reset/;
      $file          =~ m/(.+)\/(.+)$/;
      $COG2014_dir   = $1;
      $COG2014_file  = $2;
   }

   unless (-e $COG2014_dir.'/'.$COG2014_file) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No COG2014 genomes file selected, aborting.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error selecting COG2014 genomes file:".
                     "\ndir...$COG2014_dir ...file:...$COG2014_file ...\n\n";
      close ERRORLOG;
      return (0);
   }

   #read COG2014 genomes
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $COG2014_dir,
                      filename     => $COG2014_file
                     );
   ${$file_ref} = "\n".${$file_ref};
   @COG2014     = split /\n/, ${$file_ref};
   #remove header line
   shift @COG2014;

   undef ${$file_ref};
   foreach my $entry (@COG2014) {
      next unless ($entry =~ m/\w+/);
      next if ($entry =~ m/^\#/);
      my ($COG_abb, $COG_genome);
      'reset' =~ m/reset/;
      $entry  =~ m/^([^\t]+)\t[^\t]+\t([^\t]+)\n?\r?/;
      ($COG_abb, $COG_genome) = ($1, $2);

      unless (defined $COG_genome) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing COG2014 genomes for entry $entry.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing COG2008 genomes for file $COG2014_file".
                        "\nError parsing COG2008 genomes for entry $entry.\n\n";
         close ERRORLOG;
      }
      #remove uid extension
      $COG_genome =~ s/_uid\d+\s*$//;

      ${$args{COG2014genomes}}->{$COG_abb} = $COG_genome;
   }
   undef @COG2014;

   #link COG2014 ID to classes and annotation: cognames2003-2014.tab
   #expected default location is in the COG2014 directory
   ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
   $COG2014_dir  = $1;
   $COG2014_file = 'cognames2003-2014.tab';

   unless (-e $COG2014_dir.'/'.$COG2014_file) {
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $COG2014_dir,
                                                    -title      => 'Select COGdef file (cognames2003-2014.tab)',
                                                    -defaultextension => '',
                                                   );
      unless (defined $file) {return (0)};
      'reset'          =~ /reset/;
      $file            =~ m/(.+)\/(.+)$/;
      $COG2014_dir     = $1;
      $COG2014_file    = $2;
   }

   unless (-e $COG2014_dir.'/'.$COG2014_file) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No COG2014def file selected, aborting.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error selecting COG2014 def file:".
                     "\ndir...$COG2014_dir ...file:...$COG2014_file ...\n\n";
      close ERRORLOG;
      return (0);
   }

   #read COG2014 IDs and classes
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $COG2014_dir,
                      filename     => $COG2014_file
                     );
   ${$file_ref} = "\n".${$file_ref};
   @COG2014     = split /\n/, ${$file_ref};

   #remove header line
   shift @COG2014;
   undef ${$file_ref};

   foreach my $entry (@COG2014) {
      next if ($entry =~ m/^\#/);
      next unless ($entry =~ m/\w+/);
      my ($COG2014, $class, $name);
      'reset' =~ m/reset/;
      $entry  =~ m/^(COG\d+)\t([^\t]+)\t(.+)\n?\r?/;
      ($COG2014, $class, $name) = ($1, $2, $3);

      unless (defined $name) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing COG2014 COG to ID for entry $entry.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing COG2014 COG to ID for file $COG2014_file".
                        "\nError parsing COG2014 definition for entry $entry.\n\n";
         close ERRORLOG;
      }
      $name  =~ s/\s+$//gs;
      $name  =~ s/\n//gs;
      $class =~ s/\s+$//g;
      $class =~ s/\W//g;
      ${$args{COG2014_to_def}}->{$COG2014} = {(class      => $class,
                                               annotation => $name
                                             )};
   }
   undef @COG2014;

   #link accs to COG IDs
   #expected default location is in the COG2014 directory
   ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
   $COG2014_dir = $1;
   $COG2014_file  = 'cog2003-2014.csv';

   unless (-e $COG2014_dir.'/'.$COG2014_file) {
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $COG2014_dir,
                                                    -title      => 'Select COG2014 acc link file (cog2003-2014.csv)',
                                                    -defaultextension => '',
                                                   );
      unless (defined $file) {return (0)};
      'reset'        =~ /reset/;
      $file          =~ m/(.+)\/(.+)$/;
      $COG2014_dir     = $1;
      $COG2014_file    = $2;
   }

   unless (-e $COG2014_dir.'/'.$COG2014_file) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No COG2014 acc link file selected, aborting.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error selecting COG2014 domdat file:".
                     "\ndir...$COG2014_dir ...file:...$COG2014_file ...\n\n";
      close ERRORLOG;
      return (0);
   }

   #read COG acc link IDs
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $COG2014_dir,
                      filename     => $COG2014_file
                     );
   ${$file_ref} = "\n".${$file_ref};
   @COG2014 = split /\n/, ${$file_ref};
   undef ${$file_ref};
   foreach my $entry (@COG2014) {
      next unless ($entry =~ m/\w+/);
      my ($acc, $genome_abbr, $COG2014);
      'reset' =~ m/reset/;
      $entry =~ m/^(\d+)\,([^\,]+?)\,\d+\,\d+\,\d+\,\d+\,([^\,]+?)\,/;
      ($acc, $genome_abbr, $COG2014) = ($1, $2, $3);

      unless (defined $COG2014) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing COG2014 acc to ID for entry $entry.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing COG2014 acc to ID for file $COG2014_file".
                        "\nError parsing COG2014 definition for entry $entry.\n\n";
         close ERRORLOG;
      }
      $genome_abbr =~ s/_uid\d+\s*//;
      if ($COG2014 !~ m/COG/) {$COG2014 = ''};

      ${$args{COG2014_to_acc}}->{$acc} = {(genome_abbr => $genome_abbr,
                                           COG2014code => $COG2014
                                          )};
   }
   undef @COG2014;

   #link accs to protein IDs
   #expected default location is in the COG2014 directory
   ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
   $COG2014_dir = $1;
   $COG2014_file  = 'prot2003-2014.tab';

   unless (-e $COG2014_dir.'/'.$COG2014_file) {
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $COG2014_dir,
                                                    -title      => 'Select COG2014 ref-acc link file (prot2003-2014.tab)',
                                                    -defaultextension => '',
                                                   );
      unless (defined $file) {return (0)};
      'reset'        =~ /reset/;
      $file          =~ m/(.+)\/(.+)$/;
      $COG2014_dir     = $1;
      $COG2014_file    = $2;
   }

   unless (-e $COG2014_dir.'/'.$COG2014_file) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No COG2014 ref-acc link file selected, aborting.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error selecting COG2014 ref-acc file:".
                     "\ndir...$COG2014_dir ...file:...$COG2014_file ...\n\n";
      close ERRORLOG;
      return (0);
   }

   #read COG acc prot IDs
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $COG2014_dir,
                      filename     => $COG2014_file
                     );
   ${$file_ref} = "\n".${$file_ref};
   @COG2014 = split /\n/, ${$file_ref};
   undef ${$file_ref};
   foreach my $entry (@COG2014) {
      next unless ($entry =~ m/\w+/);
      my ($prot_id, $acc);
      'reset' =~ m/reset/;
      $entry =~ m/^(\d+)\t(\w+)\n?\r?/;
      ($acc, $prot_id) = ($1, $2);

      unless (defined $acc) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing COG2014 acc to protID for entry $entry.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing COG2014 acc to protID for file $COG2014_file".
                        "\nError parsing COG2014 definition for entry $entry.\n\n";
         close ERRORLOG;
      }

      ${$args{COG2014_to_refseq}}->{$prot_id} = $acc;
   }
   undef @COG2014;

   #read COG2014 assignment (fun2003-2014.tab)
   #expected default location is in the COG2014 directory
   ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
   $COG2014_dir = $1;
   $COG2014_file  = 'fun2003-2014.tab';

   unless (-e $COG2014_dir.'/'.$COG2014_file) {
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $COG2014_dir,
                                                    -title      => 'Select COG2014 fun file (fun2003-2014.tab)',
                                                    -defaultextension => '',
                                                   );
      unless (defined $file) {return (0)};
      'reset'        =~ /reset/;
      $file          =~ m/(.+)\/(.+)$/;
      $COG2014_dir     = $1;
      $COG2014_file    = $2;
   }

   unless (-e $COG2014_dir.'/'.$COG2014_file) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No COG2014 fun file selected, aborting.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error selecting COG2014 fun file:".
                     "\ndir...$COG2014_dir ...file:...$COG2014_file ...\n\n";
      close ERRORLOG;
      return (0);
   }

   #read COG function file
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $COG2014_dir,
                      filename     => $COG2014_file
                     );
   ${$file_ref} = "\n".${$file_ref};
   @COG2014 = split /\n/, ${$file_ref};

   #remove header line
   shift @COG2014;

   undef ${$file_ref};

   foreach my $entry (@COG2014) {
      my ($letter, $family);
      'reset' =~ m/reset/;
      $entry  =~ m/(\w)\s(.+)\n?\r?/;
      $letter = $1;
      $family = $2;
      ${$args{COGletter_to_family}}->{$letter} = $family;
   }
   undef @COG2014;
   return (1);
}

sub translate_POG2013 {
   my %args = @_;

   my ($file_ref, $POG2013_dir, $POG2013_file, @POG2013);

   #pogs.txt contains all information available for POGs
   #expected default location is in the COG2014 directory
   ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
   $POG2013_dir = $1;
   $POG2013_file  = 'pogs.txt';

   unless (-e $POG2013_dir.'/'.$POG2013_file) {
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $POG2013_dir,
                                                    -title      => 'Select POG2013 info file (pogs.txt)',
                                                    -defaultextension => '',
                                                   );
      unless (defined $file) {return (0)};
      'reset'        =~ /reset/;
      $file          =~ m/(.+)\/(.+)$/;
      $POG2013_dir   = $1;
      $POG2013_file  = $2;
   }

   unless (-e $POG2013_dir.'/'.$POG2013_file) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No POG2013 info file selected, aborting.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error selecting POG2013 info file:".
                     "\ndir...$POG2013_dir ...file:...$POG2013_file ...\n\n";
      close ERRORLOG;
      return (0);
   }

   #read POG2013 info
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $POG2013_dir,
                      filename     => $POG2013_file
                     );
   ${$file_ref} = "\n".${$file_ref};
   @POG2013 = split /\n/, ${$file_ref};

   undef ${$file_ref};
   foreach my $entry (@POG2013) {
      next unless ($entry =~ m/\w+/);
      my ($POG_class, $POG_ntacc, $POG_phage, $POG_gi, $POG_ref, $POG_annotation, $POG_phylogeny);
      my @tmp = split /\:/, $entry;

      $POG_class      = $tmp[0];
      $POG_ntacc      = $tmp[1];
      $POG_phage      = $tmp[2];
      $POG_gi         = $tmp[3];
      if (defined $tmp[5] && $tmp[5] =~ m/\w+/) {
         $POG_annotation = $tmp[5];
      } else {
         $POG_annotation = '';
         next; #ignore hits that can't be traced back.
      }
      if (defined $tmp[-1] && $tmp[-1] =~ m/\w+/) {
         $POG_phylogeny  = $tmp[-1];
      } else {
         $POG_phylogeny  = '';
      }

      unless (defined $POG_gi) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing POG2013 info for entry $entry.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing POG2013 info for file $POG2013_file".
                        "\nError parsing POG2013 info for entry $entry.\n\n";
         close ERRORLOG;
      }
      #clean up annotation
      $POG_ref = $POG_annotation;
      if ($POG_annotation =~ m/\|/) {
         $POG_annotation  =~ s/^\>[^ ]+\s//;
      } else {
         $POG_annotation  = '';
      }
      $POG_phylogeny =~ s/^\[\d+\]//;
      $POG_phylogeny =~ s/\"//g;

      #get ref acc
      $POG_ref        =~ s/^\>.*?ref\|([^\|]+)\|.*/$1/;
      $POG_annotation =~ s/\n//gs;

      ${$args{POG2013_to_gi}}->{$POG_ref} = {(class      => $POG_class,
                                              ntacc      => $POG_ntacc,
                                              phage      => $POG_phage,
                                              gi         => $POG_gi,
                                              annotation => $POG_annotation,
                                              phylogeny  => $POG_phylogeny
                                            )};

   }
   undef @POG2013;

   return (1);
}

sub translate_COG2003 {
   my %args = @_;
   my ($file_ref, @COGs, @ORGs, %ORGs, $COG_whog_dir, $COG_whog, $COG_fun_dir, $COG_fun, $COG_org_dir, $COG_org);

   #test for files
   ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
   $COG_whog_dir  = ${$args{ini_ref}}{COG_db_path};
   $COG_fun_dir   = ${$args{ini_ref}}{COG_db_path};
   $COG_org_dir   = ${$args{ini_ref}}{COG_db_path};
   $COG_whog      = 'whog';
   $COG_fun       = 'fun.txt';
   $COG_org       = 'org.txt';

   #try and find COG2003 location; assume initially that all parsing files are in the same directory
   if (-e ${$args{ini_ref}}{COG_db_path}.'/'.$COG_whog) {
      $COG_whog_dir = ${$args{ini_ref}}{COG_db_path};
      $COG_fun_dir  = ${$args{ini_ref}}{COG_db_path};
      $COG_org_dir  = ${$args{ini_ref}}{COG_db_path};
   } elsif (-e ${$args{ini_ref}}{COG_db_path}.'/COG2003/'.$COG_whog) {
      $COG_whog_dir = ${$args{ini_ref}}{COG_db_path}.'/COG2003';
      $COG_fun_dir  = ${$args{ini_ref}}{COG_db_path}.'/COG2003';
      $COG_org_dir  = ${$args{ini_ref}}{COG_db_path}.'/COG2003';
   } else {
      #create starting position
      $COG_whog_dir  = ${$args{ini_ref}}{COG_db_path};
      $COG_fun_dir   = ${$args{ini_ref}}{COG_db_path};
      $COG_org_dir   = ${$args{ini_ref}}{COG_db_path};
      #ask for whog file
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $COG_whog_dir,
                                                    -title      => 'Select COG whog file',
                                                    -defaultextension => '',
                                                   );
      unless (defined $file) {return (0)};
      'reset'        =~ /reset/;
      $file          =~ m/(.+)\/(.+)$/;
      $COG_whog_dir  = $1;
      $COG_whog      = $2;

      #ask for fun file
      $file = ${$args{main_window}}->getOpenFile(-initialdir => $COG_fun_dir,
                                                 -title      => 'Select COG fun file',
                                                 -defaultextension => '',
                                                );
      unless (defined $file) {return (0)};
      'reset'   =~ /reset/;
      $file     =~ m/(.+)\/(.+)$/;
      $COG_fun_dir  = $1;
      $COG_fun = $2;

      #ask for org file
      $file = ${$args{main_window}}->getOpenFile(-initialdir => $COG_org_dir,
                                                 -title      => 'Select COG org.txt file',
                                                 -defaultextension => '',
                                                );
      unless (defined $file) {return (0)};
      'reset'        =~ /reset/;
      $file          =~ m/(.+)\/(.+)$/;
      $COG_org_dir  = $1;
      $COG_org      = $2;
   }

   #read COG codes
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $COG_whog_dir,
                      filename     => $COG_whog
                     );
   ${$file_ref} = "\n".${$file_ref};
   (@COGs) = (${$file_ref} =~ m/\n\r?(.*?)_______\n\r?/gs);
   undef ${$file_ref}; #free up memory
   #read organism file
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $COG_org_dir,
                      filename     => $COG_org
                     );
   ${$file_ref} = "\n".${$file_ref};
   @ORGs = split /\n/, ${$file_ref};
   foreach my $org (@ORGs) {
      'reset' =~ m/reset/;
      $org =~ m/^(\S+)\s+\d+\s+\S+\s+(.+)/;
      $ORGs{$1} = $2;
   }
   undef ${$file_ref};

   undef ${$file_ref}; #free up memory
   foreach my $entry (@COGs) {
      my ($header, $number, $letter, $COG,  @code);
      'reset' =~ m/reset/;
      $entry  =~ m/\[(\w+)\]\s+(COG\d+)\s+([^\n\r]*?)\n\r?/;
      $letter = $1;
      $number = $2;
      $header = $3;
      $entry =~ s/^.*?\n\r?//s;
      $header =~ s/\n//gs;

      #@code = ($entry =~ m/(\S+)\s+/gs);
      @code = split /\n/, $entry;
      #iterate through all COG codes
      my $phylo   = '';
      foreach my $code (@code) {
         next if ($code !~ /\w+/);
         my @entries = ();
         'reset'     =~ m/reset/;
         #assume multi-line entry; carries over $phylo from original entry line
         if ($code       =~ m/^\s+(\w\w\w)\:/) {
            $code       =~ m/^\s+(\w\w\w)\:/;
            $phylo      = $1;
         }
         unless (defined $phylo) {
            print "\nError parsing COG2003 phylogeny information in whog: $code\n";
         }
         $code       =~ s/^\s+\w\w\w\:\s+//;

         @entries = split /\s+/, $code;
         foreach my $entry (@entries) {
            next unless ($entry =~ m/\w+/);
            if (defined ${$args{COGcode_to_number}}->{$entry}) {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Duplicate COG code assignment for $entry.",
                                                             -buttons => ['OK'],
                                                             -bitmap  => 'info');
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error parsing COG codes for file $COG_whog".
                              "\nDuplicate COG code assignment for $entry.\n\n";
               close ERRORLOG;
            }

            ${$args{COGcode_to_number}}->{$entry} = $number;
            ${$args{COGcode_to_header}}->{$entry} = $header;
            ${$args{COGcode_to_letter}}->{$entry} = $letter;
            ${$args{COGcode_to_phylo}}->{$entry}  = $ORGs{$phylo};
         }
      }
      undef @code;
   }
   undef (@COGs);
   undef %ORGs;

   #read COG assignment
   $file_ref = &slurp(main_window  => $args{main_window},
                      auto_ini_ref => $args{auto_ini_ref},
                      directory    => $COG_fun_dir,
                      filename     => $COG_fun
                     );
   @COGs = (${$file_ref} =~ m/\n\r?( +[^\n\r]+)/gs);
   foreach my $entry (@COGs) {
      my ($letter, $family);
      'reset' =~ m/reset/;
      $entry  =~ m/\[(\w+)\]\s+(.+)/;
      $letter = $1;
      $family = $2;
      ${$args{COGletter_to_family}}->{$letter} = $family;
   }
   undef @COGs;
   return (1);
}


1;

