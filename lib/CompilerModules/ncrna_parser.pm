#!/opt/ActivePerl-5.8/bin/perl

#RFam result parser
#input arguments: main_window, progress_bar, gb_file, gene_model

package CompilerModules::ncrna_parser;
use strict;
use initialise::read_me  qw(:DEFAULT);
use Basics::progress_bar qw(:DEFAULT);
use File::Copy;

use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&ncrna_parser &translate_rfam);
use vars qw();


#local variables
my (%args, $value);

format ENTRY =
~~                   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$value
.

sub ncrna_parser {
   my %args = @_;
   my (@res_files, $max_count, $count, $file, $file_ref, $local_threshold, $domains,
       $global_counter, $domain_count);

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling non-coding RNA results",
                   label        => ''
                  );
   &show_pbar_2;
   #define local score threshold
   #$file_ref = slurp(main_window => $args{main_window},
   #                  directory   => ${$args{ini_ref}}{input_files},
   #                  filename    => $args{filename}
   #                 );
   #if ($file_ref eq '0') {return (0)};
   #remove header
   #$$file_ref =~ s/^\>[^\n]*?\n//;
   #undef $file_ref;

   #$local_threshold = log(length(${$file_ref}) * 2);
   #$local_threshold += (${$args{auto_ini_ref}}{ncrna_threshold} / 100) * $local_threshold; # add to baseline threshold to make hits more conservative
   #if ($local_threshold =~ /\d+\.(\d){2,}/) {
   #    $local_threshold=sprintf("%.".(2)."f", $local_threshold);
   #}

   #get all res files for input file
   opendir SEQINPUT, ${$args{ini_ref}}{rfam_results};
   @res_files = grep /$args{filename}_chunk_\d+\.Rfam$/, readdir(SEQINPUT);
   closedir SEQINPUT;
   $max_count = $#res_files + 2;
   $global_counter = 1;
   $count = 1;
   $domain_count = 1;

   #parse through ncRNA results
   foreach my $res_file (@res_files) {
      my ($acc, $offset, @blocks);

      &update_pbar_2(title        => "Compiling non-coding RNA results",
                     label        => "Parsing result files $count of $max_count",
                     progress     => ($count / $max_count) * 100,
                    );
      $count++;

      #get offset
      'reset' =~ m/reset/;
      $res_file =~ m/$args{filename}_chunk_(\d+)\.Rfam$/;
      $offset = $1;
      unless (defined $offset) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse ncRNA filename $res_file for offset value.",
                                                       -bitmap  => 'error',
                                                       -buttons => ['ok']);
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing ncRNAs".
                        "\nCould not parse ncRNA filename $res_file for offset value.\n\n";
         close ERRORLOG;
         next;
      };

      $file_ref = slurp(main_window => $args{main_window},
                        directory   => ${$args{ini_ref}}{rfam_results},
                        filename    => $res_file
                       );
      if ($file_ref eq '0') {
         &hide_pbar_2;
         return (0);
      }

      #capture individual RFam models
      @blocks = split /\n\/\//, $$file_ref;

      #iterate through model
      foreach my $block (@blocks) {
         #skip no hits
         next if ($block =~ m/No hits detected that satisfy reporting thresholds/);
         next if ($block =~ m/\n?\[ok\]/is);
         #skip rRNA and tRNA is selected
         next if (${$args{auto_ini_ref}}{ncrna_ignore_rrna} == 1 && $block =~ /Description\:\s+.+?S\s+ribosomal RNA/s);
         next if (${$args{auto_ini_ref}}{ncrna_ignore_trna} == 1 && $block =~ /Description\:\s+tRNA/s );

         #get acc
         'reset' =~ m/reset/;
         $block =~ m/Accession\:\s+(RF\d+)/s;
         my $acc = $1;
         unless (defined $acc) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse ncRNA $res_file for acc value.",
                                                          -bitmap  => 'error',
                                                          -buttons => ['ok']);
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing ncRNAs".
                           "\nCould not parse ncRNA $res_file for acc value.\nBlock:\n$block\n\n";
            close ERRORLOG;
            next;
         };
         #get individual entries
         $block =~ s/^.+?\nHit scores\:(.+)\nHit alignments.+/$1/s;
         my @entries = ($block =~ m/\s+(\(\d+\)[^\n]+?)\n/gs);

         #iterate through individual ncRNA entries
         foreach my $entry (@entries) {
            $entry =~ m/^\s*\(\d+\)\s*[\?\!]?\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+([\+\-])\s/;
            my ($evalue, $score, $start, $end, $orientation) = ($1, $2, $3, $4, $5);
            unless (defined $orientation) {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Could not parse ncRNA $res_file for domain values.",
                                                             -bitmap  => 'error',
                                                             -buttons => ['ok']);
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error parsing ncRNAs".
                              "\nCould not parse ncRNA $res_file for domain values.\nBlock:\n$block\nEntry: $entry\n\n";
               close ERRORLOG;
               next;
            };

            #evalue below threshold?
            next if ($evalue > ${$args{auto_ini_ref}}{ncrna_threshold});

            #put domains into hash
            if ($orientation eq '+') {
               #get absolute boundaries
               my ($abs_left, $abs_right);
               $abs_left  = $offset + $start;
               $abs_right = $offset + $end;
               $domains->{$acc}->{$domain_count} = { 'orientation' => 'sense',
                                                     'score'       => $score,
                                                     'evalue'      => $evalue,
                                                     'left_bd'     => $abs_left,
                                                     'right_bd'    => $abs_right,
                                                   };
            } else {
               #get absolute boundaries
               my ($abs_left, $abs_right);
               $abs_left  = $offset + $end;
               $abs_right = $offset + $start;
               $domains->{$acc}->{$domain_count} = { 'orientation' => 'antisense',
                                                     'score'       => $score,
                                                     'evalue'      => $evalue,
                                                     'left_bd'     => $abs_left,
                                                     'right_bd'    => $abs_right,
                                                   };
            }
            $global_counter++;
            $domain_count++;
         }
      }

      #remove any overlaps and keep hits with highest score for any given acc
      #&remove_overlaps(domains        => \$domains,
      #                 acc            => $acc,
      #                 global_counter => \$global_counter
      #                );

   }

   #remove Rfams within ORFs if selected
   if (${$args{auto_ini_ref}}{ncrna_IG_only} == 1) {
      &update_pbar_2(title        => "Non-coding RNA analysis",
                     label        => "Removing Rfams within ORFs",
                     progress     => 1,
                    );
      &remove_ORF_Rfams(domains          => \$domains,
                        global_counter   => \$global_counter,
                        feature_list_ref => $args{feature_list_ref}
                        );
   }

   #remove overlapping non-identical RFams if selected
   if (${$args{auto_ini_ref}}{ncrna_overlap} == 0) {
      &update_pbar_2(title        => "Non-coding RNA analysis",
                     label        => "Removing global overlaps",
                     progress     => 1,
                    );
      &remove_global_overlaps(domains        => \$domains,
                              global_counter => \$global_counter
                              );
   }

   #generate new entries
   $count = 1;
   foreach my $acc ( %{ $domains } ) {
      foreach my $counter ( %{ $domains->{$acc} } ) {
         #skip hash references
         next if ($counter =~ /\D+/);

         #update status bar
         if (($count % 10) == 0) {
            &update_pbar_2(title        => "Compiling non-coding RNA results",
                           label        => "Compiling non-coding RNA results for $args{filename}",
                           progress     => ($count / $global_counter) * 100,
                          );
         }
         $count++;

         my ($boundary_ncrna, $label, $ncrna_function, $score, $colour, $note,
             $feature_ncRNA, $key_ncRNA, $ncrna_product, $product_label);

         #make sure entry is defined
         if (defined $domains->{$acc}->{$counter}->{'orientation'}) {
            my $key;
            #increase ID
            $args{counter}++;

            #set key
            if (${$args{rfam_model}}->{$acc}->{'ID'} =~ /_rRNA/ || ${$args{rfam_model}}->{$acc}->{'DE'} =~ /S rRNA\)/) {
               $key = 'rRNA            ';
            } elsif (${$args{rfam_model}}->{$acc}->{'ID'} =~ /tRNA/ || ${$args{rfam_model}}->{$acc}->{'DE'} =~ /tRNA/) {
               $key = 'tRNA            ';
            } else {
               $key = 'ncRNA           ';
            }

            #create boundary entry
            if ($domains->{$acc}->{$counter}->{'orientation'} eq 'sense') {
               $boundary_ncrna  = '     '.$key.$domains->{$acc}->{$counter}->{'left_bd'}.'..'.$domains->{$acc}->{$counter}->{'right_bd'}."\n";
            } elsif ($domains->{$acc}->{$counter}->{'orientation'} eq 'antisense') {
               $boundary_ncrna  = '     '.$key.'complement('.$domains->{$acc}->{$counter}->{'left_bd'}.'..'.$domains->{$acc}->{$counter}->{'right_bd'}."\)\n";
            } else {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Could not set boundaries for ACC $acc counter $counter.",
                                                             -bitmap  => 'error',
                                                             -buttons => ['ok']);
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error parsing ncRNAs".
                              "\nCould not set boundaries for ACC $acc counter $counter.\n\n";
               close ERRORLOG;
               &hide_pbar_2;
               return (0);
            };

            #create ncRNA label
            my $label;
            if (${$args{rfam_model}}->{$acc}->{'ID'} =~ /\S+/) {
               $value         = '/label='.${$args{rfam_model}}->{$acc}->{'ID'}."\n";
               $product_label = ${$args{rfam_model}}->{$acc}->{'ID'};
            }
            #elsif (${$args{rfam_model}}->{$acc}->{'DE'} =~ /\S+/) {
            #   $value         = '/label='.${$args{rfam_model}}->{$acc}->{'DE'}."\n";
            #   $product_label = ${$args{rfam_model}}->{$acc}->{'DE'};
            #}
            else {
               $value         = '/label=non-coding RNA structure';
               $product_label = 'non-coding RNA structure';
            }
            open  ENTRY, '>', \$label;
            write ENTRY;
            close ENTRY;

            #add product
            $value = '/product="'.$acc.', Function: '.$product_label."\"\n";
            open  ENTRY, '>', \$ncrna_product;
            write ENTRY;
            close ENTRY;

            #add verbose description if selected
            if (${$args{auto_ini_ref}}{ncrna_verbose} == 1) {
               if (${$args{rfam_model}}->{$acc}->{'DE'} eq '\N') {
                  $value = '/function="'.${$args{rfam_model}}->{$acc}->{'OR'}.'"'."\n";
               } else {
                  ${$args{rfam_model}}->{$acc}->{'DE'} =~ s/\"/\'/gs; #subsitute quotation marks
                  $value = '/function="'.${$args{rfam_model}}->{$acc}->{'DE'}.'"'."\n";
               }
            }
            open  ENTRY, '>', \$ncrna_function;
            write ENTRY;
            close ENTRY;

            #create score
            $value = '/score='.$domains->{$acc}->{$counter}->{'score'}.'; E-value='.$domains->{$acc}->{$counter}->{'evalue'}."\n";
            open  ENTRY, '>', \$score;
            write ENTRY;
            close ENTRY;

            #set colour
            $value = '/colour=7'."\n";
            open  ENTRY, '>', \$colour;
            write ENTRY;
            close ENTRY;

            #combine feature
            $feature_ncRNA = $boundary_ncrna.
                             $label.
                             $score.
                             $ncrna_product.
                             $ncrna_function.
                             $colour;
            $key_ncRNA     = $domains->{$acc}->{$counter}->{'left_bd'}.'_'.$domains->{$acc}->{$counter}->{'right_bd'}.'_ncRNA_'.$args{counter}; #create  uniqe ID instead of ORFnumber

            ${$args{genbank_ref}}{$key_ncRNA} = $feature_ncRNA;
            push (@{$args{feature_list_ref}}, $key_ncRNA);
         }
      }
   }
   &hide_pbar_2;
   undef $domains;
   return;
}

sub remove_ORF_Rfams {
   my %args = @_;
   my $orf_model_hash;
   #create left_bd-right_bd has from ORF model
   foreach my $entry (@{$args{feature_list_ref}}) {
      next unless ($entry =~ m/_CDS_/);
      #get rigth and left bd
      "reset" =~ m/reset/;
      $entry =~ m/^(\d+)_(\d+)_.+_(\d+)$/;
      my ($left_bd, $right_bd, $counter) = ($1, $2, $3);
      unless (defined $counter) {
         next;
      }
      $orf_model_hash->{$counter} = {left_bd  => $left_bd,
                                     right_bd => $right_bd
                                    };
   }

   #iterate each domain and check if in ORF
   foreach my $first_acc (keys %{ ${$args{domains}} }) {
      next if ($first_acc =~ /\W+/);
      foreach my $first_counter (keys %{ ${$args{domains}} -> { $first_acc } }) {
         next if ($first_counter =~ /\D+/);
         next if (${$args{domains}} -> { $first_acc } -> { $first_counter } -> { 'score' } !~ /\w+/); #ignore deleted
         foreach my $orf_counter (keys %{$orf_model_hash}) {
            if (${$args{domains}} -> { $first_acc } -> { $first_counter } -> { 'left_bd' }  >= $orf_model_hash->{ $orf_counter}->{ 'left_bd'} &&
                ${$args{domains}} -> { $first_acc } -> { $first_counter } -> { 'right_bd' } <= $orf_model_hash->{ $orf_counter}->{ 'right_bd'}) {
               delete ${$args{domains}} -> { $first_acc } -> { $first_counter };
               last;
            }
         }
      }
   }
}

sub remove_global_overlaps {
   my %args = @_;
   my $iterate;
   foreach my $first_acc (keys %{ ${$args{domains}} }) {
      next if ($first_acc =~ /\W+/);
      foreach my $first_counter (keys %{ ${$args{domains}} -> { $first_acc } }) {
         next if ($first_counter =~ /\D+/);
         next if (${$args{domains}} -> { $first_acc } -> { $first_counter } -> { 'score' } !~ /\w+/); #ignore deleted
         my $first_score    = ${$args{domains}} -> { $first_acc } -> { $first_counter } -> { 'score' };
         my $first_left_bd  = ${$args{domains}} -> { $first_acc } -> { $first_counter } -> { 'left_bd' };
         my $first_right_bd = ${$args{domains}} -> { $first_acc } -> { $first_counter } -> { 'right_bd' };
         foreach my $second_acc (keys %{ ${$args{domains}} }) {
            next if ($second_acc =~ /\W+/);
            foreach my $second_counter (keys %{ ${$args{domains}} -> { $second_acc } }) {
               next if ($second_counter =~ /\D+/);
               next if ($first_counter == $second_counter && $first_acc eq $second_acc); #ignore self
               next if (${$args{domains}} -> { $second_acc } -> { $second_counter } -> { 'score' } !~ /\w+/); #ignore deleted
               my $second_score    = ${$args{domains}} -> { $second_acc } -> { $second_counter } -> { 'score' };
               my $second_left_bd  = ${$args{domains}} -> { $second_acc } -> { $second_counter } -> { 'left_bd' };
               my $second_right_bd = ${$args{domains}} -> { $second_acc } -> { $second_counter } -> { 'right_bd' };

               if (($second_left_bd  >= $first_left_bd  && $second_left_bd  <= $first_right_bd) ||
                   ($second_right_bd >= $first_left_bd  && $second_right_bd <= $first_right_bd) ||
                   ($second_left_bd  <= $first_left_bd  && $second_right_bd >= $first_right_bd)) {
                  $iterate = 1;
                  ${$args{global_counter}}--;
                  if ($first_score >= $second_score) {
                     delete ${$args{domains}} -> { $second_acc } -> { $second_counter };
                  } else {
                     delete ${$args{domains}} -> { $first_acc } -> { $first_counter };
                  }
               }
            }
         }
      }
   }
   if ($iterate == 1) {
      #restart iteration
      &remove_global_overlaps(domains        => $args{domains},
                              global_counter => $args{global_counter}
                             );
   }
   return;
}

sub remove_overlaps {
   my %args = @_;
   my $iterate;

   foreach my $first_counter (keys %{ ${$args{domains}} -> { $args{acc} } }) {
      next if (${$args{domains}} -> { $args{acc} } -> { $first_counter } -> { 'score' } !~ /\w+/); #ignore deleted
      my $first_score    = ${$args{domains}} -> { $args{acc} } -> { $first_counter } -> { 'score' };
      my $first_left_bd  = ${$args{domains}} -> { $args{acc} } -> { $first_counter } -> { 'left_bd' };
      my $first_right_bd = ${$args{domains}} -> { $args{acc} } -> { $first_counter } -> { 'right_bd' };
      foreach my $second_counter (keys %{ ${$args{domains}} -> { $args{acc} } }) {
         next if ($first_counter == $second_counter); #ignore self
         next if (${$args{domains}} -> { $args{acc} } -> { $second_counter } -> { 'score' } !~ /\w+/); #ignore deleted
         my $second_score    = ${$args{domains}} -> { $args{acc} } -> { $second_counter } -> { 'score' };
         my $second_left_bd  = ${$args{domains}} -> { $args{acc} } -> { $second_counter } -> { 'left_bd' };
         my $second_right_bd = ${$args{domains}} -> { $args{acc} } -> { $second_counter } -> { 'right_bd' };

         if (($second_left_bd  >= $first_left_bd && $second_left_bd  <= $first_right_bd) ||
             ($second_right_bd >= $first_left_bd && $second_right_bd <= $first_right_bd) ||
             ($second_left_bd  <= $first_left_bd && $second_right_bd >= $first_right_bd)) {
            $iterate = 1;
            ${$args{global_counter}}--;
            if ($first_score >= $second_score) {
               delete ${$args{domains}} -> { $args{acc} } -> { $second_counter };
               last;
            } else {
               delete ${$args{domains}} -> { $args{acc} } -> { $first_counter };
               last;
            }
         }
      }
   }
   if ($iterate == 1) {
      #restart iteration
      &remove_overlaps(domains        => $args{domains},
                       acc            => $args{acc},
                       global_counter => $args{global_counter}
                      );
   }
   return;
}


sub translate_rfam {
   my %args = (retry_counter   => 0,
               RFam_descriptor => '',
               @_
              );
   my ($seen, $file, $block, %descriptor, $rfam_model, $acc);

   #Rfam descriptor file to open
   if ($args{RFam_descriptor} =~ m/\w+/) {
      $file = $args{RFam_descriptor};
   } else {
      $file = ${$args{ini_ref}}{rfam_db_path}.'/rfam.txt';
   }

   #test if exists
   unless (-e $file) {
      #try to uncompress if archive is present
      if (-e ${$args{ini_ref}}{rfam_db_path}.'/rfam.txt.gz') {
         ${$args{progress_bar}}->configure(-label=>" Uncompressing rfam.txt ");
         ${$args{main_window}}->update;
         `gzip -d $file\.gz`;
      }
      #try gzip if tar fails
      unless (-e $file) {
         #try file in  archive folder
         copy (${$args{auto_ini_ref}}{work_dir}.'/lib/Archives/rfam.txt.gz', ${$args{ini_ref}}{rfam_db_path}.'/rfam.txt.gz');
         `gzip -d $file\.gz`;
      }
      #still not there? point to archive an re-try (only once though)
      unless (-e $file) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not find Rfam descriptor archive \n$file",
                                           -type    => 'OK',
                                           -icon    => 'info');
         my $types = [ ['gz',  ['.gz']],
                       ['All', ['.*' ]]
                     ];
         my $RFam_descriptor = ${$args{main_window}}->getOpenFile(-initialdir => ${$args{ini_ref}}{work_dir},
                                                                  -title      => "Select RFam descriptor archive $file",
                                                                  -filetypes  => $types
                                                                 );
         unless (defined $RFam_descriptor) {
            return (0);
         }
         if (defined $RFam_descriptor && $args{retry_counter} <= 1) {
            copy ($RFam_descriptor, ${$args{ini_ref}}{rfam_db_path}.'/rfam.txt.gz');
            &translate_rfam(main_window            => $args{main_window},
                            progress_bar           => $args{progress_bar},
                            auto_ini_ref           => $args{auto_ini_ref},
                            ini_ref                => $args{ini_ref},
                            rfam_model             => \$rfam_model,
                            retry_counter          => $args{retry_counter}++,
                            Rfam_descriptor        => ${$args{ini_ref}}{rfam_db_path}.'/rfam.txt'
                           );
         } elsif ($args{retry_counter} > 1) {
            return (0);
         }
      }
   }

   #read RFam descriptions and  keys
   open READ, '<'.$file or do {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Could not open Rfam.full in directory ${$args{ini_ref}}{rfam_db_path}.",
                                                    -bitmap  => 'error',
                                                    -buttons => ['ok']);
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error parsing ncRNAs".
                     "\nCould not open Rfam.full in directory ${$args{ini_ref}}{rfam_db_path}.\n\n";
      close ERRORLOG;
      return (0);
   };

   #define entry blocks
   while (<READ>) {
      next unless ($_ =~ /^\d/);
      my @fields = split /\t/, $_;
      #assign fields
      ${$args{rfam_model}}->{$fields[2]} = {'AC' => $fields[2],
                                            'ID' => $fields[3],
                                            'OR' => $fields[4],
                                            'DE' => $fields[11],
                                           };
   }
   close READ;
   return (1);
}



1;

