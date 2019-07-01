#!/opt/ActivePerl-5.8/bin/perl

#PFam result parser
#input arguments: main_window, progress_bar, gb_file, gene_model

package CompilerModules::pfam_parser;
use strict;
use initialise::read_me  qw(:DEFAULT);
use Basics::progress_bar qw(:DEFAULT);

use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&pfam_parser &pfam3_parser &translate_pfam &translate_interpro);
use vars qw();


#local variables
my (%args, %seen, $value);

format ENTRY =
~~                   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$value
.

sub pfam_parser {
   my %args = @_;

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling PFam results",
                   label        => ''
                  );
   &show_pbar_2;

   #parse through PFam results
   foreach (my $i = 1; $i < $args{counter}; $i++) {
      my ($filename, $id, $left_bd, $right_bd, $orientation, @pfam_db,
          $pfam_to_domain, $model_to_domain, @pfam_key, $counter);

      #get filename, left_bd, right_bd, orientation
      'reset' =~ m/reset/;
      $args{orf_to_id}->{$i} =~ m/^(.*?)___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
      $filename    = $1;
      $id          = $2;
      $left_bd     = $3;
      $right_bd    = $4;
      $orientation = $5;
      unless (defined $orientation) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not parse from entry $args{orf_to_id}->{$i}\.",
                                           -icon    => 'error',
                                           -type    => 'OK');
         &hide_pbar_2;
         return (0);
      }

      #update status bar
      if (($i % 20) == 0) {
         &update_pbar_2(title        => "Compiling PFam results",
                        label        => "Compiling PFam results for $filename",
                        progress     => ($i / $args{counter}) * 100,
                       );
      }

      #iterate through all selected databases
      @pfam_db = split / ; /,${$args{auto_ini_ref}}{Pfam_db};

      $counter = 1;
      foreach my $pfam_db (@pfam_db) {
         my ($family, $domain, @family);

         #test if file exists
         next unless (-e ${$args{ini_ref}}{pfam_results}.'/'.$filename.'_'.$pfam_db.'_'.$id);

         #read input file
         (my $file_ref) = slurp(main_window => $args{main_window},
                                directory   => ${$args{ini_ref}}{pfam_results},
                                filename    => $filename.'_'.$pfam_db.'_'.$id
                               );
         if ($file_ref eq '0') {return (0)};

         #skip if no PFam results
         next if (${$file_ref} =~ /\[no hits above thresholds\]/s);

         #grab family and domain box
         'reset'      =~ /reset/;
         ${$file_ref} =~ /Scores for sequence family classification(.*?)Parsed for domains(.*?)Alignments of top-scoring domains/s;
         $family      = $1;
         $domain      = $2;

         unless (defined $domain) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not parse PFam results from entry $filename\_$pfam_db\_$id.",
                                              -icon    => 'error',
                                              -type    => 'OK');
            &hide_pbar_2;
            return (0);
         }

         #clean up
         $family =~ s/.*?------- ---\n\r?//s;
         $domain =~ s/.*?-------\n\r?//s;

         #save domain block away
         $pfam_to_domain->{$pfam_db} = $domain;

         #create arrays
         @family = split /\n\r?/, $family;

         foreach my $family (@family) {
            next if ($family !~ /\w+/);
            my ($model, $score, $evalue, $key);
            'reset'      =~ m/reset/;
            $family      =~ m/(\S+)\s+.{41}\s*(\S+)\s+(\S+)\s+\d/;
            $model       = $1;
            $score       = $2;
            $evalue      = $3;

            unless (defined $evalue) {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Could not parse PFam results from entry $family in file $filename\_$pfam_db\_$id.",
                                                 -icon    => 'error',
                                                 -type    => 'OK');
               &hide_pbar_2;
               return (0);
            }

            #check evalue
            if ($evalue =~ /^e/) {$evalue = '1'.$evalue};

            $key = $counter.'___'.$pfam_db.'___'.$model.'___'.$evalue;

            push (@pfam_key, $key);
            undef $key;
            $counter++;
         }
      }

      #sort Pfam array and then get the No. of best Pfam family models indicated in the current selection

      my @sorted =
         map  $_->[0] =>
         sort { $a->[1] <=> $b->[1] }
         map  [ $_, m/^\d+___\S+___\S+___(.+)$/ ]
         => @pfam_key;
      undef @pfam_key;

      @pfam_key = splice(@sorted,0,${$args{auto_ini_ref}}{Pfam_family_max});
      undef @sorted;

      #now process top models and create features
      foreach my $best_model (@pfam_key) {
         my ($db, $model, $evalue, @domains);
         'reset' =~ /reset/;
         $best_model =~ m/\d+___(.+)___(.+)___(.+)/;
         $db     = $1;
         $model  = $2;
         $evalue = $3;

         unless (defined $evalue) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not parse PFam key $best_model in file $filename\_$db\_$id.",
                                              -icon    => 'error',
                                              -type    => 'OK');
            &hide_pbar_2;
            return (0);
         }

         #skip if above threshold
         next if ($evalue > ${$args{ini_ref}}{PFam_cut});

         #grab domains from domain block
         @domains = ($pfam_to_domain->{$db} =~ m/($model[^\n\r]+)\n\r?/gs);

         #iterate individual domains
         foreach my $domain (@domains) {
            my ($pfam_lb, $pfam_rb, $pfam_lba, $pfam_rba, $boundary_pfam, $score, $domain_evalue,
                $label, $pfam_match, $interpro, $note, $product,
                $feature_pfam, $key_pfam);
            'reset' =~ /reset/;
            $domain =~ m/\S+\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+\d+\s+\d+\s+\S+\s+(\S+)\s+(\S+)/;
            $pfam_lb       = $1;
            $pfam_rb       = $2;
            $score         = $3;
            $domain_evalue = $4;

            unless (defined $domain_evalue) {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Could not parse domain entry $domain in file $filename\_$db\_$id.",
                                                 -icon    => 'error',
                                                 -type    => 'OK');
                  &hide_pbar_2;
                  return (0);
            }

            #test for domain threshold if selected
            if (${$args{auto_ini_ref}}{Pfam_restrict_domains} == 1) {
               next if ($domain_evalue > ${$args{ini_ref}}{PFam_cut});
            }

            #create boundary entry
            if ($orientation eq 'sense') {
               $pfam_rba = $left_bd + ($pfam_rb * 3) - 1;
               $pfam_lba = $left_bd + ($pfam_lb * 3) - 3;

               #modify boundaries is necessary
               if ($pfam_lba == $left_bd) {
                  $pfam_lba += 3;  #add one codon to avoid equal to start position
               }
               if ($pfam_rba == $right_bd) {
                  $pfam_rba -= 3;  #remove one codon to avoid equal to stop position
               }

               $boundary_pfam  = '     PFAM_match      '.$pfam_lba.'..'.$pfam_rba."\n";
            } elsif ($orientation eq 'antisense') {
               my $length = $pfam_rb - $pfam_lb;
               $pfam_rba = $right_bd - ($pfam_lb * 3) + 3;
               $pfam_lba = $pfam_rba - ($length * 3) - 2;

               #modify boundaries is necessary
               if ($pfam_lba == $left_bd) {
                  $pfam_lba += 3;  #add one codon to avoid equal to stop position
               }
               if ($pfam_rba == $right_bd) {
                  $pfam_rba -= 3;  #remove one codon to avoid equal to start position
               }

               $boundary_pfam  = '     PFAM_match      complement('.$pfam_lba.'..'.$pfam_rba."\)\n";
            }


            #create Pfam label
            $value = '/label='.$model."\n";
            open  ENTRY, '>', \$label;
            write ENTRY;
            close ENTRY;

            #create Pfam match
            $value = '/pfam_match="'.${$args{Pfamcode_to_number}}->{$model}.': '.
               ${$args{Pfamcode_to_family}}->{$model}.': ';

            #add verbose description if selected
            if (${$args{auto_ini_ref}}{pfam_verbose} =~ /(verbose|both)/) {
               $value .= ${$args{Pfamcode_to_descriptor}}->{$model};
            }
            $value =~ s/\:?\s+$//g;
            $value =~ s/pfam_match\=\s*\: ?$//;
            $value .= '"'."\n";
            open  ENTRY, '>', \$pfam_match;
            write ENTRY;
            close ENTRY;

            #create InterPro descriptor if selected
            if (${$args{auto_ini_ref}}{pfam_verbose} =~ /(InterPro|both)/ && ${$args{IPaccess_to_descriptor}}->{${$args{Pfamcode_to_access}}->{$model}} =~ /\w+/) {
               $value = '/interpro="'.${$args{IPaccess_to_code}}->{${$args{Pfamcode_to_access}}->{$model}}.': '.
                  ${$args{IPaccess_to_descriptor}}->{${$args{Pfamcode_to_access}}->{$model}};
               $value =~ s/\:?\s+$//g;
               $value .= '"'."\n";
               open  ENTRY, '>', \$interpro;
               write ENTRY;
               close ENTRY;
            } else {
               $interpro = "";
            }

            #create Pfam note
            $value = '/note="Cumulative evalue for family '.$model.': '.$evalue.'"'."\n";
            open  ENTRY, '>', \$note;
            write ENTRY;
            close ENTRY;

            #create product entry
            $value = '/product="'.$model.', '.${$args{Pfamcode_to_family}}->{$model}.
                     ' Score='.$score.' Expect='.$domain_evalue.'"'."\n";
            open  ENTRY, '>', \$product;
            write ENTRY;
            close ENTRY;

            #combine feature
            $feature_pfam = $boundary_pfam.
                            $label.
                            $pfam_match.
                            $interpro.
                            $note.
                            $product;
            $key_pfam     = $pfam_lba.'_'.$pfam_rba.'_PFAM_match_'.$id; #maintain original uniqe ID instead of ORFnumber

            ${$args{genbank_ref}}{$key_pfam} = $feature_pfam;
            push (@{$args{feature_list_ref}}, $key_pfam);
         }
      }
   }
   &hide_pbar_2;
   return;
}

sub pfam3_parser {
   my %args = @_;

   #define alignment symbology
   my %alignment_symbols = ();
   $alignment_symbols{'..'} = 'both ends of the alignment ended internally';
   $alignment_symbols{'[]'} = 'both ends of the alignment were full-length flush to the ends of the query or target';
   $alignment_symbols{'[.'} = 'only the left end was flush or full length';
   $alignment_symbols{'.]'} = 'only the right end was flush or full length';


   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling PFam results",
                   label        => ''
                  );
   &show_pbar_2;

   #parse through PFam results
   foreach (my $i = 1; $i < $args{counter}; $i++) {
      my ($filename, $id, $left_bd, $right_bd, $orientation, @pfam_db,
          $pfam_to_domain, $model_to_domain, @pfam_key, $counter);

      #get filename, left_bd, right_bd, orientation
      'reset' =~ m/reset/;
      $args{orf_to_id}->{$i} =~ m/^(.*?)___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
      $filename    = $1;
      $id          = $2;
      $left_bd     = $3;
      $right_bd    = $4;
      $orientation = $5;
      unless (defined $orientation) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not parse from entry $args{orf_to_id}->{$i}\.",
                                           -icon    => 'error',
                                           -type    => 'OK');
         &hide_pbar_2;
         return (0);
      }

      #update status bar
      if (($i % 20) == 0) {
         &update_pbar_2(title        => "Compiling PFam3 results",
                        label        => "Compiling PFam3 results for $filename",
                        progress     => ($i / $args{counter}) * 100,
                       );
      }

      #iterate through all selected databases
      @pfam_db = split / ; /,${$args{auto_ini_ref}}{Pfam_db};

      $counter = 1;
      foreach my $pfam_db (@pfam_db) {
         my ($model, $position, @model);

         #test if file exists
         next unless (-e ${$args{ini_ref}}{pfam_results}.'/'.$filename.'_'.$pfam_db.'_'.$id);

         #read input file
         (my $file_ref) = slurp(main_window => $args{main_window},
                                directory   => ${$args{ini_ref}}{pfam_results},
                                filename    => $filename.'_'.$pfam_db.'_'.$id
                               );
         if ($file_ref eq '0') {
            &hide_pbar_2;
            return (0)
         };

         #skip if no PFam results
         next if (${$file_ref} =~ /\[No hits detected that satisfy reporting thresholds\]/s);

         #grab family and domain box
         'reset'      =~ /reset/;
         ${$file_ref} =~ /Scores for complete sequence \(score includes all domains\)\:(.*?)Domain annotation for each model \(and alignments\)\:(.*?)Internal pipeline statistics summary/s;
         $model       = $1;
         $position    = $2;

         unless (defined $position) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not parse PFam3 results from entry $filename\_$pfam_db\_$id.",
                                              -icon    => 'error',
                                              -type    => 'OK');
            &hide_pbar_2;
            return (0);
         }

         #clean up
         $model =~ s/^.*?--------      -----------\n\r?//s;
         $model =~ s/\s*------ inclusion threshold ------\s*\n//s;

         #save domain block away
         $pfam_to_domain->{$pfam_db} = $position;

         #create arrays
         @model = split /\n\r?/, $model;

         foreach my $entry (@model) {
            next if ($entry !~ /\w+/);
            my ($local_model, $score, $evalue, $key);
            'reset'      =~ m/reset/;
            $entry       =~ m/^\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s/;
            $evalue      = $1;
            $score       = $2;
            $local_model = $3;

            unless (defined $local_model) {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Could not parse PFam3 results from entry $entry in file $filename\_$pfam_db\_$id.",
                                                 -icon    => 'error',
                                                 -type    => 'OK');
               &hide_pbar_2;
               return (0);
            }

            #check evalue
            if ($evalue =~ /^e/) {$evalue = '1'.$evalue};

            $key = $counter.'___'.$pfam_db.'___'.$local_model.'___'.$evalue;
            push (@pfam_key, $key);
            undef $key;
            $counter++;
         }
      }

      #sort Pfam array and then get the No. of best Pfam family models indicated in the current selection
      my @sorted =
         map  $_->[0] =>
         sort { $a->[1] <=> $b->[1] }
         map  [ $_, m/^\d+___\S+___\S+___(.+)$/ ]
         => @pfam_key;
      undef @pfam_key;

      @pfam_key = splice(@sorted,0,${$args{auto_ini_ref}}{Pfam_family_max});
      undef @sorted;

      #now process top models and create of features
      foreach my $best_model (@pfam_key) {
         my ($db, $model, $evalue, @domains);
         'reset' =~ /reset/;
         $best_model =~ m/\d+___(.+)___(.+)___(.+)/;
         $db     = $1;
         $model  = $2;
         $evalue = $3;

         unless (defined $evalue) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not parse PFam key $best_model in file $filename\_$db\_$id.",
                                              -icon    => 'error',
                                              -type    => 'OK');
            &hide_pbar_2;
            return (0);
         }

         #skip if above threshold
         next if ($evalue > ${$args{ini_ref}}{PFam_cut});

         #grab domains from domain block
         @domains = ($pfam_to_domain->{$db} =~ m/(>>\s*$model.+?)Alignments for each domain/gs);

         #iterate individual domains
         foreach my $list (@domains) {
            my (@domain);
            #clean up
            $list =~ s/^.*?------- -------    ----\n\r?//s;
            $list =~ s/\n\n//gs;
            'reset' =~ /reset/;
            #capture individual domains
            @domain = split /\n/, $list;

            foreach my $entry (@domain) {
               my ($pfam_lb, $pfam_rb, $pfam_lba, $pfam_rba, $boundary_pfam, $score, $i_evalue, $symbol,
                   $satisfy_threshold, $label, $pfam_match, $interpro, $note, $alignment, $product,
                   $feature_pfam, $key_pfam);

               'reset' =~ /reset/;
               $entry =~ m/^\s+\d+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+(\S+)\s/;
               ($satisfy_threshold, $score, $i_evalue, $pfam_lb, $pfam_rb, $symbol) = ($1, $2, $3, $4, $5, $6);

               unless (defined $symbol) {
                  ${$args{main_window}}->messageBox(-title   => 'Error',
                                                    -message => "Could not parse domain entry $entry in file $filename\_$db\_$id.",
                                                    -icon    => 'error',
                                                    -type    => 'OK');
                     &hide_pbar_2;
                     return (0);
               }

               #test for domain threshold if selected
               if (${$args{auto_ini_ref}}{Pfam_restrict_domains} == 1) {
                  next if ($i_evalue > ${$args{ini_ref}}{PFam_cut});
               }

               #create boundary entry
               if ($orientation eq 'sense') {
                  $pfam_rba = $left_bd + ($pfam_rb * 3) - 1;
                  $pfam_lba = $left_bd + ($pfam_lb * 3) - 3;

                  #modify boundaries is necessary
                  if ($pfam_lba == $left_bd) {
                     $pfam_lba += 3;  #add one codon to avoid equal to start position
                  }
                  if ($pfam_rba == $right_bd) {
                     $pfam_rba -= 3;  #remove one codon to avoid equal to stop position
                  }

                  $boundary_pfam  = '     PFAM_match      '.$pfam_lba.'..'.$pfam_rba."\n";
               } elsif ($orientation eq 'antisense') {
                  my $length = $pfam_rb - $pfam_lb;
                  $pfam_rba = $right_bd - ($pfam_lb * 3) + 3;
                  $pfam_lba = $pfam_rba - ($length * 3) - 2;

                  #modify boundaries is necessary
                  if ($pfam_lba == $left_bd) {
                     $pfam_lba += 3;  #add one codon to avoid equal to stop position
                  }
                  if ($pfam_rba == $right_bd) {
                     $pfam_rba -= 3;  #remove one codon to avoid equal to start position
                  }

                  $boundary_pfam  = '     PFAM_match      complement('.$pfam_lba.'..'.$pfam_rba."\)\n";
               }

               #create Pfam label
               $value = '/label='.$model."\n";
               open  ENTRY, '>', \$label;
               write ENTRY;
               close ENTRY;

               #create Pfam alignment note
               $value = '/note='.$alignment_symbols{$symbol}."\n";
               open  ENTRY, '>', \$alignment;
               write ENTRY;
               close ENTRY;

               #create Pfam match
               $value = '/pfam_match="'.${$args{Pfamcode_to_number}}->{$model}.': '.
                  ${$args{Pfamcode_to_family}}->{$model}.': ';

               #add verbose description if selected
               if (${$args{auto_ini_ref}}{pfam_verbose} =~ /(verbose|both)/) {
                  $value .= ${$args{Pfamcode_to_descriptor}}->{$model};
               }
               $value =~ s/\:?\s+$//g;
               $value .= '"'."\n";
               open  ENTRY, '>', \$pfam_match;
               write ENTRY;
               close ENTRY;

               #create InterPro descriptor if selected
               if (${$args{auto_ini_ref}}{pfam_verbose} =~ /(InterPro|both)/ &&
                   ${$args{IPaccess_to_descriptor}}->{${$args{Pfamcode_to_access}}->{$model}} =~ /\w+/) {
                  $value = '/interpro="'.${$args{IPaccess_to_code}}->{${$args{Pfamcode_to_access}}->{$model}}.': '.
                     ${$args{IPaccess_to_descriptor}}->{${$args{Pfamcode_to_access}}->{$model}};
                  $value =~ s/\:?\s+$//g;
                  $value .= '"'."\n";
                  open  ENTRY, '>', \$interpro;
                  write ENTRY;
                  close ENTRY;
               } else {
                  $interpro = "";
               }

               #create Pfam note
               $value = '/note="Cumulative evalue for family '.$model.': '.$evalue.'"'."\n";
               open  ENTRY, '>', \$note;
               write ENTRY;
               close ENTRY;

               #create product entry
               if ($satisfy_threshold eq '!') {
                  $value = '/product="'.$model.', '.${$args{Pfamcode_to_family}}->{$model}.
                           ' Score='.$score.' Expect='.$i_evalue.'"'."\n";
               } else {
                  $value = '/product="['.$model.'], '.${$args{Pfamcode_to_family}}->{$model}.
                           ' Score='.$score.' Expect='.$i_evalue.' (putative)"'."\n";
               }
               open  ENTRY, '>', \$product;
               write ENTRY;
               close ENTRY;

               #combine feature
               $feature_pfam = $boundary_pfam.
                               $label.
                               $alignment.
                               $pfam_match.
                               $interpro.
                               $note.
                               $product;
               $key_pfam     = $pfam_lba.'_'.$pfam_rba.'_PFAM_match_'.$id; #maintain original uniqe ID instead of ORFnumber

               ${$args{genbank_ref}}{$key_pfam} = $feature_pfam;
               push (@{$args{feature_list_ref}}, $key_pfam);
            }
         }
      }
   }
   &hide_pbar_2;
   return;
}


sub translate_pfam {
   my %args = @_;
   my ($file_ref, $max_length);
   my $release_number = 0;

   #determine database release version
   open READ, "\<${$args{ini_ref}}{pfam_db_path}".'/relnotes.txt' or do {
      ${$args{main_window}}->Dialog(-title   => 'Error',
                                    -text    => "Error opening PFam relnote file. Aborting",
                                    -buttons => ['OK'],
                                    -bitmap  => 'info');
      return (0);
   };

   while (<READ>) {
      next unless ($_ =~ m/^\s+RELEASE\s/);
      m/^\s+RELEASE\s+(\d+)/;
      $release_number = $1;
      unless (defined $release_number && $release_number =~ m/\d+/) {
         ${$args{main_window}}->Dialog(-title   => 'Error',
                                       -text    => "Error processing PFam Release information.\n$_\n\n",
                                       -buttons => ['OK'],
                                       -bitmap  => 'info');
         close READ;
         return (0);
      }
      last;
   }
   close READ;

   #go through descriptor files for all selected PFam databases
   my @selected_dbs = split /;/, ${$args{auto_ini_ref}}{full_Pfam_db};
   foreach my $pfam_db (@selected_dbs) {
      next unless ($pfam_db =~ m/\w+/);
      my $file                   = '';
      my $choose_descriptor_file = '';
      #set default values for PfamA and PfamB or pfam fs/ls
      if ($pfam_db =~ m/^.+\/pfam.a/i) {
         #Pfam descriptor file to open
         $file = ${$args{ini_ref}}{Pfam_descriptor}.'/pfamA.txt';
         $choose_descriptor_file = 'pfamA.txt'
      } elsif ($pfam_db =~ m/^.+\/pfam.b/i) {
         #Pfam descriptor file to open
         $file = ${$args{ini_ref}}{Pfam_descriptor}.'/pfamB.txt';
         $choose_descriptor_file = 'pfamB.txt'
      } else {
         #user needs to define file
         $file = '';
         $choose_descriptor_file = $pfam_db;
      }

      #ignore HMMER3 PfamB descriptor for now, no useful information stored
      next if ($file =~ m/pfamB.txt/ && $args{HMMERver} == 3);

      #read Pfam descriptions and domain keys
      my $types = [ ['txt', ['.txt']],
                    ['All', ['.*'  ]]
                  ];
      unless (-e $file) {
         $file = ${$args{main_window}}->getOpenFile(-initialdir  => ${$args{ini_ref}}{Pfam_descriptor},
                                                    -title       => "Select Pfam descriptor file \($choose_descriptor_file\)",
                                                    -filetypes   => $types,
                                                    -initialfile => ''
                                                   );
         unless (defined $file) {return (0)};
      }

      #change path to Pfam_descriptor to allow to save changes
      ${$args{ini_ref}}{Pfam_descriptor} = $file;
      ${$args{ini_ref}}{Pfam_descriptor} =~ s/^(.+)\/.+?$/$1/;

      #read PFam descriptions and domain keys
      open READ, '<'.$file or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not open Pfam descriptor in directory ${$args{ini_ref}}{Pfam_descriptor}.",
                                           -icon    => 'error',
                                           -type    => 'OK');
         return (0);
      };

      while (<READ>) {
         my @tmp;
         next if ($_ !~ /\t/);
         @tmp = split/\t/, $_;
         my $PFam_access     = '';
         my $PFam_number     = '';
         my $PFam_code       = '';
         my $PFam_family     = '';
         my $PFam_descriptor = '';

         if ($release_number < 29) {
            $PFam_access     = $tmp[0];
            $PFam_number     = $tmp[1];
            $PFam_code       = $tmp[2];
            $PFam_family     = $tmp[4];
            $PFam_descriptor = $tmp[9];
         } elsif ($release_number >= 29) {
            $PFam_access     = $tmp[0]; #new format in interpro/Pfam does not have access numbers anymore
            $PFam_number     = $tmp[0];
            $PFam_code       = $tmp[1];
            $PFam_family     = $tmp[3];
            $PFam_descriptor = $tmp[8];
            $PFam_access     =~ s/[\"\']//g;
            $PFam_access     =~ s/PF//ig; #same foramt as for interpro, no 'PF' header
         } else {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Unknown Version Error in parsing PFam descriptors.",
                                              -icon    => 'error',
                                              -type    => 'OK');
            return (0);
         };
         undef @tmp;
         unless (defined $PFam_descriptor) {
             ${$args{main_window}}->messageBox(-title   => 'Error',
                                               -message => "Error parsing Pfam descriptor line $_.\nin file $file",
                                               -icon    => 'error',
                                               -type    => 'OK');
             return (0);
         }

         if ($PFam_descriptor eq '\N' || $PFam_descriptor eq 'NULL') {$PFam_descriptor = ''};
         $PFam_access     =~ s/[\"\']//g;
         $PFam_code       =~ s/[\"\']//g;
         $PFam_code       =~ s/\s+//g;
         $PFam_family     =~ s/[\"\']//g;
         $PFam_number     =~ s/[\"\']//g;
         $PFam_number     =~ s/PF//ig;
         $PFam_descriptor =~ s/^[\"\']//;
         $PFam_descriptor =~ s/[\"\']$//;
         $PFam_descriptor =~ s/^\s+//;
         $PFam_descriptor =~ s/[\;\,]* ?$//;
         $PFam_descriptor =~ s/\"/\'/gs;
         $PFam_descriptor =~ s/([ \.\,\;]?)\[\d+\]([ \.\,\;])/$1$2/g;
         $PFam_descriptor =~ s/\s+(\.\s*)/$1/g;

         ${$args{Pfamcode_to_descriptor}}->{$PFam_code} = $PFam_descriptor;
         ${$args{Pfamcode_to_number}}->{$PFam_code}     = $PFam_number;
         ${$args{Pfamcode_to_family}}->{$PFam_code}     = $PFam_family;
         ${$args{Pfamcode_to_access}}->{$PFam_code}     = $PFam_access;
      }
      close READ;
   }
}

sub translate_interpro {
   my %args = @_;
   my (@inter, $interpro_dir, $interpro_txt, $file_ref);
   #Interpro file to open
   $interpro_dir = ${$args{ini_ref}}{Interpro_descriptor};
   $interpro_txt = 'interpro.txt';

   #read InterPro descriptions and domain keys
   unless (-e $interpro_dir.'/'.$interpro_txt) {
      my $file = ${$args{main_window}}->getOpenFile(-initialdir => $interpro_dir,
                                                    -title      => 'Select Interpro descriptor file'
                                                   );
      unless (defined $file) {return (0)};
      'reset'        =~ /reset/;
      $file          =~ m/(.+)\/(.+)$/;
      $interpro_dir  = $1;
      $interpro_txt  = $2;
   }

   open READ, '<'.$interpro_dir.'/'.$interpro_txt or do {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Could not open Interpro descripor in directory $interpro_dir.",
                                        -icon    => 'error',
                                        -type    => 'OK');

      return (0);
   };

   #read  codes
   $file_ref = &slurp(main_window  => $args{main_window},
                      directory    => $interpro_dir,
                      filename     => $interpro_txt
                     );

   if ($$file_ref =~ /[\"\']\n[\"\']/) {
      @inter = split/[\"\']\n+[\"\']/,${$file_ref}; #"
   } elsif ($$file_ref =~ /\\\n\n+/) {
      @inter = split/\\\n\n/,${$file_ref};
   } elsif ($$file_ref =~ /\nPF\d+\s+IPR\d+/) {
      $$file_ref = "\n".$$file_ref;
      #release 29 onwards
      @inter = split/\nPF/,${$file_ref};
   } else {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Unknown format in interpro.txt.",
                                        -icon    => 'error',
                                        -type    => 'OK');

      return (0);
   };

   foreach (@inter) {
      next unless ($_ =~ m/\S+/);
      s/[\"\']//gs;
      'reset' =~ /reset/;
      m/(\d+)\t(IPR\d+)\t\\(.+)/s;
      my $Interpro_access = $1;
      my $Interpro_code = $2;
      my $Interpro_descriptor = $3;

      #clean up HTML
      $Interpro_descriptor =~ s/^n\<p\>//;
      $Interpro_descriptor =~ s/^n//;
      'reset' =~ m/reset/;
      $Interpro_descriptor =~ s/\(\<a\s+class[^\>]+?www.uniprot.org\/uniprot[^\>]+?(\d+).?\>SWISSPROT\<\/a\>\)/\(SWISSPROT $1\)/gs; #retain this
      'reset' =~ m/reset/;
      $Interpro_descriptor =~ s/\<a\s+class[^\>]+?www.uniprot.org\/uniprot[^\>]+?\d+.?\>SWISSPROT\<\/a\>//gs; #remove this, internal link only
      'reset' =~ m/reset/;
      $Interpro_descriptor =~ s/\<a\s+class[^\>]+?www.ncbi.nlm.nih.gov[^\>]+?PubMed[^\>]+?\>(PUBMED.\d+)\<\/a\>/$1/gs;
      'reset' =~ m/reset/;
      $Interpro_descriptor =~ s/\<a\s+class[^\>]+?www.ebi.ac.uk\/interpro\/entry\/(IPR\d+)\>INTERPRO\<\/a\>/INTERPRO $1/gs;
      'reset' =~ m/reset/;
      $Interpro_descriptor =~ s/\<a\s+[^\>]+?www.ncbi.nlm.nih.gov\/Taxonomy\/Browser[^\>]+?id\=\d+\>[\-\_\w\s]+?\<\/a\>\s//gs;
      'reset' =~ m/reset/;
      $Interpro_descriptor =~ s/\<a\s+class[^\>]+?expasy.org\/prosite\/(PDOC\d+)\>PROSITEDOC\<\/a\>/PROSITE $1/gs;
      $Interpro_descriptor =~ s/\<a.+?\/a\>//gs;
      $Interpro_descriptor =~ s/\\n\<\/p\>\\n//gs;
      $Interpro_descriptor =~ s/\<\/p\>\\n\<p\>//gs;
      $Interpro_descriptor =~ s/\<\/p\>\\n//gs;
      $Interpro_descriptor =~ s/\\n/ /gs;
      $Interpro_descriptor =~ s/ {3,20}/\n/gs;
      $Interpro_descriptor =~ s/\(e\.g\.\s+\)//gs;
      $Interpro_descriptor =~ s/\(\)//gs;
      $Interpro_descriptor =~ s/\<sup\>//gs;
      $Interpro_descriptor =~ s/\<\/?p\>/\n/gs;
      $Interpro_descriptor =~ s/\n+/ /gs;
      'reset' =~ m/reset/;
      $Interpro_descriptor =~ s/ ([\.\,])/$1/gs;
      $Interpro_descriptor =~ s/\<\/?ul\>/ /gs;
      $Interpro_descriptor =~ s/\<\/?li\>/ /gs;
      $Interpro_descriptor =~ s/\<\/?ol\>/ /gs;
      $Interpro_descriptor =~ s/\<\/?pre\>/ /gs;
      $Interpro_descriptor =~ s/\<\/?su[pb]\>/ /gs;
      $Interpro_descriptor =~ s/ +/ /gs;
      $Interpro_descriptor =~ s/\(\)//sg;
      'reset' =~ m/reset/;
      $Interpro_descriptor =~ s/\s+([\.\,]) /$1 /gs;
      $Interpro_descriptor =~ s/\"/\'/gs;

      #from release 29 onwards, Interpro_access equals PFam domain number
      ${$args{IPaccess_to_code}}->{$Interpro_access}       = $Interpro_code;
      ${$args{IPaccess_to_descriptor}}->{$Interpro_access} = $Interpro_descriptor;
   }
}

1;

