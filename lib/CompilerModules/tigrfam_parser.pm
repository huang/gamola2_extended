#!/opt/ActivePerl-5.8/bin/perl

#TIGRfam result parser
#input arguments: main_window, progress_bar, gb_file, gene_model

package CompilerModules::tigrfam_parser;
use strict;
use ProgrammeModules::genbank_parser qw(:DEFAULT);
use initialise::read_me              qw(:DEFAULT);
use Basics::progress_bar             qw(:DEFAULT);

use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&tigrfam_parser &tigrfam3_parser &translate_tigrfam &tigr_annotation);
use vars qw();


#local variables
my (%args, %seen, $value);

format ENTRY =
~~                   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$value
.

sub tigrfam_parser {
   my %args = @_;
   my @feature_list;
   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling TIGRfam results",
                   label        => ''
                  );
   &show_pbar_2;

   #check if feature_list has CDS/gene features. if not, recompile gb file and try getting CDS/gene from there
   #this is needed for the optional TIGRfam gene annotation
   {
      #get respective gene and CDS features
      my @features = grep (/\_(gene|CDS)\_/, @{$args{feature_list_ref}});
      #try and get feature model from original Genbank input file if none found
      if ($#features < 0) {
         my $gb_file = $args{input_file};
         $gb_file =~ s/_GAMOLAdna$//;
         my (undef, undef, undef, undef, undef,
             $feature_list_ref, $genbank_ref, undef) = &gb_parser (main_window   => $args{main_window},
                                                                   progress_bar  => $args{progress_bar},
                                                                   auto_ini_ref  => $args{auto_ini_ref},
                                                                   ini_ref       => $args{ini_ref},
                                                                   gb_file       => ${$args{ini_ref}}{move_gb}.'/'.$gb_file
                                                                  );
         @features = grep (/_(gene|CDS)\_/, @{$feature_list_ref});
         #add gene and CDS features to global genbank ref if transfer annotation is selected
         if (${$args{auto_ini_ref}}{TIGR_annotation} == 1) {
            foreach my $entry (@features) {
               $args{genbank_ref}{$entry} = ${$genbank_ref}{$entry};
            }
         }
         undef $feature_list_ref;
         undef $genbank_ref;
      }
      @feature_list = @features;
      undef @features;
   }

   #parse through TIGRfam results
   foreach (my $i = 1; $i < $args{counter}; $i++) {
      my ($filename, $id, $left_bd, $right_bd, $orientation, @tigrfam_db,
          $tigrfam_to_domain, $model_to_domain, @tigrfam_key, $counter,
          $best_evalue, $annotation);

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
                                           -message => "Could not parse boundaries from entry $args{orf_to_id}->{$i}\.",
                                           -icon    => 'error',
                                           -type    => 'OK');
         &hide_pbar_2;
         return (0);
      }

      #update status bar
      if (($i % 20) == 0) {
         &update_pbar_2(title        => "Compiling TIGRfam results",
                        label        => "Compiling TIGRfam results for $filename",
                        progress     => ($i / $args{counter}) * 100,
                       );
      }

      #iterate through all selected databases
      @tigrfam_db = split / ; /,${$args{auto_ini_ref}}{TIGRfam_db};

      $counter = 1;
      foreach my $tigrfam_db (@tigrfam_db) {
         my ($family, $domain, @family);

         #test if file exists
         next unless (-e ${$args{ini_ref}}{TIGRfam_results}.'/'.$filename.'_'.$tigrfam_db.'_'.$id);

         #read input file
         (my $file_ref) = slurp(main_window => $args{main_window},
                                directory   => ${$args{ini_ref}}{TIGRfam_results},
                                filename    => $filename.'_'.$tigrfam_db.'_'.$id
                               );
         if ($file_ref eq '0') {&hide_pbar_2; return (0)};

         #skip if no TIGRFam results
         next if (${$file_ref} =~ /\[no hits above thresholds\]/s || ${$file_ref} =~ /\[No hits detected that satisfy reporting thresholds\]/);

         #grab family and domain box
         'reset'      =~ /reset/;
         ${$file_ref} =~ /Scores for sequence family classification(.*?)Parsed for domains(.*?)Alignments of top-scoring domains/s;
         $family      = $1;
         $domain      = $2;

         unless (defined $domain) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not parse TIGRFam results for family or domain from entry $filename\_$tigrfam_db\_$id.",
                                              -icon    => 'error',
                                              -type    => 'OK');
            &hide_pbar_2;
            return (0);
         }

         #clean up
         $family =~ s/.*?------- ---\n\r?//s;
         $domain =~ s/.*?-------\n\r?//s;

         #save domain block away
         $tigrfam_to_domain->{$tigrfam_db} = $domain;

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
                                                 -message => "Could not parse TIGRFam results for score and evalue from entry $family in file $filename\_$tigrfam_db\_$id.",
                                                 -icon    => 'error',
                                                 -type    => 'OK');
               &hide_pbar_2;
               return (0);
            }

            #check evalue
            if ($evalue =~ /^e/) {$evalue = '1'.$evalue};

            $key = $counter.'___'.$tigrfam_db.'___'.$model.'___'.$evalue;
            push (@tigrfam_key, $key);
            undef $key;
            $counter++;

         }
      }

      #sort TIGRfam array and then get the best TIGRfam models as indicated in current selection
      my @sorted =
         map  $_->[0] =>
         sort { $a->[1] <=> $b->[1] }
         map  [ $_, m/^\d+___\S+___\S+___(.+)$/,  ]
         => @tigrfam_key;
      undef @tigrfam_key;

      #delete duplicates
      my %seen = ();
      @sorted = grep { ! $seen{$_} ++ } @sorted;

      @tigrfam_key = splice(@sorted,0,${$args{auto_ini_ref}}{TIGRfam_family_max});
      undef @sorted;

      #now process top  models and create features
      $best_evalue = 100;
      foreach my $best_model (@tigrfam_key) {
         my ($db, $model, $evalue, @domains,
             $AC, $DE, $EN, $EC, $CC, $GS);
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
         next if ($evalue > ${$args{ini_ref}}{TIGRfam_cut});

         #grab domains from domain block
         @domains = ($tigrfam_to_domain->{$db} =~ m/($model [^\n\r]+)\n\r?/gs);

         #iterate individual domains
         foreach my $domain (@domains) {
            my ($tigrfam_lb, $tigrfam_rb, $tigrfam_lba, $tigrfam_rba, $boundary_tigrfam, $score, $domain_evalue,
                $label, $tigrfam_match, $go_interpro, $note, $product, $EC_number, $colour,
                $feature_tigrfam, $key_tigrfam);
            'reset' =~ /reset/;
            $domain =~ m/\S+\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+\d+\s+\d+\s+\S+\s+(\S+)\s+(\S+)/;
            $tigrfam_lb    = $1;
            $tigrfam_rb    = $2;
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
            #parse TIGRfam info
            {
               if (${$args{TIGRdomain}}->{$model} =~ /AC\s+/) {
                  'reset' =~ /reset/;
                  ${$args{TIGRdomain}}->{$model} =~ /AC\s+(.*?)\n\r?/s;
                  $AC = $1;
               } else {
                  $AC = "";
               }

               if (${$args{TIGRdomain}}->{$model} =~ /DE\s+/) {
                  'reset' =~ /reset/;
                  my @tmp = (${$args{TIGRdomain}}->{$model} =~ /DE\s+(.*?)\n\r?/gs);
                  $DE = join (' ',@tmp);
                  undef @tmp;
               } else {
                  $DE = "";
               }

               if (${$args{TIGRdomain}}->{$model} =~ /EN\s+/) {
                  'reset' =~ /reset/;
                  my @tmp = (${$args{TIGRdomain}}->{$model} =~ /EN\s+(.*?)\n\r?/gs);
                  $EN = join (' ',@tmp);
                  undef @tmp;
               } else {
                  $EN = "";
               }

               if (${$args{TIGRdomain}}->{$model} =~ /GS\s+/) {
                  'reset' =~ /reset/;
                  my @tmp = (${$args{TIGRdomain}}->{$model} =~ /GS\s+(.*?)\n\r?/gs);
                  $GS = join (' ',@tmp);
                  undef @tmp;
               } else {
                  $GS = "";
               }

               if (${$args{TIGRdomain}}->{$model} =~ /EC\s+/) {
                  'reset' =~ /reset/;
                  ${$args{TIGRdomain}}->{$model} =~ /\n\r?EC\s+(.*?)\n\r?/s;
                  $EC = $1;
               } else {
                  $EC = "";
               }

               if (${$args{TIGRdomain}}->{$model} =~ /CC\s+/) {
                  'reset' =~ /reset/;
                  my @tmp = (${$args{TIGRdomain}}->{$model} =~ /CC\s+(.*?)\n\r?/gs);
                  $CC = join (' ',@tmp);
                  undef @tmp;
               } else {
                  $CC = "";
               }
            }

            # catch quotatiom marks
            $CC =~ s/\"/\'/gs;

            #create boundary entry
            if ($orientation eq 'sense') {
               $tigrfam_rba = $left_bd + ($tigrfam_rb * 3) - 1;
               $tigrfam_lba = $left_bd + ($tigrfam_lb * 3) - 3;

               #modify boundaries is necessary
               if ($tigrfam_lba == $left_bd) {
                  $tigrfam_lba += 3;  #add one codon to avoid equal to start position
               }
               if ($tigrfam_rba == $right_bd) {
                  $tigrfam_rba -= 3;  #remove one codon to avoid equal to stop position
               }

               $boundary_tigrfam  = '     TIGR_match      '.$tigrfam_lba.'..'.$tigrfam_rba."\n";
            } elsif ($orientation eq 'antisense') {
               my $length = $tigrfam_rb - $tigrfam_lb;
               $tigrfam_rba = $right_bd - ($tigrfam_lb * 3) + 3;
               $tigrfam_lba = $tigrfam_rba - ($length * 3) - 2;

               #modify boundaries is necessary
               if ($tigrfam_lba == $left_bd) {
                  $tigrfam_lba += 3;  #add one codon to avoid equal to stop position
               }
               if ($tigrfam_rba == $right_bd) {
                  $tigrfam_rba -= 3;  #remove one codon to avoid equal to start position
               }

               $boundary_tigrfam  = '     TIGR_match      complement('.$tigrfam_lba.'..'.$tigrfam_rba."\)\n";
            }

            #create TIGRfam label
            unless (defined ${$args{TIGRrole_GO}}->{$AC}) {
               ${$args{TIGRrole_GO}}->{$AC} = ' ';
            }
            if ($GS =~ /\S+/) {
               $value = '/label='.$GS.'; '.${$args{TIGRrole_GO}}->{$AC};
            } else {
               $value = '/label='.$model.'; '.${$args{TIGRrole_GO}}->{$AC};
            }
            $value =~ s/GO\:\s*$//;
            $value =~ s/;\s+$//;
            $value =~ s/\//_/g;
            $value =~ s/^_label/\/label/;
            $value .= "\n";

            open  ENTRY, '>', \$label;
            write ENTRY;
            close ENTRY;

            #create TIGRfam product
            $value = '/product="'.$AC.' (Role:'.${$args{TIGRrole_number}}->{$AC}.'): ';
            if ($DE =~ /\S+/) {
               $value .= $DE.' ';
            } else {
               $value .= $EN.' ';
            }

            #no role defined?
            if ($value =~ /(Role:)/) {
               $value =~ s/\(Role:\)//;
            }

            $value .= '; Score= '.$score.' Expect= '.$domain_evalue.'"';
            open  ENTRY, '>', \$product;
            write ENTRY;
            close ENTRY;

            #create TIGRfam note
            $value = '/note="Role:'.${$args{TIGRrole_number}}->{$AC}.' = '.
                     ${$args{TIGRrole_name}}->{${$args{TIGRrole_number}}->{$AC}};

            #add verbose description if selected
            if (${$args{auto_ini_ref}}{tigrfam_verbose} =~ /verbose/) {
               $value .= '; '.$CC;
            }
            $value .= '"'."\n";

            #nothing defined?
            if ($value =~ /note\=\"Role: \= ; \"/) {
               $value = '/note="Role: No role defined"'."\n";
            } elsif ($value =~ /note\=\"Role: \= ;\s*\w+/) {
               $value =~ s/note\=\"Role: \= ;/note\=\"Role: /s;
            }
            open  ENTRY, '>', \$note;
            write ENTRY;
            close ENTRY;

            #create TIGRfam EC number if defined
            if ($EC =~ /\S+/) {
               $EC_number = '                     /EC_number="'.$EC.'"'."\n";
            } else {
               $EC_number = '';
            }

            #create TIGRfam colour
            $colour = '                     /colour=4'."\n";

            #create GO from interpro entry
            if (defined ${$args{TIGRrole_GO}}->{$AC} && ${$args{TIGRrole_GO}}->{$AC} =~ /\S+/) {
               my ($GO_name, $GO_def);
               #parse GO entry
               'reset' =~ /reset/;
               ${$args{TIGRrole_GoClass}}->{${$args{TIGRrole_GO}}->{$AC}} =~ m/\<name\>(.*?)\<\/name\>/s;
               $GO_name = $1;
               $GO_name =~ s/\n\r?//gs;
               $GO_name =~ s/\s+/ /gs;

               ${$args{TIGRrole_GoClass}}->{${$args{TIGRrole_GO}}->{$AC}} =~ m/\<defstr\>(.*?)\<\/defstr\>/s;
               $GO_def = $1;
               $GO_def =~ s/\n\r?//gs;
               $GO_def =~ s/\s+/ /gs;
               $GO_def =~ s/\"/\'/gs;

               $value = '/go_from_interpro="'.${$args{TIGRrole_GO}}->{$AC}.', '.$GO_name.'; '.$GO_def.'"';
               open  ENTRY, '>', \$go_interpro;
               write ENTRY;
               close ENTRY;
            } else {
               $go_interpro = '';
            }


            #combine feature
            $feature_tigrfam = $boundary_tigrfam.
                               $label.
                               $EC_number.
                               $product.
                               $note.
                               $go_interpro.
                               $colour;
            $key_tigrfam     = $tigrfam_lba.'_'.$tigrfam_rba.'_TIGR_match_'.$id; #maintain original uniqe ID instead of ORFnumber

            ${$args{genbank_ref}}{$key_tigrfam} = $feature_tigrfam;
            push (@{$args{feature_list_ref}}, $key_tigrfam);

            if ($evalue < $best_evalue) {
               $best_evalue = $evalue;
               $annotation = $feature_tigrfam;
            }
         }
      }

      #add gene name, CDS name and EC number if so selected, if TIGRfam has been found
      if (${$args{auto_ini_ref}}{TIGR_annotation} == 1 && $annotation =~ /\w+/) {
         #get respective gene and CDS features
         my @features = grep (/$left_bd\_$right_bd\_(gene|CDS)\_/, @feature_list);

         #get TIGRfam info
         my ($label, $EC, $descriptor) = &tigr_annotation(main_window  => $args{main_window},
                                                          progress_bar => $args{progress_bar},
                                                          auto_ini_ref => $args{auto_ini_ref},
                                                          ini_ref      => $args{ini_ref},
                                                          annotation   => $annotation,
                                                         );
         #no results? next
         next if ($label !~ /\w+/ && $EC !~ /\w+/ && $descriptor !~ /\w+/);

         #update annotation
         foreach my $feature (@features) {
            #update gene/CDS annotation?
            if (${$args{auto_ini_ref}}{TIGR_ECCDS} == 1) {
               #first gene
               if ($feature =~ /gene/ && $label =~ /\w+/ && $args{genbank_ref}{$feature} =~ /\/gene/) {
                  my ($id, $annotation);
                  'reset' =~ m/reset/;
                  $args{genbank_ref}{$feature} =~ m/(\/gene\=\")[^\"]*?(_?\s*\d*\s*\")/s;
                  $id = $2;
                  unless (defined $id) {$id = '"'};
                  $value = '/gene="'.$label.$id;
                  open  ENTRY, '>', \$annotation;
                  write ENTRY;
                  close ENTRY;
                  $args{genbank_ref}{$feature} =~ s/\n\s*?\/gene\=\"[^\"]*?_?\s*\d*\s*\"\s*\n/\n$annotation/s;
               } elsif ($feature =~ /gene/ && $label =~ /\w+/ && $args{genbank_ref}{$feature} =~ /\/note/) {
                  my ($id, $annotation);
                  'reset' =~ m/reset/;
                  $args{genbank_ref}{$feature} =~ m/(\/note\=\")[^\"]*?(_?\s*\d*\s*\")/s;
                  $id = $2;
                  unless (defined $id) {$id = '"'};
                  $value = '/note="'.$label.$id;
                  open  ENTRY, '>', \$annotation;
                  write ENTRY;
                  close ENTRY;
                  $args{genbank_ref}{$feature} =~ s/\n\s*\/note\=\"[^\"].*?_?\s*\d*\s*\"\s*\n/\n$annotation/s;
               }
               #then CDS
               if ($feature =~ /CDS/ && $descriptor =~ /\w+/ && $args{genbank_ref}{$feature} =~ m/\/gene/) {
                  if ($label =~ /\w+/) {$label = ', '.$label};
                  my ($id, $annotation);
                  'reset' =~ m/reset/;
                  $args{genbank_ref}{$feature} =~ m/(\/gene\=\")[^\"]*?(_?\s*\d*\s*\")/s;
                  $id = $2;
                  unless (defined $id) {$id = '"'};
                  $value = '/gene="'.$descriptor.$label.$id;
                  $value =~ s/ , /, /gs;
                  open  ENTRY, '>', \$annotation;
                  write ENTRY;
                  close ENTRY;
                  $args{genbank_ref}{$feature} =~ s/\n\s*\/gene\=\"[^\"]*?_?\s*\d*\s*\"\s*\n/\n$annotation/s;
               } elsif ($feature =~ /CDS/ && $descriptor =~ /\w+/ && $args{genbank_ref}{$feature} =~ m/\/product/) {
                  if ($label =~ /\w+/) {$label = ', '.$label};
                  my ($id, $annotation);
                  'reset' =~ m/reset/;
                  $args{genbank_ref}{$feature} =~ m/(\/product\=\")[^\"]*?(_?\s*\d*\s*\")/s;
                  $id = $2;
                  unless (defined $id) {$id = '"'};
                  $value = '/product="'.$descriptor.$label.$id;
                  $value =~ s/ , /, /gs;
                  open  ENTRY, '>', \$annotation;
                  write ENTRY;
                  close ENTRY;
                  $args{genbank_ref}{$feature} =~ s/\n\s*\/product\=\"[^\"]*?_?\s*\d*\s*\"\s*\n/\n$annotation/s;
               } elsif ($feature =~ /CDS/ && $descriptor =~ /\w+/ && $args{genbank_ref}{$feature} =~ m/\/note/) {
                  if ($label =~ /\w+/) {$label = ', '.$label};
                  my ($id, $annotation);
                  'reset' =~ m/reset/;
                  $args{genbank_ref}{$feature} =~ m/(\/note\=\")[^\"]*?(_?\s*\d*\s*\")/s;
                  $id = $2;
                  unless (defined $id) {$id = '"'};
                  $value = '/note="'.$descriptor.$label.$id;
                  $value =~ s/ , /, /gs;
                  open  ENTRY, '>', \$annotation;
                  write ENTRY;
                  close ENTRY;
                  $args{genbank_ref}{$feature} =~ s/\n\s*\/note\=\"[^\"]*?_?\s*\d*\s*\"\s*\n/\n$annotation/s;
               }
            }

            #generate ECnumber entries
            my @EC = split /\s+/, $EC; #if there are multiple EC numbers, put each in a separate entry
            $EC    = '';
            foreach my $entry (@EC) {
               next unless ($entry =~ /\S+/);
               $EC .= '                     /EC_number='.$entry."\n";
            }
            $EC =~ s/^                     \/EC_number\=//;
            $EC =~ s/\n$//s;

            #delete existing EC number entries if selected
            #if ($args{genbank_ref}{$feature} =~ /\/EC_number/ && $EC =~ /\w+/ && ${$args{auto_ini_ref}}{TIGR_ECdiscard} == 1) {
            #   $args{genbank_ref}{$feature} =~ s/\/EC_number.*?(\s+?\/)/$1/gs;
            #}

            #add EC number
            if ($feature                     =~ /CDS/ &&
                $EC                          =~ /\w+/ &&
                $args{genbank_ref}{$feature} =~ /\/EC_number\=/ &&
                ${$args{auto_ini_ref}}{TIGR_ECdiscard} == 0) {
               next;
            } elsif ($feature                     =~ /CDS/ &&
                     $EC                          =~ /\w+/ &&
                     $args{genbank_ref}{$feature} =~ /\/(gene|product|note)/ &&
                     $args{genbank_ref}{$feature} !~ /\/EC_number\=/ &&
                     ${$args{auto_ini_ref}}{TIGR_ECdiscard} == 0) {
               if ($args{genbank_ref}{$feature} =~ /\/gene/) {
                  $args{genbank_ref}{$feature} =~ s/(\/gene\=\"[^\"]*?\")/$1\n                     \/EC_number\=$EC/s;
               } elsif ($args{genbank_ref}{$feature} =~ /\/product/) {
                  $args{genbank_ref}{$feature} =~ s/(\/product\=\"[^\"]*?\")/$1\n                     \/EC_number\=$EC/s;
               } elsif ($args{genbank_ref}{$feature} =~ /\/note/) {
                  $args{genbank_ref}{$feature} =~ s/(\/note\=\"[^\"]*?\")/$1\n                     \/EC_number\=$EC/s;
               }
            } elsif ($feature                     =~ /CDS/ &&
                     $EC                          =~ /\w+/ &&
                     $args{genbank_ref}{$feature} !~ /\/(gene|product|note)/ &&
                     $args{genbank_ref}{$feature} !~ /\/EC_number\=/ &&
                     ${$args{auto_ini_ref}}{TIGR_ECdiscard} == 0) {
               $args{genbank_ref}{$feature} =~ s/(                     \/)/                     \/EC_number\=$EC\n$1/s;
            } elsif ($feature                     =~ /CDS/ &&
                     $EC                          =~ /\w+/ &&
                     $args{genbank_ref}{$feature} =~ /\/EC_number/ &&
                     ${$args{auto_ini_ref}}{TIGR_ECdiscard} == 1) {
               $args{genbank_ref}{$feature} =~ s/(\/EC_number\=\"?.*?\"?)\n/\/EC_number\=$EC/s;
            } elsif ($feature                     =~ /CDS/ &&
                     $EC                          =~ /\w+/ &&
                     $args{genbank_ref}{$feature} !~ /\/EC_number/ &&
                     $args{genbank_ref}{$feature} =~ /\/(gene|product|note)/ &&
                     ${$args{auto_ini_ref}}{TIGR_ECdiscard} == 1) {
               if ($args{genbank_ref}{$feature} =~ /\/gene/) {
                  $args{genbank_ref}{$feature} =~ s/(\/gene\=\"[^\"]*?\")/$1\n                     \/EC_number\=$EC/s;
               } elsif ($args{genbank_ref}{$feature} =~ /\/product/) {
                  $args{genbank_ref}{$feature} =~ s/(\/product\=\"[^\"]*?\")/$1\n                     \/EC_number\=$EC/s;
               } elsif ($args{genbank_ref}{$feature} =~ /\/note/) {
                  $args{genbank_ref}{$feature} =~ s/(\/note\=\"[^\"]*?\")/$1\n                     \/EC_number\=$EC/s;
               }
            }  elsif ($feature                     =~ /CDS/ &&
                     $EC                          =~ /\w+/ &&
                     $args{genbank_ref}{$feature} !~ /\/EC_number/ &&
                     $args{genbank_ref}{$feature} !~ /\/(gene|product|note)/ &&
                     ${$args{auto_ini_ref}}{TIGR_ECdiscard} == 1) {
               $args{genbank_ref}{$feature} =~ s/(                     \/)/                     \/EC_number\=$EC\n$1/s;
            }
         }
      }
   }
   undef @feature_list;
   &hide_pbar_2;
   return;
}

sub tigrfam3_parser {
   my %args = @_;
   my @feature_list;

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
                   title        => "Compiling TIGRfam results",
                   label        => ''
                  );
   &show_pbar_2;

   #check if feature_list has CDS/gene features. if not, recompile gb file and try getting CDS/gene from there
   #this is needed for the optional TIGRfam gene annotation
   {
      #get respective gene and CDS features
      my @features = grep (/\_(gene|CDS)\_/, @{$args{feature_list_ref}});
      #try and get feature model from original Genbank input file if none found
      if ($#features < 0) {
         my $gb_file = $args{input_file};
         $gb_file =~ s/_GAMOLAdna$//;
         my (undef, undef, undef, undef, undef,
             $feature_list_ref, $genbank_ref, undef) = &gb_parser (main_window   => $args{main_window},
                                                                   progress_bar  => $args{progress_bar},
                                                                   auto_ini_ref  => $args{auto_ini_ref},
                                                                   ini_ref       => $args{ini_ref},
                                                                   gb_file       => ${$args{ini_ref}}{move_gb}.'/'.$gb_file
                                                                  );
         @features = grep (/_(gene|CDS)\_/, @{$feature_list_ref});
         #add gene and CDS features to global genbank ref if transfer annotation is selected
         if (${$args{auto_ini_ref}}{TIGR_annotation} == 1) {
            foreach my $entry (@features) {
               $args{genbank_ref}{$entry} = ${$genbank_ref}{$entry};
            }
         }
         undef $feature_list_ref;
         undef $genbank_ref;
      }
      @feature_list = @features;
      undef @features;
   }

   #parse through TIGRfam results
   foreach (my $i = 1; $i < $args{counter}; $i++) {
      my ($filename, $id, $left_bd, $right_bd, $orientation, @tigrfam_db,
          $tigrfam_to_domain, $model_to_domain, @tigrfam_key, $counter,
          $best_evalue, $annotation);

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
         &update_pbar_2(title        => "Compiling TIGRfam results",
                        label        => "Compiling TIGRfam results for $filename",
                        progress     => ($i / $args{counter}) * 100,
                       );
      }

      #iterate through all selected databases
      @tigrfam_db = split / ; /,${$args{auto_ini_ref}}{TIGRfam_db};

      $counter = 1;
      foreach my $tigrfam_db (@tigrfam_db) {
         my ($model, $position, @model);

         #test if file exists
         next unless (-e ${$args{ini_ref}}{TIGRfam_results}.'/'.$filename.'_'.$tigrfam_db.'_'.$id);

         #read input file
         (my $file_ref) = slurp(main_window => $args{main_window},
                                directory   => ${$args{ini_ref}}{TIGRfam_results},
                                filename    => $filename.'_'.$tigrfam_db.'_'.$id
                               );
         if ($file_ref eq '0') {&hide_pbar_2; return (0)};

         #skip if no TIGRFam results
         next if (${$file_ref} =~ /\[No hits detected that satisfy reporting thresholds\]/);

         #grab family and domain box
         'reset'      =~ /reset/;
         ${$file_ref} =~ /Scores for complete sequence \(score includes all domains\)\:(.*?)Domain annotation for each model \(and alignments\)\:(.*?)Internal pipeline statistics summary/s;
         $model       = $1;
         $position    = $2;

         unless (defined $position) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not parse TIGRFam results from entry $filename\_$tigrfam_db\_$id.",
                                              -icon    => 'error',
                                              -type    => 'OK');
            &hide_pbar_2;
            return (0);
         }

         #clean up
         $model =~ s/^.*?--------\s+-----------\s*\n\r?//s;
         $model =~ s/\s*------ inclusion threshold ------\s*\n//s;

         #save domain block away
         $tigrfam_to_domain->{$tigrfam_db} = $position;

         #create arrays
         @model = split /\n\r?/, $model;

         foreach my $entry (@model) {
            next if ($entry !~ /\w+/);
            my ($local_model, $local_description, $score, $evalue, $key);
            'reset'            =~ m/reset/;
            $entry             =~ m/^\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s/;
            $evalue            = $1;
            $score             = $2;
            $local_model       = $3;
            $local_description = $4;

            unless (defined $local_description) {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Could not parse TIGRFam results from entry $entry in file $filename\_$tigrfam_db\_$id.",
                                                 -icon    => 'error',
                                                 -type    => 'OK');
               &hide_pbar_2;
               return (0);
            }

            #check evalue
            if ($evalue =~ /^e/) {$evalue = '1'.$evalue};
            #cleanup descriptor
            $local_description =~ s/\:$//;

            $key = $counter.'___'.$tigrfam_db.'___'.$local_model.'___'.$evalue.'___'.$local_description;
            push (@tigrfam_key, $key);
            undef $key;
            $counter++;
         }
      }

      #sort TIGRfam array and then get the best TIGRfam models as indicated in current selection
      my @sorted =
         map  $_->[0] =>
         sort { $a->[1] <=> $b->[1] }
         map  [ $_, m/^\d+___\S+___\S+___(.+)___/,  ]
         => @tigrfam_key;
      undef @tigrfam_key;

      #delete duplicates
      my %seen = ();
      @sorted = grep { ! $seen{$_} ++ } @sorted;

      @tigrfam_key = splice(@sorted,0,${$args{auto_ini_ref}}{TIGRfam_family_max});
      undef @sorted;

      #now process top  models and create features
      $best_evalue = 100;
      foreach my $best_model (@tigrfam_key) {
         my ($db, $model, $evalue, $descriptor, @domains,
             $AC, $DE, $EN, $EC, $CC, $GS);
         'reset'     =~ /reset/;
         $best_model =~ m/\d+___(.+)___(.+)___(.+)___(.+)/;
         $db         = $1;
         $model      = $2;
         $evalue     = $3;
         $descriptor = $4;

         unless (defined $descriptor) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not parse TIGRFam key $best_model in file $filename\_$db\_$id.",
                                              -icon    => 'error',
                                              -type    => 'OK');
            &hide_pbar_2;
            return (0);
         }

         #skip if above threshold
         next if ($evalue > ${$args{ini_ref}}{TIGRfam_cut});

         #grab domains from domain block
         @domains = ($tigrfam_to_domain->{$db} =~ m/(>>\s*$model.+?)Alignments for each domain/gs);

         #iterate individual domains
         foreach my $list (@domains) {
            my (@domain);
            #clean up
            $list =~ s/^.*?-------\s+-------\s+----\n\r?//s;
            $list =~ s/\n\n//gs;
            'reset' =~ /reset/;
            #capture individual domains
            @domain = split /\n/, $list;

            foreach my $entry (@domain) {
               my ($tigrfam_lb, $tigrfam_rb, $tigrfam_lba, $tigrfam_rba, $boundary_tigrfam, $score, $i_evalue, $symbol,
                   $label, $tigrfam_match, $go_interpro, $note, $note1, $product, $EC_number, $colour,
                   $feature_tigrfam, $key_tigrfam, $tmp);
               'reset' =~ /reset/;
               my @tmp = split /\s+/, $entry;
               $score      = $tmp[3];
               $i_evalue   = $tmp[6];
               $tigrfam_lb = $tmp[10];
               $tigrfam_rb = $tmp[11];
               $symbol     = $tmp[12];
               @tmp = ();

               unless (defined $symbol) {
                  ${$args{main_window}}->messageBox(-title   => 'Error TIGRFam3',
                                                    -message => "Could not parse domain entry $entry in file $filename\_$db\_$id.",
                                                    -icon    => 'error',
                                                    -type    => 'OK');
                     &hide_pbar_2;
                     return (0);
               }
               #parse TIGRfam info
               if (defined ${$args{TIGRdomain}}->{$model}) {
                  if (${$args{TIGRdomain}}->{$model} =~ /AC\s+/) {
                     'reset' =~ /reset/;
                     ${$args{TIGRdomain}}->{$model} =~ /AC\s+(.*?)\n\r?/s;
                     $AC = $1;
                  } else {
                     $AC = "";
                  }
                  if (${$args{TIGRdomain}}->{$model} =~ /DE\s+/) {
                     'reset' =~ /reset/;
                     my @tmp = (${$args{TIGRdomain}}->{$model} =~ /DE\s+(.*?)\n\r?/gs);
                     $DE = join (' ',@tmp);
                     undef @tmp;
                  } else {
                     $DE = "";
                  }

                  if (${$args{TIGRdomain}}->{$model} =~ /EN\s+/) {
                     'reset' =~ /reset/;
                     my @tmp = (${$args{TIGRdomain}}->{$model} =~ /EN\s+(.*?)\n\r?/gs);
                     $EN = join (' ',@tmp);
                     undef @tmp;
                  } else {
                     $EN = "";
                  }

                  if (${$args{TIGRdomain}}->{$model} =~ /GS\s+/) {
                     'reset' =~ /reset/;
                     my @tmp = (${$args{TIGRdomain}}->{$model} =~ /GS\s+(.*?)\n\r?/gs);
                     $GS = join (' ',@tmp);
                     undef @tmp;
                  } else {
                     $GS = "";
                  }

                  if (${$args{TIGRdomain}}->{$model} =~ /EC\s+/) {
                     'reset' =~ /reset/;
                     ${$args{TIGRdomain}}->{$model} =~ /\n\r?EC\s+(.*?)\n\r?/s;
                     $EC = $1;
                  } else {
                     $EC = "";
                  }

                  if (${$args{TIGRdomain}}->{$model} =~ /CC\s+/) {
                     'reset' =~ /reset/;
                     my @tmp = (${$args{TIGRdomain}}->{$model} =~ /CC\s+(.*?)\n\r?/gs);
                     $CC = join (' ',@tmp);
                     undef @tmp;
                  } else {
                     $CC = "";
                  }
               } else {
                  if (${$args{TIGRdomain}}->{$descriptor} =~ /AC\s+/) {
                     'reset' =~ /reset/;
                     ${$args{TIGRdomain}}->{$descriptor} =~ /AC\s+(.*?)\n\r?/s;
                     $AC = $1;
                  } else {
                     $AC = "";
                  }
                  if (${$args{TIGRdomain}}->{$descriptor} =~ /DE\s+/) {
                     'reset' =~ /reset/;
                     my @tmp = (${$args{TIGRdomain}}->{$descriptor} =~ /DE\s+(.*?)\n\r?/gs);
                     $DE = join (' ',@tmp);
                     undef @tmp;
                  } else {
                     $DE = "";
                  }

                  if (${$args{TIGRdomain}}->{$descriptor} =~ /EN\s+/) {
                     'reset' =~ /reset/;
                     my @tmp = (${$args{TIGRdomain}}->{$descriptor} =~ /EN\s+(.*?)\n\r?/gs);
                     $EN = join (' ',@tmp);
                     undef @tmp;
                  } else {
                     $EN = "";
                  }

                  if (${$args{TIGRdomain}}->{$descriptor} =~ /GS\s+/) {
                     'reset' =~ /reset/;
                     my @tmp = (${$args{TIGRdomain}}->{$descriptor} =~ /GS\s+(.*?)\n\r?/gs);
                     $GS = join (' ',@tmp);
                     undef @tmp;
                  } else {
                     $GS = "";
                  }

                  if (${$args{TIGRdomain}}->{$descriptor} =~ /EC\s+/) {
                     'reset' =~ /reset/;
                     ${$args{TIGRdomain}}->{$descriptor} =~ /\n\r?EC\s+(.*?)\n\r?/s;
                     $EC = $1;
                  } else {
                     $EC = "";
                  }

                  if (${$args{TIGRdomain}}->{$descriptor} =~ /CC\s+/) {
                     'reset' =~ /reset/;
                     my @tmp = (${$args{TIGRdomain}}->{$descriptor} =~ /CC\s+(.*?)\n\r?/gs);
                     $CC = join (' ',@tmp);
                     undef @tmp;
                  } else {
                     $CC = "";
                  }
               }

               # catch quotatiom marks and '//'
               $CC =~ s/\"/\'/gs;
               $CC =~ s/\n\/\// /gs;
               $CC =~ s/\/\// /gs;

               #create boundary entry
               if ($orientation eq 'sense') {
                  $tigrfam_rba = $left_bd + ($tigrfam_rb * 3) - 1;
                  $tigrfam_lba = $left_bd + ($tigrfam_lb * 3) - 3;

                  #modify boundaries is necessary
                  if ($tigrfam_lba == $left_bd) {
                     $tigrfam_lba += 3;  #add one codon to avoid equal to start position
                  }
                  if ($tigrfam_rba == $right_bd) {
                     $tigrfam_rba -= 3;  #remove one codon to avoid equal to stop position
                  }

                  $boundary_tigrfam  = '     TIGR_match      '.$tigrfam_lba.'..'.$tigrfam_rba."\n";
               } elsif ($orientation eq 'antisense') {
                  my $length = $tigrfam_rb - $tigrfam_lb;
                  $tigrfam_rba = $right_bd - ($tigrfam_lb * 3) + 3;
                  $tigrfam_lba = $tigrfam_rba - ($length * 3) - 2;

                  #modify boundaries is necessary
                  if ($tigrfam_lba == $left_bd) {
                     $tigrfam_lba += 3;  #add one codon to avoid equal to stop position
                  }
                  if ($tigrfam_rba == $right_bd) {
                     $tigrfam_rba -= 3;  #remove one codon to avoid equal to start position
                  }

                  $boundary_tigrfam  = '     TIGR_match      complement('.$tigrfam_lba.'..'.$tigrfam_rba."\)\n";
               }

               #create TIGRfam label
               #iteration over all GO role links;
               #at this point only GO role numbers will be incorporated into the label
               $value = '';
               if ($GS =~ /\S+/) {
                  $value = '/label='.$GS.'; ';#.${$args{TIGRrole_GO}}->{$AC}->{$entry};
               } else {
                  $value = '/label='.$model.'; ';#.${$args{TIGRrole_GO}}->{$AC}->{$entry};
               }
               $tmp = '';
               foreach my $entry (keys %{${$args{TIGRrole_GO}}->{$AC}}) {
                  unless (defined ${$args{TIGRrole_GO}}->{$AC}->{$entry}) {
                     ${$args{TIGRrole_GO}}->{$AC}->{$entry} = ' ';
                  }
                  $tmp = ${$args{TIGRrole_GO}}->{$AC}->{$entry};
                  $tmp =~ s/GO\:\s*$//;
                  $tmp =~ s/;\s+$//;
                  $tmp =~ s/\//_/g;
                  #$tmp =~ s/^_label/\/label/;
                  $value .= $tmp.' ';
               }
               $value =~ s/\;$//;

               open  ENTRY, '>', \$label;
               write ENTRY;
               close ENTRY;

               #create TIGRfam product
               #keep iteration in - will slow downa bit, but if required in future it is already nearly fully build
               #only last role will be parsed
               $value = '/product="'.$AC;
               foreach my $entry (keys %{${$args{TIGRrole_number}}->{$AC}}) {
                  my $tmp = ' (Role:'.${$args{TIGRrole_number}}->{$AC}->{$entry}.'): ';

                  #no role defined?
                  if ($tmp =~ /\(Role:\)/) {
                     $tmp =~ s/\(Role:\)//;
                  }
                  $value .= $tmp;
               }
               #add separator if no Roles defined
               unless ($value =~ /Role\:/) {
                  $value .= ': ';
               }

               if ($DE =~ /\S+/) {
                  $value .= $DE.' ';
               } else {
                  $value .= $EN.' ';
               }
               $value   .= '; Score= '.$score.' Expect= '.$i_evalue.'"'."\n";
               open  ENTRY, '>', \$product;
               write ENTRY;
               close ENTRY;

               #create cumulativ evalue note
               $value = '/note="Cumulative evalue:  '.$evalue.'"'."\n";
               open  ENTRY, '>', \$note1;
               write ENTRY;
               close ENTRY;

               #create TIGRfam note
               $value = '/note="Role:';
               foreach my $entry (keys %{${$args{TIGRrole_number}}->{$AC}}) {
                  $value .= ${$args{TIGRrole_number}}->{$AC}->{$entry}.' = '.
                            ${$args{TIGRrole_name}}->{${$args{TIGRrole_number}}->{$AC}->{$entry}};
               }

               #add verbose description if selected
               if (${$args{auto_ini_ref}}{tigrfam_verbose} =~ /verbose/ && $value =~ /Role\:\s*\w+/) {
                  $value .= '; '.$CC;
               } else {
                  $value .= $CC;
               }
               $value .= '"'."\n";

               #nothing defined?
               if ($value =~ /note\=\"Role: \= ; \"/ || $value !~ /[A-Za-z]+/) {
                  $value = '/note="Role: No role defined"'."\n";
               } elsif ($value =~ /note\=\"Role: \= ;\s*\w+/) {
                  $value =~ s/note\=\"Role: \= ;/note\=\"Role: /s;
               }
               open  ENTRY, '>', \$note;
               write ENTRY;
               close ENTRY;

               #create TIGRfam EC number if defined
               if ($EC =~ /\S+/) {
                  $EC_number = '                     /EC_number="'.$EC.'"'."\n";
               } else {
                  $EC_number = '';
               }

               #create TIGRfam colour
               $colour = '                     /colour=4'."\n";

               #create GO from interpro entry
               foreach my $entry (keys %{${$args{TIGRrole_GO}}->{$AC}}) {
                  if (defined ${$args{TIGRrole_GO}}->{$AC}->{$entry} && ${$args{TIGRrole_GO}}->{$AC}->{$entry} =~ /\S+/) {
                     my ($GO_name, $GO_def);
                     #parse GO entry
                     'reset' =~ /reset/;
                     ${$args{TIGRrole_GoClass}}->{${$args{TIGRrole_GO}}->{$AC}->{$entry}} =~ m/\<name\>(.*?)\<\/name\>/s;
                     $GO_name = $1;
                     $GO_name =~ s/\n\r?//gs;
                     $GO_name =~ s/\s+/ /gs;

                     ${$args{TIGRrole_GoClass}}->{${$args{TIGRrole_GO}}->{$AC}->{$entry}} =~ m/\<defstr\>(.*?)\<\/defstr\>/s;
                     $GO_def = $1;
                     $GO_def =~ s/\n\r?//gs;
                     $GO_def =~ s/\s+/ /gs;
                     $GO_def =~ s/\"/\'/gs;

                     $value = '/go_from_interpro="'.${$args{TIGRrole_GO}}->{$AC}->{$entry}.', '.$GO_name.'; '.$GO_def.'"';
                     open  ENTRY, '>', \$go_interpro;
                     write ENTRY;
                     close ENTRY;
                  } else {
                     $go_interpro = '';
                  }
               }

               #combine feature
               $feature_tigrfam = $boundary_tigrfam.
                                  $label.
                                  $EC_number.
                                  $product.
                                  $note1.
                                  $note.
                                  $go_interpro.
                                  $colour;
               $key_tigrfam     = $tigrfam_lba.'_'.$tigrfam_rba.'_TIGR_match_'.$id; #maintain original uniqe ID instead of ORFnumber

               ${$args{genbank_ref}}{$key_tigrfam} = $feature_tigrfam;
               push (@{$args{feature_list_ref}}, $key_tigrfam);

               if ($evalue < $best_evalue) {
                  $best_evalue = $evalue;
                  $annotation = $feature_tigrfam;
               }
            }
         }
      }

      #add gene name, CDS name and EC number if so selected, if TIGRfam has been found
      if (${$args{auto_ini_ref}}{TIGR_annotation} == 1 && $annotation =~ /\w+/) {
         #get respective gene and CDS features
         my @features = grep (/$left_bd\_$right_bd\_(gene|CDS)\_/, @feature_list);

         #get TIGRfam info
         my ($label, $EC, $descriptor) = &tigr_annotation(main_window  => $args{main_window},
                                                          progress_bar => $args{progress_bar},
                                                          auto_ini_ref => $args{auto_ini_ref},
                                                          ini_ref      => $args{ini_ref},
                                                          annotation   => $annotation,
                                                         );
         #no results? next
         next if ($label !~ /\w+/ && $EC !~ /\w+/ && $descriptor !~ /\w+/);

         #update annotation
         foreach my $feature (@features) {
            #update gene/CDS annotation?
            if (${$args{auto_ini_ref}}{TIGR_ECCDS} == 1) {
               #first gene
               if ($feature =~ /gene/ && $label =~ /\w+/ && $args{genbank_ref}{$feature} =~ /\/gene/) {
                  my ($id, $annotation);
                  'reset' =~ m/reset/;
                  $args{genbank_ref}{$feature} =~ m/(\/gene\=\")[^\"]*?(_?\s*\d*\s*\")/s;
                  $id = $2;
                  unless (defined $id) {$id = '"'};
                  $value = '/gene="'.$label.$id;
                  open  ENTRY, '>', \$annotation;
                  write ENTRY;
                  close ENTRY;
                  $args{genbank_ref}{$feature} =~ s/\n\s*?\/gene\=\"[^\"]*?_?\s*\d*\s*\"\s*\n/\n$annotation/s;
               } elsif ($feature =~ /gene/ && $label =~ /\w+/ && $args{genbank_ref}{$feature} =~ /\/note/) {
                  my ($id, $annotation);
                  'reset' =~ m/reset/;
                  $args{genbank_ref}{$feature} =~ m/(\/note\=\")[^\"]*?(_?\s*\d*\s*\")/s;
                  $id = $2;
                  unless (defined $id) {$id = '"'};
                  $value = '/note="'.$label.$id;
                  open  ENTRY, '>', \$annotation;
                  write ENTRY;
                  close ENTRY;
                  $args{genbank_ref}{$feature} =~ s/\n\s*\/note\=\"[^\"].*?_?\s*\d*\s*\"\s*\n/\n$annotation/s;
               }
               #then CDS
               if ($feature =~ /CDS/ && $descriptor =~ /\w+/ && $args{genbank_ref}{$feature} =~ m/\/gene/) {
                  if ($label =~ /\w+/) {$label = ', '.$label};
                  my ($id, $annotation);
                  'reset' =~ m/reset/;
                  $args{genbank_ref}{$feature} =~ m/(\/gene\=\")[^\"]*?(_?\s*\d*\s*\")/s;
                  $id = $2;
                  unless (defined $id) {$id = '"'};
                  $value = '/gene="'.$descriptor.$label.$id;
                  $value =~ s/ , /, /gs;
                  open  ENTRY, '>', \$annotation;
                  write ENTRY;
                  close ENTRY;
                  $args{genbank_ref}{$feature} =~ s/\n\s*\/gene\=\"[^\"]*?_?\s*\d*\s*\"\s*\n/\n$annotation/s;
               } elsif ($feature =~ /CDS/ && $descriptor =~ /\w+/ && $args{genbank_ref}{$feature} =~ m/\/product/) {
                  if ($label =~ /\w+/) {$label = ', '.$label};
                  my ($id, $annotation);
                  'reset' =~ m/reset/;
                  $args{genbank_ref}{$feature} =~ m/(\/product\=\")[^\"]*?(_?\s*\d*\s*\")/s;
                  $id = $2;
                  unless (defined $id) {$id = '"'};
                  $value = '/product="'.$descriptor.$label.$id;
                  $value =~ s/ , /, /gs;
                  open  ENTRY, '>', \$annotation;
                  write ENTRY;
                  close ENTRY;
                  $args{genbank_ref}{$feature} =~ s/\n\s*\/product\=\"[^\"]*?_?\s*\d*\s*\"\s*\n/\n$annotation/s;
               } elsif ($feature =~ /CDS/ && $descriptor =~ /\w+/ && $args{genbank_ref}{$feature} =~ m/\/note/) {
                  if ($label =~ /\w+/) {$label = ', '.$label};
                  my ($id, $annotation);
                  'reset' =~ m/reset/;
                  $args{genbank_ref}{$feature} =~ m/(\/note\=\")[^\"]*?(_?\s*\d*\s*\")/s;
                  $id = $2;
                  unless (defined $id) {$id = '"'};
                  $value = '/note="'.$descriptor.$label.$id;
                  $value =~ s/ , /, /gs;
                  open  ENTRY, '>', \$annotation;
                  write ENTRY;
                  close ENTRY;
                  $args{genbank_ref}{$feature} =~ s/\n\s*\/note\=\"[^\"]*?_?\s*\d*\s*\"\s*\n/\n$annotation/s;
               }
            }

            #generate ECnumber entries
            my @EC = split /\s+/, $EC; #if there are multiple EC numbers, put each in a separate entry
            $EC    = '';
            foreach my $entry (@EC) {
               next unless ($entry =~ /\S+/);
               $EC .= '                     /EC_number='.$entry."\n";
            }
            $EC =~ s/^                     \/EC_number\=//;
            $EC =~ s/\n$//s;

            #delete existing EC number entries if selected
            #if ($args{genbank_ref}{$feature} =~ /\/EC_number/ && $EC =~ /\w+/ && ${$args{auto_ini_ref}}{TIGR_ECdiscard} == 1) {
            #   $args{genbank_ref}{$feature} =~ s/\/EC_number.*?(\s+?\/)/$1/gs;
            #}

            #add EC number
            if ($feature                     =~ /CDS/ &&
                $EC                          =~ /\w+/ &&
                $args{genbank_ref}{$feature} =~ /\/EC_number\=/ &&
                ${$args{auto_ini_ref}}{TIGR_ECdiscard} == 0) {
               next;
            } elsif ($feature                     =~ /CDS/ &&
                     $EC                          =~ /\w+/ &&
                     $args{genbank_ref}{$feature} =~ /\/(gene|product|note)/ &&
                     $args{genbank_ref}{$feature} !~ /\/EC_number\=/ &&
                     ${$args{auto_ini_ref}}{TIGR_ECdiscard} == 0) {
               if ($args{genbank_ref}{$feature} =~ /\/gene/) {
                  $args{genbank_ref}{$feature} =~ s/(\/gene\=\"[^\"]*?\")/$1\n                     \/EC_number\=$EC/s;
               } elsif ($args{genbank_ref}{$feature} =~ /\/product/) {
                  $args{genbank_ref}{$feature} =~ s/(\/product\=\"[^\"]*?\")/$1\n                     \/EC_number\=$EC/s;
               } elsif ($args{genbank_ref}{$feature} =~ /\/note/) {
                  $args{genbank_ref}{$feature} =~ s/(\/note\=\"[^\"]*?\")/$1\n                     \/EC_number\=$EC/s;
               }
            } elsif ($feature                     =~ /CDS/ &&
                     $EC                          =~ /\w+/ &&
                     $args{genbank_ref}{$feature} !~ /\/(gene|product|note)/ &&
                     $args{genbank_ref}{$feature} !~ /\/EC_number\=/ &&
                     ${$args{auto_ini_ref}}{TIGR_ECdiscard} == 0) {
               $args{genbank_ref}{$feature} =~ s/(                     \/)/                     \/EC_number\=$EC\n$1/s;
            } elsif ($feature                     =~ /CDS/ &&
                     $EC                          =~ /\w+/ &&
                     $args{genbank_ref}{$feature} =~ /\/EC_number/ &&
                     ${$args{auto_ini_ref}}{TIGR_ECdiscard} == 1) {
               $args{genbank_ref}{$feature} =~ s/(\/EC_number\=\"?.*?\"?)\n/\/EC_number\=$EC/s;
            } elsif ($feature                     =~ /CDS/ &&
                     $EC                          =~ /\w+/ &&
                     $args{genbank_ref}{$feature} !~ /\/EC_number/ &&
                     $args{genbank_ref}{$feature} =~ /\/(gene|product|note)/ &&
                     ${$args{auto_ini_ref}}{TIGR_ECdiscard} == 1) {
               if ($args{genbank_ref}{$feature} =~ /\/gene/) {
                  $args{genbank_ref}{$feature} =~ s/(\/gene\=\"[^\"]*?\")/$1\n                     \/EC_number\=$EC/s;
               } elsif ($args{genbank_ref}{$feature} =~ /\/product/) {
                  $args{genbank_ref}{$feature} =~ s/(\/product\=\"[^\"]*?\")/$1\n                     \/EC_number\=$EC/s;
               } elsif ($args{genbank_ref}{$feature} =~ /\/note/) {
                  $args{genbank_ref}{$feature} =~ s/(\/note\=\"[^\"]*?\")/$1\n                     \/EC_number\=$EC/s;
               }
            }  elsif ($feature                     =~ /CDS/ &&
                     $EC                          =~ /\w+/ &&
                     $args{genbank_ref}{$feature} !~ /\/EC_number/ &&
                     $args{genbank_ref}{$feature} !~ /\/(gene|product|note)/ &&
                     ${$args{auto_ini_ref}}{TIGR_ECdiscard} == 1) {
               $args{genbank_ref}{$feature} =~ s/(                     \/)/                     \/EC_number\=$EC\n$1/s;
            }
         }
      }
   }
   undef @feature_list;
   &hide_pbar_2;
   return;
}


#read TIGRfam infos into memory
sub translate_tigrfam {
   my %args = @_;
   my (@infofiles, @temp_dir_array);
   my $counter;

   ##################
   #reading TIGRinfo
   ${$args{progress_bar}}->configure(-label=>"Reading TIGRinfo information");
   ${$args{main_window}}->update;
   opendir INFODIR, ${$args{ini_ref}}{TIGRfam_info} or do {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Could not open TIGR info directory ${$args{ini_ref}}{TIGRfam_info}.",
                                        -icon    => 'error',
                                        -type    => 'OK');

      return (0);
   };
   @infofiles = grep /\.INFO$/i, readdir(INFODIR);
   closedir (INFODIR);

   #create hashes for each TIGRfam
   foreach my $entry (@infofiles) {
      my ($file_ref, $domain);

      #read TIGRFam codes
      $file_ref = &slurp(main_window  => $args{main_window},
                         directory    => ${$args{ini_ref}}{TIGRfam_info},
                         filename     => $entry
                        );

      'reset' =~ /reset/;
      ${$file_ref} =~ m/^ID\s+(.+)\n\r?/;
      $domain = $1;

      unless (defined $domain) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not TIGR info file $entry.",
                                           -icon    => 'error',
                                           -type    => 'OK');

         return (0);
      }
      $domain =~ s/\s//g;
      ${$args{TIGRdomain}}->{$domain} = ${$file_ref};
      undef $file_ref;
   }
   undef @infofiles;

   #########################
   #reading TIGR role links
   ${$args{progress_bar}}->configure(-label=>"Reading TIGRfam role link information");
   ${$args{main_window}}->update;

   open READ, ${$args{ini_ref}}{TIGRfam_db_path}.'/TIGRFAMS_ROLE_LINK' or do {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Could not open TIGRfam role link file in directory ${$args{ini_ref}}{TIGRfam_db_path}.",
                                        -icon    => 'error',
                                        -type    => 'OK');

      return (0);
   };
   #create hashes for each TIGRfam-role
   $counter = 0;
   while (<READ>) {
      if ($_ =~ /^TIGR/) {
         $counter++;
         my @tmp = ();
         @tmp = split /\s+/;
         ${$args{TIGRrole_number}}->{$tmp[0]}->{$counter} = $tmp[1];
      }
   }
   close READ;

   #########################
   #reading TIGR role names
   ${$args{progress_bar}}->configure(-label=>"Reading TIGRfam role name information");
   ${$args{main_window}}->update;
   open READ, ${$args{ini_ref}}{TIGRfam_db_path}.'/TIGR_ROLE_NAMES' or do {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Could not open TIGRfam role names file in directory ${$args{ini_ref}}{TIGRfam_db_path}.",
                                        -icon    => 'error',
                                        -type    => 'OK');

      return (0);
   };
   while (<READ>) {
      if ($_ =~ /^role_id\:/) {
         my @tmp = ();
         @tmp = split /\t+/;
         chomp $tmp[3];
         unless (defined (${$args{TIGRrole_name}}->{$tmp[1]})) {${$args{TIGRrole_name}}->{$tmp[1]} = ""};
         ${$args{TIGRrole_name}}->{$tmp[1]} = ${$args{TIGRrole_name}}->{$tmp[1]}."\, Subrole\:".$tmp[3];
         ${$args{TIGRrole_name}}->{$tmp[1]} =~ s/^\, Subrole\://;
      }
   }
   close READ;

   ########################
   #read TGIRfam-GO classification
   ${$args{progress_bar}}->configure(-label=>"Reading TIGRfam-GO link information");
   ${$args{main_window}}->update;
   open READ, ${$args{ini_ref}}{TIGRfam_db_path}.'/TIGRFAMS_GO_LINK' or do {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Could not open TIGRfam-GO link file in directory ${$args{ini_ref}}{TIGRfam_db_path}.",
                                        -icon    => 'error',
                                        -type    => 'OK');

      return (0);
   };
   $counter = 0;
   while (<READ>) {
      if ($_ =~ /^TIGR/) {
         $counter++;
         my @tmp = ();
         @tmp = split /\s+/;
         $tmp[1] =~ s/\"/\'/gs;
         ${$args{TIGRrole_GO}}->{$tmp[0]}->{$counter} = $tmp[1];
      }
   }
   close READ;

   ########################
   #read GO classification
   ${$args{progress_bar}}->configure(-label=>"Reading GO classification");
   ${$args{main_window}}->update;
   opendir INFODIR, ${$args{ini_ref}}{TIGRfam_db_path} or do {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Could not open directory ${$args{ini_ref}}{TIGRfam_db_path}.",
                                        -icon    => 'error',
                                        -type    => 'OK');

      return (0);
   };
   #check if it is a full build, NOT a daily one
   if (-e ${$args{ini_ref}}{TIGRfam_db_path}.'/go_daily-termdb.obo-xml') {
      ${$args{main_window}}->messageBox(-title   => 'Warn',
                                        -message => "At least one daily GO build was detected. Make sure to use a full build instead.",
                                        -icon    => 'error',
                                        -type    => 'OK');
   }

   @infofiles = ();
   @infofiles = grep /^go_\d+\-termdb\.obo\-xml$/, readdir(INFODIR);
   #were there any files?
   if (scalar @infofiles < 1) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Could not find a valid GO classification file in directory ${$args{ini_ref}}{TIGRfam_db_path}.\nGO terms will not be available",
                                        -icon    => 'error',
                                        -type    => 'OK');
      return;
   }

   #sort array
   my @sorted =
         map  $_->[0] =>
         sort { $a->[1] <=> $b->[1] }
         map  [ $_, m/^go_(\d+)\-termdb/ ]
         => @infofiles;

   #use latest file
   my $file_ref = &slurp(main_window  => $args{main_window},
                         directory    => ${$args{ini_ref}}{TIGRfam_db_path},
                         filename     => $sorted[-1]
                        );
   undef @infofiles;
   undef @sorted;

   my @tmp = ();
   @tmp = split /\<term\>/,${$file_ref};
   foreach (@tmp) {
      my $GoClass;
      next if ($_ =~ /\n\<obo\>\n/);
      'reset' =~ /reset/;
      if ($_ =~ m/\<id\>(GO\:\d+)\<\/id\>/i) {
         $GoClass = $1;
         $_ =~ s/\"/\'/gs;
         ${$args{TIGRrole_GoClass}}->{$GoClass} = $_;
      }
   }
   return;
}

sub tigr_annotation {
   my %args = @_;
   my ($label, $EC, $descriptor);

   #get label
   'reset' =~ m/reset/;
   $args{annotation} =~ m/\/label\=(.*?);/s;
   $label = $1;
   unless (defined $label) {$label = ''};

   #parse label
   'reset' =~ m/reset/;
   if    ($label =~ m/^([a-z]{3,4}[A-Z])(_|\b)/)  {$label= $1} #regular gene name
   elsif ($label =~ m/^([A-Z][a-z]{2,3}[A-Z])(_|\b)/)  {$label= $1}
   elsif ($label =~ m/^([a-z]{3,5})(_|\b)/)  {$label= $1}
   elsif ($label =~ m/(_|\b)([a-z]{3,4}[A-Z])$/) {$label= $2} #gene name at end of expression
   elsif ($label =~ m/(_|\b)([A-Z][a-z]{2,3}[A-Z])$/) {$label= $2}
   elsif ($label =~ m/(_|\b)([a-z]{3,5})$/) {$label= $2}
   else  {$label = ''}; #no valid name

   #parse EC number
   'reset' =~ m/reset/;
   $args{annotation} =~ m/\/EC_number\=\"?(.*?)\"?\n\r?/s;
   $EC = $1;
   unless (defined $EC) {
      'reset' =~ m/reset/;
      $args{annotation} =~ m/\/label\=.+?;\s*EC\s*([^\/]+?)\//s;
      $EC = $1;
      $EC =~ s/\s+/ /gs;
      $EC =~ s/\; GO.+$//;
   }
   unless (defined $EC) {$EC = ''};

   #parse descriptor
   'reset' =~ m/reset/;
   $args{annotation} =~ m/\/product\=\".+?:\s*(.*?);\s+(Score|e-value|Expect)/s;
   $descriptor = $1;
   unless (defined $descriptor) {$descriptor = ''};
   $descriptor =~ s/^\d+\)\:\s*//;
   $descriptor =~ s/\s+/ /gs;
   $descriptor =~ s/\-\s+/\-/g;
   $descriptor =~ s/\Q$label\E\s*$//;

   return ($label, $EC, $descriptor);
}






1;

