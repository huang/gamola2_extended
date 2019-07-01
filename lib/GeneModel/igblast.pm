#!/opt/ActivePerl-5.8/bin/perl

#run intergenic Blast on exisiting gene model
#input arguments: progress_bar, directory, filename, ini_ref, auto_ini_ref


package GeneModel::igblast;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&run_igblast);
use vars qw(@additional_orfs);

use Basics::MesgBox;
use File::Copy;
use Cwd;
use Basics::progress_bar         qw(:DEFAULT);
use initialise::read_me          qw(:DEFAULT);
use GeneModel::read_genemodel    qw(:DEFAULT);
use ProgrammeModules::blast      qw(:DEFAULT);
use ProgrammeModules::sequence   qw(:DEFAULT);


#local variables
my (%args, $tl, @additional_orfs);
#our @additional_orfs : shared = ();


sub run_igblast {
   my %args = @_;
   my ($orf_list_ref, $sense_ref, $antisense_ref, $status, $header, $dna_seq,
       );
   my (@sense_orfs, @antisense_orfs, @remove_small_orfs);

   #clear previous IG results if not re-used
   if (${$args{auto_ini_ref}}{reuseIGresults} == 0) {
      ${$args{progress_bar}}->configure(-label=>"Clearing IG Blast Result Directory");
      ${$args{main_window}}->update;
      chdir ${$args{ini_ref}}{ig_results};
      unlink <*>;
      chdir ${$args{auto_ini_ref}}{work_dir};
   }

   #set up progress bar
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'IG Blast analysis',
                   label        => 'IG Blast analysis'
                  );
   &show_pbar_2;

   #read existing gene model
   ${$args{progress_bar}}->configure(-label=>"Reading gene model");
   ${$args{main_window}}->update;
   ### add option if .combined fle does not exist -> file selector, then rename file to .combined to make compatible
   my @types = (["Gene Model Files", '.mod', 'TEXT'], ["All Files", "*"] );

   unless (-e ${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.combined') {
      my $model = ${$args{main_window}}->getOpenFile(-title      => 'Select gene model file',
                                                     -filetypes  => \@types,
                                                     -initialdir => ${$args{ini_ref}}{genemodel_output}
                                                    );
      unless (defined $model && -e $model) {return(0)};
      copy($model, ${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.combined');
   }

   ($orf_list_ref) = &read_internal_gene_model (main_window  => $args{main_window},
                                                progress_bar => $args{progress_bar},
                                                auto_ini_ref => $args{auto_ini_ref},
                                                ini_ref      => $args{ini_ref},
                                                directory    => $args{directory},
                                                input_file   => $args{input_file},
                                               );
   if ($orf_list_ref eq 0) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Error while reading internal gene model for $args{input_file}. Aborting analysis",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error generating Intergenic Blast gene model for file $args{input_file}".
                     "\nError while reading internal gene model for $args{input_file}. Aborting analysis\n\n";
      close ERRORLOG;
      &hide_pbar_2;
      return (0); #return to main
   }

   #reading input file
   ${$args{progress_bar}}->configure(-label=>"Reading source file $args{filename}");
   ${$args{main_window}}->update;

   my ($file_ref) = &slurp(main_window => \$args{progress_bar},
                           directory   => $args{directory},
                           filename    => $args{input_file}
                          );
   #parsing input file
   if (${$file_ref} =~ /^\>/) {
      ${$file_ref}  =~ m/^(\>[^\n]+)\n(.+)/s;
      $header       = $1;
      $dna_seq      = $2;
   } else {
      $header  = '>'.$args{filename};
      $dna_seq = ${$file_ref};
   }
   $dna_seq =~ s/\s//gs;

   #setup for smart or full ig regions
   $status = &igmodel (main_window           => $args{main_window},
                       progress_bar          => $args{progress_bar},
                       auto_ini_ref          => $args{auto_ini_ref},
                       ini_ref               => $args{ini_ref},
                       orf_list_ref          => $orf_list_ref,
                       remove_small_orfs_ref => \@remove_small_orfs,
                       input_file            => $args{input_file},
                       dna_length            => length($dna_seq),
                       dna_seq_ref           => \$dna_seq,
                       header_ref            => \$header,
                      );
   if ($status eq '0') {
      &hide_pbar_2;
      return (0);
   }

   #parse additional ORFs into a proper gene model
   ${$args{progress_bar}}->configure(-label=>"Adding ". ($#additional_orfs + 1) . " new ORFs to gene model");
   ${$args{main_window}}->update;
   $status = &add_additional_orfs(main_window           => $args{main_window},
                                  auto_ini_ref          => $args{auto_ini_ref},
                                  ini_ref               => $args{ini_ref},
                                  progress_bar          => $args{progress_bar},
                                  input_file            => $args{input_file},
                                  orf_list              => $orf_list_ref,
                                  remove_small_orfs_ref => \@remove_small_orfs,
                                  dna_length            => length($dna_seq),
                                  dna_ref               => \$dna_seq,
                                  header                => \$header
                                 );
   if ($status eq '0') {
      &hide_pbar_2;
      return (0);
   }

   #clear results if checked
   if (${$args{auto_ini_ref}}{reuseIGresults} == 0 && ${$args{auto_ini_ref}}{reuse_results} == 0) {
      my @inputfile = ();
      my $curdir = getcwd();
      ${$args{progress_bar}}->configure(-label=>"Cleaning up IG Blast result folder");
      ${$args{main_window}}->update;
      chdir ${$args{ini_ref}}{ig_results};
      unlink <*>;
      chdir $curdir;
   }
   undef @additional_orfs;
   undef @sense_orfs;
   undef @antisense_orfs;
   return (1);

}


sub igmodel {
   my %args = @_;
   my ($ig_ref, $sense_ref, $antisense_ref);
   my (@sense_orfs, @antisense_orfs);

   #setup for smart or full
   if (${$args{auto_ini_ref}}{igblastsmart} == 1) {
      ${$args{progress_bar}}->configure(-label=>"Grabbing intergenic regions");
      ${$args{main_window}} ->update;

      $ig_ref     = &determine_igs (main_window           => $args{main_window},
                                    auto_ini_ref          => $args{auto_ini_ref},
                                    ini_ref               => $args{ini_ref},
                                    progress_bar          => $args{progress_bar},
                                    input_file            => $args{input_file},
                                    dna_length            => $args{dna_length},
                                    orf_list_ref          => $args{orf_list_ref},
                                    remove_small_orfs_ref => $args{remove_small_orfs_ref},
                                   );
      #multithreading blast
      ${$args{progress_bar}}->configure(-label=>"Identifying potential IG ORFs for $args{input_file}");
      ${$args{main_window}} ->update;

      &multi_blast_ig (main_window   => $args{main_window},
                       progress_bar  => $args{progress_bar},
                       auto_ini_ref  => $args{auto_ini_ref},
                       ini_ref       => $args{ini_ref},
                       input_file    => $args{input_file},
                       ig            => $ig_ref,
                       dna_ref       => $args{dna_seq_ref},
                       dna_length    => $args{dna_length},
                       header        => $args{header_ref},
                       orientation   => 'sense',
                       orf_list_ref  => $args{orf_list_ref},
                      );
      &multi_blast_ig (main_window   => $args{main_window},
                       progress_bar  => $args{progress_bar},
                       auto_ini_ref  => $args{auto_ini_ref},
                       ini_ref       => $args{ini_ref},
                       input_file    => $args{input_file},
                       ig            => $ig_ref,
                       dna_ref       => $args{dna_seq_ref},
                       dna_length    => $args{dna_length},
                       header        => $args{header_ref},
                       orientation   => 'antisense',
                       orf_list_ref  => $args{orf_list_ref},
                      );

   } elsif (${$args{auto_ini_ref}}{igblastsmart} == 0) {
      ${$args{progress_bar}}->configure(-label=>"Grabbing sense and antisene models");
      ${$args{main_window}} ->update;

      #grab all sense ORFs
      @sense_orfs = grep (/___sense$/, @{$args{orf_list_ref}});

      #grab all antisense ORFs
      @antisense_orfs = grep (/___antisense$/, @{$args{orf_list_ref}});

      #determine start-stop of intergenic regions
      $sense_ref     = &determine_igs (main_window           => $args{main_window},
                                       progress_bar          => $args{progress_bar},
                                       auto_ini_ref          => $args{auto_ini_ref},
                                       ini_ref               => $args{ini_ref},
                                       input_file            => $args{input_file},
                                       dna_length            => $args{dna_length},
                                       orf_list_ref          => \@sense_orfs,
                                       remove_small_orfs_ref => $args{remove_small_orfs_ref},
                                      );
      $antisense_ref = &determine_igs (main_window           => $args{main_window},
                                       progress_bar          => $args{progress_bar},
                                       auto_ini_ref          => $args{auto_ini_ref},
                                       ini_ref               => $args{ini_ref},
                                       input_file            => $args{input_file},
                                       dna_length            => $args{dna_length},
                                       orf_list_ref          => \@antisense_orfs,
                                       remove_small_orfs_ref => $args{remove_small_orfs_ref},
                                      );
      undef @sense_orfs;
      undef @antisense_orfs;

      #multithreading blast
      ${$args{progress_bar}}->configure(-label=>"Identifying potential strand specific IG ORFs for $args{input_file}");
      ${$args{main_window}}->update;
      &multi_blast_ig (main_window   => $args{main_window},
                       progress_bar  => $args{progress_bar},
                       auto_ini_ref  => $args{auto_ini_ref},
                       ini_ref       => $args{ini_ref},
                       input_file    => $args{input_file},
                       ig            => $sense_ref,
                       dna_ref       => $args{dna_seq_ref},
                       dna_length    => $args{dna_length},
                       header        => $args{header_ref},
                       orientation   => 'sense',
                       orf_list_ref  => $args{orf_list_ref},
                      );
      &multi_blast_ig (main_window   => $args{main_window},
                       progress_bar  => $args{progress_bar},
                       auto_ini_ref  => $args{auto_ini_ref},
                       ini_ref       => $args{ini_ref},
                       input_file    => $args{input_file},
                       ig            => $antisense_ref,
                       dna_ref       => $args{dna_seq_ref},
                       dna_length    => $args{dna_length},
                       header        => $args{header_ref},
                       orientation   => 'antisense',
                       orf_list_ref  => $args{orf_list_ref},
                      );
   }

   #parse IG results
   if (${$args{auto_ini_ref}}{igcluster} > 0) {
      &parse_clustered_IG(main_window   => $args{main_window},
                          progress_bar  => $args{progress_bar},
                          auto_ini_ref  => $args{auto_ini_ref},
                          ini_ref       => $args{ini_ref},
                          input_file    => $args{input_file},
                          dna_length    => $args{dna_length},
                         );
   } else {
      &parse_IG(main_window   => $args{main_window},
                progress_bar  => $args{progress_bar},
                auto_ini_ref  => $args{auto_ini_ref},
                ini_ref       => $args{ini_ref},
                input_file    => $args{input_file},
                dna_length    => $args{dna_length},
               );
   }
   return(1);
}

sub multi_blast_ig {
   my %args = @_;
   my ($original_db, $igresultlist, $strand_start);
   my (@ig_results, @orf_subset);
   my $i = 1;

   undef $strand_start;

   #grab orf subset, depending on sense, antisense or smart selection
   if (${$args{auto_ini_ref}}{igblastsmart} == 0) {
      if ($args{orientation} eq 'sense') {
         if ($#{$args{orf_list_ref}} >= 0) {
            @orf_subset = grep (/___sense$/, @{$args{orf_list_ref}});
         } else {
            push (@orf_subset, '___1___1___3___sense'); #add artificial ORF if there are none
         }
      } elsif ($args{orientation} eq 'antisense') {
         if ($#{$args{orf_list_ref}} >= 0) {
            @orf_subset = grep (/___antisense$/, @{$args{orf_list_ref}});
         } else {
            push (@orf_subset, '___1___1___3___antisense'); #add artificial ORF if there are none
         }
      }

      #delete internal ORFs
      {
         #convert to hash
         my %orfs = ();
         foreach my $entry (@orf_subset) {
            my ($id, $left_bd, $right_bd, $orientation, $frame);
            'reset' =~ /reset/;
            $entry =~ m/(\d+)___(\d+)___(\d+)___((anti)?sense)$/;
            $id          = $1;
            $left_bd     = $2;
            $right_bd    = $3;
            $orientation = $4;
            $orfs{$id} = $left_bd.'_'.$right_bd.'_'.$orientation;
         }
         #delete internal ORFs
         &delete_internal_orfs(main_window   => $args{main_window},
                               progress_bar  => $args{progress_bar},
                               auto_ini_ref  => $args{auto_ini_ref},
                               ini_ref       => $args{ini_ref},
                               orfs          => \%orfs
                              );
         #reconsitute ORF model
         foreach my $id (sort {$a <=> $b} keys %orfs) {
            my ($left_bd, $right_bd, $orientation);
            'reset'      =~ m/reset/;
            $orfs{$id}   =~ m/(\d+)_(\d+)_(.+)/;
            $left_bd     = $1;
            $right_bd    = $2;
            $orientation = $3;
            push (@orf_subset, '___'.$id.'___'.$left_bd.'___'.$right_bd.'___'.$orientation);
         }
      }
   }

   #flatten ORF list for smart search
   if (${$args{auto_ini_ref}}{igblastsmart} == 1) {
      my %orfs = ();
      undef (@orf_subset);
      if ($#{$args{orf_list_ref}} >= 0) {
         foreach my $entry (@{$args{orf_list_ref}}) {
            my ($id, $left_bd, $right_bd, $orientation, $frame);
            'reset' =~ /reset/;
            $entry =~ m/(\d+)___(\d+)___(\d+)___((anti)?sense)$/;
            $id          = $1;
            $left_bd     = $2;
            $right_bd    = $3;
            $orientation = $4;

            unless (defined $left_bd && defined $right_bd) {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Problems parsing internal gene model file $args{input_file} at line\n$entry",
                                                             -buttons => ['OK'],
                                                             -bitmap  => 'info');
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error generating Intergenic Blast gene model for file $args{input_file}".
                              "\nProblems parsing internal gene model file $args{input_file} at line\n$entry\n\n";
               close ERRORLOG;
               next;
            }
            $orfs{$id} = $left_bd.'_'.$right_bd.'_'.$orientation;
         }
      } else {
         $orfs{'1'} = '1_3_sense';
      }

      #delete internal ORFs
      &delete_internal_orfs(main_window   => $args{main_window},
                            progress_bar  => $args{progress_bar},
                            auto_ini_ref  => $args{auto_ini_ref},
                            ini_ref       => $args{ini_ref},
                            orfs          => \%orfs
                           );

      #merge ORFs if overlap
      #&merge_overlap_orfs(main_window   => $args{main_window},
      #                    progress_bar  => $args{progress_bar},
      #                    auto_ini_ref  => $args{auto_ini_ref},
      #                    ini_ref       => $args{ini_ref},
      #                    orfs          => \%orfs
      #                   );

      #reconstiture flattened ORF model
      foreach my $id (sort {$a <=> $b} keys %orfs) {
         my ($left_bd, $right_bd, $orientation);
         'reset'      =~ m/reset/;
         $orfs{$id}   =~ m/(\d+)_(\d+)_(.+)/;
         $left_bd     = $1;
         $right_bd    = $2;
         $orientation = $3;
         push (@orf_subset, '___'.$id.'___'.$left_bd.'___'.$right_bd.'___'.$orientation);
      }

      undef %orfs;
   }

   #if superstition is allowed, remove ORFs < 200nt from ORIGINAL gene model
   #this has to be done again, as original ORF model is used (not related to IGs search)
   if (${$args{auto_ini_ref}}{IGsuperstition} == 1) {
      &delete_small_orfs(main_window           => $args{main_window},
                         progress_bar          => $args{progress_bar},
                         auto_ini_ref          => $args{auto_ini_ref},
                         ini_ref               => $args{ini_ref},
                         orfs                  => \@orf_subset,
                         remove_small_orfs_ref => $args{remove_small_orfs_ref},
                        );
   }

   #generating frame - start position hash
   foreach my $entry (@orf_subset) {
      my ($left_bd, $right_bd, $frame, $orientation, $orf, $N_term, $C_term);
      'reset' =~ m/reset/;
      $entry =~ m/___\d+___(\d+)___(\d+)___((anti)?sense)$/;
      $left_bd     = $1;
      $right_bd    = $2;
      $orientation = $3;

      #skip artificial ORF is added above
      next if ($left_bd == 1 && $right_bd == 3);

      #get nt seq
      $orf = substr(${$args{dna_ref}}, ($left_bd - 1), ($right_bd - $left_bd + 1));
      #antisense?

      if ($orientation eq 'antisense') {
         $orf = &orient_sequence(seq         => $orf,
                                 orientation => 'antisense'
                                );
      }
      #translate to aa sequence
      ($orf) = &nt2aa_unchecked(
                                sequence     => $orf
                               );
      #reduce size to IG flanking
      $N_term = substr($orf, 0, (int(${$args{auto_ini_ref}}{igflanking} / 3)));
      $N_term =~ s/^.//;#remove start codon
      $C_term = substr($orf, (0 - (int(${$args{auto_ini_ref}}{igflanking} / 3))));
      $C_term =~ s/.$//;#remove stop codon
      if ($orientation eq 'sense') {
         $frame = &frame(orientation => $orientation,
                         left_bd     => $left_bd
                         );
         $strand_start->{$frame}->{$left_bd}->{'N'}  = $N_term;
         $strand_start->{$frame}->{$left_bd}->{'C'} = $C_term;
      } else {
         $frame = &frame(orientation => $orientation,
                         left_bd     => $right_bd
                         );
         $strand_start->{$frame}->{$right_bd}->{'N'} = $N_term;
         $strand_start->{$frame}->{$right_bd}->{'C'}  = $C_term;
      }
   }
   undef @orf_subset;

   #total number of instances
   my $total_regions = $#{$args{ig}} + 2;

   #read existing IG Blast results for input file
   opendir SEQINPUT, ${$args{ini_ref}}{ig_results};
   @ig_results = grep /^$args{input_file}/, readdir(SEQINPUT);
   foreach (@ig_results) {
      $igresultlist->{$_} = 1;
   }
   closedir SEQINPUT;
   undef @ig_results;

   #clustered IG analysis?
   if (${$args{auto_ini_ref}}{igcluster} > 0) {
      &clustered_ig_blast ( main_window    => $args{main_window},
                            progress_bar   => $args{progress_bar},
                            instance       => $i,
                            all_ig_regions => $args{ig},
                            total_regions  => $total_regions,
                            dna_ref        => $args{dna_ref},
                            dna_length     => $args{dna_length},
                            header         => $args{header},
                            orientation    => $args{orientation},
                            auto_ini_ref   => $args{auto_ini_ref},
                            ini_ref        => $args{ini_ref},
                            input_file     => $args{input_file},
                            igresultlist   => $igresultlist,
                            strand_start   => $strand_start,
                           );
      $i++;
   }
   #run normal IG
   else {
      foreach my $ig (@{$args{ig}}) {
         ${$args{progress_bar}}->configure(-label=>"Identifying potential ORFs in $args{input_file} section $ig.");
         ${$args{main_window}}->update;
         &ig_blast ( main_window   => $args{main_window},
                     progress_bar  => $args{progress_bar},
                     instance      => $i,
                     ig_region     => $ig,
                     total_regions => $total_regions,
                     dna_ref       => $args{dna_ref},
                     dna_length    => $args{dna_length},
                     header        => $args{header},
                     orientation   => $args{orientation},
                     auto_ini_ref  => $args{auto_ini_ref},
                     ini_ref       => $args{ini_ref},
                     input_file    => $args{input_file},
                     igresultlist  => $igresultlist,
                     strand_start  => $strand_start,
                    );
         $i++;
      }
   }
   return;
}

sub ig_blast {
   my %args = @_;
   my ($orf_def, $curdir, $blastresult_ref, $left_bd, $right_bd, $dna_fragment, $min_match_length);
   my (@blast_list, @ig_results, @strand_orfs);

   #determine minimum length for string match (ORF on same frame as IG); also converts from nt to aa length
   #make min length smaller by 4 to count for required truncations in mathc later on
   if (${$args{auto_ini_ref}}{igflanking} < ${$args{auto_ini_ref}}{igblastminlength}) {
      $min_match_length = int(${$args{auto_ini_ref}}{igflanking} / 3) - 4;
   } else {
      $min_match_length = int(${$args{auto_ini_ref}}{igblastminlength} / 3) - 4;
   }
   if ($min_match_length < 0) {$min_match_length = 0;}

   #get current nt seq plus igflanking nt overlap on both sides
   'reset' =~ m/reset/;
   $args{ig_region} =~ m/(\d+)_(\d+)/;
   $left_bd = $1;
   $right_bd = $2;

   unless (defined $left_bd && defined $right_bd && $left_bd =~ /^\d+$/ && $right_bd =~ /^\d+$/) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No Left or right boundary for $args{ig_region} in file $args{input_file}",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error generating Intergenic Blast gene model for file $args{input_file}".
                     "\nNo Left or right boundary for $args{ig_region} in file $args{input_file}\n\n";
      close ERRORLOG;
      push (@additional_orfs, $args{ig_region}.'___ERROR___PARSING');
      return;
   }

   #expand IG region by chosen flanking nt size
   if ($left_bd > ${$args{auto_ini_ref}}{igflanking}) {
      $left_bd  -= ${$args{auto_ini_ref}}{igflanking};
   }
   if ($right_bd < ($args{dna_length} - ${$args{auto_ini_ref}}{igflanking})) {
      $right_bd += ${$args{auto_ini_ref}}{igflanking};
   }

   #get IG nt seq
   $dna_fragment = substr(${$args{dna_ref}}, ($left_bd - 1), ($right_bd - $left_bd + 1));

   #antisense?
   if ($args{orientation} eq 'antisense') {
      $dna_fragment = &orient_sequence(seq         => $dna_fragment,
                                       orientation => 'antisense'
                                      );
   }

   #find putative open reading frames
   my $temp = $dna_fragment;
   for (my $frame = 0; $frame <= 2; $frame++) {#frames
      my @temp = ($temp =~ m/(\w\w\w)/g);
      my $stop_counter = 0;
      my ($string, $first, $position, $nt_index, $arrayindex, $size_of_array);

      #define size of array
      $size_of_array = scalar(@temp);

      $nt_index   = 0;
      $arrayindex = 1;
      foreach my $codon (@temp) {
         #not 3 codons? skip
         next unless ($codon =~ m/\w\w\w/);

         $string .= $codon;
         $arrayindex++;

         if ($codon =~ m/^(tag|tga|taa)$/i || $arrayindex >= $size_of_array) { #NO STOP CODON AT
            #remove stop codon again
            $string =~ s/$codon$//;

            my $skip_ig = 0;
            #add IG region if string is long enough
            if ($args{orientation} eq 'sense' && length($string) >= ${$args{auto_ini_ref}}{igblastminlength}) {
               my ($aa_seq, $aa_ref);
               #determine current position
               $position = $left_bd + $nt_index - length($string) + $frame;

               #format into fasta
               my $ig_frame = &frame(orientation => $args{orientation},
                                     left_bd     => $position
                                     );

               #translate to aa sequence
               ($aa_seq) = &nt2aa_unchecked(
                                            sequence     => $string
                                           );

               #remove last aa in case it would prevent matching to ORF
               chop ($aa_seq);

               #does the IG aa contain an ORF seq? if so, same frame and same orf, skip
               foreach my $ORF_bd (keys %{$args{strand_start}->{$ig_frame} }) {
                  #truncate if IG contains ORF aa seq
                  if (exists $args{strand_start}->{$ig_frame}->{$ORF_bd}->{'N'}) {
                     #build substring matches
                     my $match = join '|' => map { substr $aa_seq, -$_} reverse ($min_match_length .. length($aa_seq));

                     if ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'N'} =~ m/.?$match$/) {
                        my ($substring) = ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'N'} =~ m/.?$match$/g);
                        $aa_seq =~ s/.?$substring.?//i;
                        $string = substr($string, 0, (length($aa_seq) * 3));
                        undef $substring;
                     }

                     if (length($string) < ${$args{auto_ini_ref}}{igblastminlength}) {
                        $string = '';
                        $skip_ig = 1;
                        last;
                     }
                  }
                  if (exists $args{strand_start}->{$ig_frame}->{$ORF_bd}->{'C'}) {
                     #build substring matches
                     my $match = join '|' => map { substr $aa_seq, 0, $_} reverse ($min_match_length .. length($aa_seq));
                     if ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'C'} =~ m/.?$match$/) {
                        my ($substring) = ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'C'} =~ m/.?$match?/g);
                        $aa_seq =~ s/.?$substring.?//i;
                        $string = substr($string, 0, (length($aa_seq) * 3));
                        undef $substring;
                     }

                     if (length($string) < ${$args{auto_ini_ref}}{igblastminlength}) {
                        $string = '';
                        $skip_ig = 1;
                        last;
                     }
                  }
               }
               next if ($skip_ig == 1);

               push (@blast_list, $position.'_'.$string.'_'.$ig_frame);
               $aa_seq   = '';
               $aa_ref   = '';
               $position = '';
            } elsif ($args{orientation} eq 'antisense' && length($string) >= ${$args{auto_ini_ref}}{igblastminlength}) {

               my ($aa_seq, $aa_ref);
               #determine current position
               $position = $left_bd + length($dna_fragment) - $nt_index - $frame + 3; #left bd, in sense and antisense direction; '+3' to remove putative stop codon

               my $ig_frame = &frame(orientation => $args{orientation},
                                     left_bd     => $position + 2 #end of codon
                                     );

               #translate to aa sequence
               ($aa_seq) = &nt2aa_unchecked(
                                            sequence     => $string
                                           );
               #remove last aa in case it would prevent matching to ORF
               chop ($aa_seq);

               #does this co-incide with same frame ORF?
               foreach my $ORF_bd (keys %{$args{strand_start}->{$ig_frame} }) {
                  #truncate if IG contains ORF aa seq
                  if (exists $args{strand_start}->{$ig_frame}->{$ORF_bd}->{'N'}) {
                     #build substring matches
                     my $match = join '|' => map { substr $aa_seq, -$_} reverse ($min_match_length .. length($aa_seq));
                     if ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'N'} =~ m/.?$match$/) {
                        my ($substring) = ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'N'} =~ m/.?$match$/g);
                        $aa_seq =~ s/.?$substring.?//i;
                        $string = substr($string, 0, (length($aa_seq) * 3));
                        undef $substring;
                     }

                     if (length($string) < ${$args{auto_ini_ref}}{igblastminlength}) {
                        $string  = '';
                        $skip_ig = 1;
                        last;
                     }
                  }
                  if (exists $args{strand_start}->{$ig_frame}->{$ORF_bd}->{'C'}) {
                     #build substring matches
                     my $match = join '|' => map { substr $aa_seq, 0, $_} reverse ($min_match_length .. length($aa_seq));
                     if ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'C'} =~ m/.?$match$/) {
                        my ($substring) = ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'C'} =~ m/.?$match?/g);
                        $aa_seq =~ s/.?$substring.?//i;
                        $string = substr($string, 0, (length($aa_seq) * 3));
                        undef $substring;
                     }

                     if (length($string) < ${$args{auto_ini_ref}}{igblastminlength}) {
                        $string  = '';
                        $skip_ig = 1;
                        last;
                     }
                  }
               }
               next if ($skip_ig == 1);
               #reset position to correct left_bd
               $position -= 3;
               push (@blast_list, $position.'_'.$string.'_'.$ig_frame);
            }
            $string = '';
         }
         $nt_index  += 3;
      }

      #get last one if long enough
      if (length($string) > ${$args{auto_ini_ref}}{igblastminlength}) {
         my ($aa_seq, $aa_ref);
         my $ig_frame = &frame(orientation => $args{orientation},
                               left_bd     => $position
                               );
         push (@blast_list, $position.'_'.$string.'_'.$ig_frame);
      }
      undef $first;
      undef $string;
      undef $nt_index;
      undef @temp;
      #shift frame
      $temp =~ s/^.//;
   }

   #return id no hits
   if ($#blast_list < 0) {return};

   #max number of fasta entries
   my $max_count = $#blast_list + 2;
   my $label = "IG Blast analysis";

   #start multithreading IG Blast if selected
   if (${$args{auto_ini_ref}}{IG_blast_threaded} == 1) {
      for (my $i=1; $i<=$#blast_list+2; $i += ${$args{auto_ini_ref}}{IG_CPU}) {
         my $count = $i + ${$args{auto_ini_ref}}{IG_CPU};
         if ($count > $#blast_list+2) {$count = $#blast_list+2};

         &update_pbar_2(title        => 'IG Blast analysis (Threaded Blast)',
                        label        => "IG Blast analysis for $args{input_file}, instance $count out of $max_count from $args{instance} out of $args{total_regions} $args{orientation} IG regions",
                        progress     => ($count / $max_count) * 100,
                       );
         ${$args{main_window}}->update;

         my @childs = ();
         for (my $j=$i; $j<$count; $j++) {
            #start forking
            my $pid = fork();
            if ($pid) {
               # parent
               push(@childs, $pid);
            } elsif ($pid == 0) {
               # child
               my ($ig_left, $ig_right, $orfseq, $frame);
               'reset' =~ m/reset/;
               $blast_list[$j-1] =~ m/(\d+)_(.+)_([\-\d]+)/;
               $ig_left          = $1;
               $orfseq           = $2;
               $frame            = $3;
               $orfseq           = lc $orfseq;

               #exit if not parsed correctly
               unless (defined $orfseq) {
                  print "\nERROR in $blast_list[$j-1] ...\n\n";
                  CORE::exit();
               }

               #rectify start position 1
               if ($ig_left < 1) {$ig_left = 1};

               #define right boundary
               $ig_right = $ig_left + length ($orfseq) - 1;

               #check for existing files is necessary
               if (${$args{auto_ini_ref}}{reuseIGresults} == 1) {
                  my $file = $args{input_file}.'_'.$ig_left.'_'.$frame.'.blast';
                  if (exists $args{igresultlist}->{$file}) {
                     #exit safely
                     #print "\nFile:".$args{input_file}.'_'.$ig_left.'_'.$frame.'.blast'.' already exists';
                     CORE::exit();
                  }
               }

               #send to blast
               my $seq    = '>'.$ig_left.'_'.$ig_right."\n".$orfseq;
               my $aa_seq = '';
               ($blastresult_ref, $aa_seq) = &blast_ig (
                                               auto_ini_ref => $args{auto_ini_ref},
                                               ini_ref      => $args{ini_ref},
                                               sequence_nt  => $seq
                                              );
               unless (defined $blastresult_ref) {
                  print "\n-------------\nError in IG Blast. Likely cause is an incorrect Blast db path.\n-------------\n\n";
                  CORE::exit;
               }

               if ($$blastresult_ref eq '0') {
                  $$blastresult_ref = 'No Blast Hits found';
               }

               #write result file
               open WRITE, "+>".${$args{ini_ref}}{ig_results}.'/'.$args{input_file}.'_'.$ig_left.'_'.$frame.'.blast' or do {
                  print "\nCouldn't create IG Blast result file ". $args{input_file}.'_'.$ig_left.'_'.$frame.'.blast' ."...\n\n";
                  CORE::exit();
               };
               print WRITE "Input sequence:\t$orfseq\n";
               print WRITE "Translated aa seq:\t$aa_seq\n";
               print WRITE "Orientation:\t$args{orientation}\n";
               print WRITE "Frame:\t$frame\n";
               print WRITE "IG region:\t$ig_left".'_'.$ig_right."\n";
               print WRITE $$blastresult_ref;
               close WRITE;
               #exit safely
               CORE::exit();
            } else {
               die "Couldn\'t fork\: $!\n";
            }
         }
         #wait
         foreach (@childs) {
            waitpid($_, 0);
         }
      }
   } elsif (${$args{auto_ini_ref}}{IG_blast_cluster} == 1) {
      my $progress = 0;
      #iterate over array
      foreach my $entry (@blast_list) {
         my ($ig_left, $ig_right, $orfseq, $frame);
         ${$args{progress_bar}}->configure(-label=>" ");
         ${$args{main_window}}->update;
         #update progress bar
         &update_pbar_2(title        => 'IG Blast analysis (clustered Blast)',
                        label        => "Blasting $args{input_file}, $progress of $max_count",
                        progress     => ($progress / $max_count) * 100,
                       );
         $progress++;

         'reset'  =~ m/reset/;
         $entry   =~ m/(\d+)_(.+)_([\-\d]+)/;
         $ig_left = $1;
         $orfseq  = $2;
         $frame   = $3;
         $orfseq  = lc $orfseq;

         #ignore if not parsed correctly
         unless (defined $orfseq) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error in clustered IG-Blast',
                                                          -text    => "Could not parse from entry $entry\.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error in clustered IG-Blast".
                           "\nCould not parse from entry $entry\.\n\n";
            close ERRORLOG;
            next;
         }

         #rectify start position 1
         if ($ig_left < 1) {$ig_left = 1};

         #define right boundary
         $ig_right = $ig_left + length ($orfseq) - 1;

         #define IG region
         my $ig_region = $ig_left.'_'.$ig_right;

         #check if IG region is on same frame as an ORF and if so if overlap warrants analysis.
         foreach my $orf_left_bd (keys %{$args{strand_start}->{$frame}}) {
            my ($orf_right_bd);
            'reset' =~ m/reset/;
            $args{strand_start}->{$frame}->{$orf_left_bd} =~ m/\d+___\d+___(\d+)___/;
            $orf_right_bd = $1;

            #if ORF is embedded in other ORF, skip
            if (($ig_left >= $orf_left_bd && $ig_right <= $orf_right_bd) ||
                ($ig_left <  $orf_left_bd && $ig_right >  $orf_right_bd)) {
               #if Blast exists, delete (not needed, created by mistake)
               if (-e ${$args{ini_ref}}{ig_results}.'/'.$args{input_file}.'_'.$ig_left.'_'.$frame.'.blast') {
                  unlink ${$args{ini_ref}}{ig_results}.'/'.$args{input_file}.'_'.$ig_left.'_'.$frame.'.blast';
               }
               next;
            }

            #if ORF overlaps, check if remaining seq warrants IG blast, adjust sequence as necessary
            if (($ig_left < $orf_left_bd  && $ig_right > $orf_left_bd ) ||
                ($ig_left < $orf_right_bd && $ig_right > $orf_right_bd)) {
               #ignore IG if true IG is less than minimum IG length
               if ((($ig_left      - $orf_left_bd) < ${$args{auto_ini_ref}}{igblastminlength}) ||
                   (($orf_right_bd - $ig_right)    < ${$args{auto_ini_ref}}{igblastminlength}))  {
                  #if Blast exists, delete
                  if (-e ${$args{ini_ref}}{ig_results}.'/'.$args{input_file}.'_'.$ig_left.'_'.$frame.'.blast') {
                     unlink ${$args{ini_ref}}{ig_results}.'/'.$args{input_file}.'_'.$ig_left.'_'.$frame.'.blast';
                  }
                  next;
               }
               #re-define IG region
               if ($ig_right > $orf_left_bd) {
                  $ig_right = $orf_left_bd - 3;
                  #create subset of submitted sequence
                  $orfseq = substr(${$args{dna_ref}}, ($ig_left - 1), ($ig_right - $ig_left + 1));

                  #invert sequence if antisense
                  if ($args{orientation} eq 'antisense') {
                     $orfseq = &orient_sequence(seq         => $orfseq,
                                                orientation => 'antisense'
                                                );
                  }
                  last;
               }
            }
         }

         #check for existing files is necessary
         if (${$args{auto_ini_ref}}{reuseIGresults} == 1) {
            my $file = $args{input_file}.'_'.$ig_left.'_'.$frame.'.blast';
            if (exists $args{igresultlist}->{$file}) {
               #exit safely
               ${$args{progress_bar}}->configure(-label=>'File:'.$args{input_file}.'_'.$ig_left.'_'.$frame.'.blast'.' already exists');
               ${$args{main_window}}->update;
               next;
            }
         }

         #send to blast
         my $seq    = '>'.$ig_left.'_'.$ig_right."\n".$orfseq;
         my $aa_seq = '';
         ($blastresult_ref, $aa_seq) = &blast_ig_clustered (
                                                   auto_ini_ref => $args{auto_ini_ref},
                                                   ini_ref      => $args{ini_ref},
                                                   sequence_nt  => $seq
                                                  );
         unless (defined $blastresult_ref) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error in clustered IG-Blast',
                                                          -text    => "Could not process  entry $entry\.\nError in IG Blast. Likely cause is an incorrect Blast db path",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error in clustered IG-Blast".
                           "\nCould not process entry $entry\.\nError in IG Blast. Likely cause is an incorrect Blast db path\n\n";
            close ERRORLOG;
            next;
         }

         if ($$blastresult_ref eq '0') {
            $$blastresult_ref = 'No Blast Hits found';
         }

         #write result file
         open WRITE, "+>".${$args{ini_ref}}{ig_results}.'/'.$args{input_file}.'_'.$ig_left.'_'.$frame.'.blast' or do {
            print "\nCouldn't create IG Blast result file ". $args{input_file}.'_'.$ig_left.'_'.$frame.'.blast' ."...\n\n";
            CORE::exit();
         };
         print WRITE "Input sequence:\t$orfseq\n";
         print WRITE "Translated aa seq:\t$aa_seq\n";
         print WRITE "Orientation:\t$args{orientation}\n";
         print WRITE "Frame:\t$frame\n";
         print WRITE "IG region:\t$ig_left".'_'.$ig_right."\n";
         print WRITE $$blastresult_ref;
         close WRITE;
      }
   }

   undef @blast_list;
   return (1);
}

sub clustered_ig_blast {
   #combine individual IG regions to be combined into a single BLAST analysis
   #saves time but is potentially less accurate
   my %args = @_;
   my ($orf_def, $curdir, $blastresult_ref, $left_bd, $right_bd, $dna_fragment,
       $cluster_counter, $clustered_entries, $min_match_length);
   my (@blast_list, @put_orfs, @ig_results, @strand_orfs, @seen_IG_job);
   my (%full_job_list);
   #initialise cluster count
   $cluster_counter   = 0;
   $clustered_entries = '';

   #grab existing results if reuse selected
   if (${$args{auto_ini_ref}}{reuseIGresults} == 1) {
      #grab all existing IG result files and get query header info
      my $max_count = keys %{$args{igresultlist} };
      my $count     = 1;
      @seen_IG_job  = ();
      foreach my $ig_file (keys %{$args{igresultlist} }) {
         ${$args{progress_bar}}->configure(-label=>"Capturing existing IG results for $ig_file.");
         ${$args{main_window}}->update;
         &update_pbar_2(title        => 'IG Blast analysis (Threaded Blast)',
                        label        => "Capturing existing IG results",
                        progress     => ($count / $max_count) * 100,
                       );
         ${$args{main_window}}->update;

         my ($file_ref) = &slurp(main_window => $args{main_window},
                                 directory   => ${$args{ini_ref}}{ig_results},
                                 filename    => $ig_file
                                );
         my @temp = ($$file_ref =~ m/\# Query:\s+([^\n]+?)\n/gs);

         push (@seen_IG_job, @temp);
         undef @temp;

         $count++;
      }
      undef $max_count;
      undef $count;
      ${$args{progress_bar}}->configure(-label=>"Analysis continues.");
      ${$args{main_window}}->update;
      &hide_pbar_2;
   }

   #determine minimum length for string match (ORF on same frame as IG); also converts from nt to aa length
   #make min length smaller by 4 to count for required truncations in mathc later on
   if (${$args{auto_ini_ref}}{igflanking} < ${$args{auto_ini_ref}}{igblastminlength}) {
      $min_match_length = int(${$args{auto_ini_ref}}{igflanking} / 3) - 4;
   } else {
      $min_match_length = int(${$args{auto_ini_ref}}{igblastminlength} / 3) - 4;
   }
   if ($min_match_length < 0) {$min_match_length = 0;}

   #combine individual IG regions until threshold is reached
   foreach my $ig (@{$args{all_ig_regions}}) {
      #update info
      ${$args{progress_bar}}->configure(-label=>"Clustering IG region $ig for $args{input_file}.");
      ${$args{main_window}}->update;

      #get current nt seq plus igflanking nt overlap on both sides
      'reset' =~ m/reset/;
      $ig =~ m/(\d+)_(\d+)/;
      $left_bd = $1;
      $right_bd = $2;

      unless (defined $left_bd && defined $right_bd && $left_bd =~ /^\d+$/ && $right_bd =~ /^\d+$/) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "No Left or right boundary for $ig in file $args{input_file}",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error generating Intergenic Blast gene model for file $args{input_file}".
                        "\nNo Left or right boundary for $ig in file $args{input_file}\n\n";
         close ERRORLOG;
         push (@additional_orfs, $args{ig_region}.'___ERROR___PARSING');
         next;
      }

      #expand IG region by chosen flanking nt size
      if ($left_bd > ${$args{auto_ini_ref}}{igflanking}) {
         $left_bd  -= ${$args{auto_ini_ref}}{igflanking};
      }
      if ($right_bd < ($args{dna_length} - ${$args{auto_ini_ref}}{igflanking})) {
         $right_bd += ${$args{auto_ini_ref}}{igflanking};
      }

      #get IG nt seq
      $dna_fragment = substr(${$args{dna_ref}}, ($left_bd - 1), ($right_bd - $left_bd + 1));

      #antisense?
      if ($args{orientation} eq 'antisense') {
         $dna_fragment = &orient_sequence(seq         => $dna_fragment,
                                          orientation => 'antisense'
                                         );
      }

      #find putative open reading frames
      my $temp = $dna_fragment;
      for (my $frame = 0; $frame <= 2; $frame++) {#frames
         my @temp = ($temp =~ m/(\w\w\w)/g);
         my $stop_counter = 0;
         my ($string, $first, $position, $nt_index, $arrayindex, $size_of_array);

         #define size of array
         $size_of_array = scalar(@temp);

         $nt_index = 0;
         $arrayindex = 1;
         foreach my $codon (@temp) {
            #not 3 codons? skip
            next unless ($codon =~ m/^\w\w\w$/);

            $string .= $codon;
            $arrayindex++;

            if ($codon =~ m/^(tag|tga|taa)$/i || $arrayindex >= $size_of_array) { #NO STOP CODON AT end
               #remove stop codon again
               $string =~ s/$codon$//;

               my $skip_ig = 0;
               #add IG region if string is long enough
               if ($args{orientation} eq 'sense' && length($string) >= ${$args{auto_ini_ref}}{igblastminlength}) {
                  my ($aa_seq, $aa_ref);
                  #determine current position
                  $position = $left_bd + $nt_index - length($string) + $frame;

                  #format into fasta
                  my $ig_frame = &frame(orientation => $args{orientation},
                                        left_bd     => $position
                                        );

                  #translate to aa sequence
                  ($aa_seq) = &nt2aa_unchecked(
                                               sequence     => $string
                                              );

                  #remove last aa in case it would prevent matching to ORF
                  chop ($aa_seq);

                  #does the IG aa contain an ORF seq? if so, same frame and same orf, skip
                  foreach my $ORF_bd (keys %{$args{strand_start}->{$ig_frame} }) {
                     #truncate if IG contains ORF aa seq
                     if (exists $args{strand_start}->{$ig_frame}->{$ORF_bd}->{'N'}) {
                        #build substring matches
                        my $match = join '|' => map { substr $aa_seq, -$_} reverse ($min_match_length .. length($aa_seq));

                        if ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'N'} =~ m/.?$match/) {
                           my ($substring) = ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'N'} =~ m/.?$match/g);
                           $aa_seq =~ s/.?$substring.?//i;
                           $string = substr($string, 0, (length($aa_seq) * 3));
                           undef $substring;
                        }

                        if (length($string) < ${$args{auto_ini_ref}}{igblastminlength}) {
                           $string = '';
                           $skip_ig = 1;
                           last;
                        }
                     }
                     if (exists $args{strand_start}->{$ig_frame}->{$ORF_bd}->{'C'}) {
                        #build substring matches
                        my $match = join '|' => map { substr $aa_seq, 0, $_} reverse ($min_match_length .. length($aa_seq));
                        if ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'C'} =~ m/.?$match/) {
                           my ($substring) = ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'C'} =~ m/.?$match/g);
                           $aa_seq =~ s/.?$substring.?//i;
                           $string = substr($string, 0, (length($aa_seq) * 3));
                           undef $substring;
                        }

                        if (length($string) < ${$args{auto_ini_ref}}{igblastminlength}) {
                           $string = '';
                           $skip_ig = 1;
                           last;
                        }
                     }
                  }
                  next if ($skip_ig == 1);

                  ($aa_ref) = &seq2fasta(header   => $ig.' Position:'.$position.' Frame:'.$ig_frame.' Size [nt]:'.length($string),
                                         sequence => $aa_seq
                                        );
                  #add entry to full job list for later check
                  $full_job_list{$ig.' Position:'.$position.' Frame:'.$ig_frame.' Size [nt]:'.length($string)} = $$aa_ref;

               } elsif ($args{orientation} eq 'antisense' && length($string) >= ${$args{auto_ini_ref}}{igblastminlength}) {

                  my ($aa_seq, $aa_ref);
                  #determine current position
                  $position = $left_bd + length($dna_fragment) - $nt_index - $frame;  #left bd, in sense and antisense direction;

                  my $ig_frame = &frame(orientation => $args{orientation},
                                        left_bd     => $position + 2 #end of codon
                                        );

                  #translate to aa sequence
                  ($aa_seq) = &nt2aa_unchecked(
                                               sequence     => $string
                                              );

                  #remove last aa in case it would prevent matching to ORF
                  chop ($aa_seq);

                  #does this co-incide with same frame ORF?
                  foreach my $ORF_bd (keys %{$args{strand_start}->{$ig_frame} }) {
                     #truncate if IG contains ORF aa seq
                     if (exists $args{strand_start}->{$ig_frame}->{$ORF_bd}->{'N'}) {

                        #build substring matches
                        my $match = join '|' => map { substr $aa_seq, -$_} reverse ($min_match_length .. length($aa_seq));
                        if ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'N'} =~ m/.?$match/) {
                           my ($substring) = ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'N'} =~ m/.?$match/g);
                           $aa_seq                =~ s/.?$substring.?//i;
                           my $length_fullstring  = length($string);
                           $string                = substr($string, 0, (length($aa_seq) * 3));
                           my $length_truncstring = length($string);
                           #shift position (eq left_bd later on)
                           $position += ($length_fullstring - $length_truncstring);
                           undef $substring;
                        }

                        if (length($string) < ${$args{auto_ini_ref}}{igblastminlength}) {
                           $string  = '';
                           $skip_ig = 1;
                           last;
                        }
                     }
                     if (exists $args{strand_start}->{$ig_frame}->{$ORF_bd}->{'C'}) {
                        #build substring matches
                        my $match = join '|' => map { substr $aa_seq, 0, $_} reverse ($min_match_length .. length($aa_seq));
                        if ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'C'} =~ m/.?$match/) {
                           my ($substring) = ($args{strand_start}->{$ig_frame}->{$ORF_bd}->{'C'} =~ m/.?$match/g);
                           $aa_seq =~ s/.?$substring.?//i;
                           $string = substr($string, 0, (length($aa_seq) * 3));
                           undef $substring;
                        }

                        if (length($string) < ${$args{auto_ini_ref}}{igblastminlength}) {
                           $string  = '';
                           $skip_ig = 1;
                           last;
                        }
                     }
                  }
                  next if ($skip_ig == 1);
                  #format into fasta
                  ($aa_ref) = &seq2fasta(header   => $ig.' Position:'.$position.' Frame:'.$ig_frame.' Size [nt]:'.length($string),
                                         sequence => $aa_seq
                                        );
                  #add entry to full job list for later check
                  $full_job_list{$ig.' Position:'.$position.' Frame:'.$ig_frame.' Size [nt]:'.length($string)} = $$aa_ref;
               }
               $string = '';
            }
            $nt_index  += 3;
         }

         #get last one if long enough
         if (length($string) > ${$args{auto_ini_ref}}{igblastminlength}) {
            my ($aa_seq, $aa_ref);
            my $ig_frame = &frame(orientation => $args{orientation},
                                  left_bd     => $position
                                  );
            #translate to aa sequence
            ($aa_seq) = &nt2aa_unchecked(
                                         sequence     => $string
                                        );

            #format into fasta
            ($aa_ref) = &seq2fasta(header   => $ig.' Position:'.$position.' Frame:'.$ig_frame.' Size [nt]:'.length($string),
                                   sequence => $aa_seq
                                  );
            #add entry to full job list for later check
            $full_job_list{$ig.' Position:'.$position.' Frame:'.$ig_frame.' Size [nt]:'.length($string)} = $$aa_ref;
         }
         undef $first;
         undef $string;
         undef $nt_index;
         undef @temp;
         #shift frame
         $temp =~ s/^.//;
      }
   }

   #remove existing IG results if reuse is selected
   if (${$args{auto_ini_ref}}{reuseIGresults} == 1) {
      #iterate through found IG regions and remove from current job list
      foreach my $entry (@seen_IG_job) {
         if (exists $full_job_list{$entry}) {
            delete $full_job_list{$entry};
         }
      }
   }

   #create final job list for clustered option (cluster size of '0' equals to individual IG regions blasts)
   foreach my $entry (sort keys %full_job_list) {
      if ($cluster_counter == ${$args{auto_ini_ref}}{igcluster}) {
         push (@put_orfs, $clustered_entries);
         $clustered_entries = '';
         $cluster_counter   = 0;
      }
      $clustered_entries .= $full_job_list{$entry}."\n";
      $cluster_counter++;
   }
   #get last one
   push (@put_orfs, $clustered_entries);
   undef %full_job_list;
   undef $clustered_entries;

   #return if no hits
   if ($#put_orfs < 0 || $put_orfs[0] !~ m/\w+/) {return};

   #max number of fasta entries
   my $max_count = $#put_orfs + 2;
   my $label = "IG Blast analysis";

   #start multithreading IG Blast if selected
   if (${$args{auto_ini_ref}}{IG_blast_threaded} == 1) {
      for (my $i=1; $i<=$#put_orfs+2; $i += ${$args{auto_ini_ref}}{IG_CPU}) {
         my $count = $i + ${$args{auto_ini_ref}}{IG_CPU};
         if ($count > $#put_orfs+2) {$count = $#put_orfs+2};

         &update_pbar_2(title        => 'Clustered IG Blast analysis (Threaded Blast)',
                        label        => "IG Blast analysis for $args{input_file} $count from $max_count in $args{orientation} direction",
                        progress     => ($count / $max_count) * 100,
                       );
         ${$args{main_window}}->update;

         my @childs = ();
         for (my $j=$i; $j<$count; $j++) {
            #start forking
            my $pid = fork();
            if ($pid) {
               # parent
               push(@childs, $pid);
            } elsif ($pid == 0) {
               # child
               my $blastresult;
               my $result_file_index = $j;
               #get current working directory
               my $curdir = "";
               $curdir = getcwd();
               #change to Blast directory
               if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
                  chdir ${$args{ini_ref}}{blast_executables};
                  $blastresult = "";
                  if (open RESULT, "-|") { # original process
                     local $/;
                     $blastresult = <RESULT>;
                  } else { # child
                     if (open STDIN, "-|") { # child
                        exec "./blastpgp -d ${$args{auto_ini_ref}}{intergenicblastdb} -m 9"; #child
                        die "Cannot exec: $!";
                     } else { # grandchild
                        print $put_orfs[$j-1];
                        CORE::exit;
                     }
                  }
               } elsif (${$args{auto_ini_ref}}{blast_plus} == 1) {
                  chdir ${$args{ini_ref}}{blast_plus_executables};
                  $blastresult = "";
                  if (open RESULT, "-|") { # original process
                     local $/;
                     $blastresult = <RESULT>;
                  } else { # child
                     if (open STDIN, "-|") { # child
                        exec "./blastp -db ${$args{auto_ini_ref}}{intergenicblastdb} -outfmt '7 evalue score'"; #child
                        die "Cannot exec: $!";
                     } else { # grandchild
                        print $put_orfs[$j-1];
                        CORE::exit;
                     }
                  }
               }

               #change back to previous working dir
               chdir $curdir;

               #valid result?
               unless (defined $blastresult && $blastresult =~ /\w+/) {
                  print "Error in Blast analysis using database ${$args{auto_ini_ref}}{intergenicblastdb}\n";
                  open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
                  print WRITE "Error in Blast analysis using database ${$args{auto_ini_ref}}{intergenicblastdb}".
                              "\nCould not parse result: $blastresult\n\n";
                  close WRITE;
                  $blastresult = "\n******Error in Blast analysis******\/\/\n";
               };
               #check for existing files is necessary
               if (${$args{auto_ini_ref}}{reuseIGresults} == 1 && -e ${$args{ini_ref}}{ig_results}.'/'.$args{input_file}.'_'.$args{orientation}.'_'.$result_file_index.'.blast') {
                  while (-e ${$args{ini_ref}}{ig_results}.'/'.$args{input_file}.'_'.$args{orientation}.'_'.$result_file_index.'.blast') {
                     $result_file_index++;
                  }
               }
               #write result file
               open WRITE, "+>".${$args{ini_ref}}{ig_results}.'/'.$args{input_file}.'_'.$args{orientation}.'_'.$result_file_index.'.blast' or do {
                  print "\nCouldn't create IG Blast result file ". $args{input_file}.'_'.$args{orientation}.'_'.$result_file_index.'.blast' ."...\n\n";
                  CORE::exit();
               };
               print WRITE $blastresult;
               close WRITE;
               #exit safely
               CORE::exit();
            } else {
               die "Couldn\'t fork\: $!\n";
            }
         }
         #wait
         foreach (@childs) {
            waitpid($_, 0);
         }
      }
   } elsif (${$args{auto_ini_ref}}{IG_blast_cluster} == 1) {
      my $progress          = 0;
      my $result_file_index = 1;
      #iterate over array
      foreach my $entry (@put_orfs) {
         #update progress bar
         &update_pbar_2(title        => 'Clustered IG Blast analysis (clustered Blast)',
                        label        => "Blasting $args{input_file}, $progress of $max_count",
                        progress     => ($progress / $max_count) * 100,
                       );
         $progress++;
         ${$args{progress_bar}}->configure(-label=>'IG Blast: '.$args{input_file}.'_'.$args{orientation}.'_'.$result_file_index.'.blast');
         ${$args{main_window}}->update;

         my $blastresult;
         #get current working directory
         my $curdir = "";
         $curdir = getcwd();

         #change to Blast directory
         if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
            chdir ${$args{ini_ref}}{blast_executables};
            $blastresult = "";
            if (open RESULT, "-|") { # original process
               local $/;
               $blastresult = <RESULT>;
            } else { # child
               if (open STDIN, "-|") { # child
                  exec "./blastpgp -a ${$args{auto_ini_ref}}{IG_CPU} -d ${$args{auto_ini_ref}}{intergenicblastdb} -m 9"; #child
                  die "Cannot exec: $!";
               } else { # grandchild
                  print $entry;
                  CORE::exit;
               }
            }
         } elsif (${$args{auto_ini_ref}}{blast_plus} == 1) {
            chdir ${$args{ini_ref}}{blast_plus_executables};
            $blastresult = "";
            if (open RESULT, "-|") { # original process
               local $/;
               $blastresult = <RESULT>;
            } else { # child
               if (open STDIN, "-|") { # child
                  exec "./blastp -num_threads ${$args{auto_ini_ref}}{IG_CPU} -db ${$args{auto_ini_ref}}{intergenicblastdb} -outfmt '7 evalue score'"; #child
                  die "Cannot exec: $!";
               } else { # grandchild
                  print $entry;
                  CORE::exit;
               }
            }
         }
         #change back to previous working dir
         chdir $curdir;

         #valid result?
         unless (defined $blastresult && $blastresult =~ /\w+/) {
            print "Error in Blast analysis using database ${$args{auto_ini_ref}}{intergenicblastdb}\n";
            open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
            print WRITE "Error in Blast analysis using database ${$args{auto_ini_ref}}{intergenicblastdb}".
                        "\nCould not parse result: $blastresult\n\n";
            close WRITE;
            $blastresult = "\n******Error in Blast analysis******\/\/\n";
         };
         #check for existing files is necessary
         if (${$args{auto_ini_ref}}{reuseIGresults} == 1 && -e ${$args{ini_ref}}{ig_results}.'/'.$args{input_file}.'_'.$args{orientation}.'_'.$result_file_index.'.blast') {
            while (-e ${$args{ini_ref}}{ig_results}.'/'.$args{input_file}.'_'.$args{orientation}.'_'.$result_file_index.'.blast') {
               $result_file_index++;
            }
         }
         #write result file
         open WRITE, "+>".${$args{ini_ref}}{ig_results}.'/'.$args{input_file}.'_'.$args{orientation}.'_'.$result_file_index.'.blast' or do {
            print "\nCouldn't create IG Blast result file ". $args{input_file}.'_'.$args{orientation}.'_'.$result_file_index.'.blast' ."...\n\n";
            CORE::exit();
         };
         print WRITE $blastresult;
         close WRITE;

      }
   }

   undef @put_orfs;
   &hide_pbar_2;
   return (1);
}

sub parse_IG {
   my %args = @_;

   my ($count);
   my (@ig_results);

   #parse all results when multithreading is done
   ${$args{progress_bar}}->configure(-label=>"Parsing IG Blast results for $args{input_file}.");
   ${$args{main_window}}->update;
   #update results folder
   ${$args{progress_bar}}->configure(-label=>"Reading IG Blast results for $args{input_file}.");
   ${$args{main_window}}->update;

   #grap IG blast results for current input file
   opendir SEQINPUT, ${$args{ini_ref}}{ig_results};
   @ig_results = grep /^$args{input_file}/, readdir(SEQINPUT);
   closedir SEQINPUT;

   my $max_count = $#ig_results + 1;
   $count = 0;
   foreach my $result (@ig_results) {
      my ($best_hit, $orfseq, $orientation, $frame, $left_bd, $right_bd, $local_position,
          $blastresult_ref);
      $count++;
      #read file
      ($blastresult_ref) = &slurp (directory => ${$args{ini_ref}}{ig_results},
                                   filename  => $result
                                  );

      #check if best Blast hit is below threshold
      'reset' =~ m/reset/;
      #$$blastresult_ref =~ m/\nSequences producing significant alignments[^\n]*?\n\n.*?\.\.\.\s+\d+\s+([^\n]*?)\n/s;
      $$blastresult_ref =~ m/\nSequences producing significant alignments[^\n]*?\n\n([^\n]*?)\n/s;
      $best_hit = $1;
      'reset' =~ m/reset/;
      $best_hit =~ s/.*?\s+([^\s]+)\s*$/$1/;

      unless (defined $best_hit) {next};

      &update_pbar_2(title        => 'IG Blast analysis',
                     label        => "Processing instance $count out of ".$#ig_results,
                     progress     => (($count / $max_count) * 100),
                    );

      if ($best_hit =~ /^e/) {$best_hit = '1'.$best_hit};
      $best_hit =~ s/\s//gsi;

      #above threshold? ignore
      next if ($best_hit > ${$args{auto_ini_ref}}{igmaxevalue});

      #parse file
      'reset' =~ m/reset/;
      $$blastresult_ref =~ m/Input sequence:\t([^\n]*?)\n.+?Orientation:\t([^\n]*?)\nFrame:\t([^\n]*?)\nIG region:\t(\d+)_(\d+)\n/s;
      $orfseq         = $1;
      $orientation    = $2;
      $frame          = $3;
      $left_bd        = $4;
      $right_bd       = $5;
      push (@additional_orfs, $left_bd.'_'.$right_bd.'_'.$best_hit.'_'.$orientation.'_'.$frame.'_'.$orfseq);
      undef $blastresult_ref;
   }

   #throw out duplicates, just in case
   my %seen = ();
   @additional_orfs = grep { ! $seen{$_} ++ } @additional_orfs;
   undef %seen;
   &hide_pbar_2;
   return (1);
}

sub parse_clustered_IG {
   my %args = @_;

   my ($count);
   my (@ig_results, @temp_orfs);

   #parse all results when multithreading is done
   ${$args{progress_bar}}->configure(-label=>"Parsing IG Blast results for $args{input_file}.");
   ${$args{main_window}}->update;

   #grap IG blast results for current input file
   opendir SEQINPUT, ${$args{ini_ref}}{ig_results};
   @ig_results = grep /^$args{input_file}/, readdir(SEQINPUT);
   closedir SEQINPUT;

   my $max_count = $#ig_results + 1;
   $count = 0;
   foreach my $result (@ig_results) {
      my ($best_hit, $orfseq, $orientation, $frame, $left_bd, $right_bd, $length, $local_position,
          $blastresult_ref);
      my (@indiviudal_IG_results);
      $count++;
      #read file
      ($blastresult_ref) = &slurp_cmd (directory => ${$args{ini_ref}}{ig_results},
                                       filename  => $result
                                      );

      #define orientation
      if ($result =~ m/_sense_\d+\.blast/) {
         $orientation = 'sense';
      } elsif ($result =~ m/_antisense_\d+\.blast/) {
         $orientation = 'antisense';
      } else {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Wrong filename construct $result for IG Blast,\nskipping file",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error generating Intergenic Blast gene model for file $args{input_file}".
                        "\nWrong filename construct $result for IG Blast,skipping file\n\n";
         close ERRORLOG;
         next;
      }

      #parse differently for legacy or Blast plus
      if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
         #Split into single entries
         @indiviudal_IG_results = split /\# BlastP/, $$blastresult_ref;
         undef $blastresult_ref;
         #iterate through individual IGs regions
         foreach my $entry (@indiviudal_IG_results) {
            next unless ($entry =~ m/\w+/);
            #skip entries with no hits
            next if ($entry !~ m/\# Fields\:.+?bit score\s*\n\S+/s);

            #get evalue
            $entry =~ m/\# Fields\:.+?bit score\s*?\n.+?\s+(\S+)\s+[\-\.\d]+\s*\n/s;
            $best_hit = $1;

            unless (defined $best_hit) {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Error parsing $entry in file $args{input_file} for evalue,\nskipping file",
                                                             -buttons => ['OK'],
                                                             -bitmap  => 'info');
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error generating Intergenic Blast gene model for file $args{input_file}".
                              "\nError parsing $entry in file $args{input_file} for evalue, skipping file\n\n";
               close ERRORLOG;
               next;
            };

            &update_pbar_2(title        => 'IG Blast analysis',
                           label        => "Processing instance $count out of ".$#ig_results,
                           progress     => (($count / $max_count) * 100),
                          );

            if ($best_hit =~ /^e/) {$best_hit = '1'.$best_hit};
            $best_hit =~ s/\s//gsi;

            #above threshold? ignore
            next if ($best_hit > ${$args{auto_ini_ref}}{igmaxevalue});

            #add IG to additional ORFs
            'reset'   =~ m/reset/;
            $entry    =~ m/Position\:(\d+) Frame\:(\-?\d+) Size \[nt\]\:(\d+)/s;
            ($left_bd, $frame, $length) = ($1, $2, $3);
            $right_bd = $left_bd + $length - 1;

            push (@temp_orfs, $left_bd.'_'.$right_bd.'_'.$best_hit.'_'.$orientation.'_'.$frame.'_'.$length);
         }
      } elsif (${$args{auto_ini_ref}}{blast_plus} == 1) {
         #Split into single entries
         @indiviudal_IG_results = ($$blastresult_ref =~ m/(\# Query:.+?\n# BLAST)/gs);
         undef $blastresult_ref;

         #iterate through individual IGs regions
         foreach my $entry (@indiviudal_IG_results) {
            next unless ($entry =~ m/\w+/);
            #skip entries with no hits
            next if ($entry =~ m/\# 0 hits found/);

            #get evalue
            $entry =~ m/\n\#\s+\d+\s+hits found\n(\S+)\s/s;
            $best_hit = $1;

            unless (defined $best_hit) {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Error parsing $entry in file $args{input_file} for evalue,\nskipping file",
                                                             -buttons => ['OK'],
                                                             -bitmap  => 'info');
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error generating Intergenic Blast gene model for file $args{input_file}".
                              "\nError parsing $entry in file $args{input_file} for evalue, skipping file\n\n";
               close ERRORLOG;
               next;
            };

            &update_pbar_2(title        => 'IG Blast analysis',
                           label        => "Processing instance $count out of ". $#ig_results,
                           progress     => (($count / $max_count) * 100),
                          );

            if ($best_hit =~ /^e/) {$best_hit = '1'.$best_hit};
            $best_hit =~ s/\s//gsi;

            #above threshold? ignore
            next if ($best_hit > ${$args{auto_ini_ref}}{igmaxevalue});

            #add IG to additional ORFs
            'reset'   =~ m/reset/;
            $entry    =~ m/Position\:(\d+) Frame\:(\-?\d+) Size \[nt\]\:(\d+)/s;
            ($left_bd, $frame, $length) = ($1, $2, $3);
            $right_bd = $left_bd + $length - 1;

            push (@temp_orfs, $left_bd.'_'.$right_bd.'_'.$best_hit.'_'.$orientation.'_'.$frame.'_'.$length);
         }
      }
   }

   #throw out duplicates, just in case
   my %seen   = ();
   @temp_orfs = grep { ! $seen{$_} ++ } @temp_orfs;

   ###internal ORFs in the same frame should have been taken care with earlier already.

   #test for embedded ORFs in same frame, keep best hit
   #my %query          = ();
   #@query{@temp_orfs} = (1) x @temp_orfs;
   #my %subject        = %query;
   #my $iterate        = 1;
   #my $iterate_max    = 0;
   #while ($iterate == 1 && $iterate_max < 3) {
   #   $iterate = 0;
   #   my %delete_ig      = ();
   #   my %delete_partner = ();
   #   foreach my $query (sort keys %query) {
   #      next if (exists $delete_partner{$query});
   #      next if (exists $delete_ig{$query}); #avoid bi-directional deletion
   #      'reset'=~ m/reset/;
   #      $query =~ m/(\d+)_(\d+)_(.+?)_\w+?_(\-?\d)_/;
   #      my ($q_left_bd, $q_right_bd, $q_hit, $q_frame) = ($1, $2, $3, $4);
   #      foreach my $subject (sort keys %subject) {
   #         next if ($query eq $subject);
   #         'reset'=~ m/reset/;
   #         $subject =~ m/(\d+)_(\d+)_(.+?)_\w+?_(\-?\d)_/;
   #         my ($s_left_bd, $s_right_bd, $s_hit, $s_frame) = ($1, $2, $3, $4);
   #         #skip unless same frame
   #         next unless ($q_frame eq $s_frame);
   #         #skip unless intersecting boundaries
   #         next if ($q_left_bd  > $s_right_bd);
   #         next if ($q_right_bd < $s_left_bd);
   #         #found overlapping orf, iterate
   #         $iterate = 1;
   #         #keep best hit
   #         if ($q_hit < $s_hit) {
   #           $delete_ig{$subject}++;
   #           $delete_partner{$query}++;
   #         } else {
   #           $delete_ig{$query}++;
   #           $delete_partner{$subject}++;
   #         }
   #         last;
   #      }
   #   }
   #   #delete internal IGs
   #   foreach my $entry (keys %delete_ig) {
   #      delete $query{$entry};
   #   }
   #   %subject = %query;
   #   $iterate_max++;
   #}

   #assign remaining IG ORFs
   #@additional_orfs = keys %query;
   @additional_orfs =  @temp_orfs;

   undef %seen;
   #undef %subject;
   #undef %query;
   undef @temp_orfs;

   &hide_pbar_2;
   return (1);
}

sub determine_igs {
   my %args = @_;
   my (@orf_boundaries, @model, %hash, $orfs_ref);
   my (%orfs, %add, %ig);
   my ($left_bd, $right_bd, @delete, @add);
   #orfs present?
   if ($#{ $args{orf_list_ref}} >= 0) {

      #if superstition is allowed, remove ORFs < 200nt from ORIGINAL gene model
      if (${$args{auto_ini_ref}}{IGsuperstition} == 1) {
         &delete_small_orfs(main_window           => $args{main_window},
                            progress_bar          => $args{progress_bar},
                            auto_ini_ref          => $args{auto_ini_ref},
                            ini_ref               => $args{ini_ref},
                            orfs                  => $args{orf_list_ref},
                            remove_small_orfs_ref => $args{remove_small_orfs_ref},
                           );
      }

      foreach my $entry (@{ $args{orf_list_ref}} ) {
         my ($id, $left_bd, $right_bd, $orientation);
         'reset' =~ /reset/;
         $entry =~ m/(\d+)___(\d+)___(\d+)___((anti)?sense)$/;
         $id          = $1;
         $left_bd     = $2;
         $right_bd    = $3;
         $orientation = $4;

         unless (defined $left_bd && defined $right_bd) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Problems parsing internal gene model file $args{input_file} at line\n$entry",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error generating Intergenic Blast gene model for file $args{input_file}".
                           "\nProblems parsing internal gene model file $args{input_file} at line\n$entry\n\n";
            close ERRORLOG;
            return (0);
         }
         $orfs{$id} = $left_bd.'_'.$right_bd;
      }
   } else {
      #no ORFs present, add artificial boundaries
      $orfs{'1'} = '1_3';
   }

   #delete internal ORFs
   &delete_internal_orfs(main_window   => $args{main_window},
                         progress_bar  => $args{progress_bar},
                         auto_ini_ref  => $args{auto_ini_ref},
                         ini_ref       => $args{ini_ref},
                         orfs          => \%orfs
                        );

   #merge ORFs if overlap
   #&merge_overlap_orfs(main_window   => $args{main_window},
   #                    progress_bar  => $args{progress_bar},
   #                    auto_ini_ref  => $args{auto_ini_ref},
   #                    ini_ref       => $args{ini_ref},
   #                    orfs          => \%orfs
   #                   );

   #creating IG regions
   ${$args{progress_bar}}->configure(-label=>"Creating intergenic regions.");
   ${$args{main_window}} ->update;

   #create suitable hash of ORF boundaries
   foreach my $id (sort {$a <=> $b} keys %orfs) {
      'reset'       =~ m/reset/;
      $orfs{$id}    =~ m/(\d+)_(\d+)/;
      my $left_bd   = $1;
      my $right_bd  = $2;
      $ig{$left_bd} = $right_bd;
   }
   undef %orfs;

   my $previous_left = 1;
   foreach my $left_bd (sort {$a <=> $b} keys %ig) {
      #if left bigger than prev_left, push IG region
      if ($left_bd > $previous_left) {
         push (@model, $previous_left.'_'.$left_bd);
      }
      #assign current right bd as new previous left
      $previous_left = $ig{$left_bd};
   }

   #catch last
   if ($previous_left < $args{dna_length}) {
      push (@model, $previous_left.'_'.$args{dna_length});
   }
   undef %ig;

   return (\@model);
}

sub add_additional_orfs {
   my %args = @_;
   my ($new_orfs, $IG_orfs, $orf_check, %stop);

   #parse existing ORFs if '.combined' exists
   if (-e ${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.combined') {
      ${$args{progress_bar}}->configure(-label=>"Reading existing gene model for file $args{input_file}");
      ${$args{main_window}}->update;
      open READ, '<'.${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.combined' or die "\nCouldn't open gene model file";
      while (<READ>) {
         if ($_ !~ /\w+/) {next}; #skip empty lines; left and right_bd are inverted for antisense!
         m/^\s*\d+\s+(\d+)\s+(\d+)\s+(.+)/;
         my $left_bd = $1;
         my $right_bd = $2;
         my $rest = $3;

         #test if there are ORFs <200nt to be removed first, IF there is an IG region
         my $skip = 0;
         if (${$args{auto_ini_ref}}{IGsuperstition} == 1) {

            foreach my $small_ORF (@{$args{remove_small_orfs_ref}}) {
               my ($small_left_bd, $small_right_bd, $small_orientation);
               last if ($skip == 1);
               'reset'            =~ /reset/;
               $small_ORF         =~ m/___\d+___(\d+)___(\d+)___((anti)?sense)$/;
               $small_left_bd     = $1;
               $small_right_bd    = $2;
               $small_orientation = $3;

               #both in the same orientation
               next unless (($rest =~ m/\+/ && $small_orientation eq 'sense') ||
                            ($rest =~ m/\-/ && $small_orientation eq 'antisense'));

               #get the right ORF
               if ($small_orientation =~ m/anti/) {
                  next unless ($left_bd == $small_right_bd && $right_bd == ($small_left_bd + 3));
               } else {
                  next unless ($left_bd == $small_left_bd && $right_bd == ($small_right_bd - 3));
               }

               #test if there is a larger overlapping IG ORF
               foreach my $ig_orf (@additional_orfs) {
                  last if ($skip == 1);
                  'reset'            =~ /reset/;
                  $ig_orf            =~ /(\d+)_(\d+)_[^\_]*?_(sense|antisense)_/;
                  my $ig_left        = $1;
                  my $ig_right       = $2;
                  my $ig_orientation = $3;

                  #is IG_orf larger than regular ORF and in same orientation?
                  if ((($ig_orientation eq 'sense'     && $small_orientation eq 'sense')     && $ig_left <= $left_bd  && $ig_right >= $right_bd) ||
                      (($ig_orientation eq 'antisense' && $small_orientation eq 'antisense') && $ig_left <= $right_bd && $ig_right >= $left_bd )) {
                     $skip = 1;
                  }
               }
            }
         }

         #remove smaller regular ORF from gene model
         next if ($skip == 1);

         $orf_check->{$left_bd}->{$right_bd} = 1;
         if ($left_bd < $right_bd) {
            $stop{$right_bd} = 'sense';
         } elsif ($left_bd > $right_bd) {
            $stop{$left_bd} = 'antisense';
         }

         my $id = $left_bd;
         while (exists ($new_orfs->{$id})) {$id++};
         $new_orfs->{$id} = $left_bd.'  '.$right_bd.'  '.$rest;
      }
      close READ;
   }

   #throw out any duplicates that may have arisen from spanning contigs and IGs
   my %seen   = ();
   @additional_orfs = grep { ! $seen{$_} ++ } @additional_orfs;

   #parse through additional ORFs
   foreach my $model (@additional_orfs) {
      my ($ig_left, $ig_right, $best_hit, $orientation, $frame, $orfseq);
      if ($model =~ /___ERROR___PARSING/) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while processing $model in file $args{input_file}.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error generating Intergenic Blast gene model for file $args{input_file}".
                        "\nError while processing $model in file $args{input_file}.\n\n";
         close ERRORLOG;
         next;
      }

      'reset' =~ /reset/;
      $model =~ /(\d+)_(\d+)_([^\_]*?)_(sense|antisense)_(\-?\d)_(.+)/;
      $ig_left        = $1;
      $ig_right       = $2;
      $best_hit       = $3;
      $orientation    = $4;
      $frame          = $5;
      $orfseq         = $6;

      unless (defined $ig_left && defined $ig_right && defined $orientation && defined $frame && defined $orfseq) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while processing $model in file $args{input_file}.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error generating Intergenic Blast gene model for file $args{input_file}".
                        "\nError while processing $model in file $args{input_file}.\n\n";
         close ERRORLOG;
         next;
      }

      #remove possible artifical ORFs added if there was no existing gene model for the input sequence
      #these artificail ORFs are only 3nt long and can easily be identified
      next if (abs($ig_left - $ig_right) <= 3);

      #if add stop codons
      if ($orientation eq 'antisense') {
         $ig_left -= 3;
      } else {
         $ig_right += 3;
      }

      #orf extends over sequence end?
      while ($ig_right > ($args{dna_length} - 3)) {
         $ig_right -= 3;
      }

      #ORF start after sequence ends? Shouldn't happen but who knows? skip
      next if ($ig_left >= $args{dna_length});

      #get IG nt seq
      my $dna_fragment = substr(${$args{dna_ref}}, ($ig_left - 1), ($ig_right - $ig_left + 1));

      #antisense?
      if ($orientation eq 'antisense') {
         $dna_fragment = &orient_sequence(seq         => $dna_fragment,
                                          orientation => 'antisense'
                                         );
      }
      my ($aa_seq) = &nt2aa_unchecked(
                                      sequence     => $dna_fragment
                                      );

      #shift start/stop in case sequence goes beyond stop codon
      if ($aa_seq =~ m/\*\w$/i) {
         if ($orientation eq 'sense') {
            $ig_left  -=3;
            $ig_right -=3;
         } else {
            $ig_left  +=3;
            $ig_right +=3;
         }
      }

      my $id = $ig_left + 1;
      while (exists ($new_orfs->{$id})) {$id++};
      if ($orientation eq 'antisense') {
         next if (exists $orf_check->{$ig_right}->{$ig_left});
         #adjustment for Gamola
         if (${$args{auto_ini_ref}}{igcluster} > 0) {
            $new_orfs->{$id} = $ig_right. '  ' . ($ig_left + 3) . '  ['. $frame.' L='.$orfseq.' r='.$best_hit.'] IG';
            $IG_orfs->{$id}  = $ig_right. '  ' . ($ig_left + 3) . '  ['. $frame.' L='.$orfseq.' r='.$best_hit.'] IG';
         } else {
            $new_orfs->{$id} = $ig_right. '  ' . ($ig_left + 3) . '  ['. $frame.' L='.length($orfseq).' r='.$best_hit.'] IG';
            $IG_orfs->{$id}  = $ig_right. '  ' . ($ig_left + 3) . '  ['. $frame.' L='.length($orfseq).' r='.$best_hit.'] IG';
         }
      } elsif ($orientation eq 'sense') {
         next if (exists $orf_check->{$ig_left}->{$ig_right});
         #adjustment for Gamola
         if (${$args{auto_ini_ref}}{igcluster} > 0) {
            $new_orfs->{$id} = $ig_left . '  ' . ($ig_right - 3) . '  [+'. $frame.' L='.$orfseq.' r='.$best_hit.'] IG';
            $IG_orfs->{$id}  = $ig_left . '  ' . ($ig_right - 3) . '  [+'. $frame.' L='.$orfseq.' r='.$best_hit.'] IG';
         } else {
            $new_orfs->{$id} = $ig_left . '  ' . ($ig_right - 3) . '  [+'. $frame.' L='.length($orfseq).' r='.$best_hit.'] IG';
            $IG_orfs->{$id}  = $ig_left . '  ' . ($ig_right - 3) . '  [+'. $frame.' L='.length($orfseq).' r='.$best_hit.'] IG';
         }
      } else {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while processing $model in file $args{input_file}. No orientation defined",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error generating Intergenic Blast gene model for file $args{input_file}".
                        "\nError while processing $model in file $args{input_file}. No orientation defined\n\n";
         close ERRORLOG;
         next;
      }
   }

   #create new IG model file and add IG ORFs to combined gene model file
   open IGWRITE, '+>'.${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.IG_model' or die "\nCouldn't create new gene model file\n";
   open WRITE, '+>'.${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.combined' or die "\nCouldn't create new gene model file\n";
   my $counter = 1;
   foreach my $id (sort {$a <=> $b} keys %{$new_orfs}) {
      print WRITE   '  '.$counter.'  '.$new_orfs->{$id}."\n";
      $counter++;
   }
   $counter = 1;
   foreach my $id (sort {$a <=> $b} keys %{$IG_orfs}) {
      print IGWRITE '  '.$counter.'  '.$new_orfs->{$id}."\n";
      $counter++;
   }
   undef $counter;
   close WRITE;
   close IGWRITE;

   #update results folder
   ${$args{progress_bar}}->configure(-label=>"Included IG regions into gene model.");
   ${$args{main_window}}->update;

   #clear memory
   undef $orf_check;
   undef $new_orfs;
   undef %stop;

   return (1);
}

sub initial_position {
   my %args = @_;
   my ($ig_initial_pos, $ig_left, $ig_right, $ig_rest, $counter);
   if ($args{orientation} eq 'sense') {
      ${$args{dna_ref}} =~ m/(.*?)$args{orfseq}/;
      $ig_initial_pos = length ($1);
      #repeat until correct position has been found
      my $temp_seq = ${$args{dna_ref}};
      my $adjuster = 0;
      $counter = 0;
      while (($ig_initial_pos + 2) < ($args{left_bd} + $args{local_position} - $adjuster)) {
         $counter++;
         last if ($counter == 1000);
         $temp_seq =~ s/^(.*?$args{orfseq})(.+)/$2/;
         $adjuster += length($1);
         'reset' =~ m/reset/;
         $temp_seq =~ m/(.*?)$args{orfseq}/;
         $ig_initial_pos = length ($1);
      }
      undef $temp_seq;
      $ig_left  = $ig_initial_pos + $adjuster + 1;
      $ig_right = $ig_initial_pos + $adjuster + length($args{orfseq});
      $ig_rest  = '[+' . $args{frame} . ' L=' . length($args{orfseq}) . " r=$args{best_hit}\] ";
      #adjust initial position to true value
      $ig_initial_pos = $ig_left;
   } elsif ($args{orientation} eq 'antisense') {
      my $rev_seq = &orient_sequence (seq         => $args{orfseq},
                                      orientation => 'antisense'
                                     );

      ${$args{dna_ref}} =~ m/(.*?)$rev_seq/;
      $ig_initial_pos = length ($1);

      #repeat until correct position has been found
      my $temp_seq = ${$args{dna_ref}};
      my $adjuster = 0;
      $counter = 0;
      while ( ($ig_initial_pos + 2) < ($args{left_bd} + $args{local_position} - $adjuster)) {
         $counter++;
         last if ($counter == 1000);
         $temp_seq =~ s/(^.*?$rev_seq)(.+)/$2/;
         $adjuster += length($1);
         'reset' =~ m/reset/;
         $temp_seq =~ m/(.*?)$rev_seq/;
         $ig_initial_pos = length ($1);
      }
      undef $temp_seq;
      $ig_left  = $ig_initial_pos + $adjuster + 1;
      $ig_right = $ig_initial_pos + $adjuster + length($rev_seq);
      $ig_rest  = '[-' . $args{frame} . ' L=' . length($rev_seq) . " r=$args{best_hit}\] ";
      #adjust initial position to true value
      $ig_initial_pos = $ig_left;
   }

   if ($counter >= 1000) {
      return (-1, -1, -1, -1);
   } else {
      return ($ig_initial_pos, $ig_left, $ig_right, $ig_rest);
   }
}

sub frame {
   my %args = @_;
   my $frame;
   if ($args{orientation} eq 'sense') {
      if (($args{left_bd} + 2) % 3 == 0) { #END OF CODON
         $frame = 1;
      } elsif (($args{left_bd} + 2 - 1) % 3 == 0) {
         $frame = 2;
      } elsif (($args{left_bd} + 2 - 2) % 3 == 0) {
         $frame = 3;
      }
   } elsif ($args{orientation} eq 'antisense') {
      if (($args{left_bd}) % 3 == 0) { #START OF CODON
         $frame = -1;
      } elsif (($args{left_bd} - 1) % 3 == 0) {
         $frame = -2;
      } elsif (($args{left_bd} - 2) % 3 == 0) {
         $frame = -3;
      }
   }
   return ($frame);
}

sub delete_small_orfs {
   my %args = @_;
   my (@keep);

   ${$args{progress_bar}}->configure(-label=>"Deleting small ORFs ");
   ${$args{main_window}} ->update;

   #parse through ORFs and identify small ORFs
   foreach my $entry (@{$args{orfs}}) {
      my ($left_bd, $right_bd);
      'reset'   =~ /reset/;
      $entry    =~ m/___\d+___(\d+)___(\d+)___(anti)?sense$/;
      $left_bd  = $1;
      $right_bd = $2;

      if ($right_bd - $left_bd >= 200) {
         push (@keep, $entry);
      } else {
         push (@{$args{remove_small_orfs_ref}}, $entry);
      }
      undef $left_bd;
      undef $right_bd;
   }

   #transfer large ORF list
   @{$args{orfs}} = @keep;
   undef @keep;
   return;
}

sub delete_internal_orfs {
   my %args = @_;
   my (%temporfs, $instance, $iterate);

   %temporfs       = %{$args{orfs}};
   $iterate        = 1;
   $instance       = 1;

   while ($iterate == 1) {
      ${$args{progress_bar}}->configure(-label=>"Deleting internal ORFs iteration $instance");
      ${$args{main_window}} ->update;
      $instance++;
      $iterate = 0; #exit if no internal ORFs have been found
      my ($id_test, $test_left, $test_right);
      my ($iterate_left, $iterate_right);
      my %iterate  = %{$args{orfs}};
      my @delete   = ();
      my %seen_ORF = (); #prevent bidirectional deletion of ORFs
      foreach $id_test (sort {$a <=> $b} keys %{$args{orfs}} ) {
         next if (exists $seen_ORF{$id_test}); #skip if already flagged for deletion
         #get boundaries
         'reset' =~ m/reset/;
         ${$args{orfs}}{$id_test}  =~ m/(\d+)_(\d+)/;
         ($test_left, $test_right) = ($1, $2);

         foreach my $id_iterate (sort {$a <=> $b} keys %iterate ) {
            next if ($id_iterate eq $id_test); #ignore self
            'reset' =~ m/reset/;
            ${$args{orfs}}{$id_iterate} =~ m/(\d+)_(\d+)/;
            ($iterate_left, $iterate_right) = ($1, $2);
            #found internal orf? exit and delete
            if ($iterate_left  <= $test_left &&
                $iterate_right >= $test_right) {
               #add ORF to delete list
               push (@delete, $id_test);
               $seen_ORF{$id_iterate}++;
               $iterate  = 1; #set up for next round of internal ORF scanning
               #go to next ORF
               last;
            }
         }
      }
      undef %seen_ORF;
      #delete internal ORFs
      foreach my $delete_id (@delete) {
         delete $args{orfs}{$delete_id};
      }
   }


   # while ($iterate < 100 && $instance > 0) {
      # ${$args{progress_bar}}->configure(-label=>"Deleting internal ORFs iteration ".($iterate + 1));
      # ${$args{main_window}} ->update;
      # $instance = 0;
      # my $current = 0;
      # my $totalnumber = keys(%{$args{orfs}});
#
      # my ($ref, $left_bd, $right_bd);
      # my (@id, @delete);
#
      # while ( my ($id, $value) = each %{$args{orfs}}) {
         # $current++;
         # 'reset' =~ /reset/;
         # $value =~ m/(\d+)_(\d+)/;
         # $left_bd = $1;
         # $right_bd = $2;
         # while ( my ($temp_id, $temp_value) = each %temporfs) {
            # 'reset' =~ /reset/;
            # $temp_value =~ m/(\d+)_(\d+)/;
            # my $temp_left = $1;
            # my $temp_right = $2;
#
            # next if ($temp_right < $left_bd);
            # next if ($temp_left  > $right_bd);
#
            # if (($temp_left >= $left_bd && $temp_right < $right_bd) || ($temp_left > $left_bd && $temp_right <= $right_bd)) {
               # push (@delete, $temp_id);
               # $instance++;
            # }
         # }
      # }
      # #reset boundary arrays
      # foreach my $id (@delete) {
         # if (exists ${$args{orfs}}{$id}) {
            # delete ${$args{orfs}}{$id};
         # }
      # }
      # %temporfs = %{$args{orfs}};
      # undef @delete;
      # $iterate++;
   # }
   return;
}

sub merge_overlap_orfs {
   my %args = @_;
   my (%temporfs, %duplicates, $instance, $iterate);

   #create new hash to test for duplicates
   while (my ($id, $value) = each %{$args{orfs}}) {
      $duplicates{$value}++;
   }

   $iterate = 0;
   $instance = 1;
   %temporfs = %{$args{orfs}};
   while ($iterate < 100 && $instance > 0) {
      ${$args{progress_bar}}->configure(-label=>"Merging overlapping ORFs iteration ".($iterate + 1));
      ${$args{main_window}} ->update;
      $instance = 0;
      my $current = 0;
      my $totalnumber = keys(%{$args{orfs}});
      my (@delete, @add);

      while (my ($id, $value) = each %{$args{orfs}}) {
         $current++;
         'reset' =~ /reset/;
         $value =~ m/(\d+)_(\d+)/;
         my $left_bd = $1;
         my $right_bd = $2;

         while (my ($temp_id, $temp_value) = each %temporfs) {
            'reset' =~ /reset/;
            $temp_value =~ m/(\d+)_(\d+)/;
            my $temp_left = $1;
            my $temp_right = $2;

            next if ($temp_right < $left_bd);
            next if ($temp_left  > $right_bd);

            #check if ORFs overlap on left_bd
            if ($temp_right <= $right_bd && $temp_right > $left_bd && $temp_left < $left_bd) {
               push (@delete, $id);
               push (@delete, $temp_id);
               push (@add, $id.'_'.$temp_left.'_'.$right_bd);
               $instance++;
            }
            #check of ORFs overlap on right_bd
            if ($temp_left >= $left_bd && $temp_left < $right_bd && $temp_right > $right_bd) {
               push (@delete, $id);
               push (@delete, $temp_id);
               push (@add, $id.'_'.$left_bd.'_'.$temp_right);
               $instance++;
            }
         }
      }
      #delete entries
      foreach my $id (@delete) {
         delete ${$args{orfs}}{$id};
      }
      #add entries
      foreach my $entry (@add) {
         'reset' =~ /reset/;
         $entry =~ m/(\d+)_(\d+)_(\d+)/;
         my $id = $1;
         my $left_bd = $2;
         my $right_bd = $3;
         unless (exists ${$args{orfs}}{$id} || exists $duplicates{$left_bd.'_'.$right_bd}) {
            ${$args{orfs}}{$id} = $left_bd.'_'.$right_bd;
            $duplicates{$left_bd.'_'.$right_bd}++;
         }
      }
      undef @delete;
      undef @add;
      %temporfs = %{$args{orfs}};
      $iterate++;
   }

   undef %duplicates;
   return;
}

1;