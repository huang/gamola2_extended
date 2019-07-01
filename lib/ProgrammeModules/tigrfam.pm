#!/opt/ActivePerl-5.8/bin/perl

#TIGRfam program calls and reformatting of results
#input arguments: progress_bar, nt_seq, left_bd, right_bd, orientation, ini_ref, auto_ini_ref


package ProgrammeModules::tigrfam;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&TIGRfam_file);
use vars qw();

use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use ProgrammeModules::pfam     qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);
use Cwd;

#local variables
my (%args, @TIGRFAM_tasks);

sub TIGRfam_file {
   my %args = (retry_counter => 0,
               @_); #retry counter limits the number of re-runs to check that all runs have comeplted successfully
   my @list = ();
   my @TIGRfam_results = ();
   my $progress = 0;
   my $require_hmmer3 = 0;
   my $CPU_cores = ${$args{auto_ini_ref}}{CPU};
   my @TIGRfam_db = ();
   @TIGRFAM_tasks = ();
   my $types = [
               ["TIGRfam db", ['.h3f', '.LIB']],
               ["All Files"    , ["*"]]
               ];
   my $default_model = '.h3f';

   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "TIGRFam analysis",
                   label        => 'TIGRFam analysis'
                  );
   &show_pbar_2;

   #slurp up sequence
   my ($seq_ref) = &slurp(main_window  => $args{main_window},
                          progress_bar => $args{progress_bar},
                          auto_ini_ref => $args{auto_ini_ref},
                          ini_ref      => $args{ini_ref},
                          directory    => ${$args{ini_ref}}{input_files},
                          filename     => $args{input_file}
                         );
   if ($seq_ref eq '0') {
      &hide_pbar_2;
      return (0);
   }

   #clean up sequence
   ${$seq_ref} =~ s/^\>[^\n]*?\n//;
   ${$seq_ref} =~ s/\s//gs;

   #get boundary parameters
   my $ID = '';
   my $left_bd = '';
   my $right_bd = '';
   my $orientation = '';

   #define selected TIGRfam databases
   @TIGRfam_db = split /\s+\;\s+/, ${$args{auto_ini_ref}}{TIGRfam_db};

   #test if hmmr3 is required (TIGRfam ver > 10)
   if ($TIGRfam_db[0] =~ m/TIGRFAMs_\d+/) {
      $TIGRfam_db[0] =~ m/TIGRFAMs_(\d+)/;
      if ($1 >=10) {
         $require_hmmer3 = 1;
      }
   } else {
      my $TIGR_db = ${$args{main_window}}->getOpenFile(-filetypes        => $types,
                                                       -initialdir       => ${$args{ini_ref}}{TIGRfam_db_path},
                                                       -initialfile      => ${$args{auto_ini_ref}}{full_TIGRfam_db},
                                                       -defaultextension => \$default_model);
      unless (defined $TIGR_db && -e $TIGR_db) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => 'No db file selected, assuming hmmer3 model',
                                           -type    => 'OK',
                                           -icon    => 'error'
                                          );

         $require_hmmer3 = 1;
      };

      $TIGR_db =~ /m_(\d+).*?$/;
      if ($1 >=10) {
         $require_hmmer3 = 1;
      }
   }

   #map correct entries for respective input file
   @list = grep /$args{input_file}\___/, @{$args{combined_orf_ref}};

   #no ORFs? return to main
   if ($#list < 0) {
      print "\n\nTIGRfam:Error, no ORFs for file $args{input_file} \n\n";
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "TIGRFam: No gene model found for $args{input_file}",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "TIGRFAM: No gene model found for $args{input_file}".
                     "\nReturn to main\n\n";
      close ERRORLOG;
      &hide_pbar_2;
      #return 1 as this is a non-critical error
      return (1);
   }

   #remove existing entries if re-use results
   if (${$args{auto_ini_ref}}{reuse_results} == 1 || $args{retry_counter} > 0) {
      my (@jobs, %exist_tigrfam_file);
      &update_pbar_2(title        => 'TIGRFam analysis',
                     label        => "Testing for existing results",
                     progress     => 1,
                    );
      ${$args{progress_bar}}->configure(-label=>" TIGRFam ");
      ${$args{main_window}}->update;

      #get existing result files
      opendir SEQINPUT, ${$args{ini_ref}}{TIGRfam_results};
      my @TIGRfam_results = grep /$args{input_file}\_/, readdir(SEQINPUT);
      closedir SEQINPUT;
      #create hash
      %exist_tigrfam_file = ();
      foreach (@TIGRfam_results) {
         $exist_tigrfam_file{$_} = '1';
      }

      #iterate over full gene model
      @jobs = ();
      foreach my $entry (@list) {
         my ($ID, $seen);
         'reset' =~ m/reset/;
         $entry =~ m/$args{input_file}\___(\d+)/;
         $ID = $1;

         foreach my $selected_db (@TIGRfam_db) {
            $seen = 0;

            #filesize 0? -> rerun
            #if (exists $exist_tigrfam_file{$args{input_file}.'_'.$selected_db.'_'.$ID} && -s ${$args{ini_ref}}{TIGRfam_results}.'/'.$args{input_file}.'_'.$selected_db.'_'.$ID < 1) {
            #   last;
            #}

            #file exists and has content?-> test file
            #test if file actually contains data and is complete
            if (-e ${$args{ini_ref}}{TIGRfam_results}.'/'.$args{input_file}.'_'.$selected_db.'_'.$ID) {
               my ($file_ref) = &slurp(main_window  => $args{main_window},
                                      progress_bar => $args{progress_bar},
                                      auto_ini_ref => $args{auto_ini_ref},
                                      ini_ref      => $args{ini_ref},
                                      directory    => ${$args{ini_ref}}{TIGRfam_results},
                                      filename     => $args{input_file}.'_'.$selected_db.'_'.$ID
                                     );
               #something went wrong, rerun
               if ($file_ref eq '0') {
                  unlink ${$args{ini_ref}}{TIGRfam_results}.'/'.$args{input_file}.'_'.$selected_db.'_'.$ID;
                  last;
               }
               #incomplete file, rerun
               unless ($$file_ref =~ m/\n\/\//s) {
                  unlink ${$args{ini_ref}}{TIGRfam_results}.'/'.$args{input_file}.'_'.$selected_db.'_'.$ID;
                  last;
               }

               #test normal file
               if ($$file_ref =~ /\w+/) {
                  $seen            = 1;
                  #my $local_evalue = '';

                  #rerun if somethings gone really wrong
                  unless ($$file_ref =~ /Scores for sequence family classification/si || $$file_ref =~ /Scores for complete sequence/si) {
                     ${$args{main_window}}->messageBox(-title   => 'Error',
                                                       -message => "Could not parse from file $args{input_file}\_$selected_db\_$ID\.",
                                                       -icon    => 'error',
                                                       -type    => 'OK');
                     open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
                     print WRITE "\nError in TIGRfam analysis using database $selected_db for file $args{input_file}\n\n";
                     close WRITE;
                     last;
                  }

                  #no hits? test next file
                  if ($$file_ref =~ /\[no hits above thresholds\]/si || $$file_ref =~ /\[No hits detected that satisfy reporting thresholds\]/si) {
                     next;
                  }

                  #determine evalue of best hit
                  #'reset' =~ m/reset/;
                  #$$file_ref =~ m/\nScores for.*?inclusion threshold ------.*?\s+(\S*?)\s+/is;
                  #$local_evalue = $1;
                  #unless ($local_evalue =~ /^[\d\.\-e]+$/ || $local_evalue eq "") {
                  #   ${$args{main_window}}->messageBox(-title   => 'Error',
                  #                                     -message => "Could not parse evalue $local_evalue from file $args{input_file}\_$selected_db\_$ID\.",
                  #                                     -icon    => 'error',
                  #                                     -type    => 'OK');
                  #   open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
                  #   print WRITE "\nError in TIGRfam analysis using database $selected_db for file $args{input_file}\n\n";
                  #   close WRITE;
                  #   last;
                  #}

                  #hit below threshold? skip tests for other files
                  #if ($local_evalue <= ${$args{ini_ref}}{TIGRfam_cut}) {
                  #   undef $local_evalue;
                  #   undef $$file_ref;
                  #   last;
                  #}
                  #undef $local_evalue;
                  undef $$file_ref;
               }
            }
            #last unless ($seen == 1); #exit if one result file is missing
         }
         next if ($seen == 1);
         if ($seen == 0) {
            foreach my $selected_db (@TIGRfam_db) {
               unlink ${$args{ini_ref}}{TIGRfam_results}.'/'.$args{input_file}.'_'.$selected_db.'_'.$ID; #remove existing partial db hits
               #print "Deleting $seen ...".$args{input_file}.'_'.$selected_db.'_'.$ID." ...\n";
            }
         }
         push (@jobs, $entry);
         undef $seen;
      }
      @list = @jobs;
      undef @TIGRfam_results;
      undef @jobs;
      undef %exist_tigrfam_file;
   }

   #max number of fasta entries
   my $max_count = $#list + 1;

   #iterate over array
   foreach my $entry (@list) {
      #update progress bar
      $progress++;
      if (($progress % 100) == 0) {
         &update_pbar_2(title        => "TIGRFam analysis",
                        label        => "Reading gene model for $args{input_file}",
                        progress     => ($progress / $max_count) * 100,
                       );
      }
      ($ID, $left_bd, $right_bd, $orientation) = "";
      $entry =~ m/$args{input_file}\___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
      $ID = $1;
      $left_bd = $2;
      $right_bd = $3;
      $orientation = $4;
      my $query_seq_nt = "";

      unless (defined $orientation && $orientation =~ /(sense|antisense)/) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Error parsing gene model for $args{input_file}",
                                           -type    => 'OK',
                                           -icon    => 'info');
         &hide_pbar_2;
         return (0); #return to main
      }

      #extract nt sequence from ORF
      $query_seq_nt = substr(${$seq_ref}, ($left_bd - 1), ($right_bd - $left_bd + 1));

      #translate nt sequence to aa if selected
      my $query_seq_aa = '';

      $query_seq_aa = nt2aa(main_window  => $args{main_window},
                            progress_bar => $args{progress_bar},
                            auto_ini_ref => $args{auto_ini_ref},
                            ini_ref      => $args{ini_ref},
                            module       => 'TIGRFam',
                            orientation  => $orientation,
                            sequence     => $query_seq_nt,
                            filename     => $args{input_file},
                            left_bd      => $left_bd,
                            right_bd     => $right_bd
                           );
      #create array for new PFam tasks
      push (@TIGRFAM_tasks, $args{input_file}.'___'.$ID.'___'.$left_bd.'___'.$right_bd.'___'.$orientation.'___'.$query_seq_aa);
   }

   #limit number of CPUs
   if (${$args{auto_ini_ref}}{HMMER_limit} != 0) {
      $CPU_cores = ${$args{auto_ini_ref}}{HMMER_limit};
   }

   #start multithreading TIGRfam hmmer
   for (my $i=1; $i<=$#TIGRFAM_tasks+2; $i = $i+$CPU_cores) {
      my $count = $i + $CPU_cores;
      if ($count > $#TIGRFAM_tasks+2) {$count = $#TIGRFAM_tasks+2};

      &update_pbar_2(title        => "TIGRFam analysis",
                     label        => "TIGRfam analysis for $args{input_file}, $count of $max_count",
                     progress     => ($count / ($#TIGRFAM_tasks + 2)) * 100,
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
            #grab entries
            my $sequence = "";
            my $left_bd = "";
            my $right_bd = "";
            my $orientation = "";
            my $ORF_number = "";
            $TIGRFAM_tasks[$j-1] =~ m/$args{input_file}___(\d+)___(\d+)___(\d+)___([a-z]+)___(.+)/;
            $ORF_number = $1;
            $left_bd = $2;
            $right_bd = $3;
            $orientation = $4;
            $sequence = $5;

            #run Pfam for selected feature
            if ($require_hmmer3 == 0) {
               &pfam_run (auto_ini_ref  => $args{auto_ini_ref},
                          ini_ref       => $args{ini_ref},
                          input_file    => $args{input_file},
                          ORF_number    => $ORF_number,
                          left_bd       => $left_bd,
                          right_bd      => $right_bd,
                          orientation   => $orientation,
                          sequence      => $sequence,,
                          output_folder => ${$args{ini_ref}}{TIGRfam_results},
                          full_database => ${$args{auto_ini_ref}}{full_TIGRfam_db},
                          threshold     => ${$args{ini_ref}}{TIGRfam_cut}
                         );
            } else {
               &pfam3_run (auto_ini_ref  => $args{auto_ini_ref},
                           ini_ref       => $args{ini_ref},
                           input_file    => $args{input_file},
                           ORF_number    => $ORF_number,
                           left_bd       => $left_bd,
                           right_bd      => $right_bd,
                           orientation   => $orientation,
                           sequence      => $sequence,,
                           output_folder => ${$args{ini_ref}}{TIGRfam_results},
                           full_database => ${$args{auto_ini_ref}}{full_TIGRfam_db},
                           threshold     => ${$args{ini_ref}}{TIGRfam_cut}
                          );
            }
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

   #re-run to test if all Pfam hits have been successfully run
   while ($#TIGRFAM_tasks >= 0 && $args{retry_counter} < 3) {
      $args{retry_counter}++;
      &hide_pbar_2;
      TIGRfam_file(main_window      => $args{main_window},
                   progress_bar     => $args{progress_bar},
                   auto_ini_ref     => $args{auto_ini_ref},
                   ini_ref          => $args{ini_ref},
                   directory        => $args{directory},
                   input_file       => $args{input_file},
                   combined_orf_ref => $args{combined_orf_ref},
                   retry_counter    => $args{retry_counter}
                   );
   }


   &hide_pbar_2;
   return (1);
}

1;
