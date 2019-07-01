#!/opt/ActivePerl-5.8/bin/perl

#Pfam program calls and reformatting of results
#input arguments: selected databases, output folder, progress_bar, nt_seq, left_bd, right_bd, orientation, ini_ref, auto_ini_ref
#will be used for Pfam and TIGRfam analysis, separate parsing later on

package ProgrammeModules::pfam;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&Pfam_file &pfam_run &pfam3_run);
use vars qw();

use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);
use Cwd;

#local variables
my (%args, @PFAM_tasks);


sub Pfam_file {
   my %args     = (retry_counter => 0,
                   @_); #retry counter limits the number of re-runs to check that all runs have completed successfully
   my @list     = ();
   my $progress = 0;
   my @pfam_db  = ();
   @PFAM_tasks  = ();
   my $CPU_cores = ${$args{auto_ini_ref}}{CPU};

   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "PFam analysis",
                   label        => 'PFam analysis'
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

   #define selected Pfam databases
   @pfam_db = split /\s+\;\s+/, ${$args{auto_ini_ref}}{Pfam_db};

   #map correct entries for respective input file
   @list = grep /$args{input_file}\___/, @{$args{combined_orf_ref}};

   #no ORFs? return to main
   if ($#list < 0) {
      print "\n\nPfam:Error, no ORFs for file $args{input_file} \n\n";
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "PFam: No gene model found for $args{input_file}",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "PFAM: No gene model found for $args{input_file}".
                     "\nReturn to main\n\n";
      close ERRORLOG;
      &hide_pbar_2;
      #return 1 as this is a non-critical error
      return (1);
   }

   #remove existing entries from job list if re-use results or if incomplete results list
   if (${$args{auto_ini_ref}}{reuse_results} == 1 || $args{retry_counter} > 0) {
      #my (@jobs, %exist_pfam_file, @Pfam_results);
      my @jobs = ();
      &update_pbar_2(title        => 'PFam analysis',
                     label        => "Testing for existing results",
                     progress     => 1,
                    );
      ${$args{progress_bar}}->configure(-label=>" PFam ");
      ${$args{main_window}}->update;

      #get existing result files
      #opendir SEQINPUT, ${$args{ini_ref}}{pfam_results};
      #@Pfam_results = grep /$args{input_file}\_/, readdir(SEQINPUT);
      #closedir SEQINPUT;

      #create hash
      #%exist_pfam_file = ();
      #foreach (@Pfam_results) {
      #      $exist_pfam_file{$_} = '1';
      #}

      #iterate over full gene model
      foreach my $entry (@list) {

         my ($ID, $seen);
         'reset' =~ m/reset/;
         $entry =~ m/$args{input_file}\___(\d+)/;
         $ID = $1;

         foreach my $selected_db (@pfam_db) {
            $seen = 0;

            #file exists and has content?-> test file
            #test if file actually contains data and is complete
            if (-e ${$args{ini_ref}}{pfam_results}.'/'.$args{input_file}.'_'.$selected_db.'_'.$ID) {
               my ($file_ref) = &slurp(main_window  => $args{main_window},
                                      progress_bar => $args{progress_bar},
                                      auto_ini_ref => $args{auto_ini_ref},
                                      ini_ref      => $args{ini_ref},
                                      directory    => ${$args{ini_ref}}{pfam_results},
                                      filename     => $args{input_file}.'_'.$selected_db.'_'.$ID
                                     );

               #something went wrong, rerun
               if ($file_ref eq '0') {
                  unlink ${$args{ini_ref}}{pfam_results}.'/'.$args{input_file}.'_'.$selected_db.'_'.$ID;
                  print "\n..Error reading file ". $args{input_file}.'_'.$selected_db.'_'.$ID . "...$file_ref...\n";
                  last;
               }
               #incomplete file, rerun
               unless ($$file_ref =~ m/\n\/\//s) {
                  unlink ${$args{ini_ref}}{pfam_results}.'/'.$args{input_file}.'_'.$selected_db.'_'.$ID;
                  print "\n...Incomplete file: ". $args{input_file}.'_'.$selected_db.'_'.$ID . "...$file_ref...\n";
                  last;
               }

               #test normal file with content
               if ($$file_ref =~ /\w+/) {
                  #skip if somethings gone really wrong
                  unless ($$file_ref =~ /(Scores for sequence family classification|Scores for complete sequence \(score includes all domains\))/si) {
                     ${$args{main_window}}->messageBox(-title   => 'Error',
                                                       -message => "Could not parse from file $args{input_file}\_$selected_db\_$ID\.",
                                                       -icon    => 'error',
                                                       -type    => 'OK');
                     open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
                     print WRITE "\nError in Pfam analysis using database $selected_db for file $args{input_file}\n\n";
                     close WRITE;
                     last;
                  }
                  $seen = 1;

                  #no hits? test next file
                  if ($$file_ref =~ /\[(no hits above thresholds|No hits detected that satisfy reporting thresholds)\]/si) {
                     undef $$file_ref;
                     next;
                  } else {
                     undef $$file_ref;
                     last; #exit if found result file has actual results, task compelted
                  }
               }
            }
         }

         next if ($seen == 1);
         if ($seen == 0) {
            foreach my $selected_db (@pfam_db) {
               unlink ${$args{ini_ref}}{pfam_results}.'/'.$args{input_file}.'_'.$selected_db.'_'.$ID; #remove existing partial db hits
            }
            push (@jobs, $entry);
         }
      }
      @list = ();
      @list = @jobs;
      #undef @Pfam_results;
      undef @jobs;
      #undef %exist_pfam_file;
   }

   #max number of fasta entries
   my $max_count = $#list + 1;

   #iterate over array
   foreach my $entry (@list) {
      #update progress bar
      $progress++;
      if (($progress % 100) == 0) {
         &update_pbar_2(title        => 'PFam analysis',
                        label        => "Reading gene model for $args{input_file}",
                        progress     => ($progress / $max_count) * 100,
                       );
         ${$args{main_window}}->update;
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
                            module       => 'Pfam',
                            orientation  => $orientation,
                            sequence     => $query_seq_nt,
                            filename     => $args{input_file},
                            left_bd      => $left_bd,
                            right_bd     => $right_bd
                           );

      #create array for new PFam tasks
      push (@PFAM_tasks, $args{input_file}.'___'.$ID.'___'.$left_bd.'___'.$right_bd.'___'.$orientation.'___'.$query_seq_aa);
   }
   #skip exisiting results?
   ${$args{progress_bar}}->configure(-label=>" PFam ");
   ${$args{main_window}}->update;

   #limit number of CPUs
   if (${$args{auto_ini_ref}}{HMMER_limit} != 0) {
      $CPU_cores = ${$args{auto_ini_ref}}{HMMER_limit};
   }

   #start multithreading PFAM hmmer
   for (my $i=1; $i<=$#PFAM_tasks+2; $i = $i+$CPU_cores) {
      my $count = $i + $CPU_cores;
      if ($count > $#PFAM_tasks+2) {$count = $#PFAM_tasks+2};

      &update_pbar_2(title        => 'PFam analysis',
                     label        => "Pfam analysis for $args{input_file}, ". ($count - 1) ." of $max_count",
                     progress     => ($count / ($#PFAM_tasks + 2)) * 100,
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
            $PFAM_tasks[$j-1] =~ m/$args{input_file}___(\d+)___(\d+)___(\d+)___([a-z]+)___(.+)/;
            $ORF_number = $1;
            $left_bd = $2;
            $right_bd = $3;
            $orientation = $4;
            $sequence = $5;

            if (${$args{auto_ini_ref}}{use_Pfam2} == 1) {

               #run Pfam for selected feature
               &pfam_run (auto_ini_ref  => $args{auto_ini_ref},
                          ini_ref       => $args{ini_ref},
                          input_file    => $args{input_file},
                          ORF_number    => $ORF_number,
                          left_bd       => $left_bd,
                          right_bd      => $right_bd,
                          orientation   => $orientation,
                          sequence      => $sequence,
                          output_folder => ${$args{ini_ref}}{pfam_results},
                          full_database => ${$args{auto_ini_ref}}{full_Pfam_db},
                          threshold     => ${$args{ini_ref}}{PFam_cut}
                         );
            } else {
               #run Pfam for selected feature
               &pfam3_run (auto_ini_ref  => $args{auto_ini_ref},
                           ini_ref       => $args{ini_ref},
                           input_file    => $args{input_file},
                           ORF_number    => $ORF_number,
                           left_bd       => $left_bd,
                           right_bd      => $right_bd,
                           orientation   => $orientation,
                           sequence      => $sequence,
                           output_folder => ${$args{ini_ref}}{pfam_results},
                           full_database => ${$args{auto_ini_ref}}{full_Pfam_db},
                           threshold     => ${$args{ini_ref}}{PFam_cut}
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
   while ($#PFAM_tasks >= 0 && $args{retry_counter} < 3) {
      $args{retry_counter}++;

      &hide_pbar_2;
      Pfam_file(main_window      => $args{main_window},
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


# start PFAM
sub pfam_run { #used in Pfam and TIGRfam
   my %args = @_;
   my @local_db = ();
   my $local_PFAM = "";
   my $local_evalue = "";
   my $seqfile = "";
   my $write_result = "";

   $seqfile = $args{input_file}."_$args{ORF_number}";
   open WRITE, "+> $args{output_folder}\/$seqfile" or die;
   if ($args{orientation} eq 'sense') {
      print WRITE '>'.$args{ORF_number}.'_'.$args{left_bd}.'_'.$args{right_bd}."\n";
   } elsif($args{orientation} eq 'antisense') {
      print WRITE '>'.$args{ORF_number}.'_'.$args{right_bd}.'_'.$args{left_bd}."\n";
   }
   print WRITE $args{sequence};
   close WRITE;

   #generate local array for selected PFam databases
   @local_db = split/ ; /,$args{full_database};

   #iterate over local db's
   foreach my $pfam_db (@local_db) {
      #grab db name only
      my $db_name = "";
      $pfam_db =~ m/\/([^\/]+)$/;
      $db_name = $1;
      #define name for output file
      $write_result = $args{input_file}.'_'.$db_name.'_'.$args{ORF_number};

      ##############
      system "${$args{ini_ref}}{pfam_executable}/hmmpfam $pfam_db $args{output_folder}\/$seqfile  > $args{output_folder}\/$write_result";
      ##############
      {
         local( $/, *WRITE_TEST ) ;
         open( WRITE_TEST, $args{output_folder}.'/'.$write_result ) or do {  };
         $local_PFAM = <WRITE_TEST>;
         close WRITE_TEST;
      }

      #stop if somethings gone really wrong
      unless ($local_PFAM =~ /Scores for sequence family classification/si) {
         print "Error in PFam/TIGRfam analysis using database $pfam_db for file $args{input_file}";
         open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
         print WRITE "\nError in PFam/TIGRfam analysis using database $pfam_db for file $args{input_file}\n\n";
         close WRITE;
         next;
      }

      #no hits? next database
      if ($local_PFAM =~ /\[no hits above thresholds\]/si) {
         next;
      }

      #determine evalue of best hit
      'reset' =~ m/reset/;
      $local_PFAM =~ m/\nModel\s+Description[^\n]*?\n[\s\-]+\n[^\n]+?\s+(\S+)\s+\d+\s*\n/is;
      $local_evalue = $1;
      unless ($local_evalue =~ /^[\d\.\-e]+$/ || $local_evalue eq "") {
         print "Error in PFam/TIGRfam analysis for evalue ...$local_evalue ... using database $pfam_db for $write_result\n";
         open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
         print WRITE "\nError in PFam/TIGRfam analysis using database $pfam_db for file $args{input_file}\n".
                     "Could not process local evalue for entry $write_result\n\n";
         close WRITE;
         undef $local_PFAM;
         unlink $args{output_folder}.'/'.$seqfile;
         return (0);
      }

      #hit below threshold? skip other db's
      if ($local_evalue <= $args{threshold}) {
         undef $local_PFAM;
         unlink $args{output_folder}.'/'.$seqfile;
         return (1); #exit if good hit
      }

      #no hits above threshold? next database
      if ($local_PFAM =~ /\[no hits above thresholds\]/si || $local_evalue > $args{threshold}) {
         undef $local_PFAM;
         next;
      } else {
         undef $local_PFAM;
         unlink $args{output_folder}.'/'.$seqfile;
         return (1); #exit if good hit
      }
   }

   unlink $args{output_folder}.'/'.$seqfile;
   undef $local_PFAM;
   return (1);
}

# start PFAM
sub pfam3_run {
   my %args = @_;

   my @local_db = ();
   my $local_PFAM = "";
   my $local_evalue = "";
   my $seqfile = "";
   my $write_result = "";

   $seqfile = $args{input_file}."_$args{ORF_number}";

   #modify filename for hmmer3; does not like '.'
   $seqfile =~ s/\./_/g;

   open WRITE, "+> $args{output_folder}\/$seqfile" or do {
      print "Error in PFam3 analysis: could not create sequence file $seqfile in folder $args{output_folder}.\n";
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "\nError in PFam3 analysis:\n".
                  "Could not create sequence file $seqfile in folder $args{output_folder}\n\n";
      close WRITE;
      undef $local_PFAM;
      unlink $args{output_folder}.'/'.$seqfile;
      return (0);
   };

   if ($args{orientation} eq 'sense') {
      print WRITE '>'.$args{ORF_number}.'_'.$args{left_bd}.'_'.$args{right_bd}."\n";
   } elsif($args{orientation} eq 'antisense') {
      print WRITE '>'.$args{ORF_number}.'_'.$args{right_bd}.'_'.$args{left_bd}."\n";
   }
   print WRITE $args{sequence};
   close WRITE;

   #generate local array for selected PFam databases
   @local_db = split/ ; /,$args{full_database};

   #iterate over local db's
   foreach my $pfam_db (@local_db) {
      #grab db name only
      my $db_name = "";
      $pfam_db =~ m/\/([^\/]+)$/;
      $db_name = $1;
      #define name for output file
      $write_result = $args{input_file}.'_'.$db_name.'_'.$args{ORF_number};

      #clean up db name
      if ($pfam_db =~ m/\.h3.$/) {
         $pfam_db =~ s/\.h3.$//;
      }
      ##############
      system "/${$args{ini_ref}}{pfam3_executable}/hmmscan -o $args{output_folder}\/$write_result $pfam_db $args{output_folder}\/$seqfile ";
      ##############

      {
         local( $/, *WRITE_TEST ) ;
         open( WRITE_TEST, $args{output_folder}.'/'.$write_result );
         $local_PFAM = <WRITE_TEST>;
         close WRITE_TEST;
      }

      #stop if somethings gone really wrong
      unless ($local_PFAM =~ /Scores for complete sequence/si) {
         print "\nError in PFam3 analysis using database $pfam_db for file $args{input_file}\n";
         open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
         print WRITE "\nError in PFam3 analysis using database $pfam_db for file $args{input_file}\n\n";
         close WRITE;
         next;
      }

      #no hits? next database
      if ($local_PFAM =~ /\[No hits detected that satisfy reporting thresholds\]/si) {
         next;
      }

      #determine evalue of best hit
      'reset' =~ m/reset/;
      $local_PFAM =~ m/\nE-value\s+score[^\n]*?\n[\s\-]+?\n\s+(\S+)\s/is;
      $local_evalue = $1;
      unless ($local_evalue =~ /^[\d\.\-e]+$/ || $local_evalue eq "") {
         print "Error in PFam3 analysis for evalue ...$local_evalue ... using database $pfam_db for $write_result\n";
         open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
         print WRITE "\nError in PFam/TIGRfam analysis using database $pfam_db for file $args{input_file}\n".
                     "Could not process local evalue for entry $write_result\n\n";
         close WRITE;
         undef $local_PFAM;
         unlink $args{output_folder}.'/'.$seqfile;
         return (0);
      }

      #hit below threshold? skip other db's
      if ($local_evalue <= $args{threshold}) {
         undef $local_PFAM;
         unlink $args{output_folder}.'/'.$seqfile;
         return (1); #exit if good hit
      }

      #no hits above threshold? next database
      if ($local_PFAM =~ /\[No hits detected that satisfy reporting thresholds\]/si || $local_evalue > $args{threshold}) {
         undef $local_PFAM;
         next;
      } else {
         undef $local_PFAM;
         unlink $args{output_folder}.'/'.$seqfile;
         return (1); #exit if good hit
      }
   }

   unlink $args{output_folder}.'/'.$seqfile;
   undef $local_PFAM;
   return (1);
}

1;