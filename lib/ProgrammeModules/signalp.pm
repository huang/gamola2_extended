#!/opt/ActivePerl-5.8/bin/perl

#SignalP program calls and reformatting of results
#input arguments: main_window, progress_bar, auto_ini_ref, ini_ref,
#                 directory, input_array_ref, output_folder, input_type

package ProgrammeModules::signalp;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&signalp &signalp_summary);
use vars qw();

use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);
use Cwd;
use File::Copy;

#local variables
my (%args, $signalpresult);

sub signalp {
   my %args            = @_;
   my @list            = ();
   my @signalp_results = ();
   my @signalp_tasks   = ();
   my $progress        = 0;
   my $signalp_version = '';

   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "SignalP analysis",
                   label        => 'SignalP analysis'
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
   ${$seq_ref} =~ s/^\>[^\n\r]*?\n\r?//s;
   ${$seq_ref} =~ s/\s//gs;

   #map correct entries for respective input file
   @list = grep /$args{input_file}/, @{$args{combined_orf_ref}};
   if ($#list < 0) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No gene model for file $args{input_file}\nSkipping file",
                                                    -bitmap  => 'info',
                                                    -buttons => ['OK'],
                                                   );
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "SignalP:  No gene model for file $args{input_file}".
                     "\nSkipping file\n\n";
      close ERRORLOG;
      &hide_pbar_2;
      return (0);
   }

   #remove existing entries if re-use results
   if (${$args{auto_ini_ref}}{reuse_results} == 1) {
      my (@jobs, @signalp_results, %exist_signalp_file);
      &update_pbar_2(title        => 'SignalP analysis',
                     label        => "Testing for existing results",
                     progress     => 1,
                    );

      #read existing results
      opendir SEQINPUT, ${$args{ini_ref}}{signalp_results};
      @signalp_results = grep /^$args{input_file}\_SignalP_\d+$/, readdir(SEQINPUT);
      closedir SEQINPUT;
      #create hash
      foreach (@signalp_results) {
         $exist_signalp_file{$_} = '1';
      }
      #iterate over gene model
      foreach my $entry (@list) {
         my ($ID);
         $entry =~ m/$args{input_file}\___(\d+)___/;
         $ID    = $1;

         #skip SignalP if result file exists and re-use is activated
         next if (exists $exist_signalp_file{$args{input_file}.'_SignalP_'.$ID} && -s ${$args{ini_ref}}{signalp_results}.'/'.$args{input_file}.'_SignalP_'.$ID.'/output.txt' > 0);
         #add if required
         push (@jobs, $entry);
      }
      @list = @jobs;
      undef @jobs;
      undef @signalp_results;
      undef %exist_signalp_file;
   }

   #max number of fasta entries
   my $max_count = $#list + 1;
   if ($max_count < 1) {$max_count = 1};

   #iterate over array
   foreach my $entry (@list) {
      #update progress bar
      $progress++;
      &update_pbar_2(title        => "SignalP analysis",
                     label        => "Reading entry from $args{input_file}",
                     progress     => ($progress / $max_count) * 100,,
                    );

      #get boundary parameters
      my ($ID, $left_bd, $right_bd, $orientation);
      'reset' =~ m/reset/;
      $entry =~ m/$args{input_file}\___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
      $ID = $1;
      $left_bd = $2;
      $right_bd = $3;
      $orientation = $4;
      my $query_seq_nt = "";

      unless (defined $orientation && $orientation =~ /(sense|antisense)/) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing gene model for $args{input_file}",
                                                       -bitmap  => 'info',
                                                       -buttons => ['OK'],
                                                      );
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "SignalP:  Error parsing gene model for $args{input_file}".
                        "\nSReturn to main\n\n";
         close ERRORLOG;
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
                            module       => 'SignalP',
                            orientation  => $orientation,
                            sequence     => $query_seq_nt,
                            filename     => $args{input_file},
                            left_bd      => $left_bd,
                            right_bd     => $right_bd
                           );

      #create array for new SignalP tasks
      push (@signalp_tasks, $args{input_file}.'___'.$ID.'___'.$left_bd.'___'.$right_bd.'___'.$orientation.'___'.$query_seq_aa);
   }

   ${$args{progress_bar}}->configure(-label=>" SignalP ");
   ${$args{main_window}}->update;

   #which version of SignalP? v3 & 4 are supported
   my ($file_ref) = &slurp(main_window   => $args{main_window},
                           auto_ini_ref  => $args{auto_ini_ref},
                           ini_ref       => $args{ini_ref},
                           filename      => 'signalp',
                           directory     => ${$args{ini_ref}}{signalp_dir}
                          );
   if ($$file_ref =~ m/SignalP 4/si) {
      $signalp_version = 4;
   } elsif ($$file_ref =~ m/SignalP 3/si) {
      $signalp_version = 3;
   } else {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Could not determine SignalP version. Please contact developer for bugfixing. Sorry mate.",
                                                    -bitmap  => 'info',
                                                    -buttons => ['OK'],
                                                   );
      &hide_pbar_2;
      return (0);
   }
   undef $file_ref;

   #start multithreading SignalP
   for (my $i=1; $i<=$#signalp_tasks+2; $i = $i+${$args{auto_ini_ref}}{CPU}) {
      my $count = $i + ${$args{auto_ini_ref}}{CPU};
      if ($count > $#signalp_tasks+2) {$count = $#signalp_tasks+2};

      &update_pbar_2(title        => "SignalP analysis",
                     label        => "SignalP analysis for $args{input_file}, $count of $max_count",
                     progress     => ($count / $max_count) * 100,
                    );

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
            $signalp_tasks[$j-1] =~ m/$args{input_file}___(\d+)___(\d+)___(\d+)___([a-z]+)___(.+)/;
            $ORF_number = $1;
            $left_bd = $2;
            $right_bd = $3;
            $orientation = $4;
            $sequence = $5;

            #run signalP for selected feature
            #get current working directory
            my $curdir = "";
            $curdir = getcwd();

            #change to SignalP directory
            chdir ${$args{ini_ref}}{signalp_dir};

            #define input sequence in fasta format
            my $fasta = '>'.$args{input_file}.'___'.$ORF_number.'___'.$left_bd.'___'.$right_bd.'___'.$orientation."\n".$sequence;
            my $current_output_dir = ${$args{ini_ref}}{signalp_results}.'/'.$args{input_file}.'_SignalP_'.$ORF_number;
            my $filename = $current_output_dir.'/'.$args{input_file}.'_SignalP_'.$ORF_number.'.seq';

            #make current output directory
            unless (-e $current_output_dir) {
               mkdir ($current_output_dir, 0777)  or do {
                  print "\nCannot create dir $current_output_dir\n";
                  CORE::exit;
               };
            }

            open WRITE, "+>".$filename;
            print WRITE $fasta;
            close WRITE;

            #`./signalp -q -t ${$args{auto_ini_ref}}{signalp_code} -trunc ${$args{auto_ini_ref}}{signalp_trunc} -d $current_output_dir -f full -graphics eps -m ${$args{auto_ini_ref}}{signalp_meth} $filename`; #child
            if ($signalp_version == 3) {
               if (${$args{auto_ini_ref}}{signalp_code} eq 'both') {
                  my $n = 0;
                  while ($n <= 3) {
                     `./signalp -q -t gram- -trunc ${$args{auto_ini_ref}}{signalp_trunc} -d $current_output_dir -f full -m ${$args{auto_ini_ref}}{signalp_meth} $filename`; #child
                     move("$current_output_dir\/output.txt", "$current_output_dir\/output_gram_neg.txt");
                     open READ, '<'.$current_output_dir.'/output_gram_neg.txt';
                     while (<READ>) {
                        if ($_ =~ m/Signal peptide /) {$n = 10};
                     }
                     close READ;
                     $n++;
                  }
                  $n = 0;
                  while ($n <= 3) {
                     `./signalp -q -t gram+ -trunc ${$args{auto_ini_ref}}{signalp_trunc} -d $current_output_dir -f full -m ${$args{auto_ini_ref}}{signalp_meth} $filename`; #child
                     move("$current_output_dir\/output.txt", "$current_output_dir\/output_gram_pos.txt");
                     open READ, '<'.$current_output_dir.'/output_gram_pos.txt';
                     while (<READ>) {
                        if ($_ =~ m/Signal peptide /) {$n = 10};
                     }
                     close READ;
                     $n++;
                  }
                  `cat $current_output_dir\/output_gram_neg.txt $current_output_dir\/output_gram_pos.txt > $current_output_dir\/output.txt`;
                  unlink "$current_output_dir\/output_gram_neg.txt";
                  unlink "$current_output_dir\/output_gram_pos.txt";
               } else {
                  my $n = 0;
                  while ($n <= 3) {
                     `./signalp -q -t ${$args{auto_ini_ref}}{signalp_code} -trunc ${$args{auto_ini_ref}}{signalp_trunc} -d $current_output_dir -f full -m ${$args{auto_ini_ref}}{signalp_meth} $filename`; #child
                     open READ, '<'.$current_output_dir.'/output.txt';
                     while (<READ>) {
                        if ($_ =~ m/Signal peptide /) {$n = 10};
                     }
                     close READ;
                     $n++;
                  }
               }
            } elsif ($signalp_version == 4) {
               #define cutoff parameter
               my $signalp_cutoff;
               if (${$args{auto_ini_ref}}{signalp_meth} eq 'best') {
                  $signalp_cutoff = "-U ${$args{auto_ini_ref}}{signalp_cutoff}";
               } elsif (${$args{auto_ini_ref}}{signalp_meth} eq 'notm') {
                  $signalp_cutoff = "-u ${$args{auto_ini_ref}}{signalp_cutoff}";
               }
               chdir $current_output_dir;

               #make sure SignalP4 runs over all ORFs and produces output files;
               my $n = 0;
               while ($n <= 3) {
                  `${$args{ini_ref}}{signalp_dir}/signalp -t ${$args{auto_ini_ref}}{signalp_code} -c ${$args{auto_ini_ref}}{signalp_trunc} -M ${$args{auto_ini_ref}}{signalp_min} -f all -s ${$args{auto_ini_ref}}{signalp_meth} $signalp_cutoff $filename > output.txt`; #child
                  open READ, '<'.$current_output_dir.'/output.txt';
                  while (<READ>) {
                     if ($_ =~ m/SP=\'/) {$n = 10};
                  }
                  close READ;
                  $n++;
               }
            }

            #change back to previous working dir
            chdir $curdir;

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

   #remove leftover temp files from SignalP4
   if ($signalp_version == 4) {
      @signalp_results = ();
      opendir SEQINPUT, ${$args{ini_ref}}{signalp_results};
      @signalp_results = grep /^signalp_\d+$/, readdir(SEQINPUT);
      closedir SEQINPUT;
      foreach my $directory (@signalp_results) {
         cleanup (${$args{ini_ref}}{signalp_results}.'/'.$directory);
         unlink ${$args{ini_ref}}{signalp_results}.'/'.$directory;
      }
   }

   &hide_pbar_2;
   return (1);
}

sub signalp_summary {
   my %args = @_;
   my @signalp_results = ();
   my @summary = ();
   my $progress = 1;
   my $signalp_version;

   #which version of SignalP? v3 & 4 are supported
   my ($file_ref) = &slurp(main_window   => $args{main_window},
                           auto_ini_ref  => $args{auto_ini_ref},
                           ini_ref       => $args{ini_ref},
                           filename      => 'signalp',
                           directory     => ${$args{ini_ref}}{signalp_dir}
                          );
   if ($$file_ref =~ m/SignalP 4/si) {
      $signalp_version = 4;
   } elsif ($$file_ref =~ m/SignalP 3/si) {
      $signalp_version = 3;
   } else {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Could not determine SignalP version. Please contact developer for bugfixing. Sorry mate.",
                                                    -bitmap  => 'info',
                                                    -buttons => ['OK'],
                                                   );
      &hide_pbar_2;
      return (0);
   }
   undef $file_ref;

   opendir SEQINPUT, ${$args{ini_ref}}{signalp_results};
   @signalp_results = grep !/^\./, readdir(SEQINPUT);

   #max number of fasta entries
   my $max_count = $#signalp_results + 1;
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'SignalP summary',
                   label        => 'SignalP summary'
                  );
   &show_pbar_2;

   #create SignalP summary file
   foreach my $dir (@signalp_results) {
      $progress++;
      if (($progress % 100) == 0) {
         &update_pbar_2(title        => "SignalP summary",
                        label        => "Creating SignalP summary file",
                        progress     => ($progress / $max_count) * 100,
                       );
      }

      if (-d ${$args{ini_ref}}{signalp_results}.'/'.$dir) {
         my $result_ref = &slurp(main_window => $args{main_window},
                                 directory   => ${$args{ini_ref}}{signalp_results}.'/'.$dir,
                                 filename    => 'output.txt'
                                );

         if ($result_ref eq '0') {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not READ SignalP result file in folder $dir",
                                                          -bitmap  => 'info',
                                                          -buttons => ['OK'],
                                                         );
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "SignalP:  Could not READ SignalP result file in folder $dir".
                           "\nReturn to main\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            #return (0); #return to main
            next; #try next entry
         }

         if ($signalp_version == 3) {
            #split results file in case both gram neg and gram pos was selected
            my @signalP_result = split /\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\s+SignalP/, ${$result_ref};

            foreach my $entry (@signalP_result) {
               next unless ($entry =~ m/\w+/);
               my $name       = '';
               my $pred       = '';
               my $prob       = '';
               my $cleav_prob = '';
               my $pos        = '';
               my $source     = '';

               'reset' =~ m/reset/;
               $entry =~ m/trained on (Gram\-\w+) bacteria.*\n\>([^\n]*?)\nPrediction\: ([^\n]*?)\nSignal peptide probability\: ([^\n]*?)\nMax cleavage site probability\: ([\d\.]*?) between pos\. ([\-\d]+) /s;
               ($source, $name, $pred, $prob, $cleav_prob, $pos) = ($1, $2, $3, $4, $5, $6);

               unless (defined $source && defined $name && defined $pred && defined $prob && defined $cleav_prob && defined $pos) {
                  my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                                -text    => "Error while creating SignalP summary file.\nCould not parse results for $dir",
                                                                -bitmap  => 'info',
                                                                -buttons => ['OK'],
                                                               );
                  $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
                  $error_msg-> Show();
                  open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
                  print ERRORLOG "SignalP:  Error while creating SignalP summary file.".
                                 "\nCould not parse results for $dir. Return to main\n\n";
                  close ERRORLOG;
                  &hide_pbar_2;
                  #return (0); #return to main
                  next; #try next
               }
               if ($pred eq 'Signal peptide') {
                  push (@summary, $name."\t".$prob."\t".$cleav_prob."\t".$pos."\t".$source."\n");
               } else {
                  push (@summary, $name."\t".$prob."\t---"."\t---"."\t".$source."\n");
               }
            }
            @signalP_result = ();
         } else {
            #SignalP present?
            if (${$result_ref} =~ m/SP\=\'NO'/) {
               ${$result_ref} =~ m/\nName=(.+?)SP='NO'\s+D=(\S+)\s+D-cutoff=(\S+)\sNetworks=(\S+)/s;
               my ($name, $D_score, $D_cutoff, $networks) = ($1, $2, $3, $4);
               $name =~ s/\s+$//;
               push (@summary, $name."\t".$D_score."\t".$D_cutoff."\t".$networks."\t---\tNO\n");
            } elsif  (${$result_ref} =~ m/SP\=\'YES'/) {
               ${$result_ref} =~ m/\nName=(.+?)SP='YES'\s+Cleavage site between pos.\s+\d+\s+and\s+(\d+)\:.+D=(\S+)\s+D-cutoff=(\S+)\sNetworks=(\S+)/s;
               my ($name, $pos, $D_score, $D_cutoff, $networks) = ($1, $2, $3, $4, $5);
               $name =~ s/\s+$//;
               push (@summary, $name."\t".$D_score."\t".$D_cutoff."\t".$networks."\t".$pos."\tYES\n");
            }
         }
      }
   }

   #sort array
   @summary = sort(@summary);

   #create summary file
   open WRITE, "+>".${$args{ini_ref}}{signalp_results}.'/SignalP_summary.txt';
   if ($signalp_version == 3) {
      print WRITE "Entry name\tSignal peptide probability\tMax cleavage site probability\tCleavage position\tSource type\n";
   } else {
      print WRITE "Entry name\tD-score\tD-score threshold\tNetworks\tCleavage position\tSignalPeptide\n";
   }
   foreach my $entry (@summary) {print WRITE $entry};
   close WRITE;

   &hide_pbar_2;
   return (1);
}

sub cleanup {
   my $dir = shift;
   local *DIR;

   opendir DIR, $dir or do {
      print "Error opening directory $dir: $!";
      return (0);
   };
   my $found = 0;
   while ($_ = readdir DIR) {
      next if /^\.{1,2}$/;
      my $path = "$dir/$_";
      cleanup($path) if -d $path;
      #unlink $path if -f $path;
      unlink $path
   }
   closedir DIR;
   rmdir $dir or do {
      print "Error deleting directory $dir: $!";
      return (0);
   };
}

1;