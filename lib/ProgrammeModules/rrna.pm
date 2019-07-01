#!/opt/ActivePerl-5.8/bin/perl

#ribosomal RNA uses Blast program calls
#input arguments: main_window, progress_bar, auto_ini_ref, ini_ref,
#                 directory, input_array_ref, output_folder, input_type

package ProgrammeModules::rrna;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&rrna);
use vars qw();

use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);
use Cwd;
use File::Copy;
use File::Find;

#local variables
my (%args, $newdir);

sub rrna {
   my %args = @_;
   my (%threshold, @templates, $counter, $rRNA);
   my @db = qw(5S_rRNA 16S_rRNA 23S_rRNA);

   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'ribosomal RNA analysis',
                   label        => 'ribosomal RNA analysis'
                  );
   &show_pbar_2;

   my $status = &test_databases(main_window  => $args{main_window},
                                progress_bar => $args{progress_bar},
                                auto_ini_ref => $args{auto_ini_ref},
                                ini_ref      => $args{ini_ref},
                               );
   if ($status == 0) {
      &hide_pbar_2;
      return(0);
   }

   #generate list of jobs to do
   foreach my $entry (@{ $args{input_array_ref} }) {
      foreach my $db (@db) {
         #test if re-use
         if (${$args{auto_ini_ref}}{reuse_results} == 1) {
            next if (-e ${$args{ini_ref}}{rrna_results}.'/'.$entry.'.'.$db.'.blastn' && -s ${$args{ini_ref}}{rrna_results}.'/'.$entry.'.'.$db.'.blastn' > 0);
            #add if required
            push (@templates, $entry.'___'.$db);
         }
         #add always if no re-use
         elsif (${$args{auto_ini_ref}}{reuse_results} == 0) {
            push (@templates, $entry.'___'.$db);
         }
      }
   }

   #max number of fasta entries
   my $max_count = $#templates + 1;
   if ($max_count < 1) {$max_count = 1};

   #start multithreading BlastN round for rRNA
   for (my $i=0; $i<=$#templates+1; $i = $i+${$args{auto_ini_ref}}{CPU}) {
      my $count = $i + ${$args{auto_ini_ref}}{CPU};
      if ($count > $#templates+1) {$count = $#templates+1};

      &update_pbar_2(title        => "ribosomal RNA analysis",
                     label        => "ribosomal RNA BLAST analysis  $count of $max_count",
                     progress     => ($count / $max_count) * 100,
                    );

      my @childs = ();
      for (my $j=$i; $j<$count; $j++) {
         my $seq_ref = "";
         #start forking
         my $pid = fork();
         if ($pid) {
            # parent
            push(@childs, $pid);
         } elsif ($pid == 0) {
            # child
            my ($filename, $db, $option);
            my (@list, $blastresult, $blast_hits_ref);

            #parse filename and database
            'reset' =~ m/reset/;
            $templates[$j] =~ m/(.+)___(.+)/;
            if ($1) {($filename, $db) = ($1, $2)} else {
               open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
               print WRITE "Error in ribosomal RNA analysis".
                           "\nError parsing job entry $templates[$j]\n\n";
               close WRITE;
               print "\nribosomal RNA analysis: Error parsing job entry $templates[$j]\n";
               CORE::exit();
            };

            #read fasta input file and match header to filename
            my ($seq_ref) = &slurp_cmd(auto_ini_ref => $args{auto_ini_ref},
                                       ini_ref      => $args{ini_ref},
                                       directory    => ${$args{ini_ref}}{input_files},
                                       filename     => $filename
                                      );
            if ($seq_ref eq '0') {
              #exit safely
               CORE::exit();
            }
            ${$seq_ref} =~ s/^>.*?\n/\>$filename \n/;

            #set word option
            if ($db eq '5S_rRNA') {
               $option .= '-W 8';
            } elsif ($db eq '16S_rRNA') {
               $option .= '-W 20';
            } elsif ($db eq '23S_rRNA') {
               $option .= '-W 25';
            }

            #blastN against nt sequence
            chdir ${$args{ini_ref}}{blast_executables};
            if (open RESULT, "-|") { # original process
               local $/;
               $blastresult = <RESULT>;
            } else { # child
               if (open STDIN, "-|") { # child
                  #exec "./blastall -p blastn -d ${$args{ini_ref}}{rrna_db_path}/$db -n T -I T -g T -m 9 -W 56 -F F -b 1000000 -v 1000000";
                  exec "./blastall -p blastn -d ${$args{ini_ref}}{rrna_db_path}/$db -n T -m 9 -g T -b 1000000 -v 1000000 $option ";
                  die "Cannot exec: $!";
               } else { # grandchild
                  print ${$seq_ref};
                  CORE::exit();
               }
            }
            chdir ${$args{auto_ini_ref}}{work_dir};

            # write Blast output file
            open WRITE, "+>${$args{ini_ref}}{rrna_results}/$filename\.$db\.blastn" or do {
               open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
               print WRITE "Error in ribosomal RNA analysis".
                           "\nError writing BlastN result file $filename\.$db\.blastn to directory ${$args{ini_ref}}{rrna_results}\n\n";
               close WRITE;
               print "\nribosomal RNA analysis: Error writing Blast result file $filename\.$db\.blastn to directory ${$args{ini_ref}}{rrna_results}\n";
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

   #parse blast results
   &update_pbar_2(title        => "ribosomal RNA analysis",
                  label        => "Parsing all results",
                  progress     => 1,
                 );
   $counter = 0;
   foreach my $entry (@templates) {
      my ($filename, $db);
      &update_pbar_2(title        => "ribosomal RNA analysis",
                     label        => "Parsing all results",
                     progress     => ($counter / $max_count) * 100,
                    );
      $counter++;

      #parse filename and database
      'reset' =~ m/reset/;
      $entry =~ m/(.+)___(.+)/;
      if ($1) {($filename, $db) = ($1, $2)} else {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not parse entry $entry.\nAborting",
                                           -type    => 'OK',
                                           -icon    => 'info');
         &hide_pbar_2;
         return (0);
      }

      #read result Blast files
      open READ, "<${$args{ini_ref}}{rrna_results}/$filename\.$db\.blastn" or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not open file $filename\.$db\.blastn.\nAborting",
                                           -type    => 'OK',
                                           -icon    => 'info');
         &hide_pbar_2;
         return (0);
      };

      #create entries
      while (<READ>) {
         next if ($_ =~ /^\#/); #ignore header
         my @tmp = split/\t/,$_;

         #cleanup acc
         $tmp[1] =~ s/^gi\|(.+?)\|.+/$1/;

         #does entry with same boundary already exists?
         my $update;
         if (exists $rRNA->{$filename}->{$db}->{$tmp[6]}) {
            #longer alignment?
            if ($tmp[3] > $rRNA->{$filename}->{$db}->{$tmp[6]}->{'length'}) {
               $update = 1;
            }
            #better score?
            if ($tmp[11] > $rRNA->{$filename}->{$db}->{$tmp[6]}->{'score'}) {
               $update = 1;
            }
            next unless ($update);
         }

         #create new entry
         $rRNA->{$filename}->{$db}->{$tmp[6]} = { 'db'       => $db,
                                                  'acc'      => $tmp[1],
                                                  'identity' => $tmp[2],
                                                  'length'   => $tmp[3],
                                                  'qu_start' => $tmp[6],
                                                  'qu_stop'  => $tmp[7],
                                                  'score'    => $tmp[11]
                                                };
         undef $update;
      }
      close READ;

      #iterate through relevant entries and delete all overlaps.
      #Entries will be used as template only for subsequent bl2seq analysis
      &update_pbar_2(title        => "ribosomal RNA analysis",
                     label        => "Iterating through candidate domains",
                     progress     => 1,
                    );
      &iterate(rRNA_ref => \$rRNA,
               counter  => 1,
               db       => $db,
               filename => $filename
              );
   }
   undef @templates;

   #create job array
   foreach my $filename (keys %{ $rRNA }) {
      foreach my $rRNA_db (keys %{ $rRNA->{$filename} }) {
         foreach my $start (keys %{ $rRNA->{$filename}->{$rRNA_db} }) {
            next unless ($start =~ /\d+/);
            #skip entry if already exists and re-use selected
            next if (${$args{auto_ini_ref}}{reuse_results} == 1 &&
                     -e ${$args{ini_ref}}{rrna_results}.'/'.$filename.'-'.$rRNA_db.'-'.$start.'.bl2seq' &&
                     -s ${$args{ini_ref}}{rrna_results}.'/'.$filename.'-'.$rRNA_db.'-'.$start.'.bl2seq' > 0);
            push (@templates, $filename.'___'.$rRNA_db.'___'.$start);
         }
      }
   }

   #new max_count
   $max_count = $#templates + 2;

   #start multithreading bl2seq round for rRNA
   for (my $i=0; $i<=$#templates+1; $i = $i+${$args{auto_ini_ref}}{CPU}) {
      my $count = $i + ${$args{auto_ini_ref}}{CPU};
      if ($count > $#templates+1) {$count = $#templates+1};

      &update_pbar_2(title        => "ribosomal RNA analysis",
                     label        => "ribosomal RNA pairwise analysis  $count of $max_count",
                     progress     => ($count / $max_count) * 100,
                    );

      my @childs = ();
      for (my $j=$i; $j<$count; $j++) {
         my $seq_ref = "";
         #start forking
         my $pid = fork();
         if ($pid) {
            # parent
            push(@childs, $pid);
         } elsif ($pid == 0) {
            # child
            my ($filename, $rRNA_db, $qu_start, $blastresult, $boundary_length);
            #process ID
            my $ID = $$;

            #parse filename and database
            'reset' =~ m/reset/;
            $templates[$j] =~ m/(.+)___(.+)___(.+)/;
            if ($3) {($filename, $rRNA_db, $qu_start) = ($1, $2, $3)} else {
               open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
               print WRITE "Error in ribosomal RNA analysis".
                           "\nError parsing job entry $templates[$j]\n\n";
               close WRITE;
               print "\nribosomal RNA analysis: Error parsing job entry $templates[$j]\n";
               CORE::exit();
            };

            #set  boundary length
            if ($rRNA_db eq '5S_rRNA') {
               $boundary_length = 200;
            } elsif ($rRNA_db eq '16S_rRNA') {
               $boundary_length = 1600;
            } elsif ($rRNA_db eq '23S_rRNA') {
               $boundary_length = 3000;
            }

            #read fasta input file and match header to filename
            my ($seq_ref) = &slurp_cmd(auto_ini_ref => $args{auto_ini_ref},
                                       ini_ref      => $args{ini_ref},
                                       directory    => ${$args{ini_ref}}{input_files},
                                       filename     => $filename
                                      );
            if ($seq_ref eq '0') {
              #exit safely
               CORE::exit();
            }
            ${$seq_ref} =~ s/^>.*?\n//s;
            ${$seq_ref} =~ s/\s+//gs;

            #get query sequence
            my $qu_right_bd = $rRNA->{$filename}->{$rRNA_db}->{$qu_start}->{'qu_stop'} + $boundary_length;
            my $qu_left_bd  = $qu_start - $boundary_length;
            if ($qu_left_bd < 1) {$qu_left_bd = 1};
            my $query_seq = substr(${$seq_ref}, $qu_left_bd, ($qu_right_bd - $qu_left_bd + $boundary_length));
            open WRITE, "+>${$args{ini_ref}}{rrna_results}/$ID\.query";
            print WRITE "\>$filename - $qu_left_bd\-$qu_right_bd \n$query_seq\n";
            close WRITE;

            #retrieve subject sequence from database
            my $subject_seq = `${$args{ini_ref}}{blast_executables}/fastacmd -d ${$args{ini_ref}}{rrna_db_path}/$rRNA_db -p F -s $rRNA->{$filename}->{$rRNA_db}->{$qu_start}->{'acc'}`;
            unless ($subject_seq) {
               open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
               print WRITE "Error in ribosomal RNA analysis".
                           "\nError retrieving subject sequence for entry $templates[$j]\n\n";
               close WRITE;
               print "\nribosomal RNA analysis: Error retrieving subject sequence for entry $templates[$j]\n...$rRNA->{$filename}->{$rRNA_db}->{$qu_start}->{'acc'} ...\n";
               #remove temp files
               unlink ${$args{ini_ref}}{rrna_results}.'/'.$ID.'.query';

               CORE::exit();
            };
            open WRITE, "+>${$args{ini_ref}}{rrna_results}/$ID\.subject";
            print WRITE $subject_seq."\n";
            close WRITE;


            #bl2seq against nt sequence
            chdir ${$args{ini_ref}}{blast_executables};
            $blastresult = `./bl2seq -p blastn -F F -m F -i ${$args{ini_ref}}{rrna_results}\/$ID\.query -j ${$args{ini_ref}}{rrna_results}\/$ID\.subject`;
            chdir ${$args{auto_ini_ref}}{work_dir};

            # write Blast output file
            open WRITE, "+>${$args{ini_ref}}{rrna_results}/$filename\-$rRNA_db\-$qu_start\.bl2seq" or do {
               open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
               print WRITE "Error in ribosomal RNA analysis".
                           "\nError writing BlastN result file $filename\.$rRNA_db\.blastn to directory ${$args{ini_ref}}{rrna_results}\n\n";
               close WRITE;
               print "\nribosomal RNA analysis: Error writing Blast result file $filename\.$rRNA_db\.blastn to directory ${$args{ini_ref}}{rrna_results}\n";
               #remove temp files
               unlink ${$args{ini_ref}}{rrna_results}.'/'.$ID.'.query';
               unlink ${$args{ini_ref}}{rrna_results}.'/'.$ID.'.subject';

               CORE::exit();
            };
            print WRITE $blastresult;
            close WRITE;

            #remove temp files
            unlink ${$args{ini_ref}}{rrna_results}.'/'.$ID.'.query';
            unlink ${$args{ini_ref}}{rrna_results}.'/'.$ID.'.subject';

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

   &hide_pbar_2;
   return (1);
}

sub iterate {
   my %args = @_;
   my (@delete);
   &update_pbar_2(title        => "ribosomal RNA analysis",
                  label        => "Processing candidate domains, Iteration $args{counter} for $args{db}",
                  progress     => $args{counter},
                 );
   $args{counter} = 0 if (($args{counter} ) >= 100);
   $args{counter}++;

   foreach my $qu_start (keys %{ ${$args{rRNA_ref}}->{$args{filename}}->{$args{db}} } ) {
      my $qu_stop = ${$args{rRNA_ref}}->{$args{filename}}->{$args{db}}->{$qu_start}->{'qu_stop'};
      foreach my $sb_start (keys %{ ${$args{rRNA_ref}}->{$args{filename}}->{$args{db}} } ) {
         #ignore self
         next if ($qu_start == $sb_start);
         my $sb_stop = ${$args{rRNA_ref}}->{$args{filename}}->{$args{db}}->{$sb_start}->{'qu_stop'};

         if (($sb_start > $qu_start && $sb_start < $qu_stop) ||
             ($sb_stop  > $qu_start && $sb_stop  < $qu_stop) ||
             ($sb_start < $qu_start && $sb_start > $qu_stop)) {
            #which one to delete?
            #test score
            if (${$args{rRNA_ref}}->{$args{filename}}->{$args{db}}->{$qu_start}->{'score'} >= ${$args{rRNA_ref}}->{$args{filename}}->{$args{db}}->{$sb_start}->{'score'}) {
               push (@delete, $sb_start);
            } else {
               push (@delete, $qu_start);
            }
         }
      }
   }

   #delete tagged entries
   foreach (@delete) {
      delete ${$args{rRNA_ref}}->{$args{filename}}->{$args{db}}->{$_};
   }

   #iterate if tagged hits
   if ($#delete >= 0) {
      &iterate(rRNA_ref => $args{rRNA_ref},
               counter  => $args{counter},
               db       => $args{db},
               filename => $args{filename}
              );
   }
   undef @delete;

   return (1);
}

sub test_databases {
   my %args = @_;
   my ($recompile, @cm_dataset);
   my @db_extension = qw(nhr nin nsd nsi nsq);
   my @db           = qw(5S_rRNA 16S_rRNA 23S_rRNA);

   #test if all databases exists
   foreach my $set (@db) {
      foreach my $extension (@db_extension) {
         unless (-e ${$args{ini_ref}}{rrna_db_path}.'/'.$set.'.'.$extension) {
            $recompile = 1;
            last;
         }
      }
   }

   #recomile if necessary
   if ($recompile) {
      &update_pbar_2(title        => "ribosomal RNA analysis",
                     label        => "Recompiling rRNA databases ",
                     progress     => 1,
                    );
      my ($status) = &recompile_rrna_db(main_window  => $args{main_window},
                                        progress_bar => $args{progress_bar},
                                        auto_ini_ref => $args{auto_ini_ref},
                                        ini_ref      => $args{ini_ref},
                                        directory    => ${$args{ini_ref}}{rrna_db_path}
                                       );
      if ($status == 0) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not compile rRNA databases in directory ${$args{auto_ini_ref}}{rrna_db_path}",
                                           -type    => 'OK',
                                           -icon    => 'info');
         &hide_pbar_2;
         return (0); #exit to main
      }
   }

   return (1);
}

sub recompile_rrna_db {
   my %args = @_;
   my ($rRNA_archive, $short_name);
   my @db = qw(5S_rRNA 16S_rRNA 23S_rRNA);

   #find current Blast zip file in Archive directory
   &update_pbar_2(title        => 'Recompiling rRNA databases',
                  label        => 'Finding rRNA archive',
                  progress     => 25,
                 );
   find (sub {$rRNA_archive = $File::Find::name if (/rRNA.*?\.zip/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );

   unless (defined $rRNA_archive && $rRNA_archive =~ /\w+/) {
      return (0);
   }
   $short_name = $rRNA_archive;
   $short_name =~ s/.+\/(.+)$/$1/;

   #copy rRNA archive to rRNA db directory
   copy($rRNA_archive , $args{directory}.'/'.$short_name);

   #uncompress file
   &update_pbar_2(title        => 'Recompiling rRNA databases',
                  label        => 'Uncompressing rRNA archive',
                  progress     => 45,
                 );
   chdir $args{directory};
   `unzip $short_name`;

   #formatting databases
   my $progress = 45;
   foreach my $db (@db) {
      $progress += 15;
      &update_pbar_2(title        => 'Recompiling rRNA databases',
                     label        => "Recompiling $db database",
                     progress     => $progress,
                    );
      #format db
      `${$args{ini_ref}}{blast_executables}/formatdb -i $db -p F -o T -n $args{directory}/$db`;
      #remove source fasta file
      unlink $args{directory}.'/'.$db;
   }

   #cleanup
   unlink $args{directory}.'/'.$short_name;
   unlink $args{directory}.'/formatdb.log';

   #make sure everything is accessible and executable
   `chmod -R 777 *`;

   #set compiler flag
   open WRITE, '+>'.$args{directory}.'/recompiled';
   print WRITE $args{directory};
   close WRITE;

   return (1);
}

1;