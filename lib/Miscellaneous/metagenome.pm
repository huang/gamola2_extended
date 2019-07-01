#!/opt/ActivePerl-5.8/bin/perl
#transfer annotation from one GAMOLA genbank file
#to new one
#changes are logged in separate file

package Miscellaneous::metagenome;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&metagenome_ui);
use vars qw();

use Tk;
use Tk::Pane;
use initialise::read_me              qw(:DEFAULT);
use ProgrammeModules::sequence       qw(:DEFAULT);
use ProgrammeModules::genbank_parser qw(:DEFAULT);
use ProgrammeModules::blast          qw(:DEFAULT);
use Basics::progress_bar             qw(:DEFAULT);
use db_selection::databases          qw(:DEFAULT);

use Cwd;

#local variables
my (%args, $tl, $new_width, $old_width, $current_dir,
    @new_selected_files, %new_all_files, @new_short_files,
    $meta_collapse, $meta_blast, $meta_cluster, $blast_flavour,
    $blast_cluster, $blast_threaded, $full_blast_db, $short_blast_name,
    $blast_rank, $blast_cutoff, $blast_display_limit);
my (%meta_seqs);

#remember last directory for adding files
$current_dir = ${$args{auto_ini_ref}}{work_dir};

sub metagenome_ui {
   my %args = @_;
   my ($progress_bar, $frame_left, $frame_right, $frame_bottom, $start_frame,
       $new_anno_lb, $new_width, $new_add_files, $new_clear,
       );

   #set defaults
   %meta_seqs        = ();
   $new_width        = 20;
   @new_selected_files  = ();
   %new_all_files       = ();
   @new_short_files     = ();
   unless ($short_blast_name =~ /\S+/) {
      $short_blast_name = ${$args{auto_ini_ref}}{blast_db};
      $full_blast_db    = ${$args{auto_ini_ref}}{full_blast_db};
   }
   $meta_collapse       = 1;
   $meta_blast          = 1;
   $meta_cluster        = 1000;
   $blast_threaded      = 1;
   $blast_flavour       = 'BlastX';
   $blast_rank          = 1;
   $blast_cutoff        = '1e-10';
   $blast_display_limit = 10;

   #set up GUI
   $tl = ${$args{main_window}}->DialogBox(-title   => 'Metagenome analyses',
                                          -buttons => [ 'EXIT' ]
                                         );
   $tl->minsize(20, 10);
   $tl->maxsize(120, 20);
   $frame_bottom  = $tl->Frame(-borderwidth => 2, -relief => 'groove') -> pack(-side => 'bottom', -fill => 'x');
   $frame_left    = $tl->Frame(-relief => 'raised')                    -> pack(-side => 'left',   -fill => 'both', -expand => 1);
   $frame_right   = $tl->Frame(-relief => 'raised')                    -> pack(-side => 'left',   -fill => 'both', -expand => 1);

   $frame_bottom -> Button ( -text    => 'Start Metagenome analysis',
                             -command => sub {&meta_analysis(main_window         => \$tl,
                                                             progress_bar        => \$progress_bar,
                                                             ini_ref             => $args{ini_ref},
                                                             auto_ini_ref        => $args{auto_ini_ref},
                                                             blast_display_limit => $blast_display_limit,
                                                            );
                                              #reset file list
                                              %new_all_files = ();
                                              @new_short_files = ();
                                              $new_anno_lb->delete(0.0, 'end');
                                              $tl->update;

                                              undef @new_selected_files;
                                              undef %new_all_files;
                                              undef @new_short_files;
                                             }
                            )
                                                                               -> pack (-side => 'left');

   $progress_bar = $frame_bottom->Label(-text => " ")                          -> pack (-side => 'left', -fill => 'y'); #create dummy for space

   #left frame -> options
   {
      my $frame_note    = $frame_left->Frame() -> pack(-side => 'top',  -fill => 'x');
      my $frame_options = $frame_left->Frame() -> pack(-side => 'left', -fill => 'both', -expand => 1);

      $frame_note -> Label       (-text     => "Metagenome options\n")
         ->grid(-row => 0, -column => 0, -columnspan => 2, -sticky => 'ew');

      $frame_note -> Label       (-text     => "Note: Only Nucleotide sequences are assumed\n"
                                 )
         ->grid(-row => 1, -column => 0, -columnspan => 2, -sticky => 'w');

      $frame_options -> Checkbutton (-text     => 'Re-use previous results',
                                     -variable => \${$args{auto_ini_ref}}{reuse_results}
                                     )
         ->grid(-row => 1, -column => 0, -sticky => 'nw');
      $frame_options -> Checkbutton (-text     => 'Collapse identical reads',
                                     -variable => \$meta_collapse
                                     )
         ->grid(-row => 2, -column => 0, -sticky => 'nw');
      $frame_options -> Label       (-text     => 'Blast analysis',
                                     )
         ->grid(-row => 3, -column => 0, -sticky => 'nw');
      $frame_options -> BrowseEntry  (-label    => 'Cluster reads',
                                      -choices  => [1..1000],
                                      -width    => 5,
                                      -variable => \$meta_cluster
                                     )
         ->grid(-row => 3, -column => 1, -sticky => 'nw');
      $frame_options -> BrowseEntry  (-label    => 'Blast Flavour',
                                      -choices  => ['BlastN', 'BlastX', 'tBlastX'],
                                      -width    => 11,
                                      -variable => \$blast_flavour
                                     )
         ->grid (-row => 4, -column => 1, -sticky => 'nw');
      $frame_options -> Checkbutton(-text     => 'Cluster CPUs for Blast',
                                    -variable => \$blast_cluster,
                                    -command  => sub {
                                                      if ($blast_cluster == 1) {
                                                         $blast_threaded = 0;
                                                      } elsif ($blast_cluster == 0) {
                                                         $blast_threaded = 1;
                                                      }
                                                     }
                                   )
         ->grid (-row => 5, -column => 1, -sticky => 'nw');
      $frame_options  -> Checkbutton(-text     => 'Threaded Blast',
                                     -variable => \$blast_threaded,
                                     -command  => sub {
                                                       if ($blast_threaded == 1) {
                                                          $blast_cluster = 0;
                                                       } elsif ($blast_threaded == 0) {
                                                          $blast_cluster = 1;
                                                       }
                                                      }
                                    )
         ->grid (-row => 5, -column => 2, -sticky => 'nw');
      $frame_options     -> Button(-text    => 'Blast db',
                                   -width   => 18,
                                   -height  => 1,
                                   -command => sub {
                                                     my ($select_full_db, $select_short_name) = &select_meta_db(main_window    => $args{main_window},
                                                                                                                ini_ref        => $args{ini_ref},
                                                                                                                auto_ini_ref   => $args{auto_ini_ref},
                                                                                                                ini_dir        => $blast_flavour
                                                                                                                );
                                                     if (defined $select_full_db && $select_full_db =~ m/\w+/) {
                                                        $full_blast_db    = $select_full_db;
                                                        $short_blast_name = $select_short_name;
                                                     }
                                                      $tl->update;
                                                    }
                                   )
         ->grid (-row => 6, -column => 1, -sticky => 'nsw');
      $frame_options     -> Label(-justify      => 'left',
                                  -relief       => 'groove',
                                  -borderwidth  => 2,
                                  -wraplength   => 100,
                                  -width        => 16,
                                  -height       => 2,
                                  -textvariable => \$short_blast_name,
                                  )
         ->grid (-row => 6, -column => 2, -sticky => 'new');

      $frame_options     -> Checkbutton(-text     => 'Filter and rank results',
                                        -variable => \$blast_rank,
                                       )
         ->grid (-row => 7, -column => 1, -sticky => 'nw');

      $frame_options     -> BrowseEntry  (-label    => 'e-value      ',
                                          -choices  => [1, 1e-10, 1e-20, 1e-30, 1e-40, 1e-60, 1e-80, 1e-100, 1e-120, 0 ],
                                          -width    => 7,
                                          -variable => \$blast_cutoff
                                         )
         ->grid (-row => 7, -column => 2, -sticky => 'nw');

      $frame_options     -> BrowseEntry  (-label    => 'Display limit',
                                          -choices  => [1..100],
                                          -width    => 7,
                                          -variable => \$blast_display_limit
                                         )
         ->grid (-row => 8, -column => 2, -sticky => 'nw');
   }

   #right frame -> metagenome seq data files
   {
      $frame_right                       -> Label (-text => "Input files\n")
         ->grid(-row => 0, -column => 0, -columnspan => 2, -sticky => 'nsew');
      $new_add_files = $frame_right      ->Button(-text => 'Add files',
                                                  -command => sub {&populate_new(main_window  => \$tl,
                                                                                 auto_ini_ref => $args{auto_ini_ref},
                                                                                 textbox      => $new_anno_lb,
                                                                                )
                                                                  }
                                                 )
         ->grid(-row => 1, -column => 0);

      $frame_right                       ->Label(-text => 'Double right-click entry to remove.')
         ->grid(-row => 1, -column => 1, -sticky => 'nsew');
      $new_anno_lb         = $frame_right->Scrolled("Listbox",
                                                    -setgrid    => 1,
                                                    -selectmode => 'single',
                                                    -width      => $new_width
                                                    )
         ->grid(-row => 2, -column => 0, -columnspan => 2, -sticky => 'nsew');
      $new_clear         = $frame_right  ->Button(-text    => 'Clear current selection',
                                                  -command => sub {
                                                                    %new_all_files = ();
                                                                    @new_short_files = ();
                                                                    $new_anno_lb->delete(0.0, 'end');
                                                                    $tl->update;
                                                                   }
                                                 )
      ->grid(-row => 3, -column => 0, -columnspan => 2, -sticky => 'nsew');
      $new_anno_lb                       ->configure(-width   => $new_width,
                                                     -setgrid => 1 );
   }

   #bind left and right mouse button to selection and removal of entries
   $new_anno_lb->bind('<Double-ButtonRelease-3>' => sub {&remove_entry_new(textbox      => \$new_anno_lb,
                                                                           main_window  => \$tl,
                                                                           progress_bar => \$progress_bar,
                                                                           ini_ref      => $args{ini_ref},
                                                                           auto_ini_ref => $args{auto_ini_ref}
                                                                          )
                                                        });
   my $wait = $tl->Show();
   if ($wait eq 'EXIT') {
      $tl->state('withdrawn');
   }

}

sub cleanup {
   my $dir = shift;
   local *DIR;

   opendir DIR, $dir or die "opendir $dir: $!";
   my $found = 0;
   while ($_ = readdir DIR) {
      next if /^\.{1,2}$/;
      my $path = "$dir/$_";
      unlink $path if -f $path;
      cleanup($path) if -d $path;
   }
   closedir DIR;
   rmdir $dir or print "error - $!";
}

sub meta_analysis {
   my %args = @_;
   my ($max_count, $count, $meta_res_ref, $total_seq_count);
   my (%temp_header, %unique, %final_entries, %candidate_lookup);

   #return if no files selected
   $max_count = scalar (keys %new_all_files);
   if ( $max_count < 1) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Insufficient file selection",
                                        -icon    => 'error',
                                        -type    => 'ok');
      return;
   };

   #initialise progress bar
   ${$args{progress_bar}}->configure(-text => "Reading input data in FASTA format");
   &progress_bar_1(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'Reading input data',
                   label        => 'Input data'
                  );
   &show_pbar_1;

   #iterate through input files
   %meta_seqs = ();
   while (my ($key, $value) = each %new_all_files) {
      $count++;
      &update_pbar_1(title        => "Reading input data",
                     label        => "Reading $key",
                     progress     => ($count / $max_count) * 100,
                    );
      my ($file, $dir);
      my (@fasta_entries);
      $value =~ /^(.*)\/(.+)$/;
      ($dir, $file) = ($1, $2);

      unless (defined $dir) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse entry $value",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files in Metagenome Analysis".
                        "\nCould not parse entry $value\n\n";
         close ERRORLOG;
         return (0);
      }

      #read file into memory
      my ($file_ref) = &slurp(main_window   => $args{main_window},
                              auto_ini_ref  => $args{auto_ini_ref},
                              ini_ref       => $args{ini_ref},
                              filename      => $file,
                              directory     => $dir
                             );
      $$file_ref = "\n".$$file_ref;

      #put all fasta entries in hash
      @fasta_entries = split/\n\>/, $$file_ref;

      foreach my $entry (@fasta_entries) {
         my ($header, $seq);
         next unless ($entry =~ /\S+/);

         #count total number of seq reads
         $total_seq_count++;

         $entry =~ m/^(\S*?)[\s\n]([^\n]*?\n)?(.+)/s;
         ($header, $seq) = ($1, $3);
         unless (defined $seq) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse entry \n$entry \nin file $file",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error processing input files in Metagenome Analysis".
                           "\nCould not parse entry \n$entry \nin file $file\n\n";
            close ERRORLOG;
            &hide_pbar_1;
            return (0);
         }
         #clean seq
         $seq =~ s/\s//gs;

         #put into hash
         $meta_seqs{$total_seq_count} = {(header => $header,
                                          seq    => $seq
                                        )};
      }
   }

   #create a lookup hash for potentially identical sequences
   #this will drastcially reduce the searchspace for result parsing and filtering
   ${$args{progress_bar}}->configure(-text => "Creating lookup hash");
   my %temp_seq;
   $count = 1;
   $max_count = scalar (keys %meta_seqs);
   foreach my $ID (keys %meta_seqs) {
      $count++;
      if ($count % 1000 == 0) {
         &update_pbar_1(title        => "Reading input data",
                        label        => "Creating lookup table",
                        progress     => ($count / $max_count) * 100,
                       );
      }

      #find potential candidate duplicates
      $candidate_lookup{$meta_seqs{$ID}->{header}} = $meta_seqs{$ID}->{seq};
   }
   undef %temp_seq;

   #collapse identical entries if selected
   %final_entries = ();
   if ($meta_collapse == 1) {
      my %temp_seq = ();
      my $count = 0;
      my $max_count = scalar (keys %meta_seqs);
      ${$args{progress_bar}}->configure(-text => "Collapsing identical data reads");
      #going through dataset and concatenating header if sequences are identical
      foreach my $ID (keys %meta_seqs) {
         my $short_seq = substr($meta_seqs{$ID}->{seq},0,200); #grab first 200 nt of seq to stay within hash key limit
         $count++;
         if ($count % 1000 == 0) {
            &update_pbar_1(title        => "Reading input data",
                           label        => "Collapsing data entries",
                           progress     => ($count / $max_count) * 100,
                          );
         }

         #find potential candidate duplicates
         if ($temp_seq{$short_seq}) {
            if (defined $unique{$short_seq}) {
               $temp_header{$short_seq} .= $unique{$short_seq}.'___';
            }
            $temp_header{$short_seq} .= $ID.'___';
            delete $unique{$short_seq}; #delete false unique
         } else {
            #initialise potential unique hits
            $unique{$short_seq} = $ID;
         }
         $temp_seq{$short_seq}++;
      }
      undef %temp_seq;

      #investigate candidate duplicates
      $count = 1;
      my $counter = 0;
      $max_count = scalar (keys %temp_header);
      foreach my $short_seq (keys %temp_header) {
         my (@individual_ID, @combined_header, %compare, $final_header, %seen);
         $counter++;
         if ($counter % 1000 == 0) {
            &update_pbar_1(title        => "Duplicate data entries",
                           label        => "Testing for duplicates",
                           progress     => ($counter / $max_count) * 100,
                          );
         }

         #get individual headers
         @individual_ID = ();
         @individual_ID = split/___/, $temp_header{$short_seq};

         #if only 2 entries make direct comparison instead of iterating
         if ($#individual_ID == 1) {
            if ($meta_seqs{$individual_ID[0]}->{seq} eq $meta_seqs{$individual_ID[1]}->{seq}) {
               my @temp = ($meta_seqs{$individual_ID[0]}->{header}, $meta_seqs{$individual_ID[1]}->{header});
               @temp = sort @temp;
               $final_header = $temp[0].'___'.$temp[1];
               undef @temp;
               while (exists $final_entries{$count}->{header}) {$count++};

               $final_entries{$count} = {(header => $final_header,
                                          seq    => $meta_seqs{$individual_ID[0]}->{seq}
                                        )}; #add to final hash
               $count++;
               next;
            } else {
               foreach my $ID (@individual_ID) {
                  while (exists $final_entries{$count}->{header}) {$count++};
                  $final_entries{$count} = {(header => $meta_seqs{$ID}->{header},
                                             seq    => $meta_seqs{$ID}->{seq}
                                           )};
                  $count++;
               }
               next;
            }
         }

         #put into hash
         %compare = ();
         foreach my $ID (@individual_ID) {
            $compare{$ID}++
         }

         #find duplicates
         %seen = ();
         while (scalar (keys %compare) >= 1) {
            #only one hit left?
            if (scalar (keys %compare) == 1) {
               foreach my $ID (keys %compare) {
                  while (exists $final_entries{$count}->{header}) {$count++};
                  $final_entries{$count} = {(header => $meta_seqs{$ID}->{header},
                                             seq    => $meta_seqs{$ID}->{seq}
                                           )}; #add to final hash
               }
               last;
            }

            #iterate through duplicate entries
            foreach my $ID (keys %compare) {
               next unless ($ID =~ /\S+/);
               #initialise final header
               @combined_header = ();
               push (@combined_header, $meta_seqs{$ID}->{header});
               $seen{$ID}++;

               #start comparison to other candidates
               foreach my $comp (keys %compare) {
                  next unless ($comp =~ /\S+/);
                  next if ($ID eq $comp); #skip idential headers
                  if ($meta_seqs{$ID}->{seq} eq $meta_seqs{$comp}->{seq}) {
                     #remove identical hits from search list
                     $seen{$comp}++;
                     #add to combined header
                     push (@combined_header, $meta_seqs{$comp}->{header});
                  }
               }

               #sort final header
               @combined_header = sort @combined_header;
               $final_header = '';
               $final_header = join('___', @combined_header);
               $final_header =~ s/___$//;
               $final_header =~ s/^___//;

               #find unique identifier
               while (exists $final_entries{$count}->{header}) {$count++};
               $final_entries{$count} = {(header => $final_header,
                                          seq    => $meta_seqs{$ID}->{seq}
                                        )}; #add to final hash

               #exit to delete entries
               last;
            }

            #delete duplicate entries and current entry before next iteration
            foreach my $ID (keys %seen) {
               delete $compare{$ID};
            }
            undef %seen;
         }
         undef @individual_ID;
         undef %compare;
      }
      undef %temp_header;

      #add initial unique entries
      $count = 1;
      $max_count = scalar (keys %unique);
      foreach my $short_seq (keys %unique) {
         $counter++;
         if ($counter % 1000 == 0) {
            &update_pbar_1(title        => "Duplicate data entries",
                           label        => "Adding unique entries",
                           progress     => ($count / $max_count) * 100,
                          );
         }
         while (exists $final_entries{$count}->{header}) {$count++};

         $final_entries{$count} = {(header => $meta_seqs{$unique{$short_seq}}->{header},
                                    seq    => $meta_seqs{$unique{$short_seq}}->{seq}
                                  )}; #add to final hash
         $count++;
      }
      undef %unique;

      #re-initialise origine sequence hash
      %meta_seqs = %final_entries;
      undef %final_entries;
   }
   &hide_pbar_1;

   #run Blast
   ${$args{progress_bar}}->configure(-text => "Performing Blast: $blast_flavour");
   &blast_hash(main_window         => $args{main_window},
               progress_bar        => $args{progress_bar},
               auto_ini_ref        => $args{auto_ini_ref},
               ini_ref             => $args{ini_ref},
               blast_display_limit => $args{blast_display_limit},
               seq_hash            => \%meta_seqs,
               candidate_lookup    => \%candidate_lookup,
               cluster_seq         => $meta_cluster,
               cluster_blast       => $blast_cluster,
               blast_db            => $full_blast_db,
               blast_flavour       => $blast_flavour,
              );

   #filter Blast results if selected
   if ($blast_rank == 1) {
      ${$args{progress_bar}}->configure(-text => "Filtering results");
      $meta_res_ref = &filter_results(main_window      => $args{main_window},
                                      progress_bar     => $args{progress_bar},
                                      auto_ini_ref     => $args{auto_ini_ref},
                                      ini_ref          => $args{ini_ref},
                                      candidate_lookup => \%candidate_lookup,
                                      cluster_seq      => $meta_cluster,
                                      cluster_blast    => $blast_cluster,
                                      blast_db         => $full_blast_db,
                                      blast_flavour    => $blast_flavour,
                                      total_seq_count  => $total_seq_count
                                    );
   }

   #done
   ${$args{progress_bar}}->configure(-text => "Analysis finished");
   ${$args{main_window}}->messageBox(-title   => 'Info',
                                     -message => "Metagenome analysis completed",
                                     -icon    => 'info',
                                     -type    => 'ok');

   &hide_pbar_1;
   return(1);
}


sub filter_results {
   my %args = @_;
   my ($counter, $max_count, $individual_result, $entry_count);
   my (@res_list, %results, @sort_list);

   &show_pbar_1;

   #read results file line by line insted
   open READ, '<'.${$args{ini_ref}}{blast_results}.'/metagenome.'.$args{blast_flavour};
   $counter     = 1;
   $entry_count = 1;
   while (<READ>) {
      if ($counter % 100000 == 0) {
         &update_pbar_1(title        => "Filtering Blast results",
                        label        => "Reading results file, entry $entry_count",
                        progress     => ($counter / 10000000) * 100,
                       );
         if ($counter >= 10000000) {$counter = 1};
      }
      $counter++;
      if ($individual_result =~ m/\n\n\n$/) {
         push (@res_list, $individual_result);
         $individual_result = '';
         $entry_count++;
      }
      $individual_result .= $_;
   }
   push (@res_list, $individual_result);
   undef $individual_result;
   undef $entry_count;
   close READ;


   #create individual entries
   #@res_list  = split/Blast Flavour:\t/,$$seq_ref;
   #undef $seq_ref;

   #create results hash based on query name, best blast evalue and list of results
   $counter = 1;
   $max_count = $#res_list + 2;
   foreach my $entry (@res_list) {
      next unless ($entry =~ /\S+/);
      if ($counter % 100 == 0) {
         &update_pbar_1(title        => "Filtering Blast results",
                        label        => "Iterating through results",
                        progress     => ($counter / $max_count) * 100,
                       );
      }

      my ($query, $blast_db, $evalue, $result, @blast_limit);
      'reset' =~ m/reset/;
      $entry  =~ m/Blast db:\t(.*?)\nQuery name:\s+(.*?)\nBlast overview\s+(>.+)/s;
      ($blast_db, $query, $result) = ($1, $2, $3);
      unless (defined $result) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse entry \n$entry \nThe entry will be ignored and analysis continues",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files in Metagenome Analysis".
                        "\nCould not parse entry \n$entry \n\n";
         close ERRORLOG;
         next;
      }
      #clean up query
      $query =~ s/\s+$//;

      #clean up blast db
      if ($blast_db =~ m/\//) {
         $blast_db =~ s/^.*\/(.+)$/$1/;
      }

      #get best evalue
      if ($result =~ m/>No Blast Hits found/) {
         $evalue  = 100;
      } else {
         'reset' =~ m/reset/;
         $result =~ m/>.*?Expect\=(\S+)/;
         $evalue = $1;
         unless (defined $evalue) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse for evalue in entry \n$entry \nThe entry will be ignored and analysis continues",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error processing input files in Metagenome Analysis".
                           "\nCould not parse for evalue in entry \n$entry \n\n";
            close ERRORLOG;
            next;
         }
      }

      #create hash
      $results{$counter} = {(evalue   => $evalue,
                             blast_db => $blast_db,
                             header   => $query,
                             result   => $result
                            )
                           };

      #skip if above threshold
      if ($evalue > $blast_cutoff) {
         $counter++;
      } else {
         #create sortable list for summary file
         push(@sort_list, $evalue.'__'.$counter);
         $counter++;
      }
   }
   undef @res_list;

   #sort list by evalue
   &update_pbar_1(title        => "Filtering Blast results",
                  label        => "Sorting results by e-value",
                  progress     => ($counter / $max_count) * 100,
                 );
   @sort_list =
      map  $_->[0] =>
      sort { $a->[1] <=> $b->[1] }
      map  [ $_, m/^(.+?)\__/ ]
      => @sort_list;

   #create file of sorted entries
   $counter   = 1;
   $max_count = $#sort_list + 2;
   open WRITE, "+>${$args{ini_ref}}{blast_results}".'/sorted_results';
   print WRITE "Entry number\tAccess key\tSequence frequency\tBest E-value\tUsed Database\tSequence\tBlast overview\n";

   foreach my $entry (@sort_list) {
      my ($header_count, @blast_limit, $result);
      if ($counter % 100 == 0) {
         &update_pbar_1(title        => "Filtering Blast results",
                        label        => "Writing results $counter",
                        progress     => ($counter / $max_count) * 100,
                       );
      }

      my ($evalue, $ID);
      'reset' =~ m/reset/;
      $entry =~ m/(.+?)__(.+)/;
      ($evalue, $ID) = ($1, $2);

      #process entries correctly
      my $print_query = $results{$ID}->{header};
      $header_count = 0;
      $print_query =~ s/___/\t/g;
      $header_count = ($print_query =~ tr/\t//);
      $header_count++;
      $print_query =~ s/\t/ /g;

      #limit Blast results
      @blast_limit = split/>/, $results{$ID}->{'result'};
      @blast_limit = splice(@blast_limit,0,($blast_display_limit + 1));
      $result = join('>', @blast_limit);
      undef @blast_limit;

      #modify blast results for excel
      #$results{$ID}->{'result'} =~ s/\n>/\n\t\t\t\t\t\t>/gs;
      $result =~ s/\n>/\n\t\t\t\t\t\t>/gs;

      #find correct header entry to catch sequence
      'reset' =~ m/reset/;
      my $single_header;
      my @temp = split /___/,$results{$ID}->{header};
      $single_header = $temp[0];
      undef @temp;

      unless (defined $single_header) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse for single header in entry \n". $results{$ID}->{header} ."\n",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files in Metagenome Analysis".
                        "\nCould not parse for single header in entry \n". $results{$ID}->{header} ." \n\n";
         close ERRORLOG;
         &hide_pbar_1;
         return (0);
      }

      unless (defined ${$args{candidate_lookup}}{$single_header}) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not retrieve sequence from lookup hash in entry \n$ID -> ". $results{$ID}->{header} ."\n",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files in Metagenome Analysis".
                        "\nCould not retrieve sequence from lookup hash in entry \n$ID -> ". $results{$ID}->{header} ."\n\n";
         close ERRORLOG;
         &hide_pbar_1;
         return (0);
      }

      #write results to file
      print WRITE $counter."\t".$print_query."\t".$header_count."\t".$evalue."\t".$results{$ID}->{'blast_db'}."\t".${$args{candidate_lookup}}{$single_header}."\t".$result."\n";

      $counter++;
   }
   close WRITE;

   &update_pbar_1(title        => "Filtering Blast results",
                  label        => "Binning results",
                  progress     => 1
                 );
   &rank_distribution(main_window     => $args{main_window},
                      progress_bar    => $args{progress_bar},
                      auto_ini_ref    => $args{auto_ini_ref},
                      ini_ref         => $args{ini_ref},
                      meta_res_ref    => \%results,
                      total_seq_count => $args{total_seq_count}
                     );
   &hide_pbar_1;
   return();
}

sub rank_distribution {
   my %args = @_;
   my ($header_count, %distribution, %distribution_1, %frequency, %frequency_1, @binning, @binning_1);

   @binning   = qw(1 1e-10 1e-20 1e-30 1e-40 1e-50 1e-60 1e-70 1e-80 1e-90 1e-100 1e-110 1e-120 1e-130 1e-140 1e-150 1e-160);
   @binning_1 = qw(1 1e-10 1e-40 1e-100);

   #iterate through
   foreach my $key (keys %{ $args{meta_res_ref} }) {
      if (${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-10) {
         $distribution{'1'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-10  && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-20)  {
         $distribution{'1e-10'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-10'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-20  && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-30)  {
         $distribution{'1e-20'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-20'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-30  && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-40)  {
         $distribution{'1e-30'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-30'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-40  && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-50)  {
         $distribution{'1e-40'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-40'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-50  && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-60)  {
         $distribution{'1e-50'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-50'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-60  && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-70)  {
         $distribution{'1e-60'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-60'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-70  && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-80)  {
         $distribution{'1e-70'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-70'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-80  && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-90)  {
         $distribution{'1e-80'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-80'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-90  && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-100) {
         $distribution{'1e-90'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-90'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-100 && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-110) {
         $distribution{'1e-100'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-100'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-110 && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-120) {
         $distribution{'1e-110'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-110'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-120 && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-130) {
         $distribution{'1e-120'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-120'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-130 && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-140) {
         $distribution{'1e-130'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-130'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-140 && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-150) {
         $distribution{'1e-140'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-140'} = $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-150 && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-160) {
         $distribution{'1e-150'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-150'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-160 && ${$args{meta_res_ref}}{$key}->{'evalue'} >= 0    ) {
         $distribution{'1e-160'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency{'1e-160'} += $header_count;
      }

      if (${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-10) {
         $distribution_1{'1'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency_1{'1'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-10  && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-40)  {
         $distribution_1{'1e-10'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency_1{'1e-10'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-40  && ${$args{meta_res_ref}}{$key}->{'evalue'} > 1e-100) {
         $distribution_1{'1e-40'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency_1{'1e-40'} += $header_count;
      } elsif (${$args{meta_res_ref}}{$key}->{'evalue'} <= 1e-100  && ${$args{meta_res_ref}}{$key}->{'evalue'} >= 0   ) {
         $distribution_1{'1e-100'}++;
         my $print_query = ${$args{meta_res_ref}}{$key}->{'header'};
         $header_count = 0;
         $print_query =~ s/___/\t/g;
         $header_count = ($print_query =~ tr/\t//);
         $header_count++;
         $frequency_1{'1e-100'} += $header_count;
      }
   }

   #write binning
   open WRITE, "+>${$args{ini_ref}}{blast_results}".'/detailed_binning.txt';
   print WRITE "e-value range\tDistribution\tUncollapsed frequency\n";
   foreach my $entry (@binning) {
      if (defined $distribution{$entry}) {
         print WRITE $entry."\t".$distribution{$entry}."\t".$frequency{$entry}."\n";
      } else {
         print WRITE $entry."\t0\t0\n";
      }
   }
   print WRITE "\n\nTotal number of sequence reads\t$args{total_seq_count}\n";
   close WRITE;

   open WRITE, "+>${$args{ini_ref}}{blast_results}".'/binning.txt';
   print WRITE "e-value range\tDistribution\tUncollapsed frequency\n";
   foreach my $entry (@binning_1) {
      if (defined $distribution_1{$entry}) {
         print WRITE $entry."\t".$distribution_1{$entry}."\t".$frequency{$entry}."\n";
      } else {
         print WRITE $entry."\t0\t0\n";
      }
   }
   print WRITE "\n\nTotal number of sequence reads\t$args{total_seq_count}\n";
   close WRITE;

   return;
}

sub populate_new {
   my %args = @_;
   my $new_selected_files = '';
   @new_short_files       = ();

   ($new_selected_files) = ${$args{main_window}}->getOpenFile(-initialdir => $current_dir,
                                                              -title      => 'Select files for custom database',
                                                              -multiple   => 1);

   if (defined(@{$new_selected_files})) {
      #generate short_list
      foreach my $file (@{$new_selected_files}) {
         'reset' =~ m/reset/;
         $file =~ m/(.+)\/([^\/]+)$/; #display only file_name
         $current_dir = $1;
         my $short = $2;
         push (@new_short_files, $short);
         #configure new width for textbox
         if (length($short) > $new_width) {
            $new_width = length($short);
         }
      }

      #remove duplicates
      &delete_duplicates(array => \@new_short_files);

      #populate textbox with short files
      $args{textbox}->delete(0, 'end');
      foreach my $entry (@new_short_files) {
         $args{textbox}->insert('end', $entry);
      }

      #merge selected files with overall selection array
      foreach (@{$new_selected_files}) {
         'reset' =~ m/reset/;
         m/\/([^\/]+)$/; #display only file_name
         my $short = $1;
         $new_all_files{$short} = $_;
      }

      #update display
      $args{textbox}->configure(-width => $new_width);
      $tl->update;
   }
}

sub remove_entry_new {
   my %args = @_;
   return if ($#new_short_files < 0);
   my ($current_index) = ${$args{textbox}}->curselection();
   return unless ($current_index =~ /\S+/);
   #remove from shortlist array
   splice(@new_short_files, $current_index, 1);
   #remove from all files hash
   my $sel_key = ${$args{textbox}}->get($current_index);
   foreach my $key (keys %new_all_files) {
      if ($key =~ m/$sel_key$/) {
         delete $new_all_files{$key};
         ${$args{progress_bar}}->Label(-text => "Deleted entry $key")-> grid (-row => 0, -column => 0, -sticky => 'esnw');
         ${$args{main_window}}->update;
         last;
      }
   }
   #remove from display
   ${$args{textbox}}->delete($current_index);
   ${$args{main_window}}->update;

}

sub delete_duplicates {
   my %args = @_;
   my %seen = ();
   @{$args{array}} = grep { ! $seen{$_} ++ } @{$args{array}};
}

1;




