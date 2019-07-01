#!/opt/ActivePerl-5.8/bin/perl

#custom_db: generate custom nt and aa databases from fasta, msfasta and Genbank files
#input arguments: main_window, ini_ref, auto_ini_ref,


package Miscellaneous::custom_db;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&create_custom_db &create_custom_COG_db);
use vars qw();

use Tk;
use Basics::DirSelect;
use initialise::read_me              qw(:DEFAULT);
use initialise::recompile            qw(:DEFAULT);
use ProgrammeModules::sequence       qw(:DEFAULT);
use ProgrammeModules::genbank_parser qw(:DEFAULT);
use Basics::progress_bar             qw(:DEFAULT);
use Cwd;
use File::Copy;


#local variables
my (%args, $status, $tl, @selected_files, %all_files, @short_files, $width, $aa_db, $nt_db,
    $database_name, $keep_source, $keep_multi, $target_dir, $keep_orphans, $min_number_of_orthologs);

#initialise
@selected_files          = ();
@short_files             = ();
%all_files               = ();
$min_number_of_orthologs = 0;
$keep_orphans            = 1;

sub create_custom_db {
   my %args = @_;
   $width = 20;
   my ($file_lb, $aa, $nt, $frame_left, $frame_right, $progress_bar, $adjuster_frame);
   $target_dir = ${$args{ini_ref}}{blast_db_path};

   #$tl = ${$args{main_window}}->Toplevel(-title => 'Create custom databases');
   if (defined $tl && Tk::Exists($tl)) {
      $tl->state('normal');
      $progress_bar = $tl->Frame(-borderwidth => 2, -relief => 'groove');
      $frame_left   = $tl->Frame()                                      ;
      $frame_right  = $tl->Frame(-relief => 'groove')                   ;
      $tl->update;
      ${$args{main_window}}->update;
   } else {
      $tl = ${$args{main_window}}->DialogBox(-title => 'Create custom Blast databases',
                                             -buttons => [ 'Exit'  ]
                                            );
      $tl->minsize(20, 10);
      #$tl->maxsize(120, 20);
      $progress_bar = $tl->Frame(-borderwidth => 2, -relief => 'groove') -> pack(-side => 'bottom', -fill => 'x');
      $frame_left   = $tl->Frame()                                       -> pack(-side => 'left',   -fill => 'y');
      $frame_right  = $tl->Frame(-relief => 'groove')                    -> pack(-side => 'left',   -fill => 'both', -expand => 1);
   }
   #initialise
   @selected_files = ();
   @short_files = ();
   %all_files = ();

   $progress_bar->configure(-label => " "); #create dummy for space
   #buttons into left frame
   $frame_left                 ->Label(-text => 'Options for customised databases')->grid(-row => 0, -column => 0, -sticky => 'nsew'); #dummy
   $aa = $frame_left           ->Checkbutton(-text     => 'Create aa database',
                                             -variable => \$aa_db,
                                             -command  => sub {
                                                                if ($aa_db == 1) {
                                                                   $nt_db = 0;
                                                                } elsif ($aa_db == 0) {
                                                                   $nt_db = 1;
                                                                   $args{process_gb} = &extract_genbank_seq(main_window  => \$tl,
                                                                                                            auto_ini_ref => $args{auto_ini_ref});
                                                                }
                                                              }
                                    )->grid(-row => 1, -column => 0, -sticky => 'w');
   $nt = $frame_left           ->Checkbutton(-text     => 'Create nt database',
                                             -variable => \$nt_db,
                                             -command  => sub {
                                                                if ($nt_db == 1) {
                                                                   $aa_db = 0;
                                                                   $args{process_gb} = &extract_genbank_seq(main_window  => \$tl,
                                                                                                            auto_ini_ref => $args{auto_ini_ref});
                                                                } elsif ($nt_db == 0) {
                                                                   $aa_db = 1;
                                                                }
                                                              }
                                    )->grid(-row => 2, -column => 0, -sticky => 'w');
   $frame_left                 ->Label(-text => ' ')->grid(-row => 3, -column => 0); #dummy
   my $keep_fasta = $frame_left ->Checkbutton(-text     => 'Keep Source FASTA files',
                                              -variable => \$keep_source
                                       )->grid(-row => 4, -column => 0, -sticky => 'w');
   my $keep_files = $frame_left ->Checkbutton(-text     => 'Keep individual files',
                                              -variable => \$keep_multi
                                       )->grid(-row => 5, -column => 0, -sticky => 'w');
   $frame_left                 ->Label(-text => ' ')->grid(-row => 6, -column => 0); #dummy
   my $name_bt = $frame_left   ->LabEntry(-label        => 'Enter database name',
                                          -labelPack    => [-side => "left", -anchor => "w"],
                                          -justify      => 'left',
                                          -width        => 20,
                                          -textvariable => \$database_name
                                      )->grid(-row => 7, -column => 0, -sticky => 'w');
   my $dir_bt = $frame_left    ->LabEntry(-label        => 'Enter target directory',
                                          -labelPack    => [-side => "left", -anchor => "w"],
                                          -justify      => 'left',
                                          -width        => 20,
                                          -textvariable => \$target_dir
                                      )->grid(-row => 8, -column => 0, -sticky => 'w');
                 $frame_left   ->Button  (-bitmap       => '@'.${$args{auto_ini_ref}}{work_dir}.'/lib/initialise/tree.xbm',
                                          -command      => sub {
                                                                 my $ds = $tl->DirSelect(-title => 'Select location of Blast databases');
                                                                 my $dir = $ds->Show($target_dir, -popover => $tl);
                                                                 if (defined $dir && $dir =~ /\w+/) {
                                                                    $target_dir = $dir;
                                                                    $tl->update;
                                                                 }
                                                             }
                                         )->grid(-row => 8, -column => 1, -sticky => 'w');
   my $compile   = $frame_left ->Button(-text => 'Compile database',
                                        -command => sub { unless (defined $database_name && $database_name =~ /\w+/) {
                                                             $tl->messageBox(-title   => 'Error',
                                                                             -message => 'No database name entered',
                                                                             -type    => 'OK',
                                                                             -icon    => 'error');
                                                          }

                                                          unless ($#short_files >= 0) {
                                                             $tl->messageBox(-title   => 'Error',
                                                                             -message => 'No files selected',
                                                                             -type    => 'OK',
                                                                             -icon    => 'error');
                                                          }
                                                          if (defined ($nt_db) && $nt_db == 1 && $#short_files >= 0 &&
                                                              defined $database_name && $database_name =~ /\w+/) {
                                                             &parse_files(main_window   => \$tl,
                                                                          progress_bar  => \$progress_bar,
                                                                          ini_ref       => $args{ini_ref},
                                                                          auto_ini_ref  => $args{auto_ini_ref},
                                                                          file_list_ref => \%all_files,
                                                                          process_gb    => $args{process_gb},
                                                                          type          => 'nt',
                                                                          flavour       => 'blast'
                                                                         );
                                                             #return (\%all_files, 'nt', \$tl, \$progress_bar);
                                                          } elsif (defined ($aa_db) && $aa_db == 1 && $#short_files >= 0 &&
                                                                   defined $database_name && $database_name =~ /\w+/) {
                                                             &parse_files(main_window   => \$tl,
                                                                          progress_bar  => \$progress_bar,
                                                                          ini_ref       => $args{ini_ref},
                                                                          auto_ini_ref  => $args{auto_ini_ref},
                                                                          file_list_ref => \%all_files,
                                                                          process_gb    => $args{process_gb},
                                                                          type          => 'aa',
                                                                          flavour       => 'blast'
                                                                         );
                                                             #return (\%all_files, 'aa', \$tl, \$progress_bar);
                                                          } elsif ($#short_files > 0 && defined $database_name && $database_name =~ /\w+/) {
                                                             $tl->messageBox(-title   => 'Error',
                                                                             -message => 'Select database type',
                                                                             -type    => 'OK',
                                                                             -icon    => 'error');
                                                          }
                                                        }
                                       )->grid(-row => 9, -column => 0, -sticky => 'sew');

   $frame_left->gridRowconfigure    (0, -pad    => 30);
   $frame_left->gridColumnconfigure (0, -weight => 1);
   $frame_left->gridRowconfigure    (8, -weight => 1);

   #Textbox and button in right frame
   my $select_bt = $frame_right ->Button(-text => 'Add files',
                                        -command => sub {&populate(main_window  => \$tl,
                                                                   auto_ini_ref => $args{auto_ini_ref},
                                                                   textbox      => $file_lb)}
                                       )->grid(-row => 0, -column => 0, -sticky => 'w');

   $file_lb         = $frame_right    ->Scrolled("Listbox",
                                                 -setgrid    => 1,
                                                 -selectmode => 'single',
                                                 -width      => $width
                                                 )->grid(-row => 1, -column => 0, -columnspan => 2, -sticky => 'nsew');
   $file_lb                           ->configure(-width   => $width,
                                                  -setgrid => 1 );

   my $clear_selection = $frame_right ->Button(-text => 'Clear current selection',
                                               -command => sub {
                                                                 %all_files = ();
                                                                 @short_files = ();
                                                                 $file_lb->delete(0, 'end');
                                                                 $tl->update;
                                                                }
                                              )->grid(-row => 2, -column => 0, -sticky => 'w');
   $frame_right                       ->Label(-text => 'Double click entry to remove.') ->grid(-row => 2, -column => 1, -sticky => 'w');

   $frame_right->gridRowconfigure    (1, -weight => 1);
   $frame_right->gridColumnconfigure (0, -weight => 1);

   #bind right mouse button to selection and removal of entry
   $file_lb->bind('<Double-ButtonRelease-1>' => sub {&remove_entry(textbox      => \$file_lb,
                                                                   main_window  => \$tl,
                                                                   progress_bar => \$progress_bar,
                                                                   ini_ref      => $args{ini_ref},
                                                                   auto_ini_ref => $args{auto_ini_ref}
                                                                  );

                                                    });

   my $wait = $tl->Show();
   $file_lb->delete(0, 'end');
   $tl->update;
   undef @selected_files;
   undef %all_files;
   undef @short_files;
   if ($wait eq 'Exit') {
      $tl->state('withdrawn');
   }

   if (defined ($database_name) && $database_name =~ /\w+/) {
      return ($target_dir.'/'.$database_name);
   } else {
      return;
   }
}

sub create_custom_COG_db {
   my %args = @_;
   $width = 20;
   my ($file_lb, $frame_left, $frame_right, $progress_bar, $adjuster_frame, $orphans, $min_ortho);
   $target_dir = ${$args{ini_ref}}{COG_db_path};

   #$tl = ${$args{main_window}}->Toplevel(-title => 'Create custom databases');
   if (defined $tl && Tk::Exists($tl)) {
      $tl->state('normal');
      $progress_bar = $tl->Frame(-borderwidth => 2, -relief => 'groove');
      $frame_left   = $tl->Frame()                                      ;
      $frame_right  = $tl->Frame(-relief => 'groove')                   ;
   } else {
      $tl = ${$args{main_window}}->DialogBox(-title => 'Create custom COG databases',
                                             -buttons => [ 'Exit'  ]
                                            );
      $tl->minsize(20, 10);
      #$tl->maxsize(120, 20);
      $progress_bar = $tl->Frame(-borderwidth => 2, -relief => 'groove') -> pack(-side => 'bottom', -fill => 'x');
      $frame_left   = $tl->Frame()                                       -> pack(-side => 'left',   -fill => 'y');
      $frame_right  = $tl->Frame(-relief => 'groove')                    -> pack(-side => 'left',   -fill => 'both', -expand => 1);
   }
   #initialise
   @selected_files = ();
   @short_files = ();
   %all_files = ();

   $progress_bar->configure(-label => " "); #create dummy for space
   #buttons into left frame
   $frame_left                 ->Label(-text => 'Options for customised databases')->grid(-row => 0, -column => 0, -sticky => 'nsew'); #dummy

   my $keep_fasta = $frame_left ->Checkbutton(-text     => 'Keep Source FASTA files',
                                              -variable => \$keep_source
                                       )->grid(-row => 1, -column => 0, -sticky => 'w');
   my $keep_files = $frame_left ->Checkbutton(-text     => 'Keep individual files',
                                              -variable => \$keep_multi
                                       )->grid(-row => 2, -column => 0, -sticky => 'w');
   $orphans    = $frame_left ->Checkbutton(-text     => 'Allow non-clustered ORFs in database',
                                           -variable => \$keep_orphans,
                                           -command => sub {
                                                             if ($keep_orphans == 0) {
                                                                $min_number_of_orthologs = 0;
                                                                $min_ortho->configure(-state      => 'disabled',
                                                                                      -foreground => 'gray');
                                                             } elsif ($keep_orphans == 1) {
                                                                if ($#short_files >= 0) {
                                                                   $min_number_of_orthologs = $#short_files + 1;
                                                                } else {
                                                                  $min_number_of_orthologs = 1;
                                                                }
                                                                $min_ortho->configure(-state      => 'normal',
                                                                                      -foreground => 'black');
                                                             }
                                                           }
                                       )->grid(-row => 3, -column => 0, -sticky => 'w');
   $min_ortho  = $frame_left ->BrowseEntry(-label   => "Minimum number of orthologs\n\'0\' eq no limit",
                                           -choices => [0..200],
                                           -width   => 4,
                                           -variable=> \$min_number_of_orthologs,
                                           -command => sub {
                                                             if ($min_number_of_orthologs > 0) {
                                                                $keep_orphans = 0;
                                                                $orphans->configure(-state      => 'disabled',
                                                                                    -foreground => 'gray');
                                                             } elsif ($min_number_of_orthologs == 0) {
                                                                $keep_orphans = 1;
                                                                $orphans->configure(-state      => 'normal',
                                                                                    -foreground => 'black');
                                                             }
                                                           }
                                       )->grid(-row => 4, -column => 0, -sticky => 'w');
   $frame_left                 ->Label(-text => ' '
                                       )->grid(-row => 5, -column => 0); #dummy
   my $name_bt = $frame_left   ->LabEntry(-label        => 'Enter database name',
                                          -labelPack    => [-side => "left", -anchor => "w"],
                                          -justify      => 'left',
                                          -width        => 20,
                                          -textvariable => \$database_name
                                      )->grid(-row => 6, -column => 0, -sticky => 'w');
   my $dir_bt = $frame_left    ->LabEntry(-label        => 'Enter target directory',
                                          -labelPack    => [-side => "left", -anchor => "w"],
                                          -justify      => 'left',
                                          -width        => 20,
                                          -textvariable => \$target_dir
                                      )->grid(-row => 7, -column => 0, -sticky => 'w');
                 $frame_left   ->Button  (-bitmap       => '@'.${$args{auto_ini_ref}}{work_dir}.'/lib/initialise/tree.xbm',
                                          -command      => sub {
                                                                 my $ds = $tl->DirSelect(-title => 'Select location of Blast databases');
                                                                 my $dir = $ds->Show($target_dir, -popover => $tl);
                                                                 if (defined $dir && $dir =~ /\w+/) {
                                                                    $target_dir = $dir;
                                                                    $tl->update;
                                                                 }
                                                             }
                                         )->grid(-row => 7, -column => 1, -sticky => 'w');
   my $compile   = $frame_left ->Button(-text => 'Compile COG database',
                                        -command => sub { unless (defined $database_name && $database_name =~ /\w+/) {
                                                             $tl->messageBox(-title   => 'Error',
                                                                             -message => 'No database name entered',
                                                                             -type    => 'OK',
                                                                             -icon    => 'error');
                                                          }

                                                          unless ($#short_files >= 0) {
                                                             $tl->messageBox(-title   => 'Error',
                                                                             -message => 'No files selected',
                                                                             -type    => 'OK',
                                                                             -icon    => 'error');
                                                          }

                                                          #test if COGnitor is compiled
                                                          my $cog_status = &recompile(main_window   => \$tl,
                                                                                      ini_ref       => $args{ini_ref},
                                                                                      auto_ini_ref  => $args{auto_ini_ref},
                                                                                      test_cognitor => 1);
                                                          if ($cog_status == 0) {
                                                             $tl->messageBox(-title   => 'Error',
                                                                             -message => 'Could not re-compile COGNitor, aborting',
                                                                             -type    => 'OK',
                                                                             -icon    => 'error');
                                                          } else {
                                                             #if COGNitor is OK, continue with db build
                                                             if ($#short_files >= 0     &&
                                                                 defined $database_name &&
                                                                 $database_name =~ /\w+/) {
                                                                &parse_files(main_window   => \$tl,
                                                                             progress_bar  => \$progress_bar,
                                                                             ini_ref       => $args{ini_ref},
                                                                             auto_ini_ref  => $args{auto_ini_ref},
                                                                             file_list_ref => \%all_files,
                                                                             process_gb    => 'extract',
                                                                             type          => 'aa',
                                                                             flavour       => 'cog'
                                                                            );
                                                                #return (\%all_files, 'aa', \$tl, \$progress_bar);
                                                             } elsif ($#short_files > 0 && defined $database_name && $database_name =~ /\w+/) {
                                                                $tl->messageBox(-title   => 'Error',
                                                                                -message => 'Select database type',
                                                                                -type    => 'OK',
                                                                                -icon    => 'error');
                                                             }
                                                          }
                                                        }
                                       )->grid(-row => 9, -column => 0, -sticky => 'sew');

   $frame_left->gridRowconfigure    (0, -pad    => 30);
   $frame_left->gridColumnconfigure (0, -weight => 1);
   $frame_left->gridRowconfigure    (8, -weight => 1);

   #initial state
   if ($keep_orphans == 0) {
      $min_number_of_orthologs = 0;
      $min_ortho->configure(-state      => 'disabled',
                            -foreground => 'gray');
   } elsif ($keep_orphans == 1) {
      if ($#short_files >= 0) {
         $min_number_of_orthologs = $#short_files + 1;
      } else {
         $min_number_of_orthologs = 1;
      }
      $min_ortho->configure(-state      => 'normal',
                            -foreground => 'black');
   }
   if ($min_number_of_orthologs > 0) {
      $keep_orphans = 0;
      $orphans->configure(-state      => 'disabled',
                          -foreground => 'gray');
   } elsif ($min_number_of_orthologs == 0) {
      $keep_orphans = 1;
      $orphans->configure(-state      => 'normal',
                          -foreground => 'black');
   }

   #Textbox and button in right frame
   my $select_bt = $frame_right ->Button(-text => 'Add files',
                                        -command => sub {&populate(main_window  => \$tl,
                                                                   auto_ini_ref => $args{auto_ini_ref},
                                                                   textbox      => $file_lb);
                                                          if ($min_number_of_orthologs > ($#short_files + 1)) {
                                                             $min_number_of_orthologs = $#short_files + 1;
                                                          }
                                                         }
                                       )->grid(-row => 0, -column => 0, -sticky => 'w');

   $file_lb         = $frame_right    ->Scrolled("Listbox",
                                                 -setgrid    => 1,
                                                 -selectmode => 'single',
                                                 -width      => $width
                                                 )->grid(-row => 1, -column => 0, -columnspan => 2, -sticky => 'nsew');
   $file_lb                           ->configure(-width   => $width,
                                                  -setgrid => 1 );

   my $clear_selection = $frame_right ->Button(-text => 'Clear current selection',
                                               -command => sub {
                                                                 %all_files = ();
                                                                 @short_files = ();
                                                                 $file_lb->delete(0, 'end');
                                                                 $min_number_of_orthologs = 0;
                                                                 $tl->update;
                                                                }
                                              )->grid(-row => 2, -column => 0, -sticky => 'w');
   $frame_right                       ->Label(-text => 'Double click entry to remove.') ->grid(-row => 2, -column => 1, -sticky => 'w');

   $frame_right->gridRowconfigure    (1, -weight => 1);
   $frame_right->gridColumnconfigure (0, -weight => 1);

   #bind right mouse button to selection and removal of entry
   $file_lb->bind('<Double-ButtonRelease-1>' => sub {&remove_entry(textbox      => \$file_lb,
                                                                   main_window  => \$tl,
                                                                   progress_bar => \$progress_bar,
                                                                   ini_ref      => $args{ini_ref},
                                                                   auto_ini_ref => $args{auto_ini_ref}
                                                                   )
                                                    });
   my $wait = $tl->Show();
   $file_lb->delete(0, 'end');
   $tl->update;
   undef @selected_files;
   undef %all_files;
   undef @short_files;
   if ($wait eq 'Exit') {
      $tl->state('withdrawn');
   }

   if (defined ($database_name) && $database_name =~ /\w+/) {
      return ($target_dir.'/'.$database_name);
   } else {
      return;
   }
}

sub delete_duplicates {
   my %args = @_;
   my %seen = ();
   @{$args{array}} = grep { ! $seen{$_} ++ } @{$args{array}};
}

sub populate {
   my %args = @_;
   (my $selected_files) = ${$args{main_window}}->getOpenFile(-initialdir => ${$args{auto_ini_ref}}{work_dir},
                                                             -title      => 'Select files for custom database',
                                                             -multiple   => 1);

   if (defined(@{$selected_files})) {
      #generate short_list
      foreach my $file (@{$selected_files}) {
         $file =~ m/\/([^\/]+)$/; #display only file_name
         my $short = $1;
         push (@short_files, $short);
         #configure new width for textbox
         if (length($short) > $width) {
            $width = length($short);
         }
      }
      #remove duplicates
      &delete_duplicates(array => \@short_files);

      #populate textbox with short files
      $args{textbox}->delete(0, 'end');
      foreach my $entry (@short_files) {
         $args{textbox}->insert('end', $entry);
      }

      #merge selected files with overall selection array
      foreach (@{$selected_files}) {
         $all_files{$_} = 1;
      }

      #update display
      #$args{textbox}->configure(-width => $width);
      $tl->update;
   }
}

sub remove_entry {
   my %args = @_;
   return if ($#short_files < 0);
   my ($current_index) = ${$args{textbox}}->curselection();
   #remove from shortlist array
   splice(@short_files, $current_index, 1);

   #remove from all files hash
   my $sel_key = ${$args{textbox}}->get($current_index);
   foreach my $key (keys %all_files) {
      if ($key =~ m/$sel_key$/) {
         delete $all_files{$key};
         ${$args{progress_bar}}->Label(-text => "Deleted entry $key")-> grid (-row => 0, -column => 0, -sticky => 'esnw');
         ${$args{main_window}}->update;
         last;
      }
   }
   #remove from display
   ${$args{textbox}}->delete($current_index);
   ${$args{main_window}}->update;
}

sub parse_files {
   my %args = @_;
   my ($type, $counter, $max_count, @multi_files, @remove_multi);
   my $fasta_only = 1; #assume fast track, only fasta files

   #make things much faster if only msFASTA files are selected
   #test if only (ms)fasta files are present
   {
      foreach my $key (keys %{$args{file_list_ref}}) {
         my (@entries, $file_ref);
         $counter = 0;
         ${$args{progress_bar}} ->configure(-label => "Testing filetype for $key");
         ${$args{main_window}}->update;
         unless (-e $key) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "File $key cannot be accessed",
                                                          -bitmap  => 'error',
                                                          -buttons => ['ok']);
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error compiling custom databases".
                           "\nFile $key cannot be accessed\n\n";
            close ERRORLOG;
         }
         #pull filepath apart
         'reset' =~ m/reset/;
         $key =~ m/^(.+)\/([^\/]+)$/;
         my $local_path = $1;
         my $filename = $2;
         #read file and test first line. Assume consistency of format with a file
         open READ, "<$key";
         while (<READ>) {
            if ($. == 1 && $_ !~ /^\>/i) {
               $fasta_only = 0;
            }
            last;
         }
         close READ;
      }

      #only fasta?
      if ($fasta_only == 1) {
         #merge all fasta files into one
         open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/temp.fasta";
         foreach my $key (keys %{$args{file_list_ref}}) {
            ${$args{progress_bar}}->configure(-label => "Merging fasta file $key");
            ${$args{main_window}} ->update;
            open READ, "<$key";
            while (<READ>) {
               print WRITE $_;
            }
            close READ;
            print WRITE "\n";
         }
         close WRITE;

         #create index for FASTA
         ${$args{progress_bar}} ->configure(-label => "Creating msfasta index");
         ${$args{main_window}}  ->update;
         my $status = &create_index(main_window   => $args{main_window},
                                    progress_bar  => $args{progress_bar},
                                    ini_ref       => $args{ini_ref},
                                    auto_ini_ref  => $args{auto_ini_ref},
                                    flavour       => $args{flavour}
                                   );
         if ($status == 0) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "No sequences to compile.\nAborting",
                                                          -bitmap  => 'error',
                                                          -buttons => ['ok']);
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error compiling custom databases".
                           "\nNo sequences to compile.\nAborting\n\n";
            close ERRORLOG;
            unlink "${$args{auto_ini_ref}}{work_dir}/temp.fasta";
            return(0);
         }

         #compile database
         ${$args{progress_bar}} ->configure(-label => "Compiling database");
         ${$args{main_window}}  ->update;

         if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
            if ($args{type} eq 'aa') {
               `${$args{ini_ref}}{blast_executables}/formatdb -i ${$args{auto_ini_ref}}{work_dir}/temp.fasta -p T -o T -n $target_dir/$database_name`;
            } elsif ($args{type} eq 'nt') {
               `${$args{ini_ref}}{blast_executables}/formatdb -i ${$args{auto_ini_ref}}{work_dir}/temp.fasta -p F -o T -n $target_dir/$database_name`;
            }
         } else {
            if ($args{type} eq 'aa') {
               `${$args{ini_ref}}{blast_plus_executables}/makeblastdb -dbtype prot -in ${$args{auto_ini_ref}}{work_dir}/temp.fasta -parse_seqids -hash_index -out $target_dir/$database_name`;
            } elsif ($args{type} eq 'nt') {
               `${$args{ini_ref}}{blast_plus_executables}/makeblastdb -dbtype nucl -in ${$args{auto_ini_ref}}{work_dir}/temp.fasta -parse_seqids -hash_index -out $target_dir/$database_name`;
            }
         }

         ${$args{progress_bar}} ->configure(-label => "Compiled database $database_name. Task finished.");
         ${$args{main_window}}  ->update;

         #unlink temp files or copy
         if ($keep_source == 1) {
            copy(${$args{auto_ini_ref}}{work_dir}.'/temp.fasta', ${$args{auto_ini_ref}}{work_dir}.'/'.$database_name.'.msfasta');
            ${$args{main_window}}->messageBox(-title   => 'Info',
                                              -message => "Created FASTA source file $database_name\.msfasta in directory ${$args{auto_ini_ref}}{work_dir}\.",
                                              -icon    => 'info',
                                              -type    => 'ok');
         }
         unlink "${$args{auto_ini_ref}}{work_dir}/temp.fasta";

         ${$args{main_window}}  ->messageBox(-title   => 'Success',
                                             -message => "Compiled database $database_name",
                                             -type    => 'OK',
                                             -icon    => 'info');

         return (1);
      }
   }

   #test if any of the selected files contains multiple entries
   foreach my $key (keys %{$args{file_list_ref}}) {
      my (@entries, $file_ref);
      $counter = 0;
      ${$args{progress_bar}} ->configure(-label => "Processing $key");
      ${$args{main_window}}->update;
      unless (-e $key) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "File $key cannot be accessed",
                                                       -bitmap  => 'error',
                                                       -buttons => ['ok']);
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error compiling custom databases".
                        "\nFile $key cannot be accessed\n\n";
         close ERRORLOG;
      }
      #pull filepath apart
      'reset' =~ m/reset/;
      $key =~ m/^(.+)\/([^\/]+)$/;
      my $local_path = $1;
      my $filename = $2;
      #read file to memory
      $file_ref = &slurp(main_window   => $args{main_window},
                         progress_bar  => $args{progress_bar},
                         ini_ref       => $args{ini_ref},
                         auto_ini_ref  => $args{auto_ini_ref},
                         directory     => $local_path,
                         filename      => $filename
                        );

      #determine format
      if (${$file_ref} =~ /^locus/i) {
         $type = 'gb';
      } elsif (${$file_ref} =~ /^\>/) {
         $type = 'fasta';
      } else {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Unrecognised file format for $filename. Skipping filed",
                                                       -bitmap  => 'error',
                                                       -buttons => ['ok']);
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error compiling custom databases".
                        "\nUnrecognised file format for $filename. Skipping file\n\n";
         close ERRORLOG;
         next;
      }

      #multi sequence/entry files?
      if ($type eq 'gb') {
         @entries = split/\n\/\/\n/,${$file_ref};
         #provide visual feedback if number of files to process is very large
         $max_count = $#entries;
         if ($max_count > 10000) {
            &progress_bar_2(main_window  => $args{main_window},
                            progress_bar => $args{progress_bar},
                            auto_ini_ref => $args{auto_ini_ref},
                            ini_ref      => $args{ini_ref},
                            title        => 'Parsing sequences for analysis',
                            label        => 'Extracting individual sequence files'
                           );
            &show_pbar_2;
         }

         if ($#entries > 0) {
            my $progress = 0;
            foreach my $entry (@entries) {
               $progress++;
               if ($progress % 1000 == 0) {
                  ${$args{progress_bar}} ->configure(-label => "Processing");
                  &update_pbar_2(title        => 'Extracting individual sequences for analysis',
                                 label        => "Extracting, $progress of $max_count",
                                 progress     => ($progress / $max_count) * 100,
                                );
                  ${$args{main_window}}->update;
               }
               next unless ($entry =~ /\S+/);
               #get new filename from definition, organism or accession, else counter
               my $definition;
               'reset' =~ m/reset/;
               $entry =~ m/\nDEFINITION\s+(.*?)(,|:|;|\n)/;
               $definition = $1;
               unless (defined $definition && $definition =~ /\w+/) {
                  'reset' =~ m/reset/;
                  $entry =~ m/\n\s+ORGANISM\s+(.*?)(,|:|;|\n)/;
                  $definition = $1;
               }
               unless (defined $definition && $definition =~ /\w+/) {
                  'reset' =~ m/reset/;
                  $entry =~ m/\n\ACCESSION\s+(.*?)(,|:|;|\n)/;
                  $definition = $1;
               }

               if (defined $definition) {
                  #definition should not be longer than 100 char
                  $definition = substr ($definition, 0, 100);
               }

               unless (defined $definition && $definition =~ /\w+/) {
                  $counter++;
                  #definition should not be longer than 100 char
                  my $tmp = substr ($filename, 0, 100);
                  $definition = $tmp.'__'.$counter;
               }
               #clean up filename
               $definition =~ s/[^a-zA-Z0-9\.\_]/_/g;

               #test if file aready exists; if so, increase counter
               while (-e $local_path.'/'.$definition) {
                  $counter++;
                  $definition =~ s/\__\d+$//;
                  $definition .= '__'.$counter;
               }

               #clean up entry
               $entry =~ s/^\s+//;

               #write new file
               open WRITE, "+>$local_path\/$definition" or do {
                  my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                                -text    => "Could not create file $local_path\/$definition. Aborting",
                                                                -bitmap  => 'error',
                                                                -buttons => ['ok']);
                  $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
                  $error_msg-> Show();
                  open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
                  print ERRORLOG "Error compiling custom databases".
                                 "\nCould not create file $local_path\/$definition. Aborting\n\n";
                  close ERRORLOG;
                  &hide_pbar_2;
                  return (0);
               };
               print WRITE $entry.'//';
               close WRITE;
               #add entry to array; will be added later to file hash and be used to delete created files
               push (@multi_files, $local_path.'/'.$definition);
               #add key to array, will be used to remove original multi file from hash, will be added later again
               push (@remove_multi, $key);
            }
         }
         #close progress bar
         if ($max_count > 10000) {
            &hide_pbar_2;
         }

      } else {
         ${$file_ref} = "\n".${$file_ref};
         @entries = split/\n\>/,${$file_ref};
         #provide visual feedback if number of files to process is very large
         $max_count = $#entries;
         if ($max_count > 10000) {
            &progress_bar_2(main_window  => $args{main_window},
                            progress_bar => $args{progress_bar},
                            auto_ini_ref => $args{auto_ini_ref},
                            ini_ref      => $args{ini_ref},
                            title        => 'Parsing sequences for analysis',
                            label        => 'Extracting individual sequence files'
                           );
            &show_pbar_2;
         }
         if ($#entries > 0) {
            my $progress = 0;
            foreach my $entry (@entries) {
               $progress++;
               if ($progress % 1000 == 0) {
                  ${$args{progress_bar}} ->configure(-label => "Processing");
                  &update_pbar_2(title        => 'Extracting individual sequences for analysis',
                                 label        => "Extracting, $progress of $max_count",
                                 progress     => ($progress / $max_count) * 100,
                                );
                  ${$args{main_window}}->update;
               }
               next unless ($entry =~ /\S+/);
               my $definition;
               'reset' =~ m/reset/;
               $entry =~ m/^s*(.*?)(,|:|;|>|\n)/;
               $definition = $1;

               if (defined $definition) {
                  #definition should not be longer than 100 char
                  $definition = substr ($definition, 0, 100);
               }

               unless (defined $definition && $definition =~ /\w+/) {
                  $counter++;
                  #definition should not be longer than 100 char
                  my $temp = substr ($filename, 0, 100);
                  $definition = $temp.'__'.$counter;
               }
               #clean up filename
               $definition =~ s/[^a-zA-Z0-9\.\_]+/_/g;
               $definition =~ s/_$//;

               #restrict definition to 250 characters
               $definition = substr($definition, 0, 250);

               #test if file aready exists; if so, increase counter
               while (-e $local_path.'/'.$definition) {
                  $counter++;
                  $definition =~ s/\__\d+$//;
                  $definition .= '__'.$counter;
               }
               #clean up entry
               $entry =~ s/^\s+//;

               #write new file
               open WRITE, "+>$local_path\/$definition" or do {
                  my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                                -text    => "Could not create file $local_path\/$definition. Aborting",
                                                                -bitmap  => 'error',
                                                                -buttons => ['ok']);
                  $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
                  $error_msg-> Show();
                  open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
                  print ERRORLOG "Error compiling custom databases".
                                 "\nCould not create file $local_path\/$definition. Aborting\n\n";
                  close ERRORLOG;
                  #delete  temp files
                  foreach my $entry (@multi_files) {unlink $entry};
                  return (0);
               };
               print WRITE '>'.$entry;
               close WRITE;
               #add entry to array; will be added later to file hash and be used to delete created files
               push (@multi_files, $local_path.'/'.$definition);
               #add key to array, will be used to remove original multi file from hash, will be added later again
               push (@remove_multi, $key);
            }
         }
         #close progress bar
         if ($max_count > 10000) {
            &hide_pbar_2;
         }
      }
      undef $file_ref;
   }

   #add new entries if multi files
   foreach my $entry (@multi_files) {
      if (exists ${$args{file_list_ref}}{$entry}) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "File duplication event in multi-entry files, skipping",
                                                       -bitmap  => 'error',
                                                       -buttons => ['ok']);
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error compiling custom databases".
                        "\nFile duplication event in multi-entry files, skipping\n\n";
         close ERRORLOG;
         next;
      }
      ${$args{file_list_ref}}{$entry} = 1;
   }
   #remove original multi files from hash
   foreach my $entry (@remove_multi) {
      delete ${$args{file_list_ref}}{$entry};
   }

   #temp file for formatdb
   unlink "${$args{auto_ini_ref}}{work_dir}/temp.fasta";
   open WRITE, ">>${$args{auto_ini_ref}}{work_dir}/temp.fasta";

   #provide visual feedback if number of files to process is very large
   $max_count = scalar keys %{$args{file_list_ref}};
   if ($max_count > 10000) {
      &progress_bar_2(main_window  => $args{main_window},
                      progress_bar => $args{progress_bar},
                      auto_ini_ref => $args{auto_ini_ref},
                      ini_ref      => $args{ini_ref},
                      title        => 'Parsing sequences for analysis',
                      label        => 'Parsing sequence files'
                     );
      &show_pbar_2;
   }

   #iterate over selected files
   my $key_count = 0;
   foreach my $key (keys %{$args{file_list_ref}}) {
      my (@entries, $file_ref);

      #restrict length in progress bar to 60 characters
      my $display = substr($key, 0 , 60);
      $key_count++;
      if ($key_count % 1000 == 0) {
         ${$args{progress_bar}} ->configure(-label => "Processing $display");
         &update_pbar_2(title        => 'Parsing sequences for  analysis',
                        label        => "Parsing, $key_count of $max_count",
                        progress     => ($key_count / $max_count) * 100,
                       );
         ${$args{main_window}}->update;
      }

      unless (-e $key) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "File $key cannot be accessed",
                                                       -bitmap  => 'error',
                                                       -buttons => ['ok']);
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error compiling custom databases".
                        "\nFile $key cannot be accessed\n\n";
         close ERRORLOG;
         next;
      }
      #pull filepath apart
      'reset' =~ m/reset/;
      $key =~ m/^(.+)\/([^\/]+)$/;
      my $local_path = $1;
      my $filename = $2;

      open READ, $key or do {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not open file $key. Skipping file",
                                                       -bitmap  => 'error',
                                                       -buttons => ['ok']);
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error compiling custom databases".
                        "\nCould not open file $key. Skipping file\n\n";
         close ERRORLOG;
         next;
      };
      $file_ref = <READ>;
      #determine format
      if ($file_ref =~ /^locus/i) {
         $type = 'gb';
      } elsif ($file_ref =~ /^\>/) {
         $type = 'fasta';
      } else {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Unrecognised file format for $filename. Skipping file",
                                                       -bitmap  => 'error',
                                                       -buttons => ['ok']);
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error compiling custom databases".
                        "\nUnrecognised file format for $filename. Skipping file\n\n";
         close ERRORLOG;
         next;
      }
      close READ;
      undef $file_ref;

      #parse files
      if ($type eq 'gb') {
         my ($header, $source, $core, $DNAseq, $start_ref, $feature_list_ref, $genbank_ref, $counter);
         #parse gb file
         ($header, $source, $core, $DNAseq, $start_ref,
         $feature_list_ref, $genbank_ref, $counter) = &gb_parser(main_window   => $args{main_window},
                                                                 progress_bar  => $args{progress_bar},
                                                                 ini_ref       => $args{ini_ref},
                                                                 auto_ini_ref  => $args{auto_ini_ref},
                                                                 gb_file       => $key,
                                                                 process_gb    => $args{process_gb}
                                                                );
         if ($header eq '0') {next}; #ignore faulty input files
         if ($args{type} eq 'nt') {
            if ($args{process_gb} eq 'extract') {
               my $status = &extract_nt(main_window      => $args{main_window},
                                        DNAseq           => $DNAseq,
                                        genbank_ref      => $genbank_ref,
                                        feature_list_ref => $feature_list_ref,
                                        filename         => $filename,
                                       );
               if ($status eq '0') {next};
            } elsif ($args{process_gb} eq 'sequence') {
               #write into temp file
               print WRITE ">whole_sequence__filename__$key\__\n$DNAseq\n";
            }
         } elsif ($args{type} eq 'aa') {
            my $status = &extract_aa(main_window      => $args{main_window},
                                     progress_bar     => $args{progress_bar},
                                     ini_ref          => $args{ini_ref},
                                     auto_ini_ref     => $args{auto_ini_ref},
                                     DNAseq           => $DNAseq,
                                     genbank_ref      => $genbank_ref,
                                     feature_list_ref => $feature_list_ref,
                                     filename         => $filename,
                                    );
            next if ($status eq '0');
         }
      } elsif ($type eq 'fasta') {
         my @msfasta = ();
         #read file to memory
         my $file_ref = &slurp(main_window  => $args{main_window},
                               directory    => $local_path,
                               filename     => $filename,
                              );
         if ($file_ref eq '0') {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Cannot read file $filename for content. Skipping entry",
                                                          -bitmap  => 'error',
                                                          -buttons => ['ok']);
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error compiling COG databases".
                           "\nCannot read file $filename for content. Skipping entry\n\n";
            close ERRORLOG;
            next;
         }
         #split at '\n>'
         ${$file_ref} = "\n".${$file_ref};
         @msfasta = split/\n\r?\>/,${$file_ref};
         foreach my $entry (@msfasta) {
            if ($entry !~ /\w+/) {next};
            $entry =~ s/^\>//g;
            #test for correct sequence type
            my $test = "";
            my ($type, $seq, $header);
            $entry =~ m/^([^\n]+)\n(.+)/s;
            ($header, $seq) = ($1, $2);
            $header =~ s/\W//gs;
            $header =~ s/\s/_/g;
            $seq    =~ s/\s+//gs;

            if ($seq =~ /[qeilfpzjox]/is) {
               $type = 'aa';
            } else {
               $type = &determine_seq_type(seq=>$seq);
               if ($type eq '0') {
                  my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                                -text    => "No sequence in entry $entry. Skipping entry",
                                                                -bitmap  => 'error',
                                                                -buttons => ['ok']);
                  $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
                  $error_msg-> Show();
                  open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
                  print ERRORLOG "Error compiling custom databases".
                                 "\nNo sequence in entry $entry. Skipping entry\n\n";
                  close ERRORLOG;
                  next;
               }
            }

            if (($type eq 'nt' && $args{type} eq 'aa') || ($type eq 'aa' && $args{type} eq 'nt')) {
               my $conflict = ${$args{main_window}}->Dialog(-title          => 'Sequence type conflict',
                                                            -text           => "Sequence with header \>\>$header\<\< and starting with ".substr($seq,0,30)." is likely of wrong sequence type. Use anyway?",
                                                            -buttons        => ['Yes', 'No'],
                                                            -default_button => 'No'
                                                           )->Show();
               if ($conflict eq 'Yes') {
                  print WRITE '>'.$entry."\n";
                  next;
               } else {
                  next;
               }
            }
            print WRITE '>'.$entry."\n";
         }
         @msfasta = ();
      }
   }
   close WRITE;

   #close progress bar
   if ($max_count > 10000) {
      &hide_pbar_2;
   }

   #remove multi files
   foreach my $entry (@multi_files) {
      if ($keep_multi == 0) {
         unlink $entry;
      }
      delete ${$args{file_list_ref}}{$entry};
   }

   #add original multi files to hash again
   foreach my $entry (@remove_multi) {
     ${$args{file_list_ref}}{$entry} = 1;
   }

   #create indexed fasta file
   ${$args{progress_bar}} ->configure(-label => "Creating msfasta index");
   ${$args{main_window}}  ->update;
   my $status = &create_index(main_window   => $args{main_window},
                              progress_bar  => $args{progress_bar},
                              ini_ref       => $args{ini_ref},
                              auto_ini_ref  => $args{auto_ini_ref},
                              flavour       => $args{flavour}
                             );
   if ($status == 0) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No sequences to compile.\nAborting",
                                                    -bitmap  => 'error',
                                                    -buttons => ['ok']);
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error compiling custom databases".
                     "\nNo sequences to compile.\nAborting\n\n";
      close ERRORLOG;
      unlink "${$args{auto_ini_ref}}{work_dir}/temp.fasta";
      return(0);
   }

   #compile database
   ${$args{progress_bar}} ->configure(-label => "Compiling database");
   ${$args{main_window}}  ->update;

   if ($args{flavour} eq 'blast') {
      if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
         if ($args{type} eq 'aa') {
            `${$args{ini_ref}}{blast_executables}/formatdb -i ${$args{auto_ini_ref}}{work_dir}/temp.fasta -p T -o T -n $target_dir/$database_name`;
         } elsif ($args{type} eq 'nt') {
            `${$args{ini_ref}}{blast_executables}/formatdb -i ${$args{auto_ini_ref}}{work_dir}/temp.fasta -p F -o T -n $target_dir/$database_name`;
         }
      } else {
         if ($args{type} eq 'aa') {
            `${$args{ini_ref}}{blast_plus_executables}/makeblastdb -dbtype prot -in ${$args{auto_ini_ref}}{work_dir}/temp.fasta -parse_seqids -hash_index -out $target_dir/$database_name`;
         } elsif ($args{type} eq 'nt') {
            `${$args{ini_ref}}{blast_plus_executables}/makeblastdb -dbtype nucl -in ${$args{auto_ini_ref}}{work_dir}/temp.fasta -parse_seqids -hash_index -out $target_dir/$database_name`;
         }
      }
   } elsif ($args{flavour} eq 'cog') {
      #create COG temp database

      #ASSUME new BLAST+
      #if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
         ################################
      #   `${$args{ini_ref}}{blast_executables}/formatdb -i ${$args{auto_ini_ref}}{work_dir}/temp.fasta -o T -n $target_dir/COG_temp`;
      #} else {
         #################################
         my $in  = ${$args{auto_ini_ref}}{work_dir}.'/temp.fasta';
         my $out = $target_dir.'/COG_temp';
         `${$args{ini_ref}}{blast_plus_executables}/makeblastdb -dbtype prot -in $in -parse_seqids -hash_index -out $out`;
      #}
      #pass everyting to COG creation
      &make_COGs(main_window   => $args{main_window},
                 progress_bar  => $args{progress_bar},
                 ini_ref       => $args{ini_ref},
                 auto_ini_ref  => $args{auto_ini_ref},
                 target_dir    => $target_dir,
                 database_name => $database_name,
                 input_file    => ${$args{auto_ini_ref}}{work_dir}.'/temp.fasta'
                );
      &parse_COGs(main_window   => $args{main_window},
                 progress_bar  => $args{progress_bar},
                 ini_ref       => $args{ini_ref},
                 auto_ini_ref  => $args{auto_ini_ref},
                 target_dir    => $target_dir,
                 database_name => $database_name,
                 );
      unlink ${$args{auto_ini_ref}}{work_dir}.'/COG_temp.fasta'
   }

   ${$args{progress_bar}} ->configure(-label => "Compiled database $database_name. Task finished.");
   ${$args{main_window}}  ->update;

   #unlink temp files or copy
   if ($keep_source == 1) {
      copy(${$args{auto_ini_ref}}{work_dir}.'/temp.fasta', ${$args{auto_ini_ref}}{work_dir}.'/'.$database_name.'.msfasta');
      ${$args{main_window}}->messageBox(-title   => 'Info',
                                        -message => "Created FASTA source file $database_name\.msfasta in directory ${$args{auto_ini_ref}}{work_dir}\.",
                                        -icon    => 'info',
                                        -type    => 'ok');
   }
   unlink "${$args{auto_ini_ref}}{work_dir}/temp.fasta";

   ${$args{main_window}}  ->messageBox(-title   => 'Success',
                                       -message => "Compiled database $database_name",
                                       -type    => 'OK',
                                       -icon    => 'info');

   return (1);
}

sub extract_nt {
   my %args = @_;
   my @list = ();

   #get all CDS or gene entries
   @list = grep { /_CDS_/ } @{$args{feature_list_ref}};
   #no CDS feature? try gene instead
   unless (@list && $#list >= 0) {
      @list = grep { /_gene_/ } @{$args{feature_list_ref}};
   }
   #neither CDS nor gene? Error message
   unless (@list && $#list >= 0) {
      my $error = ${$args{main_window}}->Toplevel();
      $error->Label(-text => "Genbank file $args{filename} has neither CDS nor gene features.\nSkipping file")->pack();
      $error->update;
      $error->raise;
      $error->after(3000);
      $error->destroy;
      return (0);
   }

   #iterate over entries
   foreach my $entry (@list) {
      unless (exists ${$args{genbank_ref}}{$entry}) {
         my $error = ${$args{main_window}}->Toplevel();
         $error->Label(-text => "Could not access entry $entry in genbank file $args{filename}.")->pack();
         $error->update;
         $error->raise;
         $error->after(3000);
         $error->destroy;
         return (0);
      }
      #grab designation and orientation
      my ($name, $left_bd, $right_bd, $orientation, $seq);
      if (${$args{genbank_ref}}{$entry} =~ m/\/gene\=/s) {
         ${$args{genbank_ref}}{$entry} =~ m/^\w+ +(complement)?.*?   \/gene\=\"([^\"]*?)\"/s;
         $orientation = $1;
         $name = $2;
      } elsif (${$args{genbank_ref}}{$entry} =~ m/\/locus_tag\=/s) {
         ${$args{genbank_ref}}{$entry} =~ m/^\w+ +(complement)?.*?   \/locus_tag\=\"([^\"]*?)\"/s;
         $orientation = $1;
         $name = $2;
      } elsif (${$args{genbank_ref}}{$entry} =~ m/\/product\=/s) {
         ${$args{genbank_ref}}{$entry} =~ m/^[\w]+ +(complement)?.*?   \/product\=\"([^\"]*?)\"/s;
         $orientation = $1;
         $name = $2;
      }

      if (defined $orientation && $orientation =~ /complement/) {
         $orientation = 'antisense';
      } else {
         $orientation = 'sense';
      }
      if (defined $name && $name =~ /\w+/) {
         $name =~ s/\s+/ /gs;
         $name =~ s/ /_/g;
         $name =~ s/\W/_/g;
      } else {
         $name = 'no_annotation';
      }

      #grab boundaries and orientation, check for joined features
      if (${$args{genbank_ref}}{$entry} =~ m/^\s*(gene|CDS).*?join/) {
         ($seq) = &parse_joined_features (feature   => ${$args{genbank_ref}}{$entry},
                                          DNAlength => length($args{DNAseq})
                                         );
      } else {
         $left_bd = "";
         $right_bd = "";
         $entry =~ m/(\d+)_(\d+)_/;
         $left_bd = $1;
         $right_bd = $2;

         #extract nt sequence
         $seq = "";
         $seq = substr($args{DNAseq},($left_bd - 1),($right_bd-$left_bd+1));
      }

      #orient sequence
      if ($orientation eq 'antisense') {
         $seq = reverse ($seq);
         $seq =~ tr/acgtmkrybvdhACGTMKRYBVDH/tgcakmyrvbhdTGCAKMYRVBHD/;
         $seq = lc $seq;
      }

      #write into temp file
      print WRITE ">$name __filename__$args{filename}__\n$seq\n";
   }
}

sub extract_aa {
   my %args = @_;
   my @list = ();

   #get all CDS or gene entries
   @list = grep { /_CDS_/ } @{$args{feature_list_ref}};
   #no CDS feature? try gene instead
   unless (@list && $#list >= 0) {
         @list = grep { /_gene_/ } @{$args{feature_list_ref}};
   }
   #neither CDS nor gene? Error message
   unless (@list && $#list >= 0) {
      my $error = ${$args{main_window}}->Toplevel();
      $error->Label(-text => "Genbank file $args{filename} has neither CDS nor gene features.\nSkipping file")->pack();
      $error->update;
      $error->raise;
      $error->after(3000);
      $error->destroy;
      return (0);
   }

   #iterate over entries
   foreach my $entry (@list) {
      unless (exists ${$args{genbank_ref}}{$entry}) {
         my $error = ${$args{main_window}}->Toplevel();
         $error->Label(-text => "Could not access entry $entry in genbank file $args{filename}.")->pack();
         $error->update;
         $error->raise;
         $error->after(3000);
         $error->destroy;
         return (0);
      }
      #grab designation and orientation
      my ($name, $left_bd, $right_bd, $orientation, $seq);
      if (${$args{genbank_ref}}{$entry} =~ m/\/gene\=/s) {
         ${$args{genbank_ref}}{$entry} =~ m/^\w+ +(complement)?.*?   \/gene\=\"([^\"]*?)\"/s;
         $orientation = $1;
         $name = $2;
      } elsif (${$args{genbank_ref}}{$entry} =~ m/\/locus_tag\=/s) {
         ${$args{genbank_ref}}{$entry} =~ m/^\w+ +(complement)?.*?   \/locus_tag\=\"([^\"]*?)\"/s;
         $orientation = $1;
         $name = $2;
      } elsif (${$args{genbank_ref}}{$entry} =~ m/\/product\=/s) {
         ${$args{genbank_ref}}{$entry} =~ m/^[\w]+ +(complement)?.*?   \/product\=\"([^\"]*?)\"/s;
         $orientation = $1;
         $name = $2;
      }

      if (defined $orientation && $orientation =~ /complement/) {
         $orientation = 'antisense';
      } else {
         $orientation = 'sense';
      }
      if (defined $name && $name =~ /\w+/) {
         $name =~ s/\s+/ /gs;
         $name =~ s/ /_/g;
         $name =~ s/\W/_/g;
      } else {
         $name = 'no_annotation';
      }

      #grab boundaries and orientation
      $left_bd = "";
      $right_bd = "";
      $entry =~ m/(\d+)_(\d+)_/;
      $left_bd = $1;
      $right_bd = $2;

      #extract nt sequence
      $seq = "";
      $seq = substr($args{DNAseq},($left_bd - 1),($right_bd-$left_bd+1));

      #translate to aa print
      my $aa_seq  = &nt2aa(main_window      => $args{main_window},
                           progress_bar     => $args{progress_bar},
                           ini_ref          => $args{ini_ref},
                           module           => 'customDB',
                           auto_ini_ref     => $args{auto_ini_ref},
                           DNAseq           => $args{DNAseq},
                           genbank_ref      => $args{genbank_ref},
                           feature_list_ref => $args{feature_list_ref},
                           filename         => $args{filename},
                           left_bd          => $left_bd,
                           right_bd         => $right_bd,
                           orientation      => $orientation,
                           sequence         => $seq,
                          );
      if ($aa_seq eq '0') {next};

      #write into temp file
      print WRITE ">$name __filename__$args{filename}__\n$aa_seq\n";
   }
}

sub create_index {
   my %args = @_;
   my $gb_counter = 1;
   my @msfasta = ();

   #read temp file to memory
   my $fasta_ref = &slurp(main_window  => $args{main_window},
                          directory    => ${$args{auto_ini_ref}}{work_dir},
                          filename     => 'temp.fasta',
                         );
   #empty file? return
   if ($fasta_ref eq '0') {return (0)};
   unless (defined ${$fasta_ref} && ${$fasta_ref} =~ /\w+/) {
      return (0);
   }
   ${$fasta_ref} = "\n".${$fasta_ref};
   #split into array
   @msfasta = split/\n\>/,${$fasta_ref};

   if ($args{flavour} eq 'cog') {
      &make_COG_index(main_window   => $args{main_window},
                      progress_bar  => $args{progress_bar},
                      ini_ref       => $args{ini_ref},
                      auto_ini_ref  => $args{auto_ini_ref},
                      directory     => ${$args{auto_ini_ref}}{work_dir},
                      fasta_file    => $fasta_ref,
                     );
   }

   #update msfasta file to index
   unlink "${$args{auto_ini_ref}}{work_dir}/temp.fasta";
   open WRITE, "+>${$args{auto_ini_ref}}{work_dir}/temp.fasta";
   if ($args{flavour} eq 'cog') {
      open COG, "+>${$args{auto_ini_ref}}{work_dir}/COG_temp.fasta";
   }
   #add index to fasta entries
   my %seen = ();
   foreach my $entry (@msfasta) {
      if ($entry !~ /\w+/) {next};
      $entry =~ s/^>//g;

      if ($entry =~ /(.*?) __filename__(.*?)__\n/) {
         $entry =~ s/(.*?) __filename__(.*?)__\n/\>gi\|$gb_counter\|gb\|$gb_counter\| \($2\) $1\n/s;
      } elsif ($entry =~ /^\w*?\|/) {
         'reset' =~ m/reset/;
         $entry =~ m/([^\|]*?)\|/;
         my $id = $1;
         if (exists $seen{$id}) {$id = $id.'_'.$gb_counter};
         $seen{$id} = 1;
         $entry =~ s/^[^\|]*?\|/\>gi|$gb_counter\|gb\|$id\| /;

      } else {
         $entry = "\>gi\|$gb_counter\|gb\|$gb_counter\| ".$entry;
      }
      print WRITE $entry."\n";

      #modify index file for COG, real fasta header association is stored in separate file
      if ($args{flavour} eq 'cog') {
         $entry =~ s/>.+?\n//s;
         $entry = ">$gb_counter\n".$entry;
         print COG $entry."\n";
      }

      $gb_counter++;
   }
   @msfasta = ();
   close WRITE;
   if ($args{flavour} eq 'cog') {
      close COG;
   }
   undef %seen;
   return (1);
}

sub extract_genbank_seq {
   my %args = @_;

   my $answer = ${$args{main_window}}->Dialog(-title          => 'Process genbank keys?',
                                              -text           => "Extract nt sequence from keys or use whole sequence?",
                                              -default_button => 'Extract',
                                              -buttons        => [ 'Extract', 'Whole Sequence' ],
                                              -bitmap         => 'question')->Show();
   if ($answer eq 'Extract') {
      return ('extract');
   } elsif ($answer eq 'Whole Sequence') {
      return ('sequence');
   } else {
      return ('extract');
   }
}

sub determine_seq_type {
   my %args = @_;

   #get percentage of actg and if greater then 80% then assume its nucleotide, else aa
   my $count = ($args{seq} =~ tr/[AaCcGgTt]//);
   return (0) if ($count < 1);

   if ($count / length($args{seq}) > 0.8 ) {
       return ('nt');
   } else {
      return ('aa');
   }
}

sub parse_joined_features {
   my %args = @_;
   my (@boundaries, $joined_seq, $left_bd, $right_bd);
   #grab all boundaries in joined feature
   'reset' =~ m/reset/;
   $args{feature} =~ m/^([^\/]*?)\//s;
   my $boundaries = $1;
   @boundaries = ($boundaries =~ m/(\d+\D*?\d+)/g);

   #go through boundaries
   foreach my $entry (@boundaries) {
      'reset' =~ m/reset/;
      $entry =~ m/(\d+)\D*?(\d+)/;
      $left_bd = $1;
      $right_bd = $2;
      unless (defined $left_bd && defined $right_bd) {return 0};

      #extract nt sequence
      my $seq = "";
      $seq = substr($args{DNAseq},($left_bd - 1),($right_bd-$left_bd+1));
      $joined_seq .= $seq;
   }
   return ($joined_seq);
}

sub make_COGs {
   my %args             = @_;
   my $COG_db           = $args{target_dir}.'/COG_temp';
   my $Blast_no         = ${$args{ini_ref}}{COG_db_path}.'/BLASTno';
   my $Blast_ff         = ${$args{ini_ref}}{COG_db_path}.'/BLASTff';
   my $Blast_conv       = ${$args{ini_ref}}{COG_db_path}.'/BLASTconv';
   my $COG_Blast_no_res = $Blast_no.'/'.$args{database_name}.'.tab';
   my $COG_Blast_ff_res = $Blast_ff.'/'.$args{database_name}.'.tab';
   #test for directories
   unless (-d $Blast_no) {
      mkdir($Blast_no, 0777) ;
   }
   unless (-d $Blast_ff) {
      mkdir($Blast_ff, 0777) ;
   }
   unless (-d $Blast_conv) {
      mkdir($Blast_conv, 0777) ;
   }

   #run Blast
   #compile database
   #if input file >1000 individual sequences, split and submit sequentially.
   my $file_ref = &slurp_cmd(ini_ref       => $args{ini_ref},
                             auto_ini_ref  => $args{auto_ini_ref},
                             directory     => ${$args{auto_ini_ref}}{work_dir},
                             filename      => 'COG_temp.fasta'
                            );
   my $total_seq_number = ($$file_ref =~ tr/\>//s);
   undef $file_ref;
   if ($total_seq_number > 1000) {
      my $local_count    = 0;
      my $index          = 0;
      my $subset         = '';
      my %file_hash      = ();
      my $unfiltered_res = '';
      my $filtered_res   = '';

      ${$args{progress_bar}} ->configure(-label => "Splitting input file into smaller chunks");
      ${$args{main_window}}  ->update;
      open READ, "$args{input_file}";
      while (<READ>) {
         if ($_ =~ /^>/) {$local_count++};
         if ($local_count > 1000) {
            $index++;
            open WRITE, "+>$args{input_file}"."_$index";
            print WRITE $subset;
            close WRITE;
            push (@ARGV, $args{input_file}."_$index");
            $file_hash{$index}= $args{input_file}."_$index";
            $local_count = 0;
            $subset      = '';
         }
         $subset .= $_;
      }
      close READ;
      #catch last set
      $index++;
      open WRITE, "+>$args{input_file}"."_$index";
      print WRITE $subset;
      close WRITE;
      $file_hash{$index}= $args{input_file}."_$index";
      $subset = '';

      #iterate over input chunks
      foreach my $key (sort {$a<=>$b} keys %file_hash) {
         chdir ${$args{ini_ref}}{blast_plus_executables};
         ${$args{progress_bar}} ->configure(-label => "Running unfiltered COG Blast for index $key");
         ${$args{main_window}}  ->update;
         my $query = $args{input_file}.'_'.$key;
         ############################
         #`./blastpgp -a ${$args{auto_ini_ref}}{CPU} -d $COG_db -i $args{input_file}\_$key -I T -m 9 -b 1000 -v 1000 -z 100000000 -F F -t F -o $COG_Blast_no_res\.$key`; # unfiltered BLAST results in the ./BLASTno/ directory
         #`./psiblast -num_threads ${$args{auto_ini_ref}}{CPU} -query $query -db $COG_db -show_gis -outfmt 7 -num_descriptions 1000 -num_alignments 1000 -dbsize 100000000 -comp_based_stats F -seg no -out $COG_Blast_no_res`; # unfiltered BLAST results in the ./BLASTno/ directory
         `./psiblast -num_threads ${$args{auto_ini_ref}}{CPU} -query $query -db $COG_db -show_gis -outfmt 7 -max_target_seqs 1000 -dbsize 100000000 -comp_based_stats F -seg no -out $COG_Blast_no_res\.$key`; # unfiltered BLAST results in the ./BLASTno/ directory
         ${$args{progress_bar}} ->configure(-label => "Running filtered COG Blast for index $key");
         ${$args{main_window}}  ->update;
         ############################
         #`./blastpgp -a ${$args{auto_ini_ref}}{CPU} -d $COG_db -i $args{input_file}\_$key -I T -m 9 -b 1000 -v 1000 -z 100000000 -F T -t T -o $COG_Blast_ff_res\.$key`; # filtered BLAST results in the ./BLASTff/ directory
         #`./psiblast -num_threads ${$args{auto_ini_ref}}{CPU} -query $query -db $COG_db -show_gis -outfmt 7 -num_descriptions 1000 -num_alignments 1000 -dbsize 100000000 -comp_based_stats T -seg yes -out $COG_Blast_ff_res`; # filtered BLAST results in the ./BLASTff/ directory
         `./psiblast -num_threads ${$args{auto_ini_ref}}{CPU} -query $query -db $COG_db -show_gis -outfmt 7 -max_target_seqs 1000 -dbsize 100000000 -comp_based_stats T -seg yes -out $COG_Blast_ff_res\.$key`; # filtered BLAST results in the ./BLASTff/ directory
         #remove temp files
         unlink $args{input_file}."_$key";
      }

      #merge output files
      ${$args{progress_bar}} ->configure(-label => "Merging COG Blast output files");
      ${$args{main_window}}  ->update;

      foreach my $key (sort {$a<=>$b} keys %file_hash) {
         my $file_ref = &slurp_cmd(ini_ref       => $args{ini_ref},
                                   auto_ini_ref  => $args{auto_ini_ref},
                                   directory     => $Blast_no,
                                   filename      => $args{database_name}.'.tab'."\.$key"
                                  );
         $unfiltered_res .= $$file_ref;
         undef $file_ref;
         $file_ref = &slurp_cmd(ini_ref       => $args{ini_ref},
                                auto_ini_ref  => $args{auto_ini_ref},
                                directory     => $Blast_ff,
                                filename      => $args{database_name}.'.tab'."\.$key"
                              );
        $filtered_res .= $$file_ref;
         unlink $Blast_no.'/'.$args{database_name}.'.tab'."\.$key";
         unlink $Blast_ff.'/'.$args{database_name}.'.tab'."\.$key";
      }

      #write combined results
      open WRITE, "+>".$Blast_no.'/'.$args{database_name}.'.tab';
      print WRITE $unfiltered_res;
      close WRITE;
      open WRITE, "+>".$Blast_ff.'/'.$args{database_name}.'.tab';
      print WRITE $unfiltered_res;
      close WRITE;
      undef %file_hash;
      undef $unfiltered_res;
      undef $filtered_res;
   } else {
      #change to Blast directory
      chdir ${$args{ini_ref}}{blast_executables};
      ${$args{progress_bar}} ->configure(-label => "Running unfiltered COG Blast");
      ${$args{main_window}}  ->update;
      #`./blastpgp -a ${$args{auto_ini_ref}}{CPU} -d $COG_db -i $args{input_file} -I T -m 9 -b 1000 -v 1000 -z 100000000 -F F -t F -o $COG_Blast_no_res`; # unfiltered BLAST results in the ./BLASTno/ directory
      #`./psiblast -num_threads ${$args{auto_ini_ref}}{CPU} -db $COG_db -query $args{input_file} -show_gis -outfmt 7 -num_alignments 1000 -num_descriptions 1000 -dbsize 100000000 -comp_based_stats F -seg no -out $COG_Blast_no_res`; # unfiltered BLAST results in the ./BLASTno/ directory
      `./psiblast -num_threads ${$args{auto_ini_ref}}{CPU} -query $args{input_file} -db $COG_db -show_gis -outfmt 7 -max_target_seqs 1000 -dbsize 100000000 -comp_based_stats F -seg no -out $COG_Blast_no_res`; # unfiltered BLAST results in the ./BLASTno/ directory
      ${$args{progress_bar}} ->configure(-label => "Running filtered COG Blast");
      ${$args{main_window}}  ->update;
      #`./blastpgp -a ${$args{auto_ini_ref}}{CPU} -d $COG_db -i $args{input_file} -I T -m 9 -b 1000 -v 1000 -z 100000000 -F T -t T -o $COG_Blast_ff_res`; # filtered BLAST results in the ./BLASTff/ directory
      #`./psiblast -num_threads ${$args{auto_ini_ref}}{CPU} -db $COG_db -query $args{input_file} -show_gis -outfmt 7 -num_alignments 1000 -num_descriptions 1000 -dbsize 100000000 -comp_based_stats T -seg yes -out $COG_Blast_ff_res`; # filtered BLAST results in the ./BLASTff/ directory
      `./psiblast -num_threads ${$args{auto_ini_ref}}{CPU} -query $args{input_file} -db $COG_db -show_gis -outfmt 7 -max_target_seqs 1000 -dbsize 100000000 -comp_based_stats T -seg yes -out $COG_Blast_ff_res`; # filtered BLAST results in the ./BLASTff/ directory
   }

   #Preparation of the "sequence Universe"
   ${$args{progress_bar}} ->configure(-label => "Preparation of the sequence Universe: COGmakehash");
   ${$args{main_window}}  ->update;
   chdir ${$args{ini_ref}}{COGcognitor_path};
   `./COGmakehash -i=${$args{ini_ref}}{COG_db_path}\/$args{database_name}\.p2o\.csv -o=$Blast_conv -s="," -n=1`; # makes ./BLASTconv/hash.csv file

   #Processing of BLAST results
   ${$args{progress_bar}} ->configure(-label => "Processing of BLAST results: COGreadblast");
   ${$args{main_window}}  ->update;
   `./COGreadblast -d=$Blast_conv -u=$Blast_no -f=$Blast_ff -s=$Blast_no -e=0.1 -q=2 -t=2`;

   #Lineage specific expansions
   ${$args{progress_bar}} ->configure(-label => "Lineage specific expansions: COGlse");
   ${$args{main_window}}  ->update;
   `./COGlse -t=${$args{ini_ref}}{COG_db_path} -d=$Blast_conv -j=${$args{ini_ref}}{COG_db_path}\/$args{database_name}\.job\.csv -p=${$args{ini_ref}}{COG_db_path}\/$args{database_name}\.p2o\.csv -o=${$args{ini_ref}}{COG_db_path}\/$args{database_name}\.lse\.csv`;

   #Making clusters from symmetrical best hits
   ${$args{progress_bar}} ->configure(-label => "Making clusters from symmetrical best hits: COGtriangles");
   ${$args{main_window}}  ->update;
   `./COGtriangles -i=$Blast_conv -q=${$args{ini_ref}}{COG_db_path}\/$args{database_name}\.p2o\.csv -l=${$args{ini_ref}}{COG_db_path}\/$args{database_name}\.lse\.csv -o=${$args{ini_ref}}{COG_db_path}\/$args{database_name}\.cls\.csv -t=0.5 -e=0.01 -n="CLS" -s=1`;
   `./COGtriangles.reformat.pl ${$args{ini_ref}}{COG_db_path}\/$args{database_name}\.cls\.csv ${$args{ini_ref}}{COG_db_path}\/$args{database_name}\.multi_cls\.csv`;
   move (${$args{ini_ref}}{COGcognitor_path}.'/all-edges.txt', ${$args{ini_ref}}{COG_db_path}.'/all-edges.txt');
   move (${$args{ini_ref}}{COGcognitor_path}.'/cog-edges.txt', ${$args{ini_ref}}{COG_db_path}.'/cog-edges.txt');

   #cleanup temp files
   chdir ${$args{ini_ref}}{COG_db_path};
   opendir SEQINPUT, ${$args{ini_ref}}{COG_db_path};
   my @temp = grep /^_(f)?hits_.+\.tmp$/, readdir(SEQINPUT);
   closedir SEQINPUT;
   unlink @temp;
   unlink ${$args{ini_ref}}{COG_db_path}.'/_hash.tmp';
   opendir SEQINPUT, ${$args{ini_ref}}{COG_db_path};
   @temp = grep /^COG_temp/, readdir(SEQINPUT);
   closedir SEQINPUT;
   unlink @temp;
   undef @temp;
   chdir ${$args{auto_ini_ref}}{work_dir};
}


sub make_COG_index {
   my %args    = @_;
   my %genomes = ();
   my $counter = 1;
   ${$args{progress_bar}} ->configure(-label => "Creating COG index file for database $database_name");
   ${$args{main_window}}  ->update;
   #read fasta file and write into index file
   open WRITE, "+>${$args{ini_ref}}{COG_db_path}\/$database_name\.p2o\.csv";
   open COG,   "+>${$args{ini_ref}}{COG_db_path}\/$database_name\.seq_header_association";
   my @temp = split /\n/, ${$args{fasta_file}};
   my @entries = grep (/^\>/, @temp);
   undef @temp;
   foreach my $entry (@entries) {
      'reset' =~ m/reset/;
      $entry =~ m/^\>(.*?) __filename__(.*?)__/;
      my ($prot_id, $filename) = ($1, $2);
      $genomes{$filename}++;
      print WRITE $counter.','.$filename."\n";
      print COG   $counter."\t".$prot_id."\t".$filename."\n";
      $counter++;
   }
   close WRITE;
   close COG;
   #create job-description file
   ${$args{progress_bar}} ->configure(-label => "Creating COG job-description file for database $database_name");
   ${$args{main_window}}  ->update;
   open WRITE, "+>${$args{ini_ref}}{COG_db_path}\/$database_name\.job\.csv";
   foreach my $query (sort keys %genomes) {
      foreach my $subject (sort keys %genomes) {
         next if ($query eq $subject);
         print WRITE $query.','.$subject."\n";
      }
   }
   close WRITE;
}

sub parse_COGs {
   my %args     = @_;
   my $COG_cluster;
   my $COGs;

   #read seq_header_association
   ${$args{progress_bar}} ->configure(-label => "Parsing COG results: reading association");
   ${$args{main_window}}  ->update;
   open READ, "<${$args{ini_ref}}{COG_db_path}\/$args{database_name}\.seq_header_association";
   while (<READ>) {
      next unless ($_ =~ /\S+/);
      my @header                   = split (/\t/, $_);
      my ($prot_id, $ORF, $genome) = @header[0,1,2];
      unless (defined $genome && $prot_id =~ /\d+/) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error Parsing COG sequence header association",
                                                       -bitmap  => 'error',
                                                       -buttons => ['ok']);
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error creating COG database $args{database_name}".
                        "\nError Parsing COG sequence header association for entry $_\n\n";
         close ERRORLOG;
         next;
      }
      $genome =~ s/\n//;

      $COGs->{ $prot_id } = {
                             ORF        => $ORF,
                             genome_id  => $genome,
                             cluster_id => '',
                             seq        => ''
                            };
   }
   close READ;

   #read edge data
   ${$args{progress_bar}} ->configure(-label => "Parsing COG results: reading edge data");
   ${$args{main_window}}  ->update;
   open READ, "<${$args{ini_ref}}{COG_db_path}\/cog-edges.txt";
   while (<READ>) {
      next unless ($_ =~ /\S+/);
      my @edge                 = split (/,/, $_);
      my ($cluster, $prot_id ) = @edge[0,1];
      unless (defined $cluster && $prot_id =~ /\d+/) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error Parsing COG edges",
                                                       -bitmap  => 'error',
                                                       -buttons => ['ok']);
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error creating COG database $args{database_name}".
                        "\nError Parsing COG edges for entry $_\n\n";
         close ERRORLOG;
         next;
      }
      if ($COGs->{$prot_id}->{cluster_id} !~ /$cluster/) {
         $COGs->{$prot_id}->{cluster_id} = $COGs->{$prot_id}->{cluster_id}.', '.$cluster;
      }
      $COG_cluster->{$cluster}->{$prot_id}->{id_count}++;
      if ($COG_cluster->{$cluster}->{$prot_id}->{id_count} < 2) {
         $COG_cluster->{$cluster}->{count}++;
      }
   }
   close READ;

   #read sequence data
   ${$args{progress_bar}} ->configure(-label => "Parsing COG results: reading sequences");
   ${$args{main_window}}  ->update;
   my $file_ref = &slurp_cmd(ini_ref       => $args{ini_ref},
                             auto_ini_ref  => $args{auto_ini_ref},
                             directory     => ${$args{auto_ini_ref}}{work_dir},
                             filename      => 'COG_temp.fasta'
                            );
   my @seq = split (/>/, $$file_ref);
   foreach my $entry (@seq) {
      next unless ($entry =~ /\S+/);
      'reset' =~ m/reset/;
      $entry =~ m/(\d+)\n(.+)/s;
      my ($prot_id, $seq) = ($1, $2);
      unless (defined $seq && $prot_id =~ /\d+/) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error Parsing COG fasta",
                                                       -bitmap  => 'error',
                                                       -buttons => ['ok']);
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error creating COG database $args{database_name}".
                        "\nError Parsing COG fasta file for entry $entry\n\n";
         close ERRORLOG;
         next;
      }
      $seq =~ /\s/gs;

      next if ($COGs->{$prot_id}->{seq} =~ /\S+/);
      $COGs->{$prot_id}->{seq} = $seq;
   }
   undef $file_ref;
   unlink "${$args{auto_ini_ref}}{work_dir}\/COG_temp.fasta";

   #create new fasta file and foramt into Blastdb
   ${$args{progress_bar}} ->configure(-label => "Parsing COG results: writing new FASTA file");
   ${$args{main_window}}  ->update;
   open WRITE, "+>${$args{ini_ref}}{COG_db_path}\/$args{database_name}\.fasta";
   foreach my $prot_id (sort {$a<=>$b} keys %$COGs) {
      my $final_COG_clusters; #in case some do not reach min number limit
      #skip COGs without a cluster association if selected
      next if ($keep_orphans == 0 && $COGs->{$prot_id}->{cluster_id} !~ /\S+/);
      my $gi         = $prot_id + 10000;
      my $seq_length = length($COGs->{$prot_id}->{seq});
      $COGs->{$prot_id}->{cluster_id} =~ s/^\s*,\s*//;
      #skip COG if minimum number of orthologs is not reached
      my @temp = split /\s*,\s*/, $COGs->{$prot_id}->{cluster_id};
      foreach my $cluster (@temp) {
         next if ($COG_cluster->{$cluster}->{count} < $min_number_of_orthologs && $min_number_of_orthologs > 0);
         $final_COG_clusters .= $cluster.', ';
      }
      $final_COG_clusters =~ s/\,\s*$//;

      print WRITE '>gi|'. $gi . ' '. $COGs->{$prot_id}->{ORF} .','.
                  $gi .','. $seq_length . ',1,' . $seq_length .','.
                  $final_COG_clusters .','. $COGs->{$prot_id}->{genome_id}."\n".
                  $COGs->{$prot_id}->{seq}."\n";
   }
   close WRITE;

   #create COG Blast database
   ${$args{progress_bar}} ->configure(-label => "Parsing COG results: creating new COG database");
   ${$args{main_window}}  ->update;
   #ASSUME Blast+ is present
   #if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
      ################################
   #   `${$args{ini_ref}}{blast_executables}/formatdb -i ${$args{ini_ref}}{COG_db_path}\/$args{database_name}\.fasta -o T -n $target_dir/$args{database_name}`;
   #} else {
      #################################
      my $in  = ${$args{ini_ref}}{COG_db_path}.'/'.$args{database_name}.'.fasta';
      my $out = $target_dir.'/'.$args{database_name};
      `${$args{ini_ref}}{blast_plus_executables}/makeblastdb -dbtype prot -in $in -parse_seqids -hash_index -out $out`;
   #}

   #create statistics
   ${$args{progress_bar}} ->configure(-label => "Parsing COG results: writing statistics");
   ${$args{main_window}}  ->update;
   my $total_cluster_number = keys %$COG_cluster;
   open WRITE, "+>${$args{ini_ref}}{COG_db_path}\/$args{database_name}\.statistics";
   print WRITE 'Total number of COG clusters for dataset '. $args{database_name} .': '. $total_cluster_number ."\n";
   print WRITE "Number of members for each cluster\n";
   foreach my $cluster_id (sort {substr($a,3)<=>substr($b,3)} keys %$COG_cluster) {
      next if ($COG_cluster->{$cluster_id}->{count} < $min_number_of_orthologs && $min_number_of_orthologs > 0);
      print WRITE "\t$cluster_id\t$COG_cluster->{$cluster_id}->{count}\n";
   }
   close WRITE;
}

1;