#!/opt/ActivePerl-5.8/bin/perl

#update_db: update databases  via internet
#input arguments: main_window, ini_ref, auto_ini_ref,


package Miscellaneous::update_db;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&update_db);
use vars qw();

use Tk;
use initialise::read_me  qw(:DEFAULT);
use Basics::progress_bar qw(:DEFAULT);
use Cwd;
use File::Copy;
use File::stat;
use Time::localtime;
use Net::FTP;

#local variables
my (%args, $status, $tl,
    $COG_server, $COG_path, $COG,
    $COG2014_server, $COG2014_path, $COG2014,
    $arCOG2014_server, $arCOG2014_path, $arCOG2014,
    $POG2013_server, $POG2013_path, $POG2013,
    $pfam_server, $pfam_path, $pfam, $pfam_descr_path, $pfam_descr, $pfam_ipro_path, $pfam_ipro,
    $tigrfam_server, $tigrfam_path, $tigrfam, $go_server, $go_path, $go,
    $ncRNA_server, $ncRNA_path, $ncRNA);


sub update_db {
   my %args = @_;
   my ($frame_top, $frame_bottom, $progress_bar,
       %state, %update);

   if (defined $tl && Tk::Exists($tl)) {
      $tl->state('normal');
      $progress_bar = $tl->Frame(-borderwidth => 2, -relief => 'groove');
      $frame_top    = $tl->Frame()                                      ;
      $frame_bottom = $tl->Frame()                                      ;
   } else {
      $tl = ${$args{main_window}}->Toplevel(-title   => 'Update common databases');
      $progress_bar = $tl->Frame(-borderwidth => 2, -relief => 'groove') -> pack(-side => 'bottom', -fill => 'x');
      $frame_top    = $tl->Frame()                                       -> pack(-side => 'top',    -fill => 'both', -expand => 1);
      $frame_bottom = $tl->Frame()                                       -> pack(-side => 'top',    -fill => 'x');
   }


   $progress_bar->configure(-label => " "); #create dummy for space

   #options into top frame
   $frame_top->Label(-text => 'Update functional databases')->grid(-row => 0, -column => 0, -sticky => 'nsew', -columnspan => 3);
   $frame_top->gridRowconfigure(0, -pad => 30);

   #configure initial state
   unless (defined $state{'update_COG'}) {
      #set hash
      $state{'update_COG'}      = 0; $state{'update_pfam'}       = 0; $state{'update_pfam_descr'} = 0;
      $state{'update_COG2014'}  = 0; $state{'update_arCOG2014'}  = 0; $state{'update_POG2013'}    = 0;
      $state{'update_interpro'} = 0; $state{'update_tigrfam'}    = 0; $state{'update_go'}         = 0;
      $state{'update_ncRNA'}    = 0;
   }


   #COG
   {
      $COG_server   = $frame_top -> LabEntry(-label        => 'Server:',
                                              -labelPack    => [-side => "left", -anchor => "w"],
                                              -justify      => 'left',
                                              -width        => 20,
                                              -textvariable => \${$args{ini_ref}}{COG_server_URL}
                                             )
      ->grid(-row => 1, -column => 1, -sticky => 'ew');
      $COG_path   = $frame_top -> LabEntry(-label        => 'Path:',
                                              -labelPack    => [-side => "left", -anchor => "w"],
                                              -justify      => 'left',
                                              -width        => 20,
                                              -textvariable => \${$args{ini_ref}}{COG_remote_dir}
                                             )
      ->grid(-row => 1, -column => 2, -sticky => 'ew');
      $COG        = $frame_top -> Checkbutton(-text     => 'Update COG database',
                                                 -variable => \$state{'update_COG'},
                                                 -command  => sub {
                                                                   if ($state{'update_COG'} == 1) {
                                                                      $COG_server->configure(-state      => 'normal',
                                                                                             -foreground => 'black'
                                                                                             );
                                                                      $COG_path  ->configure(-state      => 'normal',
                                                                                             -foreground => 'black'
                                                                                             );
                                                                   } elsif ($state{'update_COG'} == 0) {
                                                                      $COG_server->configure(-state      => 'disabled',
                                                                                             -foreground => 'grey'
                                                                                             );
                                                                      $COG_path  ->configure(-state      => 'disabled',
                                                                                             -foreground => 'grey'
                                                                                             );
                                                                   }
                                                                 }
                                                )
      ->grid(-row => 1, -column => 0, -sticky => 'w');
   }

   #COG2014
   {
      $COG2014_server   = $frame_top -> LabEntry(-label        => 'Server:',
                                                 -labelPack    => [-side => "left", -anchor => "w"],
                                                 -justify      => 'left',
                                                 -width        => 20,
                                                 -textvariable => \${$args{ini_ref}}{COG_server_URL}
                                                )
      ->grid(-row => 2, -column => 1, -sticky => 'ew');
      $COG2014_path   = $frame_top -> LabEntry(-label        => 'Path:',
                                               -labelPack    => [-side => "left", -anchor => "w"],
                                               -justify      => 'left',
                                               -width        => 20,
                                               -textvariable => \${$args{ini_ref}}{COG2014_remote_dir}
                                              )
      ->grid(-row => 2, -column => 2, -sticky => 'ew');
      $COG2014        = $frame_top -> Checkbutton(-text     => 'Update COG2014 database',
                                                  -variable => \$state{'update_COG2014'},
                                                  -command  => sub {
                                                                    if ($state{'update_COG2014'} == 1) {
                                                                       $COG2014_server->configure(-state      => 'normal',
                                                                                                  -foreground => 'black'
                                                                                                  );
                                                                       $COG2014_path  ->configure(-state      => 'normal',
                                                                                                  -foreground => 'black'
                                                                                                  );
                                                                    } elsif ($state{'update_COG2014'} == 0) {
                                                                       $COG2014_server->configure(-state      => 'disabled',
                                                                                                  -foreground => 'grey'
                                                                                                  );
                                                                       $COG2014_path  ->configure(-state      => 'disabled',
                                                                                                  -foreground => 'grey'
                                                                                                  );
                                                                    }
                                                                 }
                                                )
      ->grid(-row => 2, -column => 0, -sticky => 'w');
   }

   #arCOG2014
   {
      $arCOG2014_server   = $frame_top -> LabEntry(-label        => 'Server:',
                                                 -labelPack    => [-side => "left", -anchor => "w"],
                                                 -justify      => 'left',
                                                 -width        => 20,
                                                 -textvariable => \${$args{ini_ref}}{COG_server_URL}
                                                )
      ->grid(-row => 3, -column => 1, -sticky => 'ew');
      $arCOG2014_path   = $frame_top -> LabEntry(-label        => 'Path:',
                                               -labelPack    => [-side => "left", -anchor => "w"],
                                               -justify      => 'left',
                                               -width        => 20,
                                               -textvariable => \${$args{ini_ref}}{arCOG2014_remote_dir}
                                              )
      ->grid(-row => 3, -column => 2, -sticky => 'ew');
      $arCOG2014        = $frame_top -> Checkbutton(-text     => 'Update arCOG2014 database',
                                                    -variable => \$state{'update_arCOG2014'},
                                                    -command  => sub {
                                                                      if ($state{'update_arCOG2014'} == 1) {
                                                                         $arCOG2014_server->configure(-state      => 'normal',
                                                                                                      -foreground => 'black'
                                                                                                      );
                                                                         $arCOG2014_path  ->configure(-state      => 'normal',
                                                                                                      -foreground => 'black'
                                                                                                      );
                                                                      } elsif ($state{'update_arCOG2014'} == 0) {
                                                                         $arCOG2014_server->configure(-state      => 'disabled',
                                                                                                      -foreground => 'grey'
                                                                                                      );
                                                                         $arCOG2014_path  ->configure(-state      => 'disabled',
                                                                                                      -foreground => 'grey'
                                                                                                      );
                                                                      }
                                                                     }
                                                    )
      ->grid(-row => 3, -column => 0, -sticky => 'w');
   }

   #POG2013
   {
      $POG2013_server   = $frame_top -> LabEntry(-label        => 'Server:',
                                                 -labelPack    => [-side => "left", -anchor => "w"],
                                                 -justify      => 'left',
                                                 -width        => 20,
                                                 -textvariable => \${$args{ini_ref}}{COG_server_URL}
                                                )
      ->grid(-row => 4, -column => 1, -sticky => 'ew');
      $POG2013_path   = $frame_top -> LabEntry(-label        => 'Path:',
                                               -labelPack    => [-side => "left", -anchor => "w"],
                                               -justify      => 'left',
                                               -width        => 20,
                                               -textvariable => \${$args{ini_ref}}{POG2013_remote_dir}
                                              )
      ->grid(-row => 4, -column => 2, -sticky => 'ew');
      $POG2013        = $frame_top -> Checkbutton(-text     => 'Update POG2013 database',
                                                  -variable => \$state{'update_POG2013'},
                                                  -command  => sub {
                                                                    if ($state{'update_POG2013'} == 1) {
                                                                       $POG2013_server->configure(-state      => 'normal',
                                                                                                  -foreground => 'black'
                                                                                                  );
                                                                       $POG2013_path  ->configure(-state      => 'normal',
                                                                                                  -foreground => 'black'
                                                                                                  );
                                                                    } elsif ($state{'update_POG2013'} == 0) {
                                                                       $POG2013_server->configure(-state      => 'disabled',
                                                                                                  -foreground => 'grey'
                                                                                                  );
                                                                       $POG2013_path  ->configure(-state      => 'disabled',
                                                                                                  -foreground => 'grey'
                                                                                                  );
                                                                    }
                                                                 }
                                                )
      ->grid(-row => 4, -column => 0, -sticky => 'w');
   }

   #PFam
   {
      $pfam_server   = $frame_top -> LabEntry(-label        => 'Server:',
                                                 -labelPack    => [-side => "left", -anchor => "w"],
                                                 -justify      => 'left',
                                                 -width        => 20,
                                                 -textvariable => \${$args{ini_ref}}{pfam_server_URL}
                                                 )
      ->grid(-row => 5, -column => 1, -sticky => 'ew');
      $pfam_path   = $frame_top -> LabEntry(-label        => 'Path:',
                                               -labelPack    => [-side => "left", -anchor => "w"],
                                               -justify      => 'left',
                                               -width        => 20,
                                               -textvariable => \${$args{ini_ref}}{pfam_remote_dir}
                                              )
      ->grid(-row => 5, -column => 2, -sticky => 'ew');
      $pfam        = $frame_top -> Checkbutton(-text     => 'Update Pfam databases',
                                                  -variable => \$state{'update_pfam'},
                                                  -command  => sub {
                                                                    if ($state{'update_pfam'} == 1) {
                                                                       $pfam_server->configure(-state      => 'normal',
                                                                                               -foreground => 'black'
                                                                                              );
                                                                       $pfam_path  ->configure(-state      => 'normal',
                                                                                               -foreground => 'black'
                                                                                              );
                                                                    } elsif ($state{'update_pfam'} == 0) {
                                                                       $pfam_server->configure(-state      => 'disabled',
                                                                                               -foreground => 'grey'
                                                                                              );
                                                                       $pfam_path  ->configure(-state      => 'disabled',
                                                                                               -foreground => 'grey'
                                                                                              );
                                                                    }
                                                                  }
                                                 )
      ->grid(-row => 5, -column => 0, -sticky => 'w');

      $pfam_descr_path   = $frame_top -> LabEntry(-label        => 'Path:',
                                                     -labelPack    => [-side => "left", -anchor => "w"],
                                                     -justify      => 'left',
                                                     -width        => 20,
                                                     -textvariable => \${$args{ini_ref}}{pfam_descriptor_dir}
                                                    )
      ->grid(-row => 6, -column => 2, -sticky => 'ew');
      $pfam_descr        = $frame_top -> Checkbutton(-text     => 'Update Pfam descriptor',
                                                        -variable => \$state{'update_pfam_descr'},
                                                        -command  => sub {
                                                                           if ($state{'update_pfam_descr'} == 1) {
                                                                              $pfam_descr_path->configure(-state      => 'normal',
                                                                                                          -foreground => 'black'
                                                                                                          );
                                                                           } elsif ($state{'update_pfam_descr'} == 0) {
                                                                              $pfam_descr_path->configure(-state      => 'disabled',
                                                                                                          -foreground => 'grey'
                                                                                                         );
                                                                           }
                                                                     }
                                                      )
      ->grid(-row => 6, -column => 0, -sticky => 'w');
      $pfam_ipro_path   = $frame_top -> LabEntry(-label        => 'Path:',
                                                    -labelPack    => [-side => "left", -anchor => "w"],
                                                    -justify      => 'left',
                                                    -width        => 20,
                                                    -textvariable => \${$args{ini_ref}}{pfam_remote_dir}
                                                    )
      ->grid(-row => 7, -column => 2, -sticky => 'ew');
      $pfam_ipro        = $frame_top -> Checkbutton(-text     => 'Update Interpro descriptor',
                                                       -variable => \$state{'update_interpro'},
                                                       -command  => sub {
                                                                         if ($state{'update_interpro'} == 1) {
                                                                            $pfam_ipro_path->configure(-state      => 'normal',
                                                                                                       -foreground => 'black'
                                                                                                       );
                                                                         } elsif ($state{'update_interpro'} == 0) {
                                                                            $pfam_ipro_path->configure(-state      => 'disabled',
                                                                                                       -foreground => 'grey'
                                                                                                       );
                                                                         }
                                                                     }
                                                      )
      ->grid(-row => 7, -column => 0, -sticky => 'w');
   }

   #TIGRFam
   {
      $tigrfam_server   = $frame_top -> LabEntry(-label        => 'Server:',
                                                    -labelPack    => [-side => "left", -anchor => "w"],
                                                    -justify      => 'left',
                                                    -width        => 20,
                                                    -textvariable => \${$args{ini_ref}}{tigrfam_server_URL}
                                                    )
      ->grid(-row => 8, -column => 1, -sticky => 'ew');
      $tigrfam_path   = $frame_top -> LabEntry(-label        => 'Path:',
                                                  -labelPack    => [-side => "left", -anchor => "w"],
                                                  -justify      => 'left',
                                                  -width        => 20,
                                                  -textvariable => \${$args{ini_ref}}{tigrfam_remote_dir}
                                                  )
      ->grid(-row => 8, -column => 2, -sticky => 'ew');
      $tigrfam        = $frame_top -> Checkbutton(-text     => 'Update TIGRfam databases',
                                                     -variable => \$state{'update_tigrfam'},
                                                     -command  => sub {
                                                                        if ($state{'update_tigrfam'} == 1) {
                                                                           $tigrfam_server->configure(-state      => 'normal',
                                                                                                      -foreground => 'black'
                                                                                                      );
                                                                           $tigrfam_path  ->configure(-state      => 'normal',
                                                                                                      -foreground => 'black'
                                                                                                      );
                                                                        } elsif ($state{'update_tigrfam'} == 0) {
                                                                           $tigrfam_server->configure(-state      => 'disabled',
                                                                                                      -foreground => 'grey'
                                                                                                      );
                                                                           $tigrfam_path  ->configure(-state      => 'disabled',
                                                                                                      -foreground => 'grey'
                                                                                                      );
                                                                        }
                                                                     }
                                                   )
      ->grid(-row => 8, -column => 0, -sticky => 'w');

      $go_server   = $frame_top -> LabEntry(-label        => 'Server:',
                                               -labelPack    => [-side => "left", -anchor => "w"],
                                               -justify      => 'left',
                                               -width        => 20,
                                               -textvariable => \${$args{ini_ref}}{go_server_URL}
                                               )
      ->grid(-row => 9, -column => 1, -sticky => 'ew');
      $go_path   = $frame_top -> LabEntry(-label        => 'Path:',
                                             -labelPack    => [-side => "left", -anchor => "w"],
                                             -justify      => 'left',
                                             -width        => 20,
                                             -textvariable => \${$args{ini_ref}}{go_remote_dir}
                                             )
      ->grid(-row => 9, -column => 2, -sticky => 'ew');
      $go        = $frame_top -> Checkbutton(-text     => 'Update GO descriptor',
                                                -variable => \$state{'update_go'},
                                                -command  => sub {
                                                                  if ($state{'update_go'} == 1) {
                                                                     $go_server->configure(-state      => 'normal',
                                                                                           -foreground => 'black'
                                                                                          );
                                                                     $go_path  ->configure(-state      => 'normal',
                                                                                           -foreground => 'black'
                                                                                          );
                                                                  } elsif ($state{'update_go'} == 0) {
                                                                     $go_server->configure(-state      => 'disabled',
                                                                                           -foreground => 'grey'
                                                                                          );
                                                                     $go_path  ->configure(-state      => 'disabled',
                                                                                           -foreground => 'grey'
                                                                                          );
                                                                  }
                                                               }
                                               )
      ->grid(-row => 9, -column => 0, -sticky => 'w');
   }

   $frame_top->Label(-text => 'Update structural databases')->grid(-row => 10, -column => 0, -sticky => 'nsew', -columnspan => 3);
   $frame_top->gridRowconfigure(10, -pad => 30);

   #non-coding RNA
   {
      $ncRNA_server   = $frame_top -> LabEntry(-label        => 'Server:',
                                                  -labelPack    => [-side => "left", -anchor => "w"],
                                                  -justify      => 'left',
                                                  -width        => 20,
                                                  -textvariable => \${$args{ini_ref}}{rfam_server_URL}
                                                  )
      ->grid(-row => 11, -column => 1, -sticky => 'ew');
      $ncRNA_path   = $frame_top -> LabEntry(-label        => 'Path:',
                                                -labelPack    => [-side => "left", -anchor => "w"],
                                                -justify      => 'left',
                                                -width        => 20,
                                                -textvariable => \${$args{ini_ref}}{rfam_remote_dir}
                                                )
      ->grid(-row => 11, -column => 2, -sticky => 'ew');
      $ncRNA        = $frame_top -> Checkbutton(-text     => 'Update Rfam databases',
                                                   -variable => \$state{'update_ncRNA'},
                                                   -command  => sub {
                                                                      if ($state{'update_ncRNA'} == 1) {
                                                                         $ncRNA_server->configure(-state      => 'normal',
                                                                                                  -foreground => 'black'
                                                                                                  );
                                                                         $ncRNA_path  ->configure(-state      => 'normal',
                                                                                                  -foreground => 'black'
                                                                                                  );
                                                                      } elsif ($state{'update_ncRNA'} == 0) {
                                                                         $ncRNA_server->configure(-state      => 'disabled',
                                                                                                  -foreground => 'grey'
                                                                                                  );
                                                                         $ncRNA_path  ->configure(-state      => 'disabled',
                                                                                                  -foreground => 'grey'
                                                                                                  );
                                                                      }
                                                                  }
                                                )
      ->grid(-row => 11, -column => 0, -sticky => 'w');
   }

   $frame_top->gridColumnconfigure(2, -weight => 1);
   $frame_top->gridColumnconfigure(1, -weight => 1);

   #initiate state
   &set_state(\%state);

   #Buttons to run
   {
      $frame_bottom -> Label(-text     => '  ')->grid(-row => 0, -column => 0);
      $frame_bottom -> Button(-text    => 'Select All',
                              -command => sub {
                                                #set state hash
                                                foreach my $key (keys %state) {
                                                   $state{$key} = 1;
                                                }
                                                &set_state(\%state);
                                              }
                            )
      ->grid(-row => 1, -column => 0, -sticky => 'w');

      $frame_bottom -> Button(-text    => 'De-select All',
                              -command => sub {
                                                #set state hash
                                                foreach my $key (keys %state) {
                                                   $state{$key} = 0;
                                                }
                                                &set_state(\%state);
                                              }
                            )
      ->grid(-row => 1, -column => 1, -sticky => 'w');

      $frame_bottom -> Button(-text    => 'Update databases',
                              -command => sub {
                                                my ($status);
                                                #test where updates are required
                                                $status = &test_for_updates(main_window     => $args{main_window},
                                                                            progress_bar    => \$progress_bar,
                                                                            auto_ini_ref    => $args{auto_ini_ref},
                                                                            ini_ref         => $args{ini_ref},
                                                                            state_ref       => \%state,
                                                                            update_ref      => \%update
                                                                           );
                                                if ($status == 0) {
                                                   ${$args{main_window}}->messageBox(-title   => 'Error',
                                                                                     -message => "Error testing for available updates",
                                                                                     -icon    => 'error',
                                                                                     -type    => 'ok');
                                                   #reset update hash
                                                   undef %update;
                                                }
                                                #download new db files
                                                $status = &update(main_window     => $args{main_window},
                                                                  progress_bar    => \$progress_bar,
                                                                  auto_ini_ref    => $args{auto_ini_ref},
                                                                  ini_ref         => $args{ini_ref},
                                                                  update_ref      => \%update
                                                                 );
                                                if ($status == 0) {
                                                   ${$args{main_window}}->messageBox(-title   => 'Error',
                                                                                     -message => "Error updating selected databases",
                                                                                     -icon    => 'error',
                                                                                     -type    => 'ok');
                                                   #reset update hash
                                                   undef %update;
                                                }
                                                #unpack downloaded files
                                                $status = &extract_db(main_window     => $args{main_window},
                                                                      progress_bar    => \$progress_bar,
                                                                      auto_ini_ref    => $args{auto_ini_ref},
                                                                      ini_ref         => $args{ini_ref},
                                                                      update_ref      => \%update
                                                                     );
                                                if ($status == 0) {
                                                   ${$args{main_window}}->messageBox(-title   => 'Error',
                                                                                     -message => "Error extracting downloaded databases",
                                                                                     -icon    => 'error',
                                                                                     -type    => 'ok');
                                                   #reset update hash
                                                   undef %update;
                                                } elsif ($status == 1) {
                                                   ${$args{main_window}}->messageBox(-title   => 'Success',
                                                                                     -message => "All eselcted databases have been successfully updated.",
                                                                                     -icon    => 'info',
                                                                                     -type    => 'ok');
                                                   #reset update hash
                                                   undef %update;
                                                }
                                              }
                              )
      ->grid(-row => 1, -column => 2, -sticky => 'w');

      $frame_bottom -> Button(-text    => 'Force database update',
                              -command => sub {
                                                my $status;
                                                #set update on selected dbs
                                                foreach my $key (keys %state) {
                                                   if ($state{$key} == 1) {
                                                      $update{$key} = 1;
                                                   }
                                                }

                                                $status = &update(main_window     => $args{main_window},
                                                                  progress_bar    => \$progress_bar,
                                                                  auto_ini_ref    => $args{auto_ini_ref},
                                                                  ini_ref         => $args{ini_ref},
                                                                  update_ref      => \%update
                                                                 );
                                                if ($status == 0) {
                                                   ${$args{main_window}}->messageBox(-title   => 'Error',
                                                                                     -message => "Error updating selected databases",
                                                                                     -icon    => 'error',
                                                                                     -type    => 'ok');
                                                   #reset update hash
                                                   undef %update;
                                                }
                                                #unpack downloaded files
                                                $status = &extract_db(main_window     => $args{main_window},
                                                                      progress_bar    => \$progress_bar,
                                                                      auto_ini_ref    => $args{auto_ini_ref},
                                                                      ini_ref         => $args{ini_ref},
                                                                      update_ref      => \%update
                                                                     );
                                                if ($status == 0) {
                                                   ${$args{main_window}}->messageBox(-title   => 'Error',
                                                                                     -message => "Error extracting downloaded databases",
                                                                                     -icon    => 'error',
                                                                                     -type    => 'ok');
                                                   #reset update hash
                                                   undef %update;
                                                } elsif ($status == 1) {
                                                   ${$args{main_window}}->messageBox(-title   => 'Success',
                                                                                     -message => "All eselcted databases have been successfully updated.",
                                                                                     -icon    => 'info',
                                                                                     -type    => 'ok');
                                                   #reset update hash
                                                   undef %update;
                                                }
                                              }
                              )
      ->grid(-row => 1, -column => 3, -sticky => 'w');

      $frame_bottom -> Button(-text    => 'Cancel',
                              -command => sub {
                                                 $tl->state('withdrawn');
                                                 return;
                                              }
                          )
      ->grid(-row => 1, -column => 4, -sticky => 'w');
   }

   return;
}

sub set_state {
   my $state_ref = shift;

   #COG
   if (${$state_ref}{'update_COG'} == 0) {
      $COG_server->configure(-state      => 'disabled',
                             -foreground => 'grey'
                             );
      $COG_path  ->configure(-state      => 'disabled',
                             -foreground => 'grey'
                             );
   } elsif (${$state_ref}{'update_COG'} == 1) {
      $COG_server->configure(-state      => 'normal',
                             -foreground => 'black'
                             );
      $COG_path  ->configure(-state      => 'normal',
                             -foreground => 'black'
                             );
   }

   #COG2014
   if (${$state_ref}{'update_COG2014'} == 0) {
      $COG2014_server->configure(-state      => 'disabled',
                                 -foreground => 'grey'
                                 );
      $COG2014_path  ->configure(-state      => 'disabled',
                                 -foreground => 'grey'
                                 );
   } elsif (${$state_ref}{'update_COG2014'} == 1) {
      $COG2014_server->configure(-state      => 'normal',
                                 -foreground => 'black'
                                 );
      $COG2014_path  ->configure(-state      => 'normal',
                                 -foreground => 'black'
                                 );
   }

   #arCOG2014
   if (${$state_ref}{'update_arCOG2014'} == 0) {
      $arCOG2014_server->configure(-state      => 'disabled',
                                   -foreground => 'grey'
                                   );
      $arCOG2014_path  ->configure(-state      => 'disabled',
                                   -foreground => 'grey'
                                   );
   } elsif (${$state_ref}{'update_arCOG2014'} == 1) {
      $arCOG2014_server->configure(-state      => 'normal',
                                   -foreground => 'black'
                                   );
      $arCOG2014_path  ->configure(-state      => 'normal',
                                   -foreground => 'black'
                                   );
   }

   #POG2013
   if (${$state_ref}{'update_POG2013'} == 0) {
      $POG2013_server->configure(-state      => 'disabled',
                                 -foreground => 'grey'
                                 );
      $POG2013_path  ->configure(-state      => 'disabled',
                                 -foreground => 'grey'
                                 );
   } elsif (${$state_ref}{'update_POG2013'} == 1) {
      $POG2013_server->configure(-state      => 'normal',
                                 -foreground => 'black'
                                 );
      $POG2013_path  ->configure(-state      => 'normal',
                                 -foreground => 'black'
                                 );
   }

   #pfam
   if (${$state_ref}{'update_pfam'} == 1) {
      $pfam_server->configure(-state      => 'normal',
                             -foreground => 'black'
                            );
     $pfam_path  ->configure(-state      => 'normal',
                             -foreground => 'black'
                            );
   } elsif (${$state_ref}{'update_pfam'} == 0) {
      $pfam_server->configure(-state      => 'disabled',
                             -foreground => 'grey'
                            );
     $pfam_path  ->configure(-state      => 'disabled',
                             -foreground => 'grey'
                            );
   }

   #pfam descriptor
   if (${$state_ref}{'update_pfam_descr'} == 1) {
      $pfam_descr_path->configure(-state      => 'normal',
                                  -foreground => 'black'
                                  );
   } elsif (${$state_ref}{'update_pfam_descr'} == 0) {
      $pfam_descr_path->configure(-state      => 'disabled',
                                  -foreground => 'grey'
                                 );
   }

   #pfam interpro
   if (${$state_ref}{'update_interpro'} == 1) {
      $pfam_ipro_path->configure(-state      => 'normal',
                                 -foreground => 'black'
                                 );
   } elsif (${$state_ref}{'update_interpro'} == 0) {
      $pfam_ipro_path->configure(-state      => 'disabled',
                                 -foreground => 'grey'
                                );
   }

   #TIGRfam
   if (${$state_ref}{'update_tigrfam'} == 1) {
      $tigrfam_server->configure(-state      => 'normal',
                                 -foreground => 'black'
                                 );
      $tigrfam_path  ->configure(-state      => 'normal',
                                 -foreground => 'black'
                                 );
   } elsif (${$state_ref}{'update_tigrfam'} == 0) {
      $tigrfam_server->configure(-state      => 'disabled',
                                 -foreground => 'grey'
                                 );
      $tigrfam_path  ->configure(-state      => 'disabled',
                                 -foreground => 'grey'
                                 );
   }

   #GO
   if (${$state_ref}{'update_go'} == 1) {
      $go_server->configure(-state      => 'normal',
                            -foreground => 'black'
                           );
      $go_path  ->configure(-state      => 'normal',
                            -foreground => 'black'
                           );
   } elsif (${$state_ref}{'update_go'} == 0) {
      $go_server->configure(-state      => 'disabled',
                            -foreground => 'grey'
                           );
      $go_path  ->configure(-state      => 'disabled',
                            -foreground => 'grey'
                           );
   }

   #Rfam
   if (${$state_ref}{'update_ncRNA'} == 1) {
      $ncRNA_server->configure(-state      => 'normal',
                               -foreground => 'black'
                               );
      $ncRNA_path  ->configure(-state      => 'normal',
                               -foreground => 'black'
                              );
   } elsif (${$state_ref}{'update_ncRNA'} == 0) {
      $ncRNA_server->configure(-state      => 'disabled',
                               -foreground => 'grey'
                               );
      $ncRNA_path  ->configure(-state      => 'disabled',
                               -foreground => 'grey'
                               );
   }
}

sub test_for_updates {
   my %args = @_;
   my ($localtime, $ftp, $remotetime, $local_release, $remote_release, @remote_dir);

   #test COG
   if (${$args{state_ref}}{'update_COG'} == 1) {
      #test local database time
      ${$args{progress_bar}}->configure(-label=>"Retrieving local COG time stamp");
      ${$args{main_window}} ->update;
      if (-e ${$args{ini_ref}}{COG_db_path}.'/COG2003/COG.phr') {
         #$localtime = ctime(stat(${$args{ini_ref}}{COG_db_path}.'/COG.phr')->mtime);
         $localtime = (stat(${$args{ini_ref}}{COG_db_path}.'/COG2003/COG.phr'));
         $localtime = @{$localtime}[9];
      } elsif (-e ${$args{ini_ref}}{COG_db_path}.'/COG2003/whog') {
         #$localtime = ctime(stat(${$args{ini_ref}}{COG_db_path}.'/whog')->mtime);
         $localtime = (stat(${$args{ini_ref}}{COG_db_path}.'/COG2003/whog'));
         $localtime = @{$localtime}[9];
      } else {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot find database or source COG2003 files in directory ${$args{ini_ref}}{COG_db_path}\nSetting timestamp to 0",
                                           -icon    => 'error',
                                           -type    => 'ok');
         $localtime = 0;
      };

      #get remote time
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote COG2003 time stamp: connect to server");
      ${$args{main_window}} ->update;
      $ftp = Net::FTP->new(${$args{ini_ref}}{COG_server_URL}, Debug => 0) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot connect to ${$args{ini_ref}}{COG_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote COG2003 time stamp: login");
      ${$args{main_window}} ->update;
      $ftp->login("anonymous",'-anonymous@') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot login to ${$args{ini_ref}}{COG_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote COG2003 time stamp: change directory");
      ${$args{main_window}} ->update;
      $ftp->cwd('/'.${$args{ini_ref}}{COG_remote_dir}) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot change to directory ${$args{ini_ref}}{COG_remote_dir}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote COG2003 time stamp: retrieve timestamp");
      ${$args{main_window}} ->update;
      $remotetime = $ftp->mdtm('org.txt') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot determine file modification date: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };

      #test for time
      if ($localtime >= $remotetime) {
         my $force_update = ${$args{main_window}}->Dialog(-title   => 'Force update',
                                                          -text    => "Local COG2003 database appears to be more recent than remote.\nUpdate anyway?",
                                                          -bitmap  => 'question',
                                                          -buttons => ['Yes', 'No']
                                                         )->Show();
         if ($force_update eq 'Yes') {
            ${$args{update_ref}}{'update_COG'} = 1;
         } else {
            ${$args{update_ref}}{'update_COG'} = 0;
         }
      } else {
         ${$args{update_ref}}{'update_COG'} = 1;
      }
      $ftp->quit;
   }

   #test COG2014
   if (${$args{state_ref}}{'update_COG2014'} == 1) {
      #test local database time
      ${$args{progress_bar}}->configure(-label=>"Retrieving local COG2014 time stamp");
      ${$args{main_window}} ->update;
      if (-e ${$args{ini_ref}}{COG_db_path}.'/COG2014/COG2014.phr') {
         #$localtime = ctime(stat(${$args{ini_ref}}{COG_db_path}.'/COG.phr')->mtime);
         $localtime = (stat(${$args{ini_ref}}{COG_db_path}.'/COG2014/COG2014.phr'));
         $localtime = @{$localtime}[9];
      } elsif (-e ${$args{ini_ref}}{COG_db_path}.'/COG2014/cognames2003-2014.tab') {
         #$localtime = ctime(stat(${$args{ini_ref}}{COG_db_path}.'/whog')->mtime);
         $localtime = (stat(${$args{ini_ref}}{COG_db_path}.'/COG2014/cognames2003-2014.tab'));
         $localtime = @{$localtime}[9];
      } else {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot find database or source COG2014 files in directory ${$args{ini_ref}}{COG_db_path}\nSetting timestamp to 0",
                                           -icon    => 'error',
                                           -type    => 'ok');
         $localtime = 0;
      };

      #get remote time
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote COG2014 time stamp: connect to server");
      ${$args{main_window}} ->update;
      $ftp = Net::FTP->new(${$args{ini_ref}}{COG_server_URL}, Debug => 0) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot connect to ${$args{ini_ref}}{COG_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote COG2014 time stamp: login");
      ${$args{main_window}} ->update;
      $ftp->login("anonymous",'-anonymous@') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot login to ${$args{ini_ref}}{COG_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote COG2014 time stamp: change directory");
      ${$args{main_window}} ->update;
      $ftp->cwd('/'.${$args{ini_ref}}{COG2014_remote_dir}) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot change to directory ${$args{ini_ref}}{COG2014_remote_dir}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote COG2014 time stamp: retrieve timestamp");
      ${$args{main_window}} ->update;
      $remotetime = $ftp->mdtm('cognames2003-2014.tab') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot determine file modification date: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };

      #test for time
      if ($localtime >= $remotetime) {
         my $force_update = ${$args{main_window}}->Dialog(-title   => 'Force update',
                                                          -text    => "Local COG2014 database appears to be more recent than remote.\nUpdate anyway?",
                                                          -bitmap  => 'question',
                                                          -buttons => ['Yes', 'No']
                                                         )->Show();
         if ($force_update eq 'Yes') {
            ${$args{update_ref}}{'update_COG2014'} = 1;
         } else {
            ${$args{update_ref}}{'update_COG2014'} = 0;
         }
      } else {
         ${$args{update_ref}}{'update_COG2014'} = 1;
      }
      $ftp->quit;
   }

   #test arCOG2014
   if (${$args{state_ref}}{'update_arCOG2014'} == 1) {
      #test local database time
      ${$args{progress_bar}}->configure(-label=>"Retrieving local arCOG2014 time stamp");
      ${$args{main_window}} ->update;
      if (-e ${$args{ini_ref}}{COG_db_path}.'/arCOG2014/arCOG2014.phr') {
         #$localtime = ctime(stat(${$args{ini_ref}}{COG_db_path}.'/COG.phr')->mtime);
         $localtime = (stat(${$args{ini_ref}}{COG_db_path}.'/arCOG2014/arCOG2014.phr'));
         $localtime = @{$localtime}[9];
      } elsif (-e ${$args{ini_ref}}{COG_db_path}.'/arCOG2014/ar14.arCOGdef.tab') {
         #$localtime = ctime(stat(${$args{ini_ref}}{COG_db_path}.'/whog')->mtime);
         $localtime = (stat(${$args{ini_ref}}{COG_db_path}.'/arCOG2014/ar14.arCOGdef.tab'));
         $localtime = @{$localtime}[9];
      } else {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot find database or source arCOG2014 files in directory ${$args{ini_ref}}{COG_db_path}\nSetting timestamp to 0",
                                           -icon    => 'error',
                                           -type    => 'ok');
         $localtime = 0;
      };

      #get remote time
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote arCOG2014 time stamp: connect to server");
      ${$args{main_window}} ->update;
      $ftp = Net::FTP->new(${$args{ini_ref}}{COG_server_URL}, Debug => 0) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot connect to ${$args{ini_ref}}{COG_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote arCOG2014 time stamp: login");
      ${$args{main_window}} ->update;
      $ftp->login("anonymous",'-anonymous@') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot login to ${$args{ini_ref}}{COG_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote arCOG2014 time stamp: change directory");
      ${$args{main_window}} ->update;
      $ftp->cwd('/'.${$args{ini_ref}}{arCOG2014_remote_dir}) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot change to directory ${$args{ini_ref}}{arCOG2014_remote_dir}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote arCOG2014 time stamp: retrieve timestamp");
      ${$args{main_window}} ->update;
      $remotetime = $ftp->mdtm('ar14.arCOGdef.tab') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot determine file modification date: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };

      #test for time
      if ($localtime >= $remotetime) {
         my $force_update = ${$args{main_window}}->Dialog(-title   => 'Force update',
                                                          -text    => "Local arCOG2014 database appears to be more recent than remote.\nUpdate anyway?",
                                                          -bitmap  => 'question',
                                                          -buttons => ['Yes', 'No']
                                                         )->Show();
         if ($force_update eq 'Yes') {
            ${$args{update_ref}}{'update_arCOG2014'} = 1;
         } else {
            ${$args{update_ref}}{'update_arCOG2014'} = 0;
         }
      } else {
         ${$args{update_ref}}{'update_arCOG2014'} = 1;
      }
      $ftp->quit;
   }

   #test POG2013
   if (${$args{state_ref}}{'update_POG2013'} == 1) {
      #test local database time
      ${$args{progress_bar}}->configure(-label=>"Retrieving local POG2013 time stamp");
      ${$args{main_window}} ->update;
      if (-e ${$args{ini_ref}}{COG_db_path}.'/POG2013/POGseqs_HighVQ.phr') {
         #$localtime = ctime(stat(${$args{ini_ref}}{COG_db_path}.'/COG.phr')->mtime);
         $localtime = (stat(${$args{ini_ref}}{COG_db_path}.'/POG2013/POGseqs_HighVQ.phr'));
         $localtime = @{$localtime}[9];
      } elsif (-e ${$args{ini_ref}}{COG_db_path}.'/POG2013/pogs.txt') {
         #$localtime = ctime(stat(${$args{ini_ref}}{COG_db_path}.'/whog')->mtime);
         $localtime = (stat(${$args{ini_ref}}{COG_db_path}.'/POG2013/pogs.txt'));
         $localtime = @{$localtime}[9];
      } else {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot find database or source POG2013 files in directory ${$args{ini_ref}}{COG_db_path}\nSetting timestamp to 0",
                                           -icon    => 'error',
                                           -type    => 'ok');
         $localtime = 0;
      };

      #get remote time
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote POG2013 time stamp: connect to server");
      ${$args{main_window}} ->update;
      $ftp = Net::FTP->new(${$args{ini_ref}}{COG_server_URL}, Debug => 0) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot connect to ${$args{ini_ref}}{COG_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote POG2013 time stamp: login");
      ${$args{main_window}} ->update;
      $ftp->login("anonymous",'-anonymous@') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot login to ${$args{ini_ref}}{COG_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote POG2013 time stamp: change directory");
      ${$args{main_window}} ->update;
      $ftp->cwd('/'.${$args{ini_ref}}{POG2013_remote_dir}) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot change to directory ${$args{ini_ref}}{POG2013_remote_dir}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote POG2013 time stamp: retrieve timestamp");
      ${$args{main_window}} ->update;
      $remotetime = $ftp->mdtm('pogs.txt') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot determine file modification date: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };

      #test for time
      if ($localtime >= $remotetime) {
         my $force_update = ${$args{main_window}}->Dialog(-title   => 'Force update',
                                                          -text    => "Local POG2013 database appears to be more recent than remote.\nUpdate anyway?",
                                                          -bitmap  => 'question',
                                                          -buttons => ['Yes', 'No']
                                                         )->Show();
         if ($force_update eq 'Yes') {
            ${$args{update_ref}}{'update_POG2013'} = 1;
         } else {
            ${$args{update_ref}}{'update_POG2013'} = 0;
         }
      } else {
         ${$args{update_ref}}{'update_POG2013'} = 1;
      }
      $ftp->quit;
   }

   #test Pfam
   if (${$args{state_ref}}{'update_pfam'} == 1) {
      #test local database time
      ${$args{progress_bar}}->configure(-label=>"Retrieving Pfam release version");
      ${$args{main_window}} ->update;
      if (-e ${$args{ini_ref}}{pfam_db_path}.'/relnotes.txt') {
         #read file into memory
         my ($file_ref) = &slurp(main_window   => $args{main_window},
                                 auto_ini_ref  => $args{auto_ini_ref},
                                 ini_ref       => $args{ini_ref},
                                 filename      => 'relnotes.txt',
                                 directory     => ${$args{ini_ref}}{pfam_db_path}
                                );
         'reset' =~ m/reset/;
         ${$file_ref} =~ m/\s+RELEASE (\S+)/s;
         if ($1) {$local_release = $1} else {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot find Pfam release version in release note\.\nSetting release to 0",
                                              -icon    => 'error',
                                              -type    => 'ok');
            $local_release = 0;
         };
         undef $file_ref;
      } else {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot find Pfam release note in directory ${$args{ini_ref}}{pfam_db_path}\nSetting release to 0",
                                           -icon    => 'error',
                                           -type    => 'ok');
         $local_release = 0;
      };

      #get remote release version
      unlink ${$args{auto_ini_ref}}{work_dir}.'/Pfam_rel_note';
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Pfam release version: connect to server");
      ${$args{main_window}} ->update;
      $ftp = Net::FTP->new(${$args{ini_ref}}{pfam_server_URL}, Debug => 0) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot connect to ${$args{ini_ref}}{pfam_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Pfam release version: login");
      ${$args{main_window}} ->update;
      $ftp->login("anonymous",'-anonymous@') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot login to ${$args{ini_ref}}{pfam_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Pfam release version: change directory");
      ${$args{main_window}} ->update;
      $ftp->cwd('/'.${$args{ini_ref}}{pfam_remote_dir}) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot change to directory ${$args{ini_ref}}{pfam_remote_dir}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Pfam release version: retrieve timestamp");
      ${$args{main_window}} ->update;
      $ftp->binary;
      $ftp->get('relnotes.txt', ${$args{auto_ini_ref}}{work_dir}.'/Pfam_rel_note') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot retrieve release note: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      #read file into memory
      my ($file_ref) = &slurp(main_window   => $args{main_window},
                              auto_ini_ref  => $args{auto_ini_ref},
                              ini_ref       => $args{ini_ref},
                              filename      => 'Pfam_rel_note',
                              directory     => ${$args{auto_ini_ref}}{work_dir}
                             );
      'reset' =~ m/reset/;
      ${$file_ref} =~ m/\s+RELEASE (\S+)/s;
      if ($1) {$remote_release = $1} else {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot find Pfam release note in remote release note\.\nSetting release to 0",
                                           -icon    => 'error',
                                           -type    => 'ok');
         $remote_release = 0;
      };
      undef $file_ref;
      unlink ${$args{auto_ini_ref}}{work_dir}.'/Pfam_rel_note';

      #test for version
      if ($local_release >= $remote_release) {
         my $force_update = ${$args{main_window}}->Dialog(-title   => 'Force update',
                                                          -text    => "Local PFam database \(version $local_release\) appears to be more recent than remote \(version $remote_release\).\nUpdate anyway?",
                                                          -bitmap  => 'question',
                                                          -buttons => ['Yes', 'No']
                                                         )->Show();
         if ($force_update eq 'Yes') {
            ${$args{update_ref}}{'update_pfam'} = 1;
         } else {
            ${$args{update_ref}}{'update_pfam'} = 0;
         }
      } else {
         ${$args{update_ref}}{'update_pfam'} = 1;
      }
      $ftp->quit;
   }

   #test Pfam descriptor
   if (${$args{state_ref}}{'update_pfam_descr'} == 1) {
      #test local database time
      ${$args{progress_bar}}->configure(-label=>"Retrieving local Pfam descriptor time stamp");
      ${$args{main_window}} ->update;
      if (-e ${$args{ini_ref}}{Pfam_descriptor}.'/pfamA.txt') {
         $localtime = (stat(${$args{ini_ref}}{Pfam_descriptor}.'/pfamA.txt'));
         $localtime = @{$localtime}[9];
      } elsif (-e ${$args{ini_ref}}{Pfam_descriptor}.'/pfamA.txt.gz') {
         $localtime = (stat(${$args{ini_ref}}{Pfam_descriptor}.'/pfamA.txt.gz'));
         $localtime = @{$localtime}[9];
      } else {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot find Pfam descriptor in directory ${$args{ini_ref}}{Pfam_descriptor}\nSetting timestamp to 0",
                                           -icon    => 'error',
                                           -type    => 'ok');
         $localtime = 0;
      };

      #get remote time
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Pfam descriptor time stamp: connect to server");
      ${$args{main_window}} ->update;
      $ftp = Net::FTP->new(${$args{ini_ref}}{pfam_server_URL}, Debug => 0) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot connect to ${$args{ini_ref}}{pfam_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Pfam descriptor time stamp: login");
      ${$args{main_window}} ->update;
      $ftp->login("anonymous",'-anonymous@') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot login to ${$args{ini_ref}}{pfam_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Pfam descriptor time stamp: change directory");
      ${$args{main_window}} ->update;
      $ftp->cwd('/'.${$args{ini_ref}}{pfam_descriptor_dir}) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot change to directory ${$args{ini_ref}}{pfam_descriptor_dir}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Pfam descriptor time stamp: retrieve timestamp");
      ${$args{main_window}} ->update;
      $remotetime = $ftp->mdtm('pfamA.txt.gz') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot determine file modification date: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };

      #test for time
      if ($localtime >= $remotetime) {
         my $force_update = ${$args{main_window}}->Dialog(-title   => 'Force update',
                                                          -text    => "Local Pfam descriptor appears to be more recent than remote.\nUpdate anyway?",
                                                          -bitmap  => 'question',
                                                          -buttons => ['Yes', 'No']
                                                         )->Show();
         if ($force_update eq 'Yes') {
            ${$args{update_ref}}{'update_pfam_descr'} = 1;
         } else {
            ${$args{update_ref}}{'update_pfam_descr'} = 0;
         }
      } else {
         ${$args{update_ref}}{'update_pfam_descr'} = 1;
      }
      $ftp->quit;
   }

   #test Pfam-InterPro descriptor
   if (${$args{state_ref}}{'update_interpro'} == 1) {
      #test local database time
      ${$args{progress_bar}}->configure(-label=>"Retrieving local InterPro descriptor time stamp");
      ${$args{main_window}} ->update;
      if (-e ${$args{ini_ref}}{Interpro_descriptor}.'/interpro.txt') {
         $localtime = (stat(${$args{ini_ref}}{Interpro_descriptor}.'/interpro.txt'));
         $localtime = @{$localtime}[9];
      } elsif (-e ${$args{ini_ref}}{Interpro_descriptor}.'/interpro.txt.gz') {
         $localtime = (stat(${$args{ini_ref}}{Interpro_descriptor}.'/interpro.txt.gz'));
         $localtime = @{$localtime}[9];
      } else {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot find InterPro descriptor in directory ${$args{ini_ref}}{Interpro_descriptor}\nSetting timestamp to 0",
                                           -icon    => 'error',
                                           -type    => 'ok');
         $localtime = 0;
      };

      #get remote time
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote InterPro descriptor time stamp: connect to server");
      ${$args{main_window}} ->update;
      $ftp = Net::FTP->new(${$args{ini_ref}}{pfam_server_URL}, Debug => 0) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot connect to ${$args{ini_ref}}{pfam_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote InterPro descriptor time stamp: login");
      ${$args{main_window}} ->update;
      $ftp->login("anonymous",'-anonymous@') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot login to ${$args{ini_ref}}{pfam_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote InterPro descriptor time stamp: change directory");
      ${$args{main_window}} ->update;
      $ftp->cwd('/'.${$args{ini_ref}}{pfam_interpro_descr_dir}) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot change to directory ${$args{ini_ref}}{pfam_interpro_descr_dir}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote InterPro descriptor time stamp: retrieve timestamp");
      ${$args{main_window}} ->update;
      $remotetime = $ftp->mdtm('interpro.txt.gz') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot determine file modification date: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };

      #test for time
      if ($localtime >= $remotetime) {
         my $force_update = ${$args{main_window}}->Dialog(-title   => 'Force update',
                                                          -text    => "Local InterPro descriptor appears to be more recent than remote.\nUpdate anyway?",
                                                          -bitmap  => 'question',
                                                          -buttons => ['Yes', 'No']
                                                         )->Show();
         if ($force_update eq 'Yes') {
            ${$args{update_ref}}{'update_interpro'} = 1;
         } else {
            ${$args{update_ref}}{'update_interpro'} = 0;
         }
      } else {
         ${$args{update_ref}}{'update_interpro'} = 1;
      }
      $ftp->quit;
   }

   #test TIGRfam
   if (${$args{state_ref}}{'update_tigrfam'} == 1) {
      #test local release version
      ${$args{progress_bar}}->configure(-label=>"Retrieving TIGRfam release version");
      ${$args{main_window}} ->update;
      #get release filenames
      opendir SEQINPUT, ${$args{ini_ref}}{TIGRfam_db_path};
      my @rel_version = grep /^TIGRFAMs_\d+/, readdir(SEQINPUT);
      closedir SEQINPUT;
      #determine highest release number
      $local_release = 0;
      foreach my $release (@rel_version) {
         'reset' =~ m/reset/;
         $release =~ m/^TIGRFAMs_(\d+\.?\d*)_/;
         if ($1 && $1 > $local_release) {$local_release = $1;}
         elsif ($1 && $1 <= $local_release) {next;}
         else {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot find local TIGRfam release version\.",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return (0);
         }
      }

      #get remote release version
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote TIGRfam release version: connect to server");
      ${$args{main_window}} ->update;
      $ftp = Net::FTP->new(${$args{ini_ref}}{tigrfam_server_URL}, Debug => 0) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot connect to ${$args{ini_ref}}{tigrfam_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote TIGRfam release version: login");
      ${$args{main_window}} ->update;
      $ftp->login("anonymous",'-anonymous@') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot login to ${$args{ini_ref}}{tigrfam_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote TIGRfam release version: change directory");
      ${$args{main_window}} ->update;
      $ftp->cwd('/'.${$args{ini_ref}}{tigrfam_remote_dir}) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot change to directory ${$args{ini_ref}}{tigrfam_remote_dir}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote TIGRfam release version: retrieve version");
      ${$args{main_window}} ->update;
      @remote_dir = $ftp->ls or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot retrieve remote directory: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      #get remote version number
      foreach my $remote (@remote_dir) {
         if ($remote =~ m/TIGRFAMs_\d/) {
            'reset' =~ m/reset/;
            $remote =~ m/TIGRFAMs_(\d+\.?\d*)_/;
            if ($1) {$remote_release = $1; last;} else {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Cannot find remote TIGRfam release version\.\nSetting release to 0",
                                                 -icon    => 'error',
                                                 -type    => 'ok');
               $remote_release = 0;
               last;
            };
         }
      }


      #test for version
      if ($local_release >= $remote_release) {
         my $force_update = ${$args{main_window}}->Dialog(-title   => 'Force update',
                                                          -text    => "Local TIGRfam database appears to be more recent than remote.\nUpdate anyway?",
                                                          -bitmap  => 'question',
                                                          -buttons => ['Yes', 'No']
                                                         )->Show();
         if ($force_update eq 'Yes') {
            ${$args{update_ref}}{'update_tigrfam'} = 1;
         } else {
            ${$args{update_ref}}{'update_tigrfam'} = 0;
         }
      } else {
         ${$args{update_ref}}{'update_tigrfam'} = 1;
      }
      $ftp->quit;
   }

   #test Gene Ontology
   if (${$args{state_ref}}{'update_go'} == 1) {
      #test local release version
      ${$args{progress_bar}}->configure(-label=>"Retrieving Gene Ontology release version");
      ${$args{main_window}} ->update;
      #get release filenames
      opendir SEQINPUT, ${$args{ini_ref}}{GO_db_path};
      my @rel_version = grep /^go_\d+\-termdb\.obo/, readdir(SEQINPUT);
      closedir SEQINPUT;
      #determine highest release number
      $local_release = 0;
      foreach my $release (@rel_version) {
         'reset' =~ m/reset/;
         $release =~ m/^go_(\d+)/;
         if ($1 && $1 > $local_release) {$local_release = $1;}
         elsif ($1 && $1 <= $local_release) {next;}
         else {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot find local Gene Ontology release version\.",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return (0);
         }
      }

      #get remote release version
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Gene Ontology release version: connect to server");
      ${$args{main_window}} ->update;
      $ftp = Net::FTP->new(${$args{ini_ref}}{go_server_URL}, Debug => 0) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot connect to ${$args{ini_ref}}{go_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Gene Ontology release version: login");
      ${$args{main_window}} ->update;
      $ftp->login("anonymous",'-anonymous@') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot login to ${$args{ini_ref}}{go_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Gene Ontology release version: change directory");
      ${$args{main_window}} ->update;
      $ftp->cwd('/'.${$args{ini_ref}}{go_remote_dir}) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot change to directory ${$args{ini_ref}}{tigrfam_remote_dir}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Gene Ontology release version: retrieve version");
      ${$args{main_window}} ->update;
      @remote_dir = $ftp->ls or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot retrieve remote directory: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      #get remote version number
      foreach my $remote (@remote_dir) {
         if ($remote =~ m/^go_\d+/) {
            'reset' =~ m/reset/;
            $remote =~ m/^go_(\d+)\-/;
            if ($1) {$remote_release = $1; last;} else {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Cannot find remote Gene Ontology release version\.\nSetting release to 0",
                                                 -icon    => 'error',
                                                 -type    => 'ok');
               $remote_release = 0;
               last;
            };
         }
      }


      #test for version
      if ($local_release >= $remote_release) {
         my $force_update = ${$args{main_window}}->Dialog(-title   => 'Force update',
                                                          -text    => "Local Gene Ontology database appears to be more recent than remote.\nUpdate anyway?",
                                                          -bitmap  => 'question',
                                                          -buttons => ['Yes', 'No']
                                                         )->Show();
         if ($force_update eq 'Yes') {
            ${$args{update_ref}}{'update_go'} = 1;
         } else {
            ${$args{update_ref}}{'update_go'} = 0;
         }
      } else {
         ${$args{update_ref}}{'update_go'} = 1;
      }
      $ftp->quit;
   }

   #test Rfam
   if (${$args{state_ref}}{'update_ncRNA'} == 1) {
      #test local database time
      ${$args{progress_bar}}->configure(-label=>"Retrieving Rfam release version");
      ${$args{main_window}} ->update;
      if (-e ${$args{ini_ref}}{rfam_db_path}.'/README') {
         #read file into memory
         my ($file_ref) = &slurp(main_window   => $args{main_window},
                                 auto_ini_ref  => $args{auto_ini_ref},
                                 ini_ref       => $args{ini_ref},
                                 filename      => 'README',
                                 directory     => ${$args{ini_ref}}{rfam_db_path}
                                );
         'reset' =~ m/reset/;
         ${$file_ref} =~ m/\sRelease\s+(\S+)\s/s;
         if ($1) {$local_release = $1} else {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot find Rfam release version in release note\.\nSetting release to 0",
                                              -icon    => 'error',
                                              -type    => 'ok');
            $local_release = 0;
         };
         undef $file_ref;
      } else {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot find Rfam release note in directory ${$args{ini_ref}}{pfam_db_path}\nSetting release to 0",
                                           -icon    => 'error',
                                           -type    => 'ok');
         $local_release = 0;
      };

      #get remote release version
      unlink ${$args{auto_ini_ref}}{work_dir}.'/Rfam_rel_note';
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Rfam release version: connect to server");
      ${$args{main_window}} ->update;
      $ftp = Net::FTP->new(${$args{ini_ref}}{rfam_server_URL}, Debug => 0) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot connect to ${$args{ini_ref}}{rfam_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Rfam release version: login");
      ${$args{main_window}} ->update;
      $ftp->login("anonymous",'-anonymous@') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot login to ${$args{ini_ref}}{rfam_server_URL}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Rfam release version: change directory");
      ${$args{main_window}} ->update;
      $ftp->cwd('/'.${$args{ini_ref}}{rfam_remote_dir}) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot change to directory ${$args{ini_ref}}{rfam_remote_dir}: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      ${$args{progress_bar}}->configure(-label=>"Retrieving remote Rfam release version: retrieve timestamp");
      ${$args{main_window}} ->update;
      $ftp->binary;
      $ftp->get('README', ${$args{auto_ini_ref}}{work_dir}.'/Rfam_rel_note') or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot retrieve release note: $@ ",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
      #read file into memory
      my ($file_ref) = &slurp(main_window   => $args{main_window},
                              auto_ini_ref  => $args{auto_ini_ref},
                              ini_ref       => $args{ini_ref},
                              filename      => 'Rfam_rel_note',
                              directory     => ${$args{auto_ini_ref}}{work_dir}
                             );
      'reset' =~ m/reset/;
      ${$file_ref} =~ m/\sRelease\s+(\S+)\s/s;
      if ($1) {$remote_release = $1} else {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot find Rfam release note in remote release note\.\nSetting release to 0",
                                           -icon    => 'error',
                                           -type    => 'ok');
         $remote_release = 0;
      };
      undef $file_ref;
      unlink ${$args{auto_ini_ref}}{work_dir}.'/Rfam_rel_note';


      #test for version
      if ($local_release >= $remote_release) {
         my $force_update = ${$args{main_window}}->Dialog(-title   => 'Force update',
                                                          -text    => "Local Rfam database appears to be more recent than remote.\nUpdate anyway?",
                                                          -bitmap  => 'question',
                                                          -buttons => ['Yes', 'No']
                                                         )->Show();
         if ($force_update eq 'Yes') {
            ${$args{update_ref}}{'update_ncRNA'} = 1;
         } else {
            ${$args{update_ref}}{'update_ncRNA'} = 0;
         }
      } else {
         ${$args{update_ref}}{'update_ncRNA'} = 1;
      }
      $ftp->quit;
   }
   return (1);
}

sub update {
   my %args = @_;
   my ($max_count, $db_count, $delete_files, $delete,
       $TIGRfam_release, $go_release, $ftp);
   $db_count = 0;
   #create status box
   &progress_bar_1(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Updating databases",
                   label        => 'Updating selected databases'
                  );
   &show_pbar_1;

   $delete_files = ${$args{main_window}}->Dialog(-title   => 'Download',
                                                 -text    => "Delete existing local files?\n If 'No', then data will be appended to existing files.",
                                                 -bitmap  => 'question',
                                                 -buttons => ['Yes', 'No']
                                                )->Show();
   if ($delete_files eq 'Yes') {
      $delete = 1;
   } else {
      $delete = 0;
   }

   #define how many dbs need to be updated
   while (my ($key, $value) = each %{ $args{update_ref} } ) {
      $max_count++ if ($value == 1);
   }

   #update COG if selected
   if (${$args{update_ref}}{'update_COG'} == 1) {
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Retrieving COG databases',
                     progress     => ($db_count / $max_count) * 100,
                    );

      my @files = qw(fun.txt myva org.txt whog);

      &retrieve(main_window  => $args{main_window},
                progress_bar => $args{progress_bar},
                auto_ini_ref => $args{auto_ini_ref},
                ini_ref      => $args{ini_ref},
                type         => 'COG',
                files        => \@files,
                delete       => $delete,
                server       => ${$args{ini_ref}}{COG_server_URL},
                path         => ${$args{ini_ref}}{COG_remote_dir},
                local        => ${$args{ini_ref}}{COG_db_path}.'/COG2003'
               );

      #test for incomplete download
      if (-e ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download') {
         unlink ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download';
         my $continue = ${$args{main_window}}->Dialog(-title   => 'Failed download',
                                                      -text    => "COG databases could not be downloaded successfully\.\nContinue with rest of updates?",
                                                      -bitmap  => 'question',
                                                      -buttons => ['Yes', 'No']
                                                     )->Show();
         if ($continue eq 'No') {
            hide_pbar_1;
            return(0);
         }
      }
   }

   #update COG2014 if selected
   if (${$args{update_ref}}{'update_COG2014'} == 1) {
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Retrieving COG2014 databases',
                     progress     => ($db_count / $max_count) * 100,
                    );

      my @files = qw(genomes2003-2014.tab cognames2003-2014.tab cog2003-2014.csv prot2003-2014.tab fun2003-2014.tab prot2003-2014.fa.gz);

      &retrieve(main_window  => $args{main_window},
                progress_bar => $args{progress_bar},
                auto_ini_ref => $args{auto_ini_ref},
                ini_ref      => $args{ini_ref},
                type         => 'COG2014',
                files        => \@files,
                delete       => $delete,
                server       => ${$args{ini_ref}}{COG_server_URL},
                path         => ${$args{ini_ref}}{COG2014_remote_dir},
                local        => ${$args{ini_ref}}{COG_db_path}.'/COG2014'
               );

      #test for incomplete download
      if (-e ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download') {
         unlink ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download';
         my $continue = ${$args{main_window}}->Dialog(-title   => 'Failed download',
                                                      -text    => "COG2014 databases could not be downloaded successfully\.\nContinue with rest of updates?",
                                                      -bitmap  => 'question',
                                                      -buttons => ['Yes', 'No']
                                                     )->Show();
         if ($continue eq 'No') {
            hide_pbar_1;
            return(0);
         }
      }
   }

   #update arCOG2014 if selected
   if (${$args{update_ref}}{'update_arCOG2014'} == 1) {
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Retrieving arCOG2014 databases',
                     progress     => ($db_count / $max_count) * 100,
                    );

      my @files = qw(ar14.arCOGdef.tab ar14.arCOG.csv funclass.tab ar14.fa.gz);

      &retrieve(main_window  => $args{main_window},
                progress_bar => $args{progress_bar},
                auto_ini_ref => $args{auto_ini_ref},
                ini_ref      => $args{ini_ref},
                type         => 'arCOG2014',
                files        => \@files,
                delete       => $delete,
                server       => ${$args{ini_ref}}{COG_server_URL},
                path         => ${$args{ini_ref}}{arCOG2014_remote_dir},
                local        => ${$args{ini_ref}}{COG_db_path}.'/arCOG2014'
               );

      #test for incomplete download
      if (-e ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download') {
         unlink ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download';
         my $continue = ${$args{main_window}}->Dialog(-title   => 'Failed download',
                                                      -text    => "arCOG2014 databases could not be downloaded successfully\.\nContinue with rest of updates?",
                                                      -bitmap  => 'question',
                                                      -buttons => ['Yes', 'No']
                                                     )->Show();
         if ($continue eq 'No') {
            hide_pbar_1;
            return(0);
         }
      }

      #rename ar14.arCOG.csv to make it compatible
      move("${$args{ini_ref}}{COG_db_path}\/arCOG2014\/ar14.arCOG.csv", "${$args{ini_ref}}{COG_db_path}\/arCOG2014\/arCOG14.txt");
   }

   #update POG2013 if selected
   if (${$args{update_ref}}{'update_POG2013'} == 1) {
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Retrieving POG2013 databases',
                     progress     => ($db_count / $max_count) * 100,
                    );

      my @files = qw(pogs.txt);

      &retrieve(main_window  => $args{main_window},
                progress_bar => $args{progress_bar},
                auto_ini_ref => $args{auto_ini_ref},
                ini_ref      => $args{ini_ref},
                type         => 'POG2013',
                files        => \@files,
                delete       => $delete,
                server       => ${$args{ini_ref}}{COG_server_URL},
                path         => ${$args{ini_ref}}{POG2013_remote_dir},
                local        => ${$args{ini_ref}}{COG_db_path}.'/POG2013'
               );

      #test for incomplete download
      if (-e ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download') {
         unlink ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download';
         my $continue = ${$args{main_window}}->Dialog(-title   => 'Failed download',
                                                      -text    => "POG2013 databases could not be downloaded successfully\.\nContinue with rest of updates?",
                                                      -bitmap  => 'question',
                                                      -buttons => ['Yes', 'No']
                                                     )->Show();
         if ($continue eq 'No') {
            hide_pbar_1;
            return(0);
         }
      }

      #retieve database files
      @files = qw(POGseqs_HighVQ.phr POGseqs_HighVQ.pin POGseqs_HighVQ.pnd POGseqs_HighVQ.pni POGseqs_HighVQ.pog POGseqs_HighVQ.psd POGseqs_HighVQ.psi POGseqs_HighVQ.psq);

      &retrieve(main_window  => $args{main_window},
                progress_bar => $args{progress_bar},
                auto_ini_ref => $args{auto_ini_ref},
                ini_ref      => $args{ini_ref},
                type         => 'POG2013',
                files        => \@files,
                delete       => $delete,
                server       => ${$args{ini_ref}}{COG_server_URL},
                path         => ${$args{ini_ref}}{POG2013_remote_dir}.'/blastdb',
                local        => ${$args{ini_ref}}{COG_db_path}.'/POG2013'
               );

      #test for incomplete download
      if (-e ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download') {
         unlink ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download';
         my $continue = ${$args{main_window}}->Dialog(-title   => 'Failed download',
                                                      -text    => "POG2013 databases could not be downloaded successfully\.\nContinue with rest of updates?",
                                                      -bitmap  => 'question',
                                                      -buttons => ['Yes', 'No']
                                                     )->Show();
         if ($continue eq 'No') {
            hide_pbar_1;
            return(0);
         }
      }
   }

   #update Pfam if selected
   if (${$args{update_ref}}{'update_pfam'} == 1) {
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Retrieving Pfam databases',
                     progress     => ($db_count / $max_count) * 100,
                    );

      my @files = qw(Pfam-A.hmm.gz Pfam-B.hmm.gz relnotes.txt);

      &retrieve(main_window  => $args{main_window},
                progress_bar => $args{progress_bar},
                auto_ini_ref => $args{auto_ini_ref},
                ini_ref      => $args{ini_ref},
                type         => 'PFam',
                files        => \@files,
                delete       => $delete,
                server       => ${$args{ini_ref}}{pfam_server_URL},
                path         => ${$args{ini_ref}}{pfam_remote_dir},
                local        => ${$args{ini_ref}}{pfam_db_path}
               );

      #test for incomplete download
      if (-e ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download') {
         unlink ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download';
         my $continue = ${$args{main_window}}->Dialog(-title   => 'Failed download',
                                                      -text    => "Pfam databases could not be downloaded successfully\.\nContinue with rest of updates?",
                                                      -bitmap  => 'question',
                                                      -buttons => ['Yes', 'No']
                                                     )->Show();
         if ($continue eq 'No') {
            hide_pbar_1;
            return(0);
         }
      }
   }

   #update Pfam descriptor if selected
   if (${$args{update_ref}}{'update_pfam_descr'} == 1) {
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Retrieving Pfam descriptor',
                     progress     => ($db_count / $max_count) * 100,
                    );

      my @files = qw(pfamA.txt.gz pfamB.txt.gz);

      &retrieve(main_window  => $args{main_window},
                progress_bar => $args{progress_bar},
                auto_ini_ref => $args{auto_ini_ref},
                ini_ref      => $args{ini_ref},
                type         => 'PFam descriptor',
                files        => \@files,
                delete       => $delete,
                server       => ${$args{ini_ref}}{pfam_server_URL},
                path         => ${$args{ini_ref}}{pfam_descriptor_dir},
                local        => ${$args{ini_ref}}{Pfam_descriptor}
               );

      #test for incomplete download
      if (-e ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download') {
         unlink ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download';
         my $continue = ${$args{main_window}}->Dialog(-title   => 'Failed download',
                                                      -text    => "Pfam descriptor could not be downloaded successfully\.\nContinue with rest of updates?",
                                                      -bitmap  => 'question',
                                                      -buttons => ['Yes', 'No']
                                                     )->Show();
         if ($continue eq 'No') {
            hide_pbar_1;
            return(0);
         }
      }
   }

   #update Interpro descriptor if selected
   if (${$args{update_ref}}{'update_interpro'} == 1) {
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Retrieving Interpro descriptor',
                     progress     => ($db_count / $max_count) * 100,
                    );

      my @files = qw(interpro.txt.gz);

      &retrieve(main_window  => $args{main_window},
                progress_bar => $args{progress_bar},
                auto_ini_ref => $args{auto_ini_ref},
                ini_ref      => $args{ini_ref},
                type         => 'Interpro descriptor',
                files        => \@files,
                delete       => $delete,
                server       => ${$args{ini_ref}}{pfam_server_URL},
                path         => ${$args{ini_ref}}{pfam_interpro_descr_dir},
                local        => ${$args{ini_ref}}{Interpro_descriptor}
               );

      #test for incomplete download
      if (-e ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download') {
         unlink ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download';
         my $continue = ${$args{main_window}}->Dialog(-title   => 'Failed download',
                                                      -text    => "Interpro descriptor could not be downloaded successfully\.\nContinue with rest of updates?",
                                                      -bitmap  => 'question',
                                                      -buttons => ['Yes', 'No']
                                                     )->Show();
         if ($continue eq 'No') {
            hide_pbar_1;
            return(0);
         }
      }
   }

   #update TIGRfam if selected
   if (${$args{update_ref}}{'update_tigrfam'} == 1) {
      my (@remote_dir);
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Retrieving TIGRfam release version',
                     progress     => ($db_count / $max_count) * 100,
                    );
      #get release version
      {
         $ftp = Net::FTP->new(${$args{ini_ref}}{tigrfam_server_URL}, Debug => 0) or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot connect to ${$args{ini_ref}}{tigrfam_server_URL}: $@ ",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return (0);
         };
         $ftp->login("anonymous",'-anonymous@') or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot login to ${$args{ini_ref}}{tigrfam_server_URL}: $@ ",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return (0);
         };
         $ftp->cwd('/'.${$args{ini_ref}}{tigrfam_remote_dir}) or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot change to directory ${$args{ini_ref}}{tigrfam_remote_dir}: $@ ",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return (0);
         };
         @remote_dir = $ftp->ls or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot retrieve remote directory: $@ ",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return (0);
         };
         #get remote version number
         foreach my $remote (@remote_dir) {
            if ($remote =~ m/TIGRFAMs_\d/) {
               'reset' =~ m/reset/;
               $remote =~ m/TIGRFAMs_(\d+\.?\d*)_/;
               if ($1) {$TIGRfam_release = $1; last;} else {
                  ${$args{main_window}}->messageBox(-title   => 'Error',
                                                    -message => "Cannot find remote TIGRfam release version\.",
                                                    -icon    => 'error',
                                                    -type    => 'ok');
                  return(0);
               };
            }
         }
      }

      my @files = ('TIGRFAMS_GO_LINK', 'TIGRFAMS_ROLE_LINK', ' TIGR_ROLE_NAMES',
                   'TIGRFAMs_'.$TIGRfam_release.'_HMM.LIB.gz',
                   'TIGRFAMs_'.$TIGRfam_release.'_INFO.tar.gz'
                   );

      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Retrieving TIGRfam databases',
                     progress     => ($db_count / $max_count) * 100,
                    );
      &retrieve(main_window  => $args{main_window},
                progress_bar => $args{progress_bar},
                auto_ini_ref => $args{auto_ini_ref},
                ini_ref      => $args{ini_ref},
                type         => 'TIGRFam',
                files        => \@files,
                delete       => $delete,
                server       => ${$args{ini_ref}}{tigrfam_server_URL},
                path         => ${$args{ini_ref}}{tigrfam_remote_dir},
                local        => ${$args{ini_ref}}{TIGRfam_db_path}
               );

      #test for incomplete download
      if (-e ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download') {
         unlink ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download';
         my $continue = ${$args{main_window}}->Dialog(-title   => 'Failed download',
                                                      -text    => "TIGRfam databases could not be downloaded successfully\.\nContinue with rest of updates?",
                                                      -bitmap  => 'question',
                                                      -buttons => ['Yes', 'No']
                                                     )->Show();
         if ($continue eq 'No') {
            hide_pbar_1;
            return(0);
         }
      }
   }

   #update GO if selected
   if (${$args{update_ref}}{'update_go'} == 1) {
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Retrieving Gene Ontology release version',
                     progress     => ($db_count / $max_count) * 100,
                    );

      #get remote release version
      {
         $ftp = Net::FTP->new(${$args{ini_ref}}{go_server_URL}, Debug => 0) or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot connect to ${$args{ini_ref}}{go_server_URL}: $@ ",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return (0);
         };
         $ftp->login("anonymous",'-anonymous@') or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot login to ${$args{ini_ref}}{go_server_URL}: $@ ",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return (0);
         };
         $ftp->cwd('/'.${$args{ini_ref}}{go_remote_dir}) or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot change to directory ${$args{ini_ref}}{tigrfam_remote_dir}: $@ ",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return (0);
         };
         my @remote_dir = $ftp->ls or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot retrieve remote directory: $@ ",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return (0);
         };
         #get remote version number
         foreach my $remote (@remote_dir) {
            if ($remote =~ m/^go_\d+/) {
               'reset' =~ m/reset/;
               $remote =~ m/^go_(\d+)\-/;
               if ($1) {$go_release = $1; last;} else {
                  ${$args{main_window}}->messageBox(-title   => 'Error',
                                                    -message => "Cannot find remote Gene Ontology release version\.",
                                                    -icon    => 'error',
                                                    -type    => 'ok');
                  return (0);
               };
            }
         }
      }
      my @files = ('go_'.$go_release.'-termdb.obo-xml.gz');

      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Retrieving Gene Ontology descriptor',
                     progress     => ($db_count / $max_count) * 100,
                    );
      &retrieve(main_window  => $args{main_window},
                progress_bar => $args{progress_bar},
                auto_ini_ref => $args{auto_ini_ref},
                ini_ref      => $args{ini_ref},
                type         => 'GO',
                files        => \@files,
                delete       => $delete,
                server       => ${$args{ini_ref}}{go_server_URL},
                path         => ${$args{ini_ref}}{go_remote_dir},
                local        => ${$args{ini_ref}}{GO_db_path}
               );

      #test for incomplete download
      if (-e ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download') {
         unlink ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download';
         my $continue = ${$args{main_window}}->Dialog(-title   => 'Failed download',
                                                      -text    => "Gene Ontology descriptor could not be downloaded successfully\.\nContinue with rest of updates?",
                                                      -bitmap  => 'question',
                                                      -buttons => ['Yes', 'No']
                                                     )->Show();
         if ($continue eq 'No') {
            hide_pbar_1;
            return(0);
         }
      }
   }

   #update Rfam if selected
   if (${$args{update_ref}}{'update_ncRNA'} == 1) {
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Retrieving Rfam databases',
                     progress     => ($db_count / $max_count) * 100,
                    );

      my @files = qw(README Rfam.fasta.gz Rfam.full.gz Rfam.tar.gz Rfam.thr.gz);

      &retrieve(main_window  => $args{main_window},
                progress_bar => $args{progress_bar},
                auto_ini_ref => $args{auto_ini_ref},
                ini_ref      => $args{ini_ref},
                type         => 'RFam',
                files        => \@files,
                delete       => $delete,
                server       => ${$args{ini_ref}}{rfam_server_URL},
                path         => ${$args{ini_ref}}{rfam_remote_dir},
                local        => ${$args{ini_ref}}{rfam_db_path}
               );

      #test for incomplete download
      if (-e ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download') {
         unlink ${$args{auto_ini_ref}}{work_dir}.'/incomplete_download';
         my $continue = ${$args{main_window}}->Dialog(-title   => 'Failed download',
                                                      -text    => "Rfam databases could not be downloaded successfully\.\nContinue with rest of updates?",
                                                      -bitmap  => 'question',
                                                      -buttons => ['Yes', 'No']
                                                     )->Show();
         if ($continue eq 'No') {
            hide_pbar_1;
            return(0);
         }
      }
   }

   &hide_pbar_1;
   return(1);
}

sub retrieve {
   my %args = (delete => 0,
               @_
              );
   my (@dir, @childs, $ftp);

   foreach my $file (@{ $args{files} }) {
      #remove existing files
      if ($args{delete} == 1) {
         unlink $args{local}.'/'.$file;
      }

      #begin forking
      my $pid = fork();
      if ($pid) {
        push(@childs, $pid);
      } elsif ($pid == 0) {
         my ($remote_filesize, $local_filesize, $retry);
         $ftp = Net::FTP->new($args{server}, Debug => 0) or do {
            print "Error in $args{type} database retrieval: Could not connect to server $args{server}";
            open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
            print WRITE "\nError in $args{type} database retrieval: Could not connect to server $args{server}";
            close WRITE;
            CORE::exit();
         };
         $ftp->login("anonymous",'-anonymous@') or do {
            print "Error in $args{type} database retrieval: Could not login to server $args{server}";
            open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
            print WRITE "\nError in $args{type} database retrieval: Could not login to server $args{server}";
            close WRITE;
            CORE::exit();
         };
         $ftp->cwd('/'.$args{path}) or do {
            print "Error in $args{type} database retrieval: Could not change to directory $args{path}";
            open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
            print WRITE "\nError in $args{type} database retrieval: Could not change to directory $args{path}";
            close WRITE;
            CORE::exit();
         };
         $ftp->binary;

         #does file exist? If not skip
         if ($ftp->size($file) == undef) {
            print "Error in $args{type} database retrieval: File $file does not exist in specified directory $args{path}";
            open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
            print WRITE "\nError in $args{type} database retrieval: File $file does not exist in specified directory $args{path}\n";
            close WRITE;
            $ftp->quit;
            CORE::exit();
         }

         #set up initial state
         $remote_filesize = $ftp->size($file);
         if (-e $args{local}.'/'.$file) {
            $local_filesize = stat($args{local}.'/'.$file)->size;
         } else {
            $local_filesize  = 0;
         }
         $retry = 0;

         #resume download when needed
         while ($local_filesize != $remote_filesize) {
            $ftp->restart($remote_filesize);
            $ftp->get($file, $args{local}.'/'.$file, $local_filesize) or do {
               print "Error in $args{type} database retrieval: Could not retrieve file $file, attempt $retry";
               open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
               print WRITE "\nError in $args{type} database retrieval: Could not retrieve file $file, attempt $retry\n";
               close WRITE;
            };
            $local_filesize = stat($args{local}.'/'.$file)->size;
            #abort after 5 tries
            if ($retry > 5) {
               open WRITE, "+>${$args{auto_ini_ref}}{work_dir}.'/incomplete_download";
               print WRITE 'dummy';
               close WRITE;
               $ftp->quit;
               CORE::exit();
            }
            $retry++;
         }
         $ftp->quit;

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
   return (1);
}

sub extract_db {
   my %args = @_;
   my ($max_count, $db_count);
   $db_count = 0;

   #create status box
   &progress_bar_1(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Updating databases",
                   label        => 'Extracting downloaded databases'
                  );
   &show_pbar_1;

   #define how many dbs need to be updated
   while (my ($key, $value) = each %{ $args{update_ref} } ) {
      $max_count++ if ($value == 1);
   }

   #compile COG if selected
   if (${$args{update_ref}}{'update_COG'} == 1) {
      unless (-e ${$args{ini_ref}}{COG_db_path}.'/COG2003/myva') {return (0)};
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Compiling COG databases',
                     progress     => ($db_count / $max_count) * 100,
                    );
      #format database
      `${$args{ini_ref}}{blast_executables}/formatdb -i ${$args{auto_ini_ref}}{COG_db_path}/myva -p T -o T -n ${$args{ini_ref}}{COG_db_path}/COG`;

      #test if all files exist
      my @test = qw(phr pin psd psi psq);
      foreach my $ext (@test) {
         return (0) unless (-e ${$args{auto_ini_ref}}{COG_db_path}.'/COG.'.$ext);
      }
   }

   #compile COG2014 if selected
   if (${$args{update_ref}}{'update_COG2014'} == 1) {
      unless (-e ${$args{ini_ref}}{COG_db_path}.'/COG2014/prot2003-2014.fa.gz') {return (0)};
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Compiling COG2014 databases',
                     progress     => ($db_count / $max_count) * 100,
                    );

      #uncompress
      `gunzip ${$args{auto_ini_ref}}{COG_db_path}/COG2014/prot2003-2014.fa.gz`;
      #format database
      `${$args{ini_ref}}{blast_plus_executables}/makeblastdb -dbtype prot -i ${$args{auto_ini_ref}}{COG_db_path}/COG2014/prot2003-2014.fa -parse_seqids -hash_index -out ${$args{auto_ini_ref}}{COG_db_path}/COG2014/COG2014`;

      #test if all files exist
      my @test = qw(phr pin psd psi psq);
      foreach my $ext (@test) {
         return (0) unless (-e ${$args{auto_ini_ref}}{COG_db_path}.'/COG2014/COG2014.'.$ext);
      }
   }

   #compile arCOG2014 if selected
   if (${$args{update_ref}}{'update_arCOG2014'} == 1) {
      unless (-e ${$args{ini_ref}}{COG_db_path}.'/arCOG2014/ar14.fa.gz') {return (0)};
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Compiling arCOG2014 databases',
                     progress     => ($db_count / $max_count) * 100,
                    );

      #uncompress
      `gunzip ${$args{auto_ini_ref}}{COG_db_path}/arCOG2014/ar14.fa.gz`;
      #format database
      `${$args{ini_ref}}{blast_plus_executables}/makeblastdb -dbtype prot -i ${$args{auto_ini_ref}}{COG_db_path}/arCOG2014/ar14.fa -parse_seqids -hash_index -out ${$args{auto_ini_ref}}{COG_db_path}/arCOG2014/arCOG2014`;

      #test if all files exist
      my @test = qw(phr pin psd psi psq);
      foreach my $ext (@test) {
         return (0) unless (-e ${$args{auto_ini_ref}}{COG_db_path}.'/arCOG2014/arCOG2014.'.$ext);
      }
   }

   #compile Pfam if selected
   if (${$args{update_ref}}{'update_pfam'} == 1) {
      unless (-e ${$args{ini_ref}}{pfam_db_path}.'/Pfam-A.hmm.gz') {return (0)};
      unless (-e ${$args{ini_ref}}{pfam_db_path}.'/Pfam-B.hmm.gz') {return (0)};
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Unpacking Pfam-A.hmm database',
                     progress     => ($db_count / $max_count) * 100,
                    );
      `gunzip ${$args{ini_ref}}{pfam_db_path}/Pfam-A.hmm.gz`;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Unpacking Pfam-B.hmm database',
                     progress     => ($db_count / $max_count) * 100,
                    );
      `gunzip ${$args{ini_ref}}{pfam_db_path}/Pfam-B.hmm.gz`;
   }

   #compile Pfam descriptor if selected
   if (${$args{update_ref}}{'update_pfam_descr'} == 1) {
      unless (-e ${$args{ini_ref}}{Pfam_descriptor}.'/pfamA.txt.gz') {return (0)};
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Unpacking PfamA descriptor',
                     progress     => ($db_count / $max_count) * 100,
                    );
      `gunzip ${$args{ini_ref}}{Pfam_descriptor}/pfamA.txt.gz`;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Unpacking PfamB descriptor',
                     progress     => ($db_count / $max_count) * 100,
                    );
      `gunzip ${$args{ini_ref}}{Pfam_descriptor}/pfamB.txt.gz`;
   }

   #compile interpro descriptor if selected
   if (${$args{update_ref}}{'update_interpro'} == 1) {
      unless (-e ${$args{ini_ref}}{Interpro_descriptor}.'/interpro.txt.gz') {return (0)};
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Unpacking Interpro descriptor',
                     progress     => ($db_count / $max_count) * 100,
                    );
      `gunzip ${$args{ini_ref}}{Interpro_descriptor}/interpro.txt.gz`;
   }

   #compile TIGRfam if selected
   if (${$args{update_ref}}{'update_tigrfam'} == 1) {
      #get release filenames
      opendir SEQINPUT, ${$args{ini_ref}}{TIGRfam_db_path};
      my @rel_version = grep /^TIGRFAMs_\d+/, readdir(SEQINPUT);
      closedir SEQINPUT;
      #determine highest release number
      my $TIGRfam_release = 0;
      foreach my $release (@rel_version) {
         'reset' =~ m/reset/;
         $release =~ m/^TIGRFAMs_(\d+\.?\d*)_/;
         if ($1 && $1 > $TIGRfam_release) {$TIGRfam_release = $1;}
         elsif ($1 && $1 <= $TIGRfam_release) {next;}
         else {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot find local TIGRfam release version\.",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return (0);
         }
      }
      return (0) if ($TIGRfam_release == 0); #return if no files found
      #test if files exist
      unless (-e ${$args{ini_ref}}{TIGRfam_db_path}.'TIGRFAMs_'.$TIGRfam_release.'_HMM.LIB.gz') {return (0)};
      unless (-e ${$args{ini_ref}}{TIGRfam_db_path}.'TIGRFAMs_'.$TIGRfam_release.'_INFO.tar.gz') {return (0)};

      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Unpacking TIGRfam HMM database',
                     progress     => ($db_count / $max_count) * 100,
                    );
      `gunzip ${$args{ini_ref}}{TIGRfam_db_path}/TIGRFAMs_$TIGRfam_release\_HMM.LIB.gz`;

      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Unpacking TIGRfam INFO database',
                     progress     => ($db_count / $max_count) * 100,
                    );
      #cleanup info direcotry if exist
      if (-d ${$args{ini_ref}}{TIGRfam_info}) {
          &cleanup(${$args{ini_ref}}{TIGRfam_info});
          mkdir ${$args{ini_ref}}{TIGRfam_info} or do {
             ${$args{main_window}}->messageBox(-title   => 'Error',
                                               -message => "Cannot create TIGRfam INFO directory ${$args{ini_ref}}{TIGRfam_info}",
                                               -icon    => 'error',
                                               -type    => 'ok');
             return (0);
          };
      }
      copy(${$args{ini_ref}}{TIGRfam_db_path}.'/TIGRFAMs_'.$TIGRfam_release.'_INFO.tar.gz', ${$args{ini_ref}}{TIGRfam_info}.'/TIGRFAMs_'.$TIGRfam_release.'_INFO.tar.gz');
      `tar -zxf ${$args{ini_ref}}{TIGRfam_db_path}/TIGRFAMs_$TIGRfam_release\_INFO.tar.gz`;
      unlink ${$args{ini_ref}}{TIGRfam_info}.'/TIGRFAMs_'.$TIGRfam_release.'_INFO.tar.gz';
   }

   #compile GO descriptor if selected
   if (${$args{update_ref}}{'update_go'} == 1) {
      #get release filenames
      opendir SEQINPUT, ${$args{ini_ref}}{GO_db_path};
      my @rel_version = grep /^go_\d+\-termdb\.obo/, readdir(SEQINPUT);
      closedir SEQINPUT;
      #determine highest release number
      my $go_release = 0;
      foreach my $release (@rel_version) {
         'reset' =~ m/reset/;
         $release =~ m/^go_(\d+)/;
         if ($1 && $1 > $go_release) {$go_release = $1;}
         elsif ($1 && $1 <= $go_release) {next;}
         else {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot find local Gene Ontology release version\.",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return (0);
         }
      }
      return (0) if ($go_release == 0); #return if no files found

      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Unpacking Gene Ontology descriptor',
                     progress     => ($db_count / $max_count) * 100,
                    );
      `gunzip ${$args{ini_ref}}{Interpro_descriptor}/go_$go_release\-termdb.obo-xml.gz`;
   }

   #compile RFam database if selected
   if (${$args{update_ref}}{'update_ncRNA'} == 1) {
      unless (-e ${$args{ini_ref}}{rfam_db_path}.'/Rfam.fasta.gz') {return (0)};
      unless (-e ${$args{ini_ref}}{rfam_db_path}.'/Rfam.full.gz') {return (0)};
      unless (-e ${$args{ini_ref}}{rfam_db_path}.'/Rfam.tar.gz') {return (0)};
      unless (-e ${$args{ini_ref}}{rfam_db_path}.'/Rfam.thr.gz') {return (0)};
      $db_count++;
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Unpacking Rfam fasta',
                     progress     => ($db_count / $max_count) * 100,
                    );
      `gunzip ${$args{ini_ref}}{rfam_db_path}/Rfam.fasta.gz`;

      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Compiling Rfam fasta',
                     progress     => ($db_count / $max_count) * 100,
                    );
      #format database
      `${$args{ini_ref}}{blast_executables}/formatdb -i ${$args{auto_ini_ref}}{rfam_db_path}/Rfam.fasta -p F -o T -n ${$args{ini_ref}}{rfam_db_path}/Rfam.fasta`;
      unlink ${$args{auto_ini_ref}}{rfam_db_path}.'/Rfam.fasta';

      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Unpacking Rfam full',
                     progress     => ($db_count / $max_count) * 100,
                    );
      `gunzip ${$args{ini_ref}}{rfam_db_path}/Rfam.full.gz`;

      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Unpacking Rfam threshold',
                     progress     => ($db_count / $max_count) * 100,
                    );
      `gunzip ${$args{ini_ref}}{rfam_db_path}/Rfam.thr.gz`;

      #cleanup info direcotry if exist
      &update_pbar_1(title        => 'Updating databases',
                     label        => 'Unpacking Rfam models',
                     progress     => ($db_count / $max_count) * 100,
                    );
      if (-d ${$args{ini_ref}}{rfam_cm_path}) {
          &cleanup(${$args{ini_ref}}{rfam_cm_path});
          mkdir ${$args{ini_ref}}{rfam_cm_path} or do {
             ${$args{main_window}}->messageBox(-title   => 'Error',
                                               -message => "Cannot create RFam model directory ${$args{ini_ref}}{rfam_cm_path}",
                                               -icon    => 'error',
                                               -type    => 'ok');
             return (0);
          };
      }
      copy(${$args{ini_ref}}{rfam_db_path}.'/Rfam.tar.gz', ${$args{ini_ref}}{rfam_cm_path}.'/Rfam.tar.gz');
      `tar -zxf ${$args{ini_ref}}{rfam_cm_path}/Rfam.tar.gz`;
      unlink ${$args{ini_ref}}{rfam_cm_path}.'/Rfam.tar.gz';
   }

   #delete downloaded source files?
   #undef $delete_files;
   #$delete_files = ${$args{main_window}}->Dialog(-title   => 'Delete Download',
   #                                              -text    => "Delete downloaded database source files?",
   #                                              -bitmap  => 'question',
   #                                              -buttons => ['Yes', 'No']
   #                                             )->Show();
   #if ($delete_files eq 'Yes') {
   #   my @delete = (${$args{ini_ref}}{COG_db_path}.'/myva',
   #                 ${$args{ini_ref}}{pfam_db_path}.'/Pfam_fs.gz',
   #                 ${$args{ini_ref}}{pfam_db_path}.'/Pfam_ls.gz',
   #                 ${$args{ini_ref}}{Pfam_descriptor}.'/pfamA.txt.gz',
   #                 ${$args{ini_ref}}{Pfam_descriptor}.'/pfamB.txt.gz',
   #                 ${$args{ini_ref}}{Pfam_descriptor}.'/interpro.txt.gz',
   #                 ${$args{ini_ref}}{TIGRfam_db_path}.'/TIGRFAMs_'.$TIGRfam_release.'_HMM.LIB.gz',
   #                 ${$args{ini_ref}}{TIGRfam_db_path}.'/TIGRFAMs_'.$TIGRfam_release.'_FRAG.LIB.gz',
   #                 ${$args{ini_ref}}{TIGRfam_db_path}.'/TIGRFAMs_'.$TIGRfam_release.'_INFO.tar.gz',
   #                 ${$args{ini_ref}}{GO_db_path}.'go_'.$go_release.'-termdb.obo-xml.gz',
   #                 ${$args{ini_ref}}{rfam_db_path}.'/Rfam.fasta.gz',
   #                 ${$args{ini_ref}}{rfam_db_path}.'/Rfam.full.gz',
   #                 ${$args{ini_ref}}{rfam_db_path}.'/Rfam.tar.gz',
   #                 ${$args{ini_ref}}{rfam_db_path}.'/Rfam.thr.gz',
   #                );
   #   foreach my $delete (@delete) {
   #      unlink $delete;
   #   }
   #}

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


1;