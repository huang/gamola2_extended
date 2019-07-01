#!/opt/ActivePerl-5.8/bin/perl

# setup advanced parameters for critica, glimmer2 and glimmer3
#input arguments: main_window, progress_bar, ini_ref, auto_ini_ref


package GeneModel::advanced_parameters;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&setup_critica &setup_glimmer2 &setup_glimmer3 &setup_prodigal);
use vars qw();

use Basics::MesgBox;
use initialise::read_me qw(:DEFAULT);


#local vars
my (%args, $tl, $file, $width);

sub setup_prodigal {
   my %args = @_;

   my ($closed_ends_cb, $translation_table_cb, $masked_sequence_cb, $bypass_cb, $procedure_cb, $override_cb);
   my ($translation_table_entry, $procedure_entry );
   my ($prodigal_tt);

   if (! Tk::Exists($tl)) {
      $tl    = ${$args{main_window}}->Toplevel();
      $tl->configure(-title => 'Advanced Prodigal Setup');
   } else {
      $tl->state('normal')
   }
   $tl -> raise();
   my $setup = $tl ->Frame(-borderwidth => 2, -relief => 'groove') ->grid (-row => 0, -column => 0, -sticky => 'nsew');
   my $ok    = $tl ->Frame(-borderwidth => 2, -relief => 'groove') ->grid (-row => 1, -column => 0, -sticky => 'nsew');

   #OK and Cancel buttons
   $ok-> Button(-text => 'OK',
                -command => sub {
                                  $tl->state('withdrawn');
                                  return;
                                 },
               )-> grid (-row => 0, -column => 0);
   #$ok-> Button(-text => 'Cancel',
   #             -command => sub { $tl->state('withdrawn');
   #                               return;},
   #            )-> grid (-row => 0, -column => 1);

   #setup buttons
   {
      $closed_ends_cb       = $setup -> Checkbutton(-text     => 'Closed ends.  Do not allow genes to run off edges',
                                                    -variable => \${$args{auto_ini_ref}}{prodigal_closed_ends},
                                                    ) -> grid (-row => 0, -column => 0, -sticky => 'w');
      $translation_table_cb = $setup -> Checkbutton (-text     => 'Specify a translation table to use (default 11)',
                                                     -variable => \$prodigal_tt,
                                                     -command  => sub {
                                                                    if (${$args{auto_ini_ref}}{prodigal_tt} == 1) {
                                                                       $translation_table_entry->configure(-state => 'normal');
                                                                    } else {
                                                                       $translation_table_entry->configure(-state => 'disabled');
                                                                    }
                                                                  }
                                                    ) -> grid (-row => 1, -column => 0, -sticky => 'w');
      $masked_sequence_cb   = $setup -> Checkbutton(-text     => "Treat runs of n\'s as masked sequence and do not build genes across them",
                                                    -variable => \${$args{auto_ini_ref}}{prodigal_masked_seq},
                                                    ) -> grid (-row => 2, -column => 0, -sticky => 'w');
      $bypass_cb            = $setup -> Checkbutton(-text     => "Bypass the Shine-Dalgarno trainer and force the program to scan for motifs",
                                                    -variable => \${$args{auto_ini_ref}}{prodigal_bypass},
                                                    ) -> grid (-row => 3, -column => 0, -sticky => 'w');
      $procedure_cb         = $setup -> Label      (-text     => "Select procedure (single or meta).",
                                                    ) -> grid (-row => 4, -column => 0, -sticky => 'w');
      $override_cb          = $setup -> Checkbutton(-text     => "Override procedure to \'Meta\' state is sequence is <100,000nt",
                                                    -variable => \${$args{auto_ini_ref}}{prodigal_override},
                                                    ) -> grid (-row => 4, -column => 1, -sticky => 'w');

   }
   #display and options
   {
      $width = 40;
      &define_initial_length(auto_ini_ref => $args{auto_ini_ref});

      $translation_table_entry = $setup       -> BrowseEntry(-label    =>'Translation table   ',
                                                             -choices  =>['1: The Standard',
                                                                          '2: The Vertebrate Mitochondrial Code',
                                                                          '3: The Yeast Mitochondrial Code',
                                                                          '4: The Mold, Protozoan, and Coelenterate Mitochondrial Code  and the Mycoplasma/Spiroplasma Code',
                                                                          '5: The Invertebrate Mitochondrial Code',
                                                                          '6: The Ciliate, Dasycladacean and Hexamita Nuclear Code',
                                                                          '9: The Echinoderm and Flatworm Mitochondrial Code',
                                                                          '10: The Euplotid Nuclear Code',
                                                                          '11: The Bacterial and Plant Plastid Code',
                                                                          '12: The Alternative Yeast Nuclear Code',
                                                                          '13: The Ascidian Mitochondrial Code',
                                                                          '14: The Alternative Flatworm Mitochondrial Code',
                                                                          '15: Blepharisma Nuclear Code',
                                                                          '16: Chlorophycean Mitochondrial Code',
                                                                          '21: Trematode Mitochondrial Code',
                                                                          '22: Scenedesmus obliquus mitochondrial Code',
                                                                          '23: Thraustochytrium Mitochondrial Code'
                                                                        ],
                                                            -width     => 40,
                                                            -variable  =>\${$args{auto_ini_ref}}{prodigal_translation_table}
                                              ) -> grid (-row => 1, -column => 1, -sticky => 'w');
      $procedure_entry = $setup              -> BrowseEntry(-label    =>'procedure   ',
                                                             -choices  =>['Single',
                                                                          'Meta'],
                                                             -variable  =>\${$args{auto_ini_ref}}{prodigal_procedure}
                                              ) -> grid (-row => 4, -column => 1, -sticky => 'w');

   }
}

sub setup_critica {
   my %args = @_;

   my ($scoring_matrix_cb, $dicodon_scores_cb, $init_scores_cb, $sd_scores_cb,
       $prom_scores_cb, $no_sdscores_cb, $prom_find_cb, $genetic_code_cb, $threshold_cb,
       $alpha_cb, $strict_threshold_cb, $frameshift_threshold_cb, $quick_stats_cb,
       $iterations_cb);
   my ($scoring_matrix_entry, $dicodon_scores_entry, $init_scores_entry, $sd_scores_entry,
       $prom_scores_entry, $no_sdscores_entry, $prom_find_entry, $genetic_code_entry, $threshold_entry,
       $alpha_entry, $strict_threshold_entry, $frameshift_threshold_entry, $quick_stats_entry,
       $iterations_entry);
   my ($scoring_matrix_bt, $dicodon_scores_bt, $init_scores_bt, $sd_scores_bt, $prom_scores_bt);

   if (! Tk::Exists($tl)) {
      $tl    = ${$args{main_window}}->Toplevel();
      $tl->configure(-title => 'Advanced Critica Setup');
   } else {
      $tl->state('normal')
   }
   $tl -> raise();
   my $setup = $tl ->Frame(-borderwidth => 2, -relief => 'groove') ->grid (-row => 0, -column => 0, -sticky => 'nsew');
   my $ok    = $tl ->Frame(-borderwidth => 2, -relief => 'groove') ->grid (-row => 1, -column => 0, -sticky => 'nsew');

   #OK and Cancel buttons
   $ok-> Button(-text => 'OK',
                -command => sub {
                                  my ($status) = &check_if_value ( main_window  => $args{main_window},
                                                    auto_ini_ref => $args{auto_ini_ref},
                                                    ini_ref      => $args{ini_ref});

                                  if ($status == 1) {
                                     $tl->state('withdrawn');
                                     return;
                                  } elsif ($status == 0) {
                                     &setup_critica(main_window  => $args{main_window},
                                                    auto_ini_ref => $args{auto_ini_ref},
                                                    ini_ref      => $args{ini_ref});
                                  }
                                 },
               )-> grid (-row => 0, -column => 0);
   #$ok-> Button(-text => 'Cancel',
   #             -command => sub { $tl->state('withdrawn');
   #                               return;},
   #            )-> grid (-row => 0, -column => 1);

   #setup buttons
   {
      $scoring_matrix_cb = $setup -> Checkbutton(-text => 'Scoring-matrix',
                                                 -variable => \${$args{auto_ini_ref}}{scoring_matrix_state},
                                                 -command  => sub {
                                                                    if (${$args{auto_ini_ref}}{scoring_matrix_state} == 1) {
                                                                       $scoring_matrix_entry->configure(-state => 'normal');
                                                                       $scoring_matrix_bt   ->configure(-state => 'normal');
                                                                    } else {
                                                                       $scoring_matrix_entry->configure(-state => 'disabled');
                                                                       $scoring_matrix_bt   ->configure(-state => 'disabled');
                                                                    }
                                                                  }
                                                 ) -> grid (-row => 0, -column => 0, -sticky => 'w');
      $dicodon_scores_cb = $setup -> Checkbutton (-text => 'Dicodon-scores',
                                                  -variable => \${$args{auto_ini_ref}}{dicodon_scores_state},
                                                  -command  => sub {
                                                                    if (${$args{auto_ini_ref}}{dicodon_scores_state} == 1) {
                                                                       $dicodon_scores_entry->configure(-state => 'normal');
                                                                       $dicodon_scores_bt   ->configure(-state => 'normal');
                                                                    } else {
                                                                       $dicodon_scores_entry->configure(-state => 'disabled');
                                                                       $dicodon_scores_bt   ->configure(-state => 'disabled');
                                                                    }
                                                                  }
                                                ) -> grid (-row => 1, -column => 0, -sticky => 'w');
      $init_scores_cb = $setup -> Checkbutton (-text => 'Init-scores',
                                               -variable => \${$args{auto_ini_ref}}{init_scores_state},
                                                -command  => sub {
                                                                    if (${$args{auto_ini_ref}}{init_scores_state} == 1) {
                                                                       $init_scores_entry->configure(-state => 'normal');
                                                                       $init_scores_bt   ->configure(-state => 'normal');
                                                                    } else {
                                                                       $init_scores_entry->configure(-state => 'disabled');
                                                                       $init_scores_bt   ->configure(-state => 'disabled');
                                                                    }
                                                                  }
                                             ) -> grid (-row => 2, -column => 0, -sticky => 'w');
      $sd_scores_cb = $setup -> Checkbutton (-text => 'SD-scores',
                                             -variable => \${$args{auto_ini_ref}}{sd_scores_state},
                                                 -command  => sub {
                                                                    if (${$args{auto_ini_ref}}{sd_scores_state} == 1) {
                                                                       $sd_scores_entry->configure(-state => 'normal');
                                                                       $sd_scores_bt   ->configure(-state => 'normal');
                                                                    } else {
                                                                       $sd_scores_entry->configure(-state => 'disabled');
                                                                       $sd_scores_bt   ->configure(-state => 'disabled');
                                                                    }
                                                                  }
                                           ) -> grid (-row => 3, -column => 0, -sticky => 'w');
      $prom_scores_cb = $setup -> Checkbutton (-text => 'Prom-scores',
                                               -variable => \${$args{auto_ini_ref}}{prom_scores_state},
                                                 -command  => sub {
                                                                    if (${$args{auto_ini_ref}}{prom_scores_state} == 1) {
                                                                       $prom_scores_entry->configure(-state => 'normal');
                                                                       $prom_scores_bt   ->configure(-state => 'normal');
                                                                    } else {
                                                                       $prom_scores_entry->configure(-state => 'disabled');
                                                                       $prom_scores_bt   ->configure(-state => 'disabled');
                                                                    }
                                                                  }
                                             ) -> grid (-row => 4, -column => 0, -sticky => 'w');
      $no_sdscores_cb = $setup -> Checkbutton (-text => 'No-sdscores',
                                               -variable => \${$args{auto_ini_ref}}{no_sdscores_state}
                                             ) -> grid (-row => 5, -column => 0, -sticky => 'w');
      $prom_find_cb = $setup -> Checkbutton (-text => 'Prom-find',
                                             -variable => \${$args{auto_ini_ref}}{prom_find_state}
                                           ) -> grid (-row => 6, -column => 0, -sticky => 'w');
      $genetic_code_cb = $setup -> Checkbutton (-text => 'Genetic code',
                                              -variable => \${$args{auto_ini_ref}}{genetic_code_state},
                                                 -command  => sub {
                                                                    if (${$args{auto_ini_ref}}{genetic_code_state} == 1) {
                                                                       $genetic_code_entry->configure(-state => 'normal');
                                                                    } else {
                                                                       $genetic_code_entry->configure(-state => 'disabled');
                                                                    }
                                                                  }
                                              ) -> grid (-row => 7, -column => 0, -sticky => 'w');
      $threshold_cb = $setup -> Checkbutton (-text => 'Threshold',
                                             -variable => \${$args{auto_ini_ref}}{threshold_state},
                                                 -command  => sub {
                                                                    if (${$args{auto_ini_ref}}{threshold_state} == 1) {
                                                                       $threshold_entry->configure(-state => 'normal');
                                                                    } else {
                                                                       $threshold_entry->configure(-state => 'disabled');
                                                                    }
                                                                  }
                                           ) -> grid (-row => 8, -column => 0, -sticky => 'w');
      $alpha_cb = $setup -> Checkbutton (-text => 'Alpha',
                                         -variable => \${$args{auto_ini_ref}}{alpha_state},
                                                 -command  => sub {
                                                                    if (${$args{auto_ini_ref}}{alpha_state} == 1) {
                                                                       $alpha_entry->configure(-state => 'normal');
                                                                    } else {
                                                                       $alpha_entry->configure(-state => 'disabled');
                                                                    }
                                                                  }
                                       ) -> grid (-row => 9, -column => 0, -sticky => 'w');
      $strict_threshold_cb = $setup -> Checkbutton (-text => 'Strict-threshold',
                                                    -variable => \${$args{auto_ini_ref}}{strict_threshold_state}
                                                  ) -> grid (-row => 10, -column => 0, -sticky => 'w');
      $frameshift_threshold_cb = $setup -> Checkbutton (-text => 'Frameshift-threshold',
                                                        -variable => \${$args{auto_ini_ref}}{frameshift_threshold_state},
                                                 -command  => sub {
                                                                    if (${$args{auto_ini_ref}}{frameshift_threshold_state} == 1) {
                                                                       $frameshift_threshold_entry->configure(-state => 'normal');
                                                                    } else {
                                                                       $frameshift_threshold_entry->configure(-state => 'disabled');
                                                                    }
                                                                  }
                                                      ) -> grid (-row => 11, -column => 0, -sticky => 'w');
      $quick_stats_cb = $setup -> Checkbutton (-text => 'Quick-stats',
                                               -variable => \${$args{auto_ini_ref}}{quick_stats_state}
                                             ) -> grid (-row => 12, -column => 0, -sticky => 'w');
      $iterations_cb = $setup -> Checkbutton (-text => 'Iterations',
                                              -variable => \${$args{auto_ini_ref}}{iterations_state},
                                                 -command  => sub {
                                                                    if (${$args{auto_ini_ref}}{iterations_state} == 1) {
                                                                       $iterations_entry->configure(-state => 'normal');
                                                                    } else {
                                                                       $iterations_entry->configure(-state => 'disabled');
                                                                    }
                                                                  }
                                             ) -> grid (-row => 13, -column => 0, -sticky => 'w');
   }
   #display and options
   {
      $width = 40;
      &define_initial_length(auto_ini_ref => $args{auto_ini_ref});

      $scoring_matrix_entry = $setup       -> Entry (-textvariable => \${$args{auto_ini_ref}}{scoring_matrix_file},
                                                     -width        => $width
                                              ) -> grid (-row => 0, -column => 1, -sticky => 'w');
      $dicodon_scores_entry = $setup       -> Entry (-textvariable => \${$args{auto_ini_ref}}{dicodon_scores_file},
                                                     -width        => $width
                                              ) -> grid (-row => 1, -column => 1, -sticky => 'w');
      $init_scores_entry = $setup          -> Entry (-textvariable => \${$args{auto_ini_ref}}{init_scores_file},
                                                     -width        => $width
                                              ) -> grid (-row => 2, -column => 1, -sticky => 'w');
      $sd_scores_entry = $setup            -> Entry (-textvariable => \${$args{auto_ini_ref}}{sd_scores_file},
                                                     -width        => $width
                                              ) -> grid (-row => 3, -column => 1, -sticky => 'w');
      $prom_scores_entry = $setup          -> Entry (-textvariable => \${$args{auto_ini_ref}}{prom_scores_file},
                                                     -width        => $width
                                              ) -> grid (-row => 4, -column => 1, -sticky => 'w');
      $no_sdscores_entry = $setup          -> Label (-text => 'Option: Select or Deselect',
                                                     -width        => $width
                                              ) -> grid (-row => 5, -column => 1, -sticky => 'w');
      $prom_find_entry = $setup            -> Label (-text => 'Option: Select or Deselect',
                                                     -width        => $width
                                              ) -> grid (-row => 6, -column => 1, -sticky => 'w');
      $genetic_code_entry = $setup         -> BrowseEntry(-choices  =>[1, 2, 3, 4, 5, 6, 9, 10, 11],
                                                          -width    => 3,
                                                          -variable =>\${$args{auto_ini_ref}}{genetic_code}
                                              ) -> grid (-row => 7, -column => 1, -sticky => 'w');
      $threshold_entry = $setup            -> Entry (-textvariable => \${$args{auto_ini_ref}}{threshold},
                                                     -width        => $width
                                              ) -> grid (-row => 8, -column => 1, -sticky => 'w');
      $alpha_entry = $setup                -> Entry (-textvariable => \${$args{auto_ini_ref}}{alpha},
                                                     -width        => $width
                                              ) -> grid (-row => 9, -column => 1, -sticky => 'w');
      $strict_threshold_entry = $setup     -> Label (-text => 'Option: Select or Deselect',
                                                     -width        => $width
                                              ) -> grid (-row => 10, -column => 1, -sticky => 'w');
      $frameshift_threshold_entry = $setup -> Entry (-textvariable => \${$args{auto_ini_ref}}{frameshift_threshold},
                                                     -width        => $width
                                              ) -> grid (-row => 11, -column => 1, -sticky => 'w');
      $quick_stats_entry = $setup          -> Label (-text => 'Option: Select or Deselect',
                                                     -width        => $width
                                              ) -> grid (-row => 12, -column => 1, -sticky => 'w');
      $iterations_entry = $setup           -> BrowseEntry(-choices  =>[1..10],
                                                          -width    => 3,
                                                          -variable =>\${$args{auto_ini_ref}}{iterations}
                                              ) -> grid (-row => 13, -column => 1, -sticky => 'w');
   }
   #add browse buttons
   {
      $scoring_matrix_bt = $setup -> Button (-text => 'Browse',
                                             -command => sub {
                                                                my ($file) = $tl->getOpenFile(-initialdir => ${$args{ini_ref}}{critica_bin});
                                                                if (defined $file && -e $file) {
                                                                   ${$args{ini_ref}}{scoring_matrix_file} = $file;
                                                                   if (length(${$args{ini_ref}}{scoring_matrix_file}) > $width) {
                                                                      $width = length(${$args{ini_ref}}{scoring_matrix_file});
                                                                      $scoring_matrix_entry->configure(-width => $width);
                                                                      $dicodon_scores_entry->configure(-width => $width);
                                                                      $init_scores_entry   ->configure(-width => $width);
                                                                      $sd_scores_entry     ->configure(-width => $width);
                                                                      $prom_scores_entry   ->configure(-width => $width);
                                                                   }
                                                                } else {
                                                                  $tl->messageBox(-title   => 'No selection',
                                                                                  -message => 'No scoring matrix file selected',
                                                                                  -type    => 'OK',
                                                                                  -icon    => 'info');
                                                                };
                                                                $tl->raise;

                                                              }
                                             ) -> grid (-row => 0, -column => 2, -sticky => 'w');
      $dicodon_scores_bt = $setup -> Button (-text => 'Browse',
                                             -command => sub {
                                                                my ($file) = $tl->getOpenFile(-initialdir => ${$args{ini_ref}}{critica_bin});
                                                                if (defined $file && -e $file) {
                                                                   ${$args{ini_ref}}{scoring_matrix_file} = $file;
                                                                   if (length(${$args{ini_ref}}{dicodon_scores_file}) > $width) {
                                                                      $width = length(${$args{ini_ref}}{dicodon_scores_file});
                                                                      $scoring_matrix_entry->configure(-width => $width);
                                                                      $dicodon_scores_entry->configure(-width => $width);
                                                                      $init_scores_entry   ->configure(-width => $width);
                                                                      $sd_scores_entry     ->configure(-width => $width);
                                                                      $prom_scores_entry   ->configure(-width => $width);
                                                                   }
                                                                } else {
                                                                  $tl->messageBox(-title   => 'No selection',
                                                                                  -message => 'No dicodon scoring file selected',
                                                                                  -type    => 'OK',
                                                                                  -icon    => 'info');
                                                                };
                                                                $tl->raise;

                                                              }
                                             ) -> grid (-row => 1, -column => 2, -sticky => 'w');
      $init_scores_bt = $setup -> Button (-text => 'Browse',
                                             -command => sub {
                                                                my ($file) = $tl->getOpenFile(-initialdir => ${$args{ini_ref}}{critica_bin});
                                                                if (defined $file && -e $file) {
                                                                   ${$args{ini_ref}}{init_scores_file} = $file;
                                                                   if (length(${$args{ini_ref}}{init_scores_file}) > $width) {
                                                                      $width = length(${$args{ini_ref}}{init_scores_file});
                                                                      $scoring_matrix_entry->configure(-width => $width);
                                                                      $dicodon_scores_entry->configure(-width => $width);
                                                                      $init_scores_entry   ->configure(-width => $width);
                                                                      $sd_scores_entry     ->configure(-width => $width);
                                                                      $prom_scores_entry   ->configure(-width => $width);
                                                                   }
                                                                } else {
                                                                  $tl->messageBox(-title   => 'No selection',
                                                                                  -message => 'No init scoring file selected',
                                                                                  -type    => 'OK',
                                                                                  -icon    => 'info');
                                                                };
                                                                $tl->raise;

                                                              }
                                             ) -> grid (-row => 2, -column => 2, -sticky => 'w');
      $sd_scores_bt = $setup -> Button (-text => 'Browse',
                                             -command => sub {
                                                                my ($file) = $tl->getOpenFile(-initialdir => ${$args{ini_ref}}{critica_bin});
                                                                if (defined $file && -e $file) {
                                                                   ${$args{ini_ref}}{sd_scores_file} = $file;
                                                                   if (length(${$args{ini_ref}}{sd_scores_file}) > $width) {
                                                                      $width = length(${$args{ini_ref}}{sd_scores_file});
                                                                      $scoring_matrix_entry->configure(-width => $width);
                                                                      $dicodon_scores_entry->configure(-width => $width);
                                                                      $init_scores_entry   ->configure(-width => $width);
                                                                      $sd_scores_entry     ->configure(-width => $width);
                                                                      $prom_scores_entry   ->configure(-width => $width);
                                                                   }
                                                                } else {
                                                                  $tl->messageBox(-title   => 'No selection',
                                                                                  -message => 'No SD scoring file selected',
                                                                                  -type    => 'OK',
                                                                                  -icon    => 'info');
                                                                };
                                                                $tl->raise;

                                                              }
                                             ) -> grid (-row => 3, -column => 2, -sticky => 'w');
      $prom_scores_bt = $setup -> Button (-text => 'Browse',
                                             -command => sub {
                                                                my ($file) = $tl->getOpenFile(-initialdir => ${$args{ini_ref}}{critica_bin});
                                                                if (defined $file && -e $file) {
                                                                   ${$args{ini_ref}}{prom_scores_file} = $file;
                                                                   if (length(${$args{ini_ref}}{prom_scores_file}) > $width) {
                                                                      $width = length(${$args{ini_ref}}{prom_scores_file});
                                                                      $scoring_matrix_entry->configure(-width => $width);
                                                                      $dicodon_scores_entry->configure(-width => $width);
                                                                      $init_scores_entry   ->configure(-width => $width);
                                                                      $sd_scores_entry     ->configure(-width => $width);
                                                                      $prom_scores_entry   ->configure(-width => $width);
                                                                   }
                                                                } else {
                                                                  $tl->messageBox(-title   => 'No selection',
                                                                                  -message => 'No prom scoring file selected',
                                                                                  -type    => 'OK',
                                                                                  -icon    => 'info');
                                                                };
                                                                $tl->raise;

                                                              }
                                             ) -> grid (-row => 4, -column => 2, -sticky => 'w');
   }
   #define initial state of browsing buttons
   {
      if (${$args{auto_ini_ref}}{scoring_matrix_state} == 0) {
         $scoring_matrix_entry->configure(-state => 'disabled');
         $scoring_matrix_bt   ->configure(-state => 'disabled');
      }
      if (${$args{auto_ini_ref}}{dicodon_scores_state} == 0) {
         $dicodon_scores_entry->configure(-state => 'disabled');
         $dicodon_scores_bt   ->configure(-state => 'disabled');
      }
      if (${$args{auto_ini_ref}}{init_scores_state} == 0) {
         $init_scores_entry->configure(-state => 'disabled');
         $init_scores_bt   ->configure(-state => 'disabled');
      }
      if (${$args{auto_ini_ref}}{sd_scores_state} == 0) {
         $sd_scores_entry->configure(-state => 'disabled');
         $sd_scores_bt   ->configure(-state => 'disabled');
      }
      if (${$args{auto_ini_ref}}{prom_scores_state} == 0) {
         $prom_scores_entry->configure(-state => 'disabled');
         $prom_scores_bt   ->configure(-state => 'disabled');
      }
      if (${$args{auto_ini_ref}}{genetic_code_state} == 0) {
         $genetic_code_entry->configure(-state => 'disabled');
      }
      if (${$args{auto_ini_ref}}{threshold_state} == 0) {
         $threshold_entry->configure(-state => 'disabled');
      }
      if (${$args{auto_ini_ref}}{alpha_state} == 0) {
         $alpha_entry->configure(-state => 'disabled');
      }
      if (${$args{auto_ini_ref}}{frameshift_threshold_state} == 0) {
         $frameshift_threshold_entry->configure(-state => 'disabled');
      }
      if (${$args{auto_ini_ref}}{iterations_state} == 0) {
         $iterations_entry->configure(-state => 'disabled');
      }

   }

}


sub define_initial_length {
   my %args = @_;
   if (length(${$args{auto_ini_ref}}{scoring_matrix_file}) > $width) { $width = length(${$args{auto_ini_ref}}{scoring_matrix_file})};
   if (length(${$args{auto_ini_ref}}{dicodon_scores_file}) > $width) { $width = length(${$args{auto_ini_ref}}{dicodon_scores_file})};
   if (length(${$args{auto_ini_ref}}{init_scores_file})    > $width) { $width = length(${$args{auto_ini_ref}}{init_scores_file})};
   if (length(${$args{auto_ini_ref}}{sd_scores_file})      > $width) { $width = length(${$args{auto_ini_ref}}{sd_scores_file})};
   if (length(${$args{auto_ini_ref}}{prom_scores_file})    > $width) { $width = length(${$args{auto_ini_ref}}{prom_scores_file})};

}

sub check_if_value {
   my %args = @_;
   if (${$args{auto_ini_ref}}{threshold} =~ /[^\de\-\.]+/) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Threshold is not a numeric value",
                                        -icon    => 'error',
                                        -type    => 'OK');
      return(0);
   }
   if (${$args{auto_ini_ref}}{alpha} =~ /[^\de\-\.]+/) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Alpha is not a numeric value",
                                        -icon    => 'error',
                                        -type    => 'OK');
      return(0);
   }
   if (${$args{auto_ini_ref}}{frameshift_threshold} =~ /[^\de\-\.]+/) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Frameshift threshold is not a numeric value",
                                        -icon    => 'error',
                                        -type    => 'OK');
      return(0);
   }
   return(1);

}
1;
