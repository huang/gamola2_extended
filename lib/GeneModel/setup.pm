#!/opt/ActivePerl-5.8/bin/perl5.8.8

package GeneModel::setup;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&setup_gene_model &check_critica_binaries &view_genemodel_options
             &forget_gm_selection $ini_ref $auto_ini_ref) ;
use vars qw(%ini %auto_ini $ini_ref $auto_ini_ref);

#local modules used
use initialise::read_me            qw(:DEFAULT);
use initialise::recompile          qw(:DEFAULT);
use Basics::MesgBox;
use Basics::CachingFind;
use GeneModel::advanced_parameters qw(:DEFAULT);
use Miscellaneous::custom_db       qw(:DEFAULT);
use File::Find;

#local variables
my (%args, $tl, $wait, $msg, $gl2_button, $gl3_button, $critica_button, $prodigal_button,
    $gl2_model_button, $gl3_model_button, $critica_database_button,
    $glimmer_widget, $prodigal_widget, $critica_widget, $RBSfinder_widget, $igblast_widget, $path_widget);

#determine what Gene Callers are there, currently supported: Glimmer 2, Glimmer 3 and Critica
sub setup_gene_model {
   my %args = @_;
   my %current_status = (); #capture current selection status for ORF callers in case of 'cancel'
   #setup parameters:
   $current_status{use_glimmer2}               = ${$args{auto_ini_ref}}{use_glimmer2};
   $current_status{use_glimmer3}               = ${$args{auto_ini_ref}}{use_glimmer3};
   $current_status{use_prodigal}               = ${$args{auto_ini_ref}}{use_prodigal};
   $current_status{prodigal_closed_ends}       = ${$args{auto_ini_ref}}{prodigal_closed_ends};
   $current_status{prodigal_translation_table} = ${$args{auto_ini_ref}}{prodigal_translation_table};
   $current_status{prodigal_masked_seq}        = ${$args{auto_ini_ref}}{prodigal_masked_seq};
   $current_status{prodigal_bypass}            = ${$args{auto_ini_ref}}{prodigal_bypass};
   $current_status{prodigal_procedure}         = ${$args{auto_ini_ref}}{prodigal_procedure};
   $current_status{prodigal_override}          = ${$args{auto_ini_ref}}{prodigal_override};
   $current_status{selected_gl_model}          = ${$args{auto_ini_ref}}{selected_gl_model};
   $current_status{gl_short_file}              = ${$args{auto_ini_ref}}{gl_short_file};
   $current_status{make_training_file}         = ${$args{auto_ini_ref}}{make_training_file};
   $current_status{use_critica}                = ${$args{auto_ini_ref}}{use_critica};
   $current_status{selected_critica_db}        = ${$args{auto_ini_ref}}{selected_critica_db};
   $current_status{full_critica_db}            = ${$args{auto_ini_ref}}{full_critica_db} ;
   $current_status{critica_bin}                = ${$args{auto_ini_ref}}{critica_bin};
   $current_status{critica_scripts}            = ${$args{auto_ini_ref}}{critica_scripts};
   $current_status{iterations}                 = ${$args{auto_ini_ref}}{iterations};
   $current_status{iterations_state}           = ${$args{auto_ini_ref}}{iterations_state};
   $current_status{scoring_matrix_state}       = ${$args{auto_ini_ref}}{scoring_matrix_state};
   $current_status{scoring_matrix_file}        = ${$args{auto_ini_ref}}{scoring_matrix_file};
   $current_status{dicodon_scores_state}       = ${$args{auto_ini_ref}}{dicodon_scores_state};
   $current_status{dicodon_scores_file}        = ${$args{auto_ini_ref}}{dicodon_scores_file};
   $current_status{init_scores_state}          = ${$args{auto_ini_ref}}{init_scores_state};
   $current_status{init_scores_file}           = ${$args{auto_ini_ref}}{init_scores_file};
   $current_status{sd_scores_state}            = ${$args{auto_ini_ref}}{sd_scores_state};
   $current_status{sd_scores_file}             = ${$args{auto_ini_ref}}{sd_scores_file};
   $current_status{no_sdscores_state}          = ${$args{auto_ini_ref}}{no_sdscores_state};
   $current_status{prom_find_state}            = ${$args{auto_ini_ref}}{prom_find_state};
   $current_status{prom_scores_state}          = ${$args{auto_ini_ref}}{prom_scores_state};
   $current_status{prom_scores_file}           = ${$args{auto_ini_ref}}{prom_scores_file};
   $current_status{genetic_code}               = ${$args{auto_ini_ref}}{genetic_code};
   $current_status{genetic_code_state}         = ${$args{auto_ini_ref}}{genetic_code_state};
   $current_status{threshold}                  = ${$args{auto_ini_ref}}{threshold};
   $current_status{threshold_state}            = ${$args{auto_ini_ref}}{threshold_state};
   $current_status{alpha}                      = ${$args{auto_ini_ref}}{alpha};
   $current_status{alpha_state}                = ${$args{auto_ini_ref}}{alpha_state} ;
   $current_status{strict_threshold_state}     = ${$args{auto_ini_ref}}{strict_threshold_state};
   $current_status{frameshift_threshold}       = ${$args{auto_ini_ref}}{frameshift_threshold};
   $current_status{frameshift_threshold_state} = ${$args{auto_ini_ref}}{frameshift_threshold_state};
   $current_status{quick_stats_state}          = ${$args{auto_ini_ref}}{quick_stats_state};
   $current_status{runintergenicblast}         = ${$args{auto_ini_ref}}{runintergenicblast};
   $current_status{intergenicblastdb}          = ${$args{auto_ini_ref}}{intergenicblastdb};
   $current_status{igblastminlength}           = ${$args{auto_ini_ref}}{igblastminlength};
   $current_status{igflanking}                 = ${$args{auto_ini_ref}}{igflanking};
   $current_status{igmaxevalue}                = ${$args{auto_ini_ref}}{igmaxevalue};
   $current_status{igblastsmart}               = ${$args{auto_ini_ref}}{igblastsmart};
   $current_status{IGsuperstition}             = ${$args{auto_ini_ref}}{IGsuperstition};
   $current_status{reuseIGresults}             = ${$args{auto_ini_ref}}{reuseIGresults};
   $current_status{runRBSfinder}               = ${$args{auto_ini_ref}}{runRBSfinder};

   #ask for ORF_callers to use
   my ($status) = &orf_caller_selection(label              => 'ORF caller supported',
                                        main_window        => $args{main_window},
                                        ini_ref            => $args{ini_ref},
                                        auto_ini_ref       => $args{auto_ini_ref},
                                        critica_button     => $args{critica_button},
                                        glimmer_button     => $args{glimmer_button},
                                        gl_training        => $args{gl_training},
                                        current_status_ref => \%current_status
                                        );

   #check for existing ORF callers
   if ($status == 1) {
      my $exe_status = &check_existing_execs(main_window  => $args{main_window},
                                             ini_ref      => $args{ini_ref},
                                             auto_ini_ref => $args{auto_ini_ref});
      if ($exe_status == 0) {
         #compile programmes
         &recompile(progress_bar  => $args{progress_bar},
                    main_window   => $args{main_window},
                    auto_ini_ref  => $args{auto_ini_ref},
                    ini_ref       => $args{ini_ref},
                    force_compile => 1
                   );
      }
   }
   &view_genemodel_options(main_window  => $args{main_window},
                           main_options => $args{main_options},
                           igm_frame    => $args{igm_frame},
                           ini_ref      => $args{ini_ref},
                           auto_ini_ref => $args{auto_ini_ref}
                          );
   ${$args{main_window}}->update;

}

sub view_genemodel_options {
   my %args = @_;

   #Glimmer 2 or 3
   if (${$args{auto_ini_ref}}{use_glimmer2} eq 1) {
      $glimmer_widget = ${$args{igm_frame}}->Label(-text => 'Glimmer 2; ')
                                           ->grid (-row => 0, -column => 0, -sticky => 'w');
   } elsif (${$args{auto_ini_ref}}{use_glimmer3} eq 1) {
      $glimmer_widget = ${$args{igm_frame}}->Label(-text => 'Glimmer 3; ')
                                           ->grid (-row => 0, -column => 0, -sticky => 'w');
   } else {
      $glimmer_widget = ${$args{igm_frame}}->Label(-text => 'Glimmer;   ',
                                                   -foreground => 'grey'
                                                      )
                                               ->grid (-row => 0, -column => 0, -sticky => 'w');
   }
   #Prodigal
   if (${$args{auto_ini_ref}}{use_prodigal} eq 1) {
      $prodigal_widget = ${$args{igm_frame}}->Label(-text => 'Prodigal; ')
                                            ->grid (-row => 0, -column => 1, -sticky => 'w');
   } else {
      $prodigal_widget = ${$args{igm_frame}}->Label(-text       => 'Prodigal; ',
                                                    -foreground => 'grey'
                                                    )
                                            ->grid (-row => 0, -column => 1, -sticky => 'w');
   }
   #Critica
   if (${$args{auto_ini_ref}}{use_critica} eq 1) {
      $critica_widget = ${$args{igm_frame}}->Label(-text => 'Critica; ')
                                              ->grid (-row => 0, -column => 2, -sticky => 'w');
   } else {
      $critica_widget = ${$args{igm_frame}}->Label(-text       => 'Critica; ',
                                                   -foreground => 'grey'
                                                     )
                                              ->grid (-row => 0, -column => 2, -sticky => 'w');
   }
   #RBSfinder
   if (${$args{auto_ini_ref}}{runRBSfinder} eq 1) {
      $RBSfinder_widget = ${$args{igm_frame}}->Label(-text => 'RBSfinder; ')
                                                ->grid (-row => 0, -column => 3, -sticky => 'w');
   } else {
      $RBSfinder_widget = ${$args{igm_frame}}->Label(-text => 'RBSfinder; ',
                                                     -foreground => 'grey'
                                                       )
                                                ->grid (-row => 0, -column => 3, -sticky => 'w');
   }
   #intergenic Blast
   if (${$args{auto_ini_ref}}{runintergenicblast} eq 1) {
      $igblast_widget = ${$args{igm_frame}}->Label(-text => 'intergenic Blast')
                                           ->grid (-row => 0, -column => 4, -sticky => 'w');
   } else {
      $igblast_widget = ${$args{igm_frame}}->Label(-text => 'intergenic Blast',
                                                   -foreground => 'grey'
                                                      )
                                           ->grid (-row => 0, -column => 4, -sticky => 'w');
   }
}

sub forget_gm_selection {
   my %args = @_;

   #Glimmer 2 or 3
   if (${$args{auto_ini_ref}}{use_glimmer2} eq 1) {
      $glimmer_widget = ${$args{igm_frame}}->Label(-text       => 'Glimmer 2; ',
                                                   -foreground => 'grey')
                                           ->grid (-row => 0, -column => 0, -sticky => 'w');
   } elsif (${$args{auto_ini_ref}}{use_glimmer3} eq 1) {
      $glimmer_widget = ${$args{igm_frame}}->Label(-text       => 'Glimmer 3; ',
                                                   -foreground => 'grey')
                                           ->grid (-row => 0, -column => 0, -sticky => 'w');
   } else {
      $glimmer_widget = ${$args{igm_frame}}->Label(-text => 'Glimmer;   ',
                                                   -foreground => 'grey'
                                                  )
                                           ->grid (-row => 0, -column => 0, -sticky => 'w');
   }
   #Prodigal
   $prodigal_widget =   ${$args{igm_frame}}->Label(-text       => 'Prodigal; ',
                                                   -foreground => 'grey'
                                                  )
                                           ->grid (-row => 0, -column => 1, -sticky => 'w');
   #Critica
   $critica_widget =    ${$args{igm_frame}}->Label(-text       => 'Critica; ',
                                                   -foreground => 'grey'
                                                  )
                                           ->grid (-row => 0, -column => 2, -sticky => 'w');
   #RBSfinder
   $RBSfinder_widget =  ${$args{igm_frame}}->Label(-text       => 'RBSfinder; ',
                                                   -foreground => 'grey'
                                                    )
                                           ->grid (-row => 0, -column => 3, -sticky => 'w');
   #intergenic Blast
   $igblast_widget =    ${$args{igm_frame}}->Label(-text       => 'intergenic Blast',
                                                   -foreground => 'grey'
                                                   )
                                           ->grid (-row => 0, -column => 4, -sticky => 'w');
}


sub check_existing_execs {
   my %args = @_;
   unless (-e ${$args{ini_ref}}{glprog2}.'/glimmer2') {
      if (${$args{auto_ini_ref}}{use_glimmer2} == 1) {
         ${$args{main_window}} -> messageBox(-title   => 'Error',
                                             -message => "Could not find Glimmer2 executable.\nRecompiling now",
                                             -type    => 'OK',
                                             -icon    => 'error');
         return (0);
      }
   }
   unless (-e ${$args{ini_ref}}{glprog3}.'/glimmer3') {
      if (${$args{auto_ini_ref}}{use_glimmer3} == 1) {
         ${$args{main_window}} -> messageBox(-title   => 'Error',
                                             -message => "Could not find Glimmer3 executable.\nRecompiling now",
                                             -type    => 'OK',
                                             -icon    =>'error');
         return (0);
      }
   }
   if (-d ${$args{ini_ref}}{prodigal}) {
      my $executable;
      find (sub {$executable = $File::Find::name if (-f && $_ =~ m/^prodigal/i) }, ${$args{ini_ref}}{prodigal});
      unless (defined ($executable)) {
         ${$args{main_window}} -> messageBox(-title   => 'Error',
                                             -message => "Could not find Prodigal executable.\nRecompiling now",
                                             -type    => 'OK',
                                             -icon    =>'error');
         return (0);
      }
   } else {
      ${$args{main_window}} -> messageBox(-title   => 'Error',
                                          -message => "Could not find Prodigal directory.\nRecompiling now",
                                          -type    => 'OK',
                                          -icon    =>'error');
      return (0);
   }

   unless (-e ${$args{auto_ini_ref}}{critica_scripts}.'/blast-contigs') {
      if (${$args{auto_ini_ref}}{use_critica} == 1) {
         ${$args{main_window}} -> messageBox(-title   => 'Error',
                                             -message => 'Verify location of Critica directory',
                                             -type    => 'OK',
                                             -icon    =>'error');
         &orf_caller_selection(main_window  => \${$args{main_window}},
                               label        => 'Verify Critica location',
                               ini_ref      => $args{ini_ref},
                               auto_ini_ref => $args{auto_ini_ref});
      }
   }
}

sub orf_caller_selection {
   my %args = @_;
   my $tl;
   my ($glimmer2_label, $glimmer3_label, $prodigal_label, $prodigal_parameters, $critica_label, $critica_parameters, $igblast, $reuseigblast, $ig_parameters);
   #set initial width for pathway label
   my $width = 40;
   if (length(${$args{ini_ref}}{glprog2})  > $width) { $width = length(${$args{ini_ref}}{glprog2})};
   if (length(${$args{ini_ref}}{glprog3})  > $width) { $width = length(${$args{ini_ref}}{glprog3})};
   if (length(${$args{ini_ref}}{prodigal}) > $width) { $width = length(${$args{ini_ref}}{prodigal})};
   if (length(${$args{ini_ref}}{critica})  > $width) { $width = length(${$args{ini_ref}}{critica})};

   if (defined $tl && Tk::Exists($tl)) {
       $tl->state('normal');
   } else {
      $tl = ${$args{main_window}}->DialogBox(-title          => $args{label},
                                             -buttons        => [ 'OK', 'Cancel' ],
                                             -default_button => 'OK'
                                            );


      $gl2_button = $tl->add ( 'Button',
                               -bitmap => '@'.${$args{auto_ini_ref}}{work_dir}.'/lib/initialise/tree.xbm',
                               -command => sub {
                                                  &dir_select(main_window  => \$tl,
                                                              label        => 'Path to Glimmer2 executable',
                                                              auto_ini_ref => $args{auto_ini_ref},
                                                              ini_ref      => $args{ini_ref},
                                                              directory    => 'glprog2');
                                                   $tl->raise;
                                                   if (length(${$args{ini_ref}}{glprog2}) > $width) {
                                                      $width = length(${$args{ini_ref}}{glprog2});
                                                      $glimmer2_label->configure(-width => $width);
                                                      $glimmer3_label->configure(-width => $width);
                                                      $critica_label ->configure(-width => $width);
                                                      $prodigal_label->configure(-width => $width);
                                                   }
                                                })      -> grid (-row => 0, -column => 3, -sticky => 'w');
      $gl3_button = $tl->add ( 'Button',
                                -bitmap => '@'.${$args{auto_ini_ref}}{work_dir}.'/lib/initialise/tree.xbm',
                                -command => sub {
                                                  &dir_select(main_window  => \$tl,
                                                              label        => 'Path to Glimmer3 executable',
                                                              auto_ini_ref => $args{auto_ini_ref},
                                                              ini_ref      => $args{ini_ref},
                                                              directory    => 'glprog3');
                                                   $tl->raise;
                                                   if (length(${$args{ini_ref}}{glprog3}) > $width) {
                                                      $width = length(${$args{ini_ref}}{glprog3});
                                                      $glimmer2_label->configure(-width => $width);
                                                      $glimmer3_label->configure(-width => $width);
                                                      $critica_label->configure(-width => $width);
                                                      $prodigal_label->configure(-width => $width);
                                                   }
                                                 })     -> grid (-row => 1, -column => 3, -sticky => 'w');
      $prodigal_button = $tl->add ( 'Button',
                                -bitmap => '@'.${$args{auto_ini_ref}}{work_dir}.'/lib/initialise/tree.xbm',
                                -command => sub {
                                                  &dir_select(main_window  => \$tl,
                                                              label        => 'Path to Prodigal executable',
                                                              auto_ini_ref => $args{auto_ini_ref},
                                                              ini_ref      => $args{ini_ref},
                                                              directory    => 'Prodigal');
                                                   $tl->raise;
                                                   if (length(${$args{ini_ref}}{prodigal}) > $width) {
                                                      $width = length(${$args{ini_ref}}{prodigal});
                                                      $glimmer2_label->configure(-width => $width);
                                                      $glimmer3_label->configure(-width => $width);
                                                      $critica_label->configure(-width => $width);
                                                      $prodigal_label->configure(-width => $width);
                                                   }
                                                 })     -> grid (-row => 2, -column => 3, -sticky => 'w');
      $prodigal_parameters = $tl->add ( 'Button',
                                       -text => 'Setup',
                                       -command => sub {
                                                         &setup_prodigal(main_window    => \$tl,
                                                                        ini_ref        => $args{ini_ref},
                                                                        auto_ini_ref   => $args{auto_ini_ref}
                                                                       );
                                                       }
                                     ) -> grid (-row => 2, -column => 4, -sticky => 'w');

      $critica_button = $tl->add ( 'Button',
                                    -bitmap => '@'.${$args{auto_ini_ref}}{work_dir}.'/lib/initialise/tree.xbm',
                                    -command => sub {
                                                      &dir_select(main_window  => \$tl,
                                                                  label        => 'Path to Critica directory',
                                                                  auto_ini_ref => $args{auto_ini_ref},
                                                                  ini_ref      => $args{ini_ref},
                                                                  directory    => 'critica');
                                                       $tl->raise;
                                                       if (length(${$args{ini_ref}}{critica}) > $width) {
                                                          $width = length(${$args{ini_ref}}{critica});
                                                          $glimmer2_label->configure(-width => $width);
                                                          $glimmer3_label->configure(-width => $width);
                                                          $critica_label->configure(-width => $width);
                                                          $prodigal_label->configure(-width => $width);
                                                      }
                                                      my ($status) = &check_critica_binaries(main_window  => \$tl,
                                                                                             ini_ref        => $args{ini_ref},
                                                                                             auto_ini_ref   => $args{auto_ini_ref}
                                                                                            );
                                                      if ($status == 0) {
                                                         #create empty paths
                                                         ${$args{auto_ini_ref}}{critica_bin} = "";
                                                         ${$args{auto_ini_ref}}{critica_scripts} = "";
                                                      }
                                                     }) -> grid (-row => 3, -column => 3, -sticky => 'w');

      $critica_parameters = $tl->add ( 'Button',
                                       -text => 'Setup',
                                       -command => sub {
                                                         &setup_critica(main_window    => \$tl,
                                                                        ini_ref        => $args{ini_ref},
                                                                        auto_ini_ref   => $args{auto_ini_ref}
                                                                       );
                                                       }
                                     ) -> grid (-row => 3, -column => 4, -sticky => 'w');


      $tl->add ( 'Checkbutton',
                 -text     => 'Use Glimmer2',
                 -variable => \${$args{auto_ini_ref}}{use_glimmer2},
                 -command  => [sub {
                                    if (${$args{auto_ini_ref}}{use_glimmer2} == 1) {
                                       ${$args{auto_ini_ref}}{use_glimmer3} = 0;
                                       $gl2_button              -> configure(-state => 'normal');
                                       $glimmer2_label          -> configure(-state => 'normal');
                                       $gl3_button              -> configure(-state => 'disabled');
                                       $glimmer3_label          -> configure(-state => 'disabled');
                                       ${$args{glimmer_button}} -> configure(-state => 'normal');
                                       ${$args{gl_training}}    -> configure(-state => 'normal');
                                    } elsif (${$args{auto_ini_ref}}{use_glimmer2} == 0) {
                                       $gl2_button              -> configure(-state => 'disabled');
                                       $glimmer2_label          -> configure(-state => 'disabled');
                                       ${$args{glimmer_button}} -> configure(-state => 'disabled');
                                       if (${$args{auto_ini_ref}}{use_glimmer3} == 0) {
                                          ${$args{gl_training}} -> configure(-state => 'disabled');
                                       }
                                    }
                                }]) ->grid (-row => 0, -column => 0, -sticky => 'w');
      $glimmer2_label = $tl->add ( 'Label',
                                   -textvariable => \${$args{ini_ref}}{glprog2},
                                   -width        => $width,
                                   -borderwidth  => 2,
                                   -relief       => 'sunken',
                                   -justify      => 'left')->grid (-row => 0, -column => 1, -columnspan => 2, -sticky => 'w');

      $tl->add( 'Checkbutton',
                 -text     => 'Use Glimmer3',
                 -variable => \${$args{auto_ini_ref}}{use_glimmer3},
                 -command  => [sub {
                                    if (${$args{auto_ini_ref}}{use_glimmer3} == 1) {
                                       ${$args{auto_ini_ref}}{use_glimmer2} = 0;
                                       $gl3_button              -> configure(-state => 'normal');
                                       $glimmer3_label          -> configure(-state => 'normal');
                                       $gl2_button              -> configure(-state => 'disabled');
                                       $glimmer2_label          -> configure(-state => 'disabled');
                                       ${$args{glimmer_button}} -> configure(-state => 'normal');
                                       ${$args{gl_training}}    -> configure(-state => 'normal');
                                    } elsif (${$args{auto_ini_ref}}{use_glimmer3} == 0) {
                                       $gl3_button              -> configure(-state => 'disabled');
                                       $glimmer3_label          -> configure(-state => 'disabled');
                                       ${$args{glimmer_button}} -> configure(-state => 'disabled');
                                       if (${$args{auto_ini_ref}}{use_glimmer2} == 0) {
                                          ${$args{gl_training}} -> configure(-state => 'disabled');
                                       }
                                    }
                                }]) ->grid (-row => 1, -column => 0, -sticky => 'w');
      $glimmer3_label = $tl->add( 'Label',
                                  -textvariable => \${$args{ini_ref}}{glprog3},
                                  -width        => $width,
                                  -borderwidth  => 2,
                                  -relief       => 'sunken',
                                  -justify      => 'left')->grid (-row => 1, -column => 1, -columnspan => 2, -sticky => 'w');

      $tl->add( 'Checkbutton',
                 -text     => 'Use Prodigal',
                 -variable => \${$args{auto_ini_ref}}{use_prodigal},
                 -command  => [sub {
                                    if (${$args{auto_ini_ref}}{use_prodigal} == 1) {
                                       $prodigal_button     -> configure(-state => 'normal');
                                       $prodigal_parameters -> configure(-state => 'normal');
                                       $prodigal_label      -> configure(-state => 'normal');
                                    } elsif (${$args{auto_ini_ref}}{use_prodigal} == 0) {
                                       $prodigal_button     -> configure(-state => 'disabled');
                                       $prodigal_parameters -> configure(-state => 'disabled');
                                       $prodigal_label      -> configure(-state => 'disabled');
                                    }
                                }]) ->grid (-row => 2, -column => 0, -sticky => 'w');
      $prodigal_label = $tl->add( 'Label',
                                  -textvariable => \${$args{ini_ref}}{prodigal},
                                  -width        => $width,
                                  -borderwidth  => 2,
                                  -relief       => 'sunken',
                                  -justify      => 'left')->grid (-row => 2, -column => 1, -columnspan => 2, -sticky => 'w');

      $tl->add( 'Checkbutton',
                 -text     => 'Use Critica',
                 -variable => \${$args{auto_ini_ref}}{use_critica},
                 -command  => [sub {
                                    if (${$args{auto_ini_ref}}{use_critica} == 1) {
                                       $critica_button          -> configure(-state => 'normal');
                                       $critica_parameters      -> configure(-state => 'normal');
                                       $critica_label           -> configure(-state => 'normal');
                                       ${$args{critica_button}} -> configure(-state => 'normal');
                                    } else {
                                       $critica_button          -> configure(-state => 'disabled');
                                       $critica_parameters      -> configure(-state => 'disabled');
                                       $critica_label           -> configure(-state => 'disabled');
                                       ${$args{critica_button}} -> configure(-state => 'disabled');
                                    }
                                }])     -> grid (-row => 3, -column => 0, -sticky => 'w');
      $critica_label = $tl->add( 'Label',
                                  -textvariable => \${$args{ini_ref}}{critica},
                                  -width        => $width,
                                  -borderwidth  => 2,
                                  -relief       => 'sunken',
                                  -justify      => 'left')-> grid (-row => 3, -column => 1, -columnspan => 2, -sticky => 'w');

      $ig_parameters = $tl->add ( 'Button',
                                 -text => 'Setup',
                                 -command => [sub {
                                                   &selectigblastdb(main_window  => \$tl,
                                                                    label        => 'Select intergentic Blast amino-acid database',
                                                                    ini_ref      => $args{ini_ref},
                                                                    auto_ini_ref => $args{auto_ini_ref}
                                                                   );
                                                 }
                                              ]
                               ) -> grid (-row => 4, -column => 4, -sticky => 'w');

      $reuseigblast = $tl->add( 'Checkbutton',
                               -text     => 'Re-use existing intergenic Blast results',
                               -variable => \${$args{auto_ini_ref}}{reuseIGresults}) ->grid (-row => 4, -column => 1, -sticky => 'w');

      $igblast = $tl->add( 'Checkbutton',
                          -text     => 'Intergenic Blast',
                          -variable => \${$args{auto_ini_ref}}{runintergenicblast},
                          -command  => [sub {
                                                if (${$args{auto_ini_ref}}{runintergenicblast} == 1) {
                                                   $ig_parameters -> configure(-state => 'normal');
                                                   $reuseigblast  -> configure(-state => 'normal');
                                                } else {
                                                   $ig_parameters -> configure(-state => 'disabled');
                                                   $reuseigblast  -> configure(-state => 'disabled');
                                                }
                                            }
                                       ]
                          )  ->grid (-row => 4, -column => 0, -sticky => 'w');

      $tl->add( 'Checkbutton',
                 -text     => 'Iterate with RBSfinder',
                 -variable => \${$args{auto_ini_ref}}{runRBSfinder}) ->grid (-row => 5, -column => 0, -sticky => 'w');


      #define initial state of browsing buttons
      if (${$args{auto_ini_ref}}{use_glimmer2} == 1) {
         ${$args{auto_ini_ref}}{use_glimmer3} = 0;
         $gl2_button     -> configure(-state => 'normal');
         $glimmer2_label -> configure(-state => 'normal');
         $gl3_button     -> configure(-state => 'disabled');
         $glimmer3_label -> configure(-state => 'disabled');
      } elsif (${$args{auto_ini_ref}}{use_glimmer3} == 1) {
         ${$args{auto_ini_ref}}{use_glimmer2} = 0;
         $gl3_button     -> configure(-state => 'normal');
         $glimmer3_label -> configure(-state => 'normal');
         $gl2_button     -> configure(-state => 'disabled');
         $glimmer2_label -> configure(-state => 'disabled');
      }
      if (${$args{auto_ini_ref}}{use_prodigal} == 1) {
         $prodigal_button     -> configure(-state => 'normal');
         $prodigal_parameters -> configure(-state => 'normal');
         $prodigal_label      -> configure(-state => 'normal');
      } elsif (${$args{auto_ini_ref}}{use_prodigal} == 0) {
         $prodigal_button     -> configure(-state => 'disabled');
         $prodigal_parameters -> configure(-state => 'disabled');
         $prodigal_label      -> configure(-state => 'disabled');
      }
      if (${$args{auto_ini_ref}}{use_critica} == 1) {
         $critica_button     -> configure(-state => 'normal');
         $critica_parameters -> configure(-state => 'normal');
         $critica_label      -> configure(-state => 'normal');
      } elsif (${$args{auto_ini_ref}}{use_critica} == 0) {
         $critica_button     -> configure(-state => 'disabled');
         $critica_parameters -> configure(-state => 'disabled');
         $critica_label      -> configure(-state => 'disabled');
      }
      if (${$args{auto_ini_ref}}{runintergenicblast} == 1) {
         $ig_parameters -> configure(-state => 'normal');
         $reuseigblast  -> configure(-state => 'normal');
      } else {
         $ig_parameters -> configure(-state => 'disabled');
         $reuseigblast  -> configure(-state => 'disabled');
      }

      #set default status for IG blast
      #if (${$args{auto_ini_ref}}{runintergenicblast} == 1) {
      #	$ig_parameters->configure(-state => 'normal');
      #} elsif (${$args{auto_ini_ref}}{runintergenicblast} == 0) {
      #	$ig_parameters->configure(-state => 'disabled');
      #}

      #set default status for re-use IG blast results
      #if (${$args{auto_ini_ref}}{reuseIGresults} == 1 && ${$args{auto_ini_ref}}{internal_gm} == 1) {
      #	$reuseigblast->configure(-state => 'normal');
      #} elsif (${$args{auto_ini_ref}}{reuseIGresults} == 0 || ${$args{auto_ini_ref}}{internal_gm} == 0) {
      #	$reuseigblast->configure(-state => 'disabled');
      #}

      $wait = $tl->Show();
   }

   if ($wait eq 'Cancel') {
      #reset ORF caller variables
      foreach my $key (%{ $args{current_status_ref} } ) {
         ${$args{auto_ini_ref}}{$key} = ${$args{current_status_ref}}{$key};
      }

      $tl->state('withdrawn');
      return(0);
   } elsif ($wait eq 'OK') {
      if (${$args{auto_ini_ref}}{use_critica} == 1) {
         my ($status) = &check_critica_binaries(main_window  => \$tl,
                                                ini_ref        => $args{ini_ref},
                                                auto_ini_ref   => $args{auto_ini_ref}
                                               );
         if ($status == 0) {
            #create empty paths
            ${$args{auto_ini_ref}}{critica_bin} = "";
            ${$args{auto_ini_ref}}{critica_scripts} = "";
         }
      }
      $tl->state('withdrawn');
      return(1);
   }
}

sub check_critica_binaries {
   my %args = @_;

   #find critica binaries and scripts
   my $critica_path = Basics::CachingFind->new(Path => [${$args{ini_ref}}{critica}]);
   ${$args{auto_ini_ref}}{critica_bin}     = $critica_path->findFirstInPath('critica');
   ${$args{auto_ini_ref}}{critica_scripts} = $critica_path->findFirstInPath('blast-contigs');

   unless (defined ${$args{auto_ini_ref}}{critica_bin} && ${$args{auto_ini_ref}}{critica_bin} =~ /\w+/) {
      ${$args{main_window}}->messageBox(-title   => 'Missing files',
                                        -message => "Critica binaries cannot be found in specified location\nTrying to recompile",
                                        -icon    => 'error',
                                        -type    => 'ok');

      my $status = &recompile(progress_bar  => $args{progress_bar},
                              main_window   => $args{main_window},
                              auto_ini_ref  => $args{auto_ini_ref},
                              ini_ref       => $args{ini_ref},
                              force_compile => 1
                             );

      if ($status eq '0') {return (0)};
   }
   unless (defined ${$args{auto_ini_ref}}{critica_scripts} && ${$args{auto_ini_ref}}{critica_scripts} =~ /\w+/) {
      ${$args{main_window}}->messageBox(-title   => 'Missing files',
                                        -message => "Critica scripts cannot be found in specified location\nTrying to recompile",
                                        -icon    => 'error',
                                        -type    => 'ok');

      my $status = &recompile(progress_bar  => $args{progress_bar},
                              main_window   => $args{main_window},
                              auto_ini_ref  => $args{auto_ini_ref},
                              ini_ref       => $args{ini_ref},
                              force_compile => 1
                             );

      if ($status eq '0') {return (0)};
   }
   #strip test names from paths
   ${$args{auto_ini_ref}}{critica_bin} =~ s/\/critica$//;
   ${$args{auto_ini_ref}}{critica_scripts} =~ s/\/blast-contigs$//;
   return (1);
}

sub selectigblastdb {
   my %args = @_;
   my ($tl_ig, $wait, $default, $default_evalue, $blastdb, $blastdb_label, $short_file,
       $igblastsmart, $igblastfull, $igclear);

   #set initial width for pathway label
   my $width = 40;
   my $types = [
               ["Amino acid db", ['.pin', '.pal']],
               ["All Files"    , ["*"]]
               ];
   my $default_model = '.pin';

   #set default number of CPUs to system value
   ${$args{auto_ini_ref}}{IG_CPU} = ${$args{auto_ini_ref}}{CPU};

   #set default cluster value
   #${$args{auto_ini_ref}}{igcluster} = 0;

   #set default blast behaviour
   #${$args{auto_ini_ref}}{IG_blast_threaded} = 1;

   #save default value
   $default = ${$args{auto_ini_ref}}{intergenicblastdb};
   $default_evalue = $args{auto_ini_ref}{igmaxevalue};

   if (defined $tl_ig && Tk::Exists($tl_ig)) {
       $tl_ig->state('normal');
   } else {
      $tl_ig = ${$args{main_window}}->DialogBox(-title          => $args{label},
                                                -buttons        => [ 'OK', 'Cancel' ],
                                                -default_button => 'OK'
                                               );

      $tl_ig->add ('Label',
                   -text => 'DB selection'
                   ) -> grid (-row => 0, -column => 0, -sticky => 'w');
      $blastdb = $tl_ig->add ( 'Button',
                              -bitmap => '@'.${$args{auto_ini_ref}}{work_dir}.'/lib/initialise/tree.xbm',
                              -command => sub {
                                                ${$args{auto_ini_ref}}{intergenicblastdb} = $tl_ig->getOpenFile(-filetypes        => $types,
                                                                                                                -initialdir       => ${$args{ini_ref}}{blast_db_path},
                                                                                                                -initialfile      => ${$args{auto_ini_ref}}{full_blast_db},
                                                                                                                -defaultextension => \$default_model);
                                                unless (defined ${$args{auto_ini_ref}}{intergenicblastdb} && -e ${$args{auto_ini_ref}}{intergenicblastdb}) {
                                                   $tl_ig->messageBox(-title   => 'Error',
                                                                      -message => 'No db file selected, maintaining original selection',
                                                                      -type    => 'OK',
                                                                      -icon    => 'error'
                                                                     );

                                                   ${$args{auto_ini_ref}}{intergenicblastdb} = $default;
                                                   return;
                                                };

                                                unless (${$args{auto_ini_ref}}{intergenicblastdb} =~ /\.(pin|pal|phr|pnd|pni|psd|psq)$/) {
                                                   $tl_ig->messageBox(-title   => 'Error',
                                                                      -message => 'Wrong db format selected, maintaining original selection',
                                                                      -type    => 'OK',
                                                                      -icon    => 'error'
                                                                     );

                                                   ${$args{auto_ini_ref}}{intergenicblastdb} = $default;
                                                   return;
                                                };

                                                ${$args{auto_ini_ref}}{intergenicblastdb} =~ m/([^\/]*?)$/;
                                                $short_file = $1;
                                                ${$args{auto_ini_ref}}{intergenicblastdb} =~ s/\.(pin|pal|phr|pnd|pni|psd|psq)$//;

                                                if (length($short_file) > $width) {
                                                    $width = length($short_file);
                                                    $blastdb_label->configure(-width => $width);
                                                }
                                              })
                    -> grid (-row => 0, -column => 3, -sticky => 'w');

      $blastdb_label = $tl_ig->add ( 'Label',
                                   -textvariable => \$short_file,
                                   -width        => $width,
                                   -borderwidth  => 2,
                                   -relief       => 'sunken',
                                   -justify      => 'left')
                    ->grid (-row => 0, -column => 1, -columnspan => 2, -sticky => 'w');

      $tl_ig->add ('Label',
                -text => 'Create custom DB'
               ) -> grid (-row => 1, -column => 0, -sticky => 'w');
      $tl_ig->add ( 'Button',
                 -bitmap => '@'.${$args{auto_ini_ref}}{work_dir}.'/lib/initialise/tree.xbm',
                 -command => sub {
                                    ${$args{auto_ini_ref}}{intergenicblastdb} = &create_custom_db(main_window    => \${$args{main_window}},
                                                                                                  ini_ref        => $args{ini_ref},
                                                                                                  auto_ini_ref   => $args{auto_ini_ref});
                                    unless (defined ${$args{auto_ini_ref}}{intergenicblastdb} && -e ${$args{auto_ini_ref}}{intergenicblastdb}) {
                                       $tl_ig->messageBox(-title   => 'Error',
                                                          -message => 'Error creating custom database, maintaining original selection',
                                                          -type    => 'OK',
                                                          -icon    => 'error');

                                       ${$args{auto_ini_ref}}{intergenicblastdb} = $default;
                                    };

                                    ${$args{auto_ini_ref}}{intergenicblastdb} =~ m/([^\/]*?)$/;
                                    $short_file = $1;

                                    if (length($short_file) > $width) {
                                        $width = length($short_file);
                                        $blastdb_label->configure(-width => $width);
                                    }
                                    $blastdb_label->configure(-textvariable => \$short_file);
                                    $tl_ig->update;
                                  })
           -> grid (-row => 1, -column => 3, -sticky => 'w');
      $tl_ig->BrowseEntry(-label   =>'Minimum intergenic ORF length [nt]',
                       -choices  =>[1..1000],
                       -width    => 6,
                       -variable =>\$args{auto_ini_ref}{igblastminlength}
                      )
           ->grid (-row => 2, -column => 0, -sticky => 'w');
      $tl_ig->BrowseEntry(-label   =>'Size of IG flanking sequence [nt]',
                       -choices  =>[0..3000],
                       -width    => 6,
                       -variable =>\$args{auto_ini_ref}{igflanking}
                      )
      ->grid (-row => 3, -column => 0, -sticky => 'w');

      $tl_ig->BrowseEntry(-label   =>'Maximum e-value cutoff',
                       -choices  =>['0', '1e-150', '1e-100', '1e-90', '1e-80', '1e-70', '1e-60', '1e-50', '1e-40',
                                    '1e-30', '1e-20', '1e-10', '1e-05', '1e-04', '1e-03', '1e-02', '1e-01', '1', '10'],
                       -width    => 6,
                       -variable =>\$args{auto_ini_ref}{igmaxevalue}
                      )
           ->grid (-row => 4, -column => 0, -sticky => 'w');
      $igblastsmart = $tl_ig->add( 'Checkbutton',
                                -text     => 'Strand independent Intergenic Blast (only use intergenic regions between both strands)',
                                -variable => \${$args{auto_ini_ref}}{igblastsmart},
                                -command  => [sub {
                                                   if (${$args{auto_ini_ref}}{igblastsmart} == 1) {
                                                      ${$args{auto_ini_ref}}{igblastfull} = 0;
                                                   } elsif (${$args{auto_ini_ref}}{igblastsmart} == 0) {
                                                      ${$args{auto_ini_ref}}{igblastfull} = 1;
                                                   }
                                                  }]
                               )
           ->grid (-row => 5, -column => 0, -sticky => 'w');

      $igblastfull = $tl_ig->add( 'Checkbutton',
                                -text     => 'Full Intergenic Blast (use strand specific intergenic regions (takes longer))',
                                -variable => \${$args{auto_ini_ref}}{igblastfull},
                                -command  => [sub {
                                                   if (${$args{auto_ini_ref}}{igblastfull} == 1) {
                                                      ${$args{auto_ini_ref}}{igblastsmart} = 0;
                                                   } elsif (${$args{auto_ini_ref}}{igblastfull} == 0) {
                                                      ${$args{auto_ini_ref}}{igblastsmart} = 1;
                                                   }
                                                  }]
                              )
           ->grid (-row => 6, -column => 0, -sticky => 'w');

      $tl_ig->add( 'Checkbutton',
                  -text     => "Remove ORFs < 200nt from ORIGINAL gene model,\nthen build IG model",
                  -variable => \${$args{auto_ini_ref}}{IGsuperstition}
                      )
            ->grid (-row => 7, -column => 0, -sticky => 'w');

      #$igclear = $tl_ig->add( 'Checkbutton',
      #                     -text     => 'Clear IG Blast results after Gene model analysis',
      #                     -variable => \${$args{auto_ini_ref}}{clearIGrun}
      #                   )
      #      ->grid (-row => 8, -column => 0, -sticky => 'w');

      $tl_ig->BrowseEntry(-label   =>'Cluster individual IG regions (0 = none)',
                       -choices  =>[0..500],
                       -width    => 6,
                       -variable =>\${$args{auto_ini_ref}}{igcluster}
                      )
           ->grid (-row => 8, -column => 0, -sticky => 'w');
      $tl_ig->BrowseEntry(-label    =>'CPUs/Cores for IG Blast',
                                      -choices  =>[1..100],
                                      -width    => 5,
                                      -variable =>\${$args{auto_ini_ref}}{IG_CPU}
                                     )
         ->grid (-row => 9, -column => 0, -sticky => 'w');
      $tl_ig->Checkbutton(-text     => 'Cluster CPUs for Blast',
                          -variable => \${$args{auto_ini_ref}}{IG_blast_cluster},
                          -command  => sub {
                                            if (${$args{auto_ini_ref}}{IG_blast_cluster} == 1) {
                                               ${$args{auto_ini_ref}}{IG_blast_threaded} = 0;
                                            } elsif (${$args{auto_ini_ref}}{IG_blast_cluster} == 0) {
                                               ${$args{auto_ini_ref}}{IG_blast_threaded} = 1;
                                            }
                                           }
                        )
         ->grid (-row => 10, -column => 0, -sticky => 'w');
      $tl_ig->Checkbutton(-text     => 'Multi-thread Blast',
                          -variable => \${$args{auto_ini_ref}}{IG_blast_threaded},
                          -command  => sub {
                                            if (${$args{auto_ini_ref}}{IG_blast_threaded} == 1) {
                                               ${$args{auto_ini_ref}}{IG_blast_cluster} = 0;
                                            } elsif (${$args{auto_ini_ref}}{IG_blast_threaded} == 0) {
                                               ${$args{auto_ini_ref}}{IG_blast_cluster} = 1;
                                            }
                                           }
                         )
         ->grid (-row => 11, -column => 0, -sticky => 'w');

      #set default status for IG blast
      if (${$args{auto_ini_ref}}{igblastsmart} == 1 ) {
         ${$args{auto_ini_ref}}{igblastfull} = 0;
      } elsif (${$args{auto_ini_ref}}{igblastsmart} == 0) {
         ${$args{auto_ini_ref}}{igblastfull} = 1;
      }

      #show default value for ig blast db'
      ${$args{auto_ini_ref}}{intergenicblastdb} =~ m/([^\/]*?)$/;
      $short_file = $1;
      $blastdb_label->configure(-textvariable => \$short_file);

      $wait = $tl_ig->Show();
   }

   if ($wait eq 'Cancel') {
      $tl_ig->state('withdrawn');
      #restore default value
      ${$args{auto_ini_ref}}{intergenicblastdb} = $default;
      $args{auto_ini_ref}{igmaxevalue} = $default_evalue;
      return(0);
   } elsif ($wait eq 'OK') {
      $tl_ig->state('withdrawn');
      return(1);
   }
}
1;