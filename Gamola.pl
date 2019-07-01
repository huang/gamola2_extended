#!/opt/ActivePerl-5.8/bin/perl -w
#GAMOLA v2.00
#GUI

#IMPORTATNT:
#export PERL5LIB="/opt/ActivePerl-5.8/lib"


use warnings;
use strict;
use 5.008;
use FindBin qw($Bin);
use lib $Bin.'/lib';
use Tk;
use Tk::ErrorDialog;
use Tk::Adjuster;
use Tk::ProgressBar;
use Tk::TList;
use Tk::ErrorDialog;
use Tk::BrowseEntry;
use Tk::NoteBook;

use initialise::init                   qw(:DEFAULT);
use initialise::read_me                qw(:DEFAULT);
use initialise::verify                 qw(:DEFAULT);
use initialise::gb_header              qw(:DEFAULT);
use GeneModel::setup                   qw(:DEFAULT);
use ProgrammeModules::CentralStation   qw(:DEFAULT);
use Basics::error_popup;
use Basics::progress_bar               qw(:DEFAULT);
use db_selection::databases            qw(:DEFAULT);
use Miscellaneous::custom_db           qw(:DEFAULT);
use Miscellaneous::rotate_gb           qw(:DEFAULT);
use Miscellaneous::sequin              qw(:DEFAULT);
use Miscellaneous::TIGRfam_annotation  qw(:DEFAULT);
use Miscellaneous::transfer_annotation qw(:DEFAULT);
use Miscellaneous::update_db           qw(:DEFAULT);
use Miscellaneous::metagenome          qw(:DEFAULT);
use File::Copy;

#require some Perl modules
require Cwd;
require File::Copy;
require File::Find;
require File::stat;
require Proc::Background;
require Time::localtime;
#require Net::FTP;

#variables for main window
use vars qw($a $b $c $welcome $text $ttt $i
            $mw $nb $mw_frame_top $mw_frame_bottom $mw_frame_left $mw_frame_middle $mw_frame_right $mw_frame_final
            $startup_options $main_options $database_options $structural_options $run_options $misc_options
            $progress_bar $state $test $top_ref $reference_text
            $select_int_gm $select_ext_gm );

#exported variables to be used later on
use vars qw($read_file %ini %ini_comment %auto_ini
            $function $value $file $comment $user $auto @types $Bin);

my ($filter_blast, $glimmer_button, $glimmer_add, $glimmer_make_training, $critica_add, $critica_clear,
    $blast_button, $blast_sel, $Blast_add, $Blast_clear, $gene_model_only, $blast_plus , $blast_legacy,
    $COG_button, $COG_sel, $Pfam_sel, $Pfam_fam_sel, $Pfam_version2, $Pfam_version3, $PFam3_thrshld, $TIGRfam_fam_sel,
    $Pfam_button, $Pfam_clear, $Pfam_add, $Rfam_button, $Rfam_clear, $Rfam_add,
    $TIGRfam_button, $TIGRfam_clear, $TIGRfam_add, $TIGRfam_sel,
    $tigr_anno, $TIGR_EConly, $TIGR_ECCDS, $TIGR_ECdiscard, @display_parameters, @auto_display_parameters, $path_widget,
    $allow_cog_button, $critica_button, @colors, $p_bar, $trna_source, $trna_sensi, $trna_save,
    $tmm_domain, $signalp_version, $signalp_source, $signalp_cutoff, $signalp_min, $signalp_trunc, $signalp_method,
    $CRISPR_minNR, $CRISPR_minRL, $CRISPR_maxRL, $CRISPR_minSL, $CRISPR_maxSL, $status,
    $select_files, $file_lb, $conc_lb, $conc_lb_max_width, $clear_selection, $clear_clusters, $read_clusters, %all_input_files, @short_input_files, $in_file_width, $current_in_dir, $refresh_files, $concatenate_files,
    $nb2, $nb3, $nb4, $System_properties, $Blast_properties, $Genbank_options, $Genbank_header_options, $Fasta_header_options,
    $msg_timeout,
    $RNA, $membrane, $DNAstructure, $others,
    $db_options, $gb_files_options, $annotation_options, $metagenome_options);
my (%concatenate_clusters, %concat_file_name);
    #$CRISPR_searchWL
our %ini;
our %auto_ini;
our $display_parameters = \@display_parameters;
our $auto_display_parameters = \@auto_display_parameters;


#start GUI
$mw = MainWindow->new;

# read last saved values if exist or default values
if (-e $Bin.'/last.saved') {
   ($status) = &read_ini (main_window             => \$mw,
                          file                    => $Bin.'/last.saved',
                          ini_ref                 => \%ini,
                          auto_ini_ref            => \%auto_ini,
                          ini_comment_ref         => \%ini_comment,
                          display_parameters      => $display_parameters,
                          auto_display_parameters => $auto_display_parameters,
                          home_dir					 => $Bin);
} else {
   ($status) = &read_ini (main_window             => \$mw,
                          file                    => $Bin.'/Gamola.default',
                          ini_ref                 => \%ini,
                          auto_ini_ref            => \%auto_ini,
                          ini_comment_ref         => \%ini_comment,
                          display_parameters      => $display_parameters,
                          auto_display_parameters => $auto_display_parameters,
                          home_dir					 => $Bin);
}

unless (defined $status && $status == 1) {
   $mw->messageBox(-title   => 'Error',
                   -message => "Can't read variables from files.\nAborting",
                   -type    => 'OK',
                   -icon    => 'error'
                  );
   exit;
}
#delete error log file
unlink $Bin.'/Error.log';

#define initial input folder setup
$current_in_dir = $auto_ini{work_dir};
$in_file_width  = 10;
$auto_ini{concatenation_counter} = 0; #this is a floating variable to keep track of number of file groupings;

#define initial message box timeout
$msg_timeout = $auto_ini{msg_timeout} / 1000;

#set COG database options
my %cog_type = ();
$cog_type{'COG2003'}   = {('info' => 'Original NCBI 2003 COG database is expected.',
                           'file' => 'COG.phr',
                           'path' => 'COG2003',
                           'db'   => 'COG'
                         )};
$cog_type{'COG2008'}   = {('info' => 'Updated NCBI 2008 COG database is expected.',
                           'file' => 'COG2008.phr',
                           'path' => 'COG2008',
                           'db'   => 'COG2008'
                         )};
$cog_type{'COG2014'}   = {('info' => 'Updated NCBI 2014 COG database is expected.',
                           'file' => 'COG2014.phr',
                           'path' => 'COG2014',
                           'db'   => 'COG2014'
                         )};
$cog_type{'arCOG'  }   = {('info' => 'Original archaeal NCBI COG database is expected.',
                           'file' => 'arCOG.phr',
                           'path' => 'arCOG',
                           'db'   => 'arCOG'
                         )};
$cog_type{'arCOG2014'} = {('info' => 'Updated archaeal NCBI 2014 COG database is expected.',
                           'file' => 'arCOG2014.phr',
                           'path' => 'arCOG2014',
                           'db'   => 'arCOG2014'
                         )};
$cog_type{'POG2013'}   = {('info' => 'Original NCBI 2013 phage COG (POG) database is expected.',
                           'file' => 'POGseqs_HighVQ.phr',
                           'path' => 'POG2013',
                           'db'   => 'POGseqs_HighVQ'
                         )};

#initially populate input file lists
($status) = &grab_input_folder_content(ini_ref      => \%ini,
                                       auto_ini_ref => \%auto_ini,
                                      );



#Main Window
{
   $mw->optionAdd('*BorderWidth' => 1);
   $mw->title("GAMOLA 2");
   #$mw->minsize(400,600);
   #$mw->maxsize(600,800);
   #$mw->resizable(0,0);

   #display welcome screen
   &welcome;
   $mw->focus;

   $mw->Label(-text => "Global Annotation of multiplexed on-site Blasted DNA-sequences")->grid (-row => 0, -column => 0, -sticky => 'nsew');
   $nb = $mw->NoteBook()                                                                ->grid (-row => 1, -column => 0, -sticky => 'nsew');
   $mw_frame_bottom = $mw->Frame(-borderwidth => 2, -relief => 'groove')                ->grid (-row => 2, -column => 0, -sticky => 'nsew');
   $progress_bar = $mw->Frame(-borderwidth => 2, -relief => 'groove')                   ->grid (-row => 3, -column => 0, -sticky => 'nsew');
   $progress_bar->configure(-label=>" ");
   #$progress_bar -> Label(-text=>" ") -> grid (-row => 0, -column => 0); #create dummy for space

   #resize proportionally
   $mw->gridColumnconfigure(0, -weight => 1);
   $mw->gridRowconfigure(1, -weight => 1);

   #page 1: startup options
   $startup_options        = $nb             ->add('page1', -label => 'Startup options');
   $nb2                    = $startup_options->NoteBook()->grid (-row => 0, -column => 0, -sticky => 'nsew');
   $System_properties      = $nb2            ->add('page7',  -label => 'System properties');
   $Blast_properties       = $nb2            ->add('page8',  -label => 'Blast properties');
   $Genbank_options        = $nb2            ->add('page9',  -label => 'Genbank options');
   $Genbank_header_options = $nb2            ->add('page10', -label => 'Genbank header options');
   $Fasta_header_options   = $nb2            ->add('page11', -label => 'Fasta header options');

   $startup_options->gridColumnconfigure(0, -weight => 1);
   $startup_options->gridRowconfigure   (0, -weight => 1);

   #page2: main options
   $main_options           = $nb             ->add('page2',  -label => 'Main options');

   #page3: structural options
   $structural_options     = $nb             ->add('page3',  -label => 'Structural analyses');
   $nb4                    = $structural_options->NoteBook()->grid (-row => 0, -column => 0, -sticky => 'nsew');
   $RNA                    = $nb4            ->add('page16', -label => 'tRNA/rRNA/ncRNA');
   $membrane               = $nb4            ->add('page17', -label => 'Translocation/Membrane');
   $DNAstructure           = $nb4            ->add('page18', -label => 'DNA structures');
   $others                 = $nb4            ->add('page19', -label => 'Other analyses');

   $structural_options->gridColumnconfigure(0, -weight => 1);
   $structural_options->gridRowconfigure   (0, -weight => 1);

   #page3: database options
   $database_options       = $nb             ->add('page4',  -label => 'Database selection');

   #page5: run
   $run_options            = $nb             ->add('page5',  -label => 'Execute');

   #page6: miscallenous options
   $misc_options           = $nb             ->add('page6',  -label => 'Miscallenous Options');
   $nb3                    = $misc_options   ->NoteBook()->grid (-row => 0, -column => 0, -sticky => 'nsew');
   $db_options             = $nb3            ->add('page12', -label => 'Databases');
   $gb_files_options       = $nb3            ->add('page13', -label => 'Genbank files');
   $annotation_options     = $nb3            ->add('page14', -label => 'Annotation');
   $metagenome_options     = $nb3            ->add('page15', -label => 'Metagenome');

   $misc_options->gridColumnconfigure(0, -weight => 1);
   $misc_options->gridRowconfigure   (0, -weight => 1);

   #are all input filenames OK? Check status from above.
   if ($status == 0) {
      $mw->DialogBox(-title          => 'Error: input filenames must not contain special characters or spaces. Please check.',
                     -buttons        => [ 'OK' ],
                     -default_button => 'OK'
                     );
   }
}

#generating needed progress bar widgets
{
   &progress_bar_1(main_window  => \$mw,
                   progress_bar => \$progress_bar,
                   auto_ini_ref => \%auto_ini,
                   ini_ref      => \%ini,
                   title        => ' ',
                   label        => ' '
                  );
   &progress_bar_2(main_window  => \$mw,
                   progress_bar => \$progress_bar,
                   auto_ini_ref => \%auto_ini,
                   ini_ref      => \%ini,
                   title        => ' ',
                   label        => ' '
                  );
   &progress_bar_3(main_window  => \$mw,
                   progress_bar => \$progress_bar,
                   auto_ini_ref => \%auto_ini,
                   ini_ref      => \%ini,
                   title        => ' ',
                   label        => ' '
                  );
}


#setup bottom frame
{
   $mw_frame_bottom->Button(-text => "Read last saved",
                            -command => sub {&read_ini(main_window             => \$mw,
                                                       file                    => $Bin.'/last.saved',
                                                       ini_ref                 => \%ini,
                                                       auto_ini_ref            => \%auto_ini,
                                                       ini_comment_ref         => \%ini_comment,
                                                       path_widget             => \$path_widget,
                                                       display_parameters      => $display_parameters,
                                                       auto_display_parameters => $auto_display_parameters,
                                                       progress_bar            => \$progress_bar);
                                             $mw->update;
                                            })    ->grid (-row => 0, -column => 0);
   $mw_frame_bottom->Button(-text => "Read default",
                            -command => sub {&read_ini(main_window             => \$mw,
                                                       file                    => $Bin.'/Gamola.default',
                                                       ini_ref                 => \%ini,
                                                       auto_ini_ref            => \%auto_ini,
                                                       ini_comment_ref         => \%ini_comment,
                                                       path_widget             => \$path_widget,
                                                       display_parameters      => $display_parameters,
                                                       auto_display_parameters => $auto_display_parameters,
                                                       progress_bar            => \$progress_bar);
                                             $mw->update;
                                            }) ->grid (-row => 0, -column => 1);
   $mw_frame_bottom->Button(-text => "Save current selection",
                            -command => sub {&save_ini(main_window             => \$mw,
                                                       file                    => $Bin.'/last.saved',
                                                       ini_ref                 => \%ini,
                                                       auto_ini_ref            => \%auto_ini,
                                                       ini_comment_ref         => \%ini_comment,
                                                       display_parameters      => $display_parameters,
                                                       auto_display_parameters => $auto_display_parameters)
                                            })     ->grid (-row => 0, -column => 2);
   $mw_frame_bottom->Button(-text => "Modify variables",
                            -command => sub {
                                             my $file;
                                             if (-e $Bin.'/last.saved') {
                                                $file = $Bin.'/last.saved';
                                             } else {
                                                $file = $Bin.'/Gamola.default';
                                             }
                                             &mod_vars(main_window             => \$mw,
                                                       file                    => $file,
                                                       ini_ref                 => \%ini,
                                                       auto_ini_ref            => \%auto_ini,
                                                       ini_comment_ref         => \%ini_comment,
                                                       display_parameters      => $display_parameters,
                                                       auto_display_parameters => $auto_display_parameters)
                                            })     ->grid (-row => 0, -column => 3);
   $mw_frame_bottom->Button(-text => "View References",
                            -command => \&view_references)     ->grid (-row => 0, -column => 4);
   $mw_frame_bottom->Button(-text => "Exit",
                            -command => sub {
                                             &restore_input_files(ini_ref      => \%ini,
                                                                  auto_ini_ref => \%auto_ini
                                                                 );
                                             &kill_all_pbar;
                                             Tk::exit(0);
                                            })                          ->grid (-row => 0, -column => 5);
}

#setup page 1: startup options
{
   {
      $System_properties->Label(-text => 'Existing Results')
         ->grid (-row => 0, -column => 0);
      $System_properties->Checkbutton(-text     => 'Re-use existing results',
                                      -variable => \$auto_ini{reuse_results},
                                      -command  => sub {
                                                        if ($auto_ini{reuse_results} == 1) {
                                                           $auto_ini{clear_results} = 0;
                                                           &read_clusters(main_window      => \$mw,
                                                                         ini_ref           => \%ini,
                                                                         auto_ini_ref      => \%auto_ini,
                                                                         textbox           => \$file_lb,
                                                                         concbox           => \$conc_lb,
                                                                         frame             => \$run_options,
                                                                         width             => \$in_file_width,
                                                                         conc_clusters_ref => \%concatenate_clusters
                                                                        );
                                                        } elsif ($auto_ini{reuse_results} == 0) {
                                                           $auto_ini{clear_results} = 1;
                                                           %concatenate_clusters            = ();
                                                           %concat_file_name                = ();
                                                           if (defined $conc_lb) {
                                                              $conc_lb->selectAll;
                                                              $conc_lb->deleteSelected;
                                                           }
                                                           #$conc_lb->delete(0.0,'end');
                                                           $mw->update;
                                                        }
                                                       }
                                      )
         ->grid (-row => 1, -column => 0, -sticky => 'w');
      $System_properties->Checkbutton(-text => 'Clear existing results',
                                    -variable => \$auto_ini{clear_results},
                                    -command => sub {
                                                     if ($auto_ini{clear_results} == 1) {
                                                        $auto_ini{reuse_results} = 0;
                                                        %concatenate_clusters            = ();
                                                        %concat_file_name                = ();
                                                        $auto_ini{concatenation_counter} = 0;
                                                        if (defined $conc_lb) {
                                                           $conc_lb->selectAll;
                                                           $conc_lb->deleteSelected;
                                                        }
                                                        $mw->update;
                                                     } elsif ($auto_ini{clear_results} == 0) {
                                                        $auto_ini{reuse_results} = 1;
                                                        &read_clusters(main_window      => \$mw,
                                                                      ini_ref           => \%ini,
                                                                      auto_ini_ref      => \%auto_ini,
                                                                      textbox           => \$file_lb,
                                                                      concbox           => \$conc_lb,
                                                                      frame             => \$run_options,
                                                                      width             => \$in_file_width,
                                                                      conc_clusters_ref => \%concatenate_clusters
                                                                     );
                                                     }
                                                    }
                                    )
         ->grid (-row => 2, -column => 0, -sticky => 'w');

      $System_properties->Label(-text => ' ')
         ->grid (-row => 4, -column => 0);

      $System_properties->Label(-text => 'GUI behaviour')
         ->grid (-row => 5, -column => 0);
      $System_properties->BrowseEntry(-label    =>'Close Message Boxes [sec]',
                                      -choices  =>[0..1000],
                                      -width    => 5,
                                      -variable =>\$msg_timeout,
                                      -command  => sub { $auto_ini{msg_timeout} = $msg_timeout * 1000 }
                                     )
         ->grid (-row => 6, -column => 0, -sticky => 'w');


      $System_properties->Label(-text => 'Distributed computing')
         ->grid (-row => 0, -column => 1, -columnspan => 2);
      $System_properties->BrowseEntry(-label    =>'CPUs/Cores in system',
                                      -choices  =>[1..100],
                                      -width    => 5,
                                      -variable =>\$auto_ini{CPU}
                                     )
         ->grid (-row => 1, -column => 1, -sticky => 'w');
      $System_properties->Checkbutton(-text     => 'Cluster CPUs for Blast',
                                     -variable => \$auto_ini{blast_cluster},
                                     -command  => sub {
                                                       if ($auto_ini{blast_cluster} == 1) {
                                                          $auto_ini{blast_threaded} = 0;
                                                       } elsif ($auto_ini{blast_cluster} == 0) {
                                                          $auto_ini{blast_threaded} = 1;
                                                       }
                                                      }
                                   )
         ->grid (-row => 2, -column => 1, -sticky => 'w');
      $System_properties->Checkbutton(-text     => 'Multi-thread Blast',
                                     -variable => \$auto_ini{blast_threaded},
                                     -command  => sub {
                                                       if ($auto_ini{blast_threaded} == 1) {
                                                          $auto_ini{blast_cluster} = 0;
                                                       } elsif ($auto_ini{blast_threaded} == 0) {
                                                          $auto_ini{blast_cluster} = 1;
                                                       }
                                                      }
                                    )
         ->grid (-row => 3, -column => 1, -sticky => 'w');
      $System_properties->Checkbutton(-text     => 'Cluster CPUs for COG analysis',
                                     -variable => \$auto_ini{COG_cluster},
                                     -command  => sub {
                                                       if ($auto_ini{COG_cluster} == 1) {
                                                          $auto_ini{COG_threaded} = 0;
                                                       } elsif ($auto_ini{COG_cluster} == 0) {
                                                          $auto_ini{COG_threaded} = 1;
                                                       }
                                                      }
                                   )
         ->grid (-row => 2, -column => 2, -sticky => 'w');
      $System_properties->Checkbutton(-text     => 'Multi-thread COG analysis',
                                      -variable => \$auto_ini{COG_threaded},
                                      -command  => sub {
                                                         if ($auto_ini{COG_threaded} == 1) {
                                                            $auto_ini{COG_cluster} = 0;
                                                         } elsif ($auto_ini{COG_threaded} == 0) {
                                                            $auto_ini{COG_cluster} = 1;
                                                         }
                                                       }
                                    )
         ->grid (-row => 3, -column => 2, -sticky => 'w');
      $System_properties->BrowseEntry(-label    => "CPU limit for HMMER analyses\n'0' no limits",
                                      -choices  => [0..100],
                                      -width    => 5,
                                      -variable => \$auto_ini{HMMER_limit},
                                   )
         ->grid (-row => 4, -column => 1, -sticky => 'w');
      $System_properties->Label(-text => 'Post analysis options')
         ->grid (-row => 0, -column => 3);
      $System_properties->Checkbutton(-text     => 'Sort results into separate folders',
                                      -variable => \$auto_ini{sort_results},
                                    )
         ->grid (-row => 1, -column => 3, -sticky => 'w');

      my $compress_data = $System_properties->LabEntry(-label        => '    Project name',
                                                       -labelPack    => [-side => "left", -anchor => "w"],
                                                       -justify      => 'left',
                                                       -relief       => 'sunken',
                                                       -width        => 20,
                                                       -borderwidth  => 2,
                                                       -takefocus    => '1',
                                                       -textvariable => \$auto_ini{project_name}
                                                      )
         ->grid (-row => 3, -column => 3, -sticky => 'w');
      $System_properties->Checkbutton(-text     => 'Compress results',
                                      -variable => \$auto_ini{compress_results},
                                      -command  => sub {
                                                         if ($auto_ini{compress_results} == 1) {
                                                            $compress_data->configure(-state => 'normal');
                                                         } elsif ($auto_ini{compress_results} == 0) {
                                                            $compress_data->configure(-state => 'disabled');
                                                         }
                                                       }
                                    )
         ->grid (-row => 2, -column => 3, -sticky => 'w');

      $System_properties->gridColumnconfigure(0, -pad => 50);
      $System_properties->gridColumnconfigure(2, -pad => 50);
      $System_properties->gridRowconfigure   (0, -pad => 20);

      #define default setup
      if ($auto_ini{compress_results} == 1) {
         $compress_data->configure(-state => 'normal');
      } elsif ($auto_ini{compress_results} == 0) {
         $compress_data->configure(-state => 'disabled');
      }
   }
   {
      $filter_blast = $Blast_properties->LabEntry(-label        => 'Filter Blast results',
                                                  -labelPack    => [-side => "left", -anchor => "w"],
                                                  -justify      => 'left',
                                                  -relief       => 'sunken',
                                                  -width        => 40,
                                                  -borderwidth  => 2,
                                                  -takefocus    => '1',
                                                  -textvariable => \$auto_ini{filter_blast}
                                                  )
         ->grid (-row => 0, -column => 0, -sticky => 'w');
      $Blast_properties->Button(-text => "Browse",
                                -command => sub {
                                                 &select_blast_filter(main_window     => \$mw,
                                                                      ini_ref         => \%ini,
                                                                      auto_ini_ref    => \%auto_ini,
                                                                      ini_comment_ref => \%ini_comment,
                                                                     );
                                                }
                               )
         ->grid (-row => 0, -column => 1, -sticky => 'w');
      $Blast_properties->Button(-text => "Clear Filter",
                                -command => sub {
                                                 $auto_ini{filter_blast} = '';
                                                }
                               )
         ->grid (-row => 0, -column => 2, -sticky => 'w');
      $Blast_properties->BrowseEntry(-label    =>'Blast entries listed',
                                     -choices  =>[1..200],
                                     -width    => 10,
                                     -variable =>\$auto_ini{blast_entry_number}
                                    )
         ->grid (-row => 1, -column => 0, -sticky => 'w', -columnspan => '2');
      $Blast_properties->BrowseEntry(-label   =>'Translation table   ',
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
                                    -width    => 40,
                                    -variable =>\$auto_ini{translation_table}
                                   )
         ->grid (-row => 2, -column => 0, -sticky => 'w', -columnspan => 2);
      $allow_cog_button = $Blast_properties->Checkbutton(-text       => 'Allow COG annotation in Blast results?',
                                                         -wraplength => 250,
                                                         -variable   => \$auto_ini{allow_COG_anno}
                                                       )
         ->grid (-row => 3, -column => 0, -sticky => 'w', -columnspan => 2);

      my $blast_choice = $Blast_properties->BrowseEntry(-label   => '    Number of hits displayed ',
                                                        -choices => [1..200],
                                                        -width   => 4,
                                                        -variable=> \$auto_ini{blast_list}
                                                        )
         ->grid (-row => 5, -column => 0, -sticky => 'w');
      $Blast_properties->Checkbutton(-text       => 'Create Blast Summary',
                                     -variable   => \$auto_ini{blast_summary},
                                     -command    => sub {
                                                         if ($auto_ini{blast_summary} == 1) {
                                                            $blast_choice->configure(-state      => 'normal',
                                                                                     -foreground => 'black');
                                                         } elsif ($auto_ini{blast_summary} == 0) {
                                                            $blast_choice->configure(-state      => 'disabled',
                                                                                     -foreground => 'gray');
                                                         }
                                                        }
                                  )
         ->grid (-row => 4, -column => 0, -sticky => 'w');

      #set initial state
      if ($auto_ini{blast_summary} == 1) {
         $blast_choice->configure(-state      => 'normal',
                                  -foreground => 'black');
      } elsif ($auto_ini{blast_summary} == 0) {
         $blast_choice->configure(-state      => 'disabled',
                                  -foreground => 'gray');
      }

   }

   {
      my $maintain_gene = $Genbank_options->Checkbutton(-text     => 'Maintain gene qualifier',
                                                        -variable => \$auto_ini{maintain_gene},
                                                        -command  => sub {
                                                                          if ($auto_ini{update_genbank} == 1) {
                                                                             $auto_ini{create_genbank} = 0;
                                                                          } elsif ($auto_ini{update_genbank} == 0) {
                                                                             $auto_ini{create_genbank} = 1;
                                                                          }
                                                                        }
                                                       )
         ->grid (-row => 2, -column => 1, -sticky => 'w');
      my $maintain_product = $Genbank_options->Checkbutton(-text    => 'Maintain product qualifier',
                                                           -variable => \$auto_ini{maintain_product},
                                                           -command => sub {
                                                                             if ($auto_ini{update_genbank} == 1) {
                                                                                $auto_ini{create_genbank} = 0;
                                                                             } elsif ($auto_ini{update_genbank} == 0) {
                                                                                $auto_ini{create_genbank} = 1;
                                                                             }
                                                                           }
                                                          )
         ->grid (-row => 3, -column => 1, -sticky => 'w');
      my $maintain_locus_tag = $Genbank_options->Checkbutton(-text    => 'Maintain locus_tag qualifier',
                                                             -variable => \$auto_ini{maintain_locus_tag},
                                                             -command => sub {
                                                                               if ($auto_ini{update_genbank} == 1) {
                                                                                  $auto_ini{create_genbank} = 0;
                                                                               } elsif ($auto_ini{update_genbank} == 0) {
                                                                                  $auto_ini{create_genbank} = 1;
                                                                               }
                                                                             }
                                                           )
         ->grid (-row => 4, -column => 1, -sticky => 'w');
      my $maintain_protein_id = $Genbank_options->Checkbutton(-text    => 'Maintain protein_id qualifier',
                                                              -variable => \$auto_ini{maintain_protein_id},
                                                              -command => sub {
                                                                                if ($auto_ini{update_genbank} == 1) {
                                                                                   $auto_ini{create_genbank} = 0;
                                                                                } elsif ($auto_ini{update_genbank} == 0) {
                                                                                   $auto_ini{create_genbank} = 1;
                                                                                }
                                                                              }
                                                          )
         ->grid (-row => 5, -column => 1, -sticky => 'w');
      my $add_qualifiers      = $Genbank_options->Checkbutton(-text    => 'Fill in missing gene or CDS qualifiers',
                                                              -variable => \$auto_ini{add_qualifiers},
                                                              -command => sub {
                                                                                if ($auto_ini{update_genbank} == 1) {
                                                                                   $auto_ini{create_genbank} = 0;
                                                                                } elsif ($auto_ini{update_genbank} == 0) {
                                                                                   $auto_ini{create_genbank} = 1;
                                                                                }
                                                                              }
                                                          )
         ->grid (-row => 6, -column => 1, -sticky => 'w');


      my $keep_annotation = $Genbank_options->Checkbutton(-text     => 'Keep existing gene and CDS annotation',
                                                          -variable => \$auto_ini{keep_annotation},
                                                          -command  => sub {
                                                                             if ($auto_ini{keep_annotation} == 1) {
                                                                                 $maintain_gene      ->configure(-state => 'normal');
                                                                                 $maintain_product   ->configure(-state => 'normal');
                                                                                 $maintain_locus_tag ->configure(-state => 'normal');
                                                                                 $maintain_protein_id->configure(-state => 'normal');
                                                                                 $add_qualifiers     ->configure(-state => 'normal');
                                                                             } elsif ($auto_ini{keep_annotation} == 0) {
                                                                                 $maintain_gene      ->configure(-state => 'disabled');
                                                                                 $maintain_product   ->configure(-state => 'disabled');
                                                                                 $maintain_locus_tag ->configure(-state => 'disabled');
                                                                                 $maintain_protein_id->configure(-state => 'disabled');
                                                                                 $add_qualifiers     ->configure(-state => 'disabled');
                                                                             }
                                                                           }
                                                         )
         ->grid (-row => 1, -column => 1, -sticky => 'w');
      my $update_date = $Genbank_options->Checkbutton(-text     => 'Update Genbank file date',
                                                      -variable => \$auto_ini{update_gb_date}
                                                      )
         ->grid (-row => 1, -column => 2, -sticky => 'w');
      $Genbank_options->Checkbutton(-text    => 'Update existing Genbank files',
                                    -variable => \$auto_ini{update_genbank},
                                    -command => sub {
                                                     if ($auto_ini{update_genbank} == 1) {
                                                        $auto_ini{create_genbank}  = 0;
                                                        $keep_annotation       ->configure(-state => 'normal');
                                                        $update_date           ->configure(-state => 'normal');
                                                        if ($auto_ini{keep_annotation} == 1) {
                                                           $maintain_gene      ->configure(-state => 'normal');
                                                           $maintain_product   ->configure(-state => 'normal');
                                                           $maintain_locus_tag ->configure(-state => 'normal');
                                                           $maintain_protein_id->configure(-state => 'normal');
                                                           $add_qualifiers     ->configure(-state => 'normal');
                                                        }
                                                     } elsif ($auto_ini{update_genbank} == 0) {
                                                        $auto_ini{create_genbank}  = 1;
                                                        $keep_annotation    ->configure(-state => 'disabled');
                                                        $maintain_gene      ->configure(-state => 'disabled');
                                                        $maintain_product   ->configure(-state => 'disabled');
                                                        $maintain_locus_tag ->configure(-state => 'disabled');
                                                        $maintain_protein_id->configure(-state => 'disabled');
                                                        $update_date        ->configure(-state => 'disabled');
                                                        $add_qualifiers     ->configure(-state => 'disabled');
                                                     }
                                                    }
                                   )
         ->grid (-row => 1, -column => 0, -sticky => 'w');
      $Genbank_options->Checkbutton(-text    => 'Create new Genbank files',
                                    -variable => \$auto_ini{create_genbank},
                                    -command => sub {
                                                     if ($auto_ini{create_genbank} == 1) {
                                                        $auto_ini{update_genbank} = 0;
                                                        $keep_annotation    ->configure(-state => 'disabled');
                                                        $maintain_gene      ->configure(-state => 'disabled');
                                                        $maintain_product   ->configure(-state => 'disabled');
                                                        $maintain_locus_tag ->configure(-state => 'disabled');
                                                        $maintain_protein_id->configure(-state => 'disabled');
                                                        $keep_annotation    ->configure(-state => 'disabled');
                                                        $update_date        ->configure(-state => 'disabled');
                                                        $add_qualifiers     ->configure(-state => 'disabled');
                                                     } elsif ($auto_ini{create_genbank} == 0) {
                                                        $auto_ini{update_genbank} = 1;
                                                        $keep_annotation       ->configure(-state => 'normal');
                                                        $update_date           ->configure(-state => 'normal');
                                                        if ($auto_ini{keep_annotation} == 1) {
                                                           $maintain_gene      ->configure(-state => 'normal');
                                                           $maintain_product   ->configure(-state => 'normal');
                                                           $maintain_locus_tag ->configure(-state => 'normal');
                                                           $maintain_protein_id->configure(-state => 'normal');
                                                           $keep_annotation    ->configure(-state => 'normal');
                                                           $add_qualifiers     ->configure(-state => 'normal');
                                                        }
                                                     }
                                                    }
                                   )
         ->grid (-row => 0, -column => 0, -sticky => 'w');

      #configure for first display
      if ($auto_ini{update_genbank} == 1) {
         $keep_annotation       ->configure(-state => 'normal');
         $update_date           ->configure(-state => 'normal');
         if ($auto_ini{keep_annotation} == 1) {
            $maintain_gene      ->configure(-state => 'normal');
            $maintain_product   ->configure(-state => 'normal');
            $maintain_locus_tag ->configure(-state => 'normal');
            $maintain_protein_id->configure(-state => 'normal');
            $add_qualifiers     ->configure(-state => 'normal');
         } elsif ($auto_ini{keep_annotation} == 0) {
            $maintain_gene      ->configure(-state => 'disabled');
            $maintain_product   ->configure(-state => 'disabled');
            $maintain_locus_tag ->configure(-state => 'disabled');
            $maintain_protein_id->configure(-state => 'disabled');
            $add_qualifiers     ->configure(-state => 'disabled');
         }
      } elsif ($auto_ini{update_genbank} == 0) {
         $keep_annotation    ->configure(-state => 'disabled');
         $maintain_gene      ->configure(-state => 'disabled');
         $maintain_product   ->configure(-state => 'disabled');
         $maintain_locus_tag ->configure(-state => 'disabled');
         $maintain_protein_id->configure(-state => 'disabled');
         $update_date        ->configure(-state => 'disabled');
         $add_qualifiers     ->configure(-state => 'disabled');
      }
   }

   {
      $Genbank_header_options->Checkbutton(-text     => 'Create new Genbank headers',
                                           -variable => \$auto_ini{gb_header_new},
                                           -command  => sub {
                                                             if ($auto_ini{gb_header_new} == 1) {
                                                                $auto_ini{gb_header_reuse} = 0;
                                                             } elsif ($auto_ini{gb_header_new} == 0) {
                                                                $auto_ini{gb_header_reuse} = 1;
                                                             }
                                                            }
                                          )
         ->grid (-row => 0, -column => 0, -sticky => 'w');

      $Genbank_header_options->Checkbutton(-text     => 'Remove DBLINK entry from existing Genbank headers',
                                           -variable => \$auto_ini{remove_dblink},
                                          )
         ->grid (-row => 0, -column => 1, -sticky => 'w');

      $Genbank_header_options->Checkbutton(-text     => 'Re-use existing Genbank headers',
                                           -variable => \$auto_ini{gb_header_reuse},
                                           -command  => sub {
                                                             if ($auto_ini{gb_header_reuse} == 1) {
                                                                $auto_ini{gb_header_new} = 0;
                                                             } elsif ($auto_ini{gb_header_reuse} == 0) {
                                                                $auto_ini{gb_header_new} = 1;
                                                             }
                                                            }
                                          )
         ->grid (-row => 1, -column => 0, -sticky => 'w');
      $Genbank_header_options->Button(-text     => 'Re-define default Genbank header',
                                      -command  => sub {
                                                         &gb_header(main_window    => \$mw,
                                                                    main_options   => \$main_options,
                                                                    ini_ref        => \%ini,
                                                                    auto_ini_ref   => \%auto_ini,
                                                                    );
                                                       }
                                     )
         ->grid (-row => 2, -column => 0, -sticky => 'w');
   }

   {
      $Fasta_header_options->Checkbutton(-text     => 'Create individual Genbank headers',
                                         -variable => \$auto_ini{fasta_header_individual},
                                         -command  => sub {
                                                           if ($auto_ini{fasta_header_individual} == 1) {
                                                              $auto_ini{fasta_header_default} = 0;
                                                           } elsif ($auto_ini{fasta_header_individual} == 0) {
                                                              $auto_ini{fasta_header_default} = 1;
                                                           }
                                                          }
                                        )
         ->grid (-row => 0, -column => 0, -sticky => 'w');

      $Fasta_header_options->Checkbutton(-text     => 'Use default Genbank header',
                                         -variable => \$auto_ini{fasta_header_default},
                                         -command  => sub {
                                                           if ($auto_ini{fasta_header_default} == 1) {
                                                              $auto_ini{fasta_header_individual} = 0;
                                                           } elsif ($auto_ini{fasta_header_default} == 0) {
                                                              $auto_ini{fasta_header_individual} = 1;
                                                           }
                                                          }
                                        )
         ->grid (-row => 1, -column => 0, -sticky => 'w');

      $Fasta_header_options->Button(-text     => 'Re-define default Genbank header',
                                    -command  => sub {
                                                       &gb_header(main_window    => \$mw,
                                                                  main_options   => \$main_options,
                                                                  ini_ref        => \%ini,
                                                                  auto_ini_ref   => \%auto_ini,
                                                                  );
                                                   }
                                   )
         ->grid (-row => 2, -column => 0, -sticky => 'w');
   }
}

#page2: main options Blast Pfam TIGRfam COG
{
   $main_options->Label(-text=>'Select options for annotation analyses')
      -> grid (-row => 0, -column => 0, -columnspan => 7);
   $main_options->gridRowconfigure(0, -pad => 40);
   my $internal_gm_frame = $main_options->Frame()->grid (-row => 1, -column => 3, -columnspan => 3, -sticky => 'w'); #put in place to keep internal GM display in line
   $path_widget = $main_options->Label(-text => ${$args{auto_ini_ref}}{ext_gm_folder})
      ->grid (-row => 2, -column => 3, -columnspan => 5, -sticky => 'w');
   $select_int_gm = $main_options->Button(-text => 'Select',
                                          -padx => 30,
                                          -command => sub {
                                                            &setup_gene_model(main_window    => \$mw,
                                                                              main_options   => \$main_options,
                                                                              igm_frame      => \$internal_gm_frame,
                                                                              ini_ref        => \%ini,
                                                                              auto_ini_ref   => \%auto_ini,
                                                                              critica_button => \$critica_button,
                                                                              glimmer_button => \$glimmer_button,
                                                                              gl_training    => \$glimmer_make_training);
                                                          }
                                         )
      -> grid (-row => 1, -column => 1, -sticky => 'w');
   $select_ext_gm = $main_options->Button(-text => 'Browse',
                                          -padx => 26,
                                          -command => sub {
                                                           &dir_select (main_window => \$mw,
                                                                        label       => 'Choose folder for external gene model',
                                                                        auto_ini_ref => \%auto_ini,
                                                                        ini_ref      => \%ini,
                                                                        directory    => 'ext_gm_folder');
                                                           #view current pathway for external gene models
                                                           $path_widget->configure(-text => $auto_ini{ext_gm_folder},
                                                                                  );
                                                          }
                                         )
      -> grid (-row => 2, -column => 1, -sticky => 'w');
   # determine initial state of button
   if ($auto_ini{internal_gm} == 1) {
      $select_int_gm         -> configure(-state => 'normal');
      $select_ext_gm         -> configure(-state => 'disabled');
      #view current choices for gene model selections
      &view_genemodel_options(main_window    => \$mw,
                              main_options   => \$main_options,
                              igm_frame      => \$internal_gm_frame,
                              ini_ref        => \%ini,
                              auto_ini_ref   => \%auto_ini
                             );
      #forget pathway display for external gene model
      $path_widget->configure(-text       => $auto_ini{ext_gm_folder},
                              -foreground => 'grey'
                             );
   } elsif ($auto_ini{external_gm} == 1) {
      $select_int_gm         -> configure(-state => 'disabled');
      $select_ext_gm         -> configure(-state => 'normal');
      #forget genemodel selection display
      &forget_gm_selection(main_window    => \$mw,
                           main_options   => \$main_options,
                           igm_frame      => \$internal_gm_frame,
                           ini_ref        => \%ini,
                           auto_ini_ref   => \%auto_ini
                          );
      #view current pathway for external gene models
      $path_widget->configure(-text       => $auto_ini{ext_gm_folder},
                              -foreground => 'black'
                             );
   }


   $main_options->Checkbutton(-text     => 'Determine gene model',
                              -variable => \$auto_ini{internal_gm},
                              -command  => [sub {
                                                 if ($auto_ini{internal_gm} == 1) {
                                                   $auto_ini{external_gm} = 0;
                                                   if ($auto_ini{use_critica} == 1) {
                                                      $critica_button     -> configure(-state => 'normal');
                                                   }
                                                   $glimmer_button        -> configure(-state => 'normal');
                                                   $select_int_gm         -> configure(-state => 'normal');
                                                   $select_ext_gm         -> configure(-state => 'disabled');
                                                   $glimmer_make_training -> configure(-state => 'normal');
                                                   #view current choices for gene model selections
                                                   &view_genemodel_options(main_window    => \$mw,
                                                                           main_options   => \$main_options,
                                                                           igm_frame      => \$internal_gm_frame,
                                                                           ini_ref        => \%ini,
                                                                           auto_ini_ref   => \%auto_ini
                                                                          );
                                                   #forget pathway display for external gene model
                                                   $path_widget->configure(-text => $auto_ini{ext_gm_folder},
                                                                           -foreground => 'grey'
                                                                          );
                                                } elsif ($auto_ini{internal_gm} == 0) {
                                                   $auto_ini{external_gm} = 1;
                                                   $critica_button        -> configure(-state => 'disabled');
                                                   $glimmer_button        -> configure(-state => 'disabled');
                                                   $select_int_gm         -> configure(-state => 'disabled');
                                                   $select_ext_gm         -> configure(-state => 'normal');
                                                   $glimmer_make_training -> configure(-state => 'disabled');
                                                   #forget genemodel selection display
                                                   &forget_gm_selection(main_window    => \$mw,
                                                                        main_options   => \$main_options,
                                                                        igm_frame      => \$internal_gm_frame,
                                                                        ini_ref        => \%ini,
                                                                        auto_ini_ref   => \%auto_ini
                                                                       );
                                                   #view current pathway for external gene models
                                                   $path_widget->configure(-text => $auto_ini{ext_gm_folder},
                                                                           -foreground => 'black'
                                                                          );
                                                }
                                               }]
                              )
      ->grid (-row => 1, -column => 0, -sticky => 'w');
   $main_options->Checkbutton(-text     => 'Use external gene model',
                               -variable => \$auto_ini{external_gm},
                               -command => [sub {
                                                 if ($auto_ini{external_gm} == 1) {
                                                    $auto_ini{internal_gm} = 0;
                                                    $critica_button        -> configure(-state => 'disabled');
                                                    $glimmer_button        -> configure(-state => 'disabled');
                                                    $select_int_gm         -> configure(-state => 'disabled');
                                                    $select_ext_gm         -> configure(-state => 'normal');
                                                    $glimmer_make_training -> configure(-state => 'disabled');
                                                    #forget genemodel selection display
                                                    &forget_gm_selection(main_window    => \$mw,
                                                                         main_options   => \$main_options,
                                                                         igm_frame      => \$internal_gm_frame,
                                                                         ini_ref        => \%ini,
                                                                         auto_ini_ref   => \%auto_ini
                                                                        );
                                                    #view current pathway for external gene models
                                                    $path_widget->configure(-text       => $auto_ini{ext_gm_folder},
                                                                            -foreground => 'black'
                                                                           );
                                                 } elsif ($auto_ini{external_gm} == 0) {
                                                    $auto_ini{internal_gm} = 1;
                                                    if ($auto_ini{use_critica} == 1) {
                                                      $critica_button      ->configure(-state => 'normal');
                                                    }
                                                    $glimmer_button        -> configure(-state => 'normal');
                                                    $select_int_gm         -> configure(-state => 'normal');
                                                    $select_ext_gm         -> configure(-state => 'disabled');
                                                    $glimmer_make_training -> configure(-state => 'normal');
                                                    #view current choices for gene model selections
                                                    &view_genemodel_options(main_window    => \$mw,
                                                                            main_options   => \$main_options,
                                                                            igm_frame      => \$internal_gm_frame,
                                                                            ini_ref        => \%ini,
                                                                            auto_ini_ref   => \%auto_ini
                                                                           );
                                                    #forget pathway display for external gene model
                                                    $path_widget->configure(-text       => $auto_ini{ext_gm_folder},
                                                                            -foreground => 'grey'
                                                                           );
                                                }
                                               }]
                              )
      ->grid (-row => 2, -column => 0, -sticky => 'w');


   $main_options->Label(-text => 'Current choices:')
      ->grid (-row => 1, -column => 2, -sticky => 'w');
   $main_options->Label(-text => 'Current folder:')
      ->grid (-row => 2, -column => 2, -sticky => 'w');



   $main_options->Checkbutton(-text     => 'Blast',
                               -variable => \$auto_ini{blast_selector},
                               -command => sub {
                                                 if ($auto_ini{blast_selector} == 1) {
                                                    $blast_button     -> configure (-state => 'normal');
                                                    $Blast_add        -> configure (-state => 'normal');
                                                    $Blast_clear      -> configure (-state => 'normal');
                                                    $blast_sel        -> configure (-state => 'normal');
                                                    $filter_blast     -> configure (-state => 'normal');
                                                    $allow_cog_button -> configure (-state => 'normal');
                                                    $gene_model_only  -> configure (-state => 'disabled');
                                                    $auto_ini{gene_model_only} = 0;
                                                 } elsif ($auto_ini{blast_selector} == 0) {
                                                    $blast_button     -> configure (-state => 'disabled');
                                                    $Blast_add        -> configure (-state => 'disabled');
                                                    $Blast_clear      -> configure (-state => 'disabled');
                                                    $blast_sel        -> configure (-state => 'disabled');
                                                    $filter_blast     -> configure (-state => 'disabled');
                                                    $allow_cog_button -> configure (-state => 'disabled');
                                                    $gene_model_only  -> configure (-state => 'normal');
                                                    $auto_ini{gene_model_only} = 1;
                                                 }
                                               })                              ->grid (-row => 3, -column => 0, -sticky => 'w');

   $blast_sel = $main_options->BrowseEntry(-choices  => ['BlastN', 'BlastX', 'tBlastX', 'BlastP', 'gappedBlastP', 'PSI-Blast'],
                                           -width    => 11,
                                           -variable => \$auto_ini{blast_type},
                                           -command  => sub {
                                                             if ($auto_ini{blast_type} =~ m/(BlastN|tBlastX)/) {
                                                                unless (-e $auto_ini{full_blast_db}        ||
                                                                        -e $auto_ini{full_blast_db}.'.nhr' ||
                                                                        -e $auto_ini{full_blast_db}.'.nal') {
                                                                   my $error_msg = $mw->Dialog(-title   => 'Error',
                                                                                               -text    => "Selected Blast db $auto_ini{blast_db} does not appear to be a nucleotide database.\n".
                                                                                                           "Please select another database.",
                                                                                               -buttons => ['OK'],
                                                                                               -bitmap  => 'info'
                                                                                              );
                                                                   $error_msg->after($auto_ini{msg_timeout}, sub {$error_msg->Exit()});
                                                                   $error_msg-> Show();
                                                                   &select_blast_db(main_window    => \$mw,
                                                                                    ini_ref        => \%ini,
                                                                                    auto_ini_ref   => \%auto_ini,
                                                                                    critica_select => 0
                                                                                   );
                                                                }
                                                             } else {
                                                                unless (-e $auto_ini{full_blast_db}        ||
                                                                        -e $auto_ini{full_blast_db}.'.phr' ||
                                                                        -e $auto_ini{full_blast_db}.'.pal') {
                                                                   my $error_msg = $mw->Dialog(-title   => 'Error',
                                                                                               -text    => "Selected Blast db $auto_ini{blast_db} does not appear to be a protein database.\n".
                                                                                                           "Please select another database.",
                                                                                               -buttons => ['OK'],
                                                                                               -bitmap  => 'info'
                                                                                              );
                                                                   $error_msg->after($auto_ini{msg_timeout}, sub {$error_msg->Exit()});
                                                                   $error_msg-> Show();
                                                                   &select_blast_db(main_window    => \$mw,
                                                                                    ini_ref        => \%ini,
                                                                                    auto_ini_ref   => \%auto_ini,
                                                                                    critica_select => 0
                                                                                   );
                                                                }
                                                             }
                                                            }
                                                ) ->grid (-row => 3, -column => 1, -sticky => 'w');

   $gene_model_only = $main_options->Checkbutton(-text     => 'Enter Gene model only',
                                                 -variable => \$auto_ini{gene_model_only},
                                                 -command => sub {
                                                                   if ($auto_ini{gene_model_only} == 0) {
                                                                      $blast_button     -> configure (-state => 'normal');
                                                                      $blast_sel        -> configure (-state => 'normal');
                                                                      $filter_blast     -> configure (-state => 'normal');
                                                                      $allow_cog_button -> configure (-state => 'normal');
                                                                      $Blast_add        -> configure (-state => 'normal');
                                                                      $Blast_clear      -> configure (-state => 'normal');
                                                                      $auto_ini{blast_selector} = 1;
                                                                   } elsif ($auto_ini{gene_model_only} == 1) {
                                                                      $blast_button     -> configure (-state => 'disabled');
                                                                      $blast_sel        -> configure (-state => 'disabled');
                                                                      $filter_blast     -> configure (-state => 'disabled');
                                                                      $allow_cog_button -> configure (-state => 'disabled');
                                                                      $Blast_add        -> configure (-state => 'disabled');
                                                                      $Blast_clear      -> configure (-state => 'disabled');
                                                                      $auto_ini{blast_selector} = 0;
                                                                   }
                                                                  }
                                                ) ->grid (-row => 3, -column => 2, -sticky => 'w');
   $blast_legacy    = $main_options->Checkbutton(-text     => 'Use Legacy Blast',
                                                 -variable => \$auto_ini{legacy_blast},
                                                 -command => sub {
                                                                  if ($auto_ini{legacy_blast} == 1) {
                                                                     $auto_ini{blast_plus} = 0;
                                                                     $blast_plus -> configure (-state => 'disabled');
                                                                  } elsif ($auto_ini{legacy_blast} == 0) {
                                                                     $auto_ini{blast_plus} = 1;
                                                                     $blast_plus -> configure (-state => 'normal');
                                                                  }
                                                                 }
                                                ) ->grid (-row => 3, -column => 3, -sticky => 'w');
   $blast_plus      = $main_options->Checkbutton(-text     => 'Use Blast Plus',
                                                 -variable => \$auto_ini{blast_plus},
                                                 -command => sub {
                                                                  if ($auto_ini{blast_plus} == 1) {
                                                                     $blast_legacy -> configure (-state => 'disabled');
                                                                     $auto_ini{legacy_blast} = 0;
                                                                  } elsif ($auto_ini{blast_plus} == 0) {
                                                                     $blast_legacy -> configure (-state => 'normal');
                                                                     $auto_ini{legacy_blast} = 1;
                                                                  }
                                                                 }
                                                ) ->grid (-row => 4, -column => 3, -sticky => 'w');

   $main_options->Checkbutton(-text     => 'COG',
                               -variable => \$auto_ini{COG_selector},
                               -command => sub {
                                                 if ($auto_ini{COG_selector} == 1) {
                                                    $COG_button   -> configure (-state => 'normal');
                                                    $COG_sel      -> configure (-state => 'normal');
                                                 } elsif ($auto_ini{COG_selector} == 0) {
                                                    $COG_button   -> configure (-state => 'disabled');
                                                    $COG_sel      -> configure (-state => 'disabled');
                                                 }
                                               })     ->grid (-row => 5, -column => 0, -sticky => 'w');
   $main_options->Checkbutton(-text     => 'PFam',
                               -variable => \$auto_ini{Pfam_selector},
                               -command => sub {
                                                 if ($auto_ini{Pfam_selector} == 1) {
                                                    $Pfam_fam_sel -> configure (-state      => 'normal',
                                                                                -foreground => 'black');
                                                    $Pfam_button  -> configure (-state      => 'normal');
                                                    $Pfam_sel     -> configure (-state      => 'normal');
                                                    $Pfam_add     -> configure (-state      => 'normal');
                                                    $Pfam_clear   -> configure (-state      => 'normal');
                                                    $Pfam_version2-> configure (-state      => 'normal');
                                                    $Pfam_version3-> configure (-state      => 'normal');
                                                    #if ($auto_ini{use_Pfam3} == 1) {
                                                    $PFam3_thrshld-> configure (-state      => 'normal');
                                                    #} elsif ($auto_ini{use_Pfam3} == 0) {
                                                    #$PFam3_thrshld-> configure (-state      => 'disabled');
                                                    #}
                                                 } elsif ($auto_ini{Pfam_selector} == 0) {
                                                    $Pfam_fam_sel -> configure (-state      => 'disabled',
                                                                                -foreground => 'gray');
                                                    $Pfam_button  -> configure (-state      => 'disabled');
                                                    $Pfam_sel     -> configure (-state      => 'disabled');
                                                    $Pfam_add     -> configure (-state      => 'disabled');
                                                    $Pfam_clear   -> configure (-state      => 'disabled');
                                                    $Pfam_version2-> configure (-state      => 'disabled');
                                                    $Pfam_version3-> configure (-state      => 'disabled');
                                                    $PFam3_thrshld-> configure (-state      => 'disabled');
                                                 }
                                                 #also test for selected Pfam db
                                                 if ($auto_ini{Pfam_selector} == 1) {
                                                    my $status = &test_pfam_db_selection(main_window    => \$mw,
                                                                                         main_options   => \$main_options,
                                                                                         ini_ref        => \%ini,
                                                                                         auto_ini_ref   => \%auto_ini
                                                                                        );
                                                    if ($status eq 'Pfam2error') {
                                                       $mw->messageBox(-title   => 'Potential database conflict',
                                                                       -message => 'It appears a Pfam3 database is selected for Hmmer2',
                                                                       -type    => 'OK',
                                                                       -icon    => 'info');
                                                    }
                                                    if ($status eq 'Pfam3error') {
                                                       $mw->messageBox(-title   => 'Potential database conflict',
                                                                       -message => 'It appears a non-Pfam3 database is selected for Hmmer3',
                                                                       -type    => 'OK',
                                                                       -icon    => 'info');
                                                    }
                                                 }


                                               })    ->grid (-row => 6, -column => 0, -sticky => 'w');
   $Pfam_sel     = $main_options->BrowseEntry(-choices  => ['brief', 'verbose', 'InterPro', 'both'],
                                              -width    => 11,
                                              -variable => \$auto_ini{pfam_verbose})    -> grid (-row => 6, -column => 1, -sticky => 'w');
   $Pfam_fam_sel = $main_options->BrowseEntry(-label    => 'No. of Families per ORF',
                                              -choices  => [1..20],
                                              -width    => 3,
                                              -variable => \$auto_ini{Pfam_family_max}) -> grid (-row => 6, -column => 2, -sticky => 'w');
   $Pfam_version2= $main_options->Checkbutton(-text     => 'Pfam2',
                                              -variable => \$auto_ini{use_Pfam2},
                                              -command  => sub{
                                                                if ($auto_ini{use_Pfam2} == 1) {
                                                                   $auto_ini{use_Pfam3} = 0;
                                                                   #$PFam3_thrshld-> configure (-state      => 'disabled');
                                                                } elsif ($auto_ini{use_Pfam2} == 0) {
                                                                   $auto_ini{use_Pfam3} = 1;
                                                                   #$PFam3_thrshld-> configure (-state      => 'normal');
                                                                }
                                                                #also test for selected Pfam db
                                                                if ($auto_ini{Pfam_selector} == 1) {
                                                                   my $status = &test_pfam_db_selection(main_window    => \$mw,
                                                                                                        main_options   => \$main_options,
                                                                                                        ini_ref        => \%ini,
                                                                                                        auto_ini_ref   => \%auto_ini
                                                                                                       );
                                                                   if ($status eq 'Pfam2error') {
                                                                      $mw->messageBox(-title   => 'Potential database conflict',
                                                                                      -message => 'It appears a Pfam3 database is selected for Hmmer2',
                                                                                      -type    => 'OK',
                                                                                      -icon    => 'info');
                                                                   }
                                                                   if ($status eq 'Pfam3error') {
                                                                      $mw->messageBox(-title   => 'Potential database conflict',
                                                                                      -message => 'It appears a non-Pfam3 database is selected for Hmmer3',
                                                                                      -type    => 'OK',
                                                                                      -icon    => 'info');
                                                                   }
                                                                }
                                                              })                        -> grid (-row => 6, -column => 3, -sticky => 'w');
   $Pfam_version3= $main_options->Checkbutton(-text     => 'Pfam3',
                                              -variable => \$auto_ini{use_Pfam3},
                                              -command  => sub{
                                                                if ($auto_ini{use_Pfam3} == 1) {
                                                                    $auto_ini{use_Pfam2} = 0;
                                                                    #$PFam3_thrshld-> configure (-state      => 'normal');
                                                                } elsif ($auto_ini{use_Pfam3} == 0) {
                                                                    $auto_ini{use_Pfam2} = 1;
                                                                    #$PFam3_thrshld-> configure (-state      => 'disabled');
                                                                }
                                                                if ($auto_ini{Pfam_selector} == 1) {
                                                                   my $status = &test_pfam_db_selection(main_window    => \$mw,
                                                                                                        main_options   => \$main_options,
                                                                                                        ini_ref        => \%ini,
                                                                                                        auto_ini_ref   => \%auto_ini
                                                                                                       );
                                                                   if ($status eq 'Pfam2error') {
                                                                      $mw->messageBox(-title   => 'Potential database conflict',
                                                                                      -message => 'It appears a Pfam3 database is selected for Hmmer2',
                                                                                      -type    => 'OK',
                                                                                      -icon    => 'info');
                                                                   }
                                                                   if ($status eq 'Pfam3error') {
                                                                      $mw->messageBox(-title   => 'Potential database conflict',
                                                                                      -message => 'It appears a non-Pfam3 database is selected for Hmmer3',
                                                                                      -type    => 'OK',
                                                                                      -icon    => 'info');
                                                                   }
                                                                }
                                                              })                        -> grid (-row => 6, -column => 4, -sticky => 'w');
   $PFam3_thrshld= $main_options->Checkbutton(-text     => 'Restrict individual domains by global threshold',
                                              -variable => \$auto_ini{Pfam_restrict_domains}
                                              )                                         -> grid (-row => 7, -column => 4, -sticky => 'w');

   $main_options->Checkbutton(-text     => 'TIGRfam',
                               -variable => \$auto_ini{TIGRfam_selector},
                               -command => sub {
                                                 if ($auto_ini{TIGRfam_selector} == 1) {
                                                    $TIGRfam_fam_sel -> configure (-state      => 'normal',
                                                                                   -foreground => 'black');
                                                    $TIGRfam_button  -> configure (-state      => 'normal');
                                                    $TIGRfam_sel     -> configure (-state      => 'normal');
                                                    $TIGRfam_add     -> configure (-state      => 'normal');
                                                    $TIGRfam_clear   -> configure (-state      => 'normal');
                                                    $tigr_anno       -> configure (-state      => 'normal');
                                                    $TIGR_EConly     -> configure (-state      => 'normal');
                                                    $TIGR_ECCDS      -> configure (-state      => 'normal');
                                                    $TIGR_ECdiscard  -> configure (-state      => 'normal');
                                                 } elsif ($auto_ini{TIGRfam_selector} == 0) {
                                                    $TIGRfam_fam_sel -> configure (-state      => 'disabled',
                                                                                   -foreground => 'gray');
                                                    $TIGRfam_button  -> configure (-state      => 'disabled');
                                                    $TIGRfam_sel     -> configure (-state      => 'disabled');
                                                    $TIGRfam_add     -> configure (-state      => 'disabled');
                                                    $TIGRfam_clear   -> configure (-state      => 'disabled');
                                                    $tigr_anno       -> configure (-state      => 'disabled');
                                                    $TIGR_EConly     -> configure (-state      => 'disabled');
                                                    $TIGR_ECCDS      -> configure (-state      => 'disabled');
                                                    $TIGR_ECdiscard  -> configure (-state      => 'disabled');
                                                 }
                                               }) ->grid (-row => 8, -column => 0, -sticky => 'w');
   $TIGRfam_sel     = $main_options->BrowseEntry(-choices  => ['verbose + GO', 'brief'],
                                                 -width    => 11,
                                                 -variable => \$auto_ini{tigrfam_verbose}) -> grid (-row => 8, -column => 1, -sticky => 'w');

   $TIGRfam_fam_sel = $main_options->BrowseEntry(-label    => 'No. of Families per ORF',
                                                 -choices  => [1..20],
                                                 -width    => 3,
                                                 -variable => \$auto_ini{TIGRfam_family_max}) -> grid (-row => 8, -column => 2, -sticky => 'w');
   $tigr_anno =   $main_options->Checkbutton(-text     => 'Use TIGRfams for annotation',
                                             -variable => \$auto_ini{TIGR_annotation},
                                             -command  => sub {
                                                                if ($auto_ini{TIGR_annotation} == 1) {
                                                                   $TIGR_EConly    -> configure (-state => 'normal');
                                                                   $TIGR_ECCDS     -> configure (-state => 'normal');
                                                                   $TIGR_ECdiscard -> configure (-state => 'normal');
                                                                } elsif ($auto_ini{TIGR_annotation} == 0) {
                                                                   $TIGR_EConly    -> configure (-state => 'disabled');
                                                                   $TIGR_ECCDS     -> configure (-state => 'disabled');
                                                                   $TIGR_ECdiscard -> configure (-state => 'disabled');
                                                                }
                                                              }
                                             )   -> grid (-row => 9, -column => 2, -sticky => 'w');
   $TIGR_EConly = $main_options->Checkbutton(-text     => 'Transfer EC numbers only',
                                             -variable => \$auto_ini{TIGR_EConly},
                                             -command  => sub {
                                                                if ($auto_ini{TIGR_EConly} == 1) {
                                                                   $auto_ini{TIGR_ECCDS} = 0;
                                                                } elsif ($auto_ini{TIGR_EConly} == 0) {
                                                                   $auto_ini{TIGR_ECCDS} = 1;
                                                                }
                                                              }
                                             )   -> grid (-row => 10, -column => 3, -sticky => 'w', -columnspan => 4);
   $TIGR_ECCDS = $main_options->Checkbutton(-text     => 'Transfer EC numbers and gene annotation',
                                            -variable => \$auto_ini{TIGR_ECCDS},
                                            -command  => sub {
                                                                if ($auto_ini{TIGR_ECCDS} == 1) {
                                                                   $auto_ini{TIGR_EConly} = 0;
                                                                } elsif ($auto_ini{TIGR_ECCDS} == 0) {
                                                                   $auto_ini{TIGR_EConly} = 1;
                                                                }
                                                              }
                                           )   -> grid (-row => 9, -column => 3, -sticky => 'w', -columnspan => 4);
   $TIGR_ECdiscard = $main_options->Checkbutton(-text     => 'Overwrite existing EC number entries',
                                                -variable => \$auto_ini{TIGR_ECdiscard},
                                               )   -> grid (-row => 11, -column => 3, -sticky => 'w', -columnspan => 4);
}

#page3: database options
{
   $glimmer_button = $database_options->Label( -text     => 'Glimmer model',
                                               -relief   => 'flat',
                                               -width    => '18',
                                             )                            -> grid (-row => 1, -column => 0, -sticky => 'w');
   $database_options->Label(-justify      => 'left',
                            -relief       => 'groove',
                            -borderwidth  => 2,
                            -wraplength   => 170,
                            -textvariable => \$auto_ini{gl_short_file},
                            -width        => 30)                          -> grid (-row => 1, -column => 1, -sticky => 'w');

   $glimmer_add = $database_options->Button(-text     => 'Add training file',
                                            -width    => '18',
                                            -command  => [sub {
                                                               &select_glimmer_model(main_window    => \$mw,
                                                                                     ini_ref        => \%ini,
                                                                                     auto_ini_ref   => \%auto_ini,);
                                                              }])                             -> grid (-row => 1, -column => 2, -sticky => 'w');

   $glimmer_make_training = $database_options -> Checkbutton (-text     => 'Make Glimmer training files on the fly',
                                                              -variable => \$auto_ini{make_training_file},
                                                              -command  => sub {
                                                                                   if ($auto_ini{make_training_file} == 1) {
                                                                                      $glimmer_add    ->configure(-state => 'disabled');
                                                                                   } elsif ($auto_ini{make_training_file} == 0) {
                                                                                      $glimmer_add    ->configure(-state => 'normal');
                                                                                   }
                                                                               }
                                                             )                                -> grid (-row => 1, -column => 3, -sticky => 'w', -columnspan => 2);
   $critica_button = $database_options->Label (-text     => 'Critica database',
                                               -width    => '18',
                                               -relief   => 'flat',
                                              )                           -> grid (-row => 2, -column => 0, -sticky => 'w');
   $critica_add = $database_options->Button(-text     => 'Add Critica db',
                                            -width    => '18',
                                            -command  => [sub {
                                                            &select_blast_db(main_window    => \$mw,
                                                                             ini_ref        => \%ini,
                                                                             auto_ini_ref   => \%auto_ini,
                                                                             critica_select => 1
                                                                             );
                                                           }])                               -> grid (-row => 2, -column => 2, -sticky => 'w');

   $critica_clear = $database_options->Button(-text     => 'Clear selection',
                                              -width    => '18',
                                              -command  => [sub {
                                                                  &clear_db_selection(main_window     => \$mw,
                                                                                      ini_ref         => \%ini,
                                                                                      auto_ini_ref    => \%auto_ini,
                                                                                      short_db        => 'selected_critica_db',
                                                                                      full_db         => 'full_critica_db',
                                                                                     );
                                                               }])                           -> grid ( -row => 2, -column => 3, -sticky => 'w');

   $database_options->Label(-justify      => 'left',
                            -relief       => 'groove',
                            -borderwidth  => 2,
                            -wraplength   => 170,
                            -textvariable => \$auto_ini{selected_critica_db},
                            -width        => 30)                                             -> grid (-row => 2, -column => 1, -sticky => 'w');

   $blast_button = $database_options->Label (-text     => 'Blast db',
                                             -relief   => 'flat',
                                             -width    => '18',
                                            )                                                -> grid (-row => 3, -column => 0, -sticky => 'w');
   $database_options->Label(-justify      => 'left',
                            -relief       => 'groove',
                            -borderwidth  => 2,
                            -wraplength   => 170,
                            -textvariable => \$auto_ini{blast_db},
                            -width        => 30)                                             -> grid (-row => 3, -column => 1, -sticky => 'w');
   $Blast_add = $database_options->Button(-text    => 'Add Blast db',
                                          -width   => '18',
                                          -command => [sub {
                                                            &select_blast_db(main_window    => \$mw,
                                                                             ini_ref        => \%ini,
                                                                             auto_ini_ref   => \%auto_ini
                                                                             );
                                                           }])                               -> grid (-row => 3, -column => 2, -sticky => 'w');

   $Blast_clear = $database_options->Button(-text    => 'Clear selection',
                                            -width   => '18',
                                            -command => [sub {
                                                               &clear_db_selection(main_window     => \$mw,
                                                                                   ini_ref         => \%ini,
                                                                                   auto_ini_ref    => \%auto_ini,
                                                                                   short_db        => 'blast_db',
                                                                                   full_db         => 'full_blast_db',
                                                                                  );
                                                            }])                             -> grid ( -row => 3, -column => 3, -sticky => 'w');

   $COG_button = $database_options->Label (-text     => 'COG db',
                                           -relief   => 'flat',
                                           -width    => '18',
                                          )                                                 -> grid (-row => 4, -column => 0, -sticky => 'w');

   $database_options->Label(-justify      => 'left',
                            -relief       => 'groove',
                            -borderwidth  => 2,
                            -wraplength   => 170,
                            -textvariable => \$auto_ini{COG_db},
                            -width        => 30)                                            -> grid (-row => 4, -column => 1, -sticky => 'w');

   $COG_sel = $database_options->BrowseEntry(-choices  => ['COG2003', 'COG2008', 'COG2014', 'arCOG', 'arCOG2014', 'POG2013'],
                                             -width    => 11,
                                             -variable => \$auto_ini{COG_type},
                                             -command  => sub {
                                                                #set notification for expectation of COGs and check for assumed db
                                                                my $msg = $mw->Dialog(-title   => 'Note',
                                                                                      -text    => $cog_type{$auto_ini{COG_type}}->{'info'},
                                                                                      -bitmap  => 'info',
                                                                                      -buttons => ['ok']);
                                                                $msg->after($auto_ini{msg_timeout}, sub {$msg->Exit()});
                                                                $msg-> Show();

                                                                #make a quick attempt at auto-selecting the correct COG database in the standard directory
                                                                my $COG_temp_dir = $ini{COG_db_path};
                                                                $COG_temp_dir =~ s/\/\w+$//; #travel up one directory to find a potential parent directory
                                                                if (-e $ini{COG_db_path}.'/'.$cog_type{$auto_ini{COG_type}}->{'file'}) {
                                                                   $auto_ini{full_COG_db} = $ini{COG_db_path}.'/'.$cog_type{$auto_ini{COG_type}}->{'db'};
                                                                   $auto_ini{COG_db}      = $cog_type{$auto_ini{COG_type}}->{'db'};
                                                                   $mw->update;
                                                                } elsif (-e $ini{COG_db_path}.'/'.$cog_type{$auto_ini{COG_type}}->{'path'}.'/'.$cog_type{$auto_ini{COG_type}}->{'file'}) {
                                                                   $auto_ini{full_COG_db} = $ini{COG_db_path}.'/'.$cog_type{$auto_ini{COG_type}}->{'path'}.'/'.$cog_type{$auto_ini{COG_type}}->{'db'};
                                                                   $auto_ini{COG_db}      = $cog_type{$auto_ini{COG_type}}->{'db'};
                                                                   $mw->update;
                                                                } elsif (-e $COG_temp_dir.'/'.$cog_type{$auto_ini{COG_type}}->{'path'}.'/'.$cog_type{$auto_ini{COG_type}}->{'file'}) {
                                                                   $ini{COG_db_path}      = $COG_temp_dir.'/'.$cog_type{$auto_ini{COG_type}}->{'path'};
                                                                   $auto_ini{full_COG_db} = $COG_temp_dir.'/'.$cog_type{$auto_ini{COG_type}}->{'path'}.'/'.$cog_type{$auto_ini{COG_type}}->{'db'};
                                                                   $auto_ini{COG_db}      = $cog_type{$auto_ini{COG_type}}->{'db'};
                                                                   $mw->update;
                                                                } else {
                                                                    &select_cog_db (main_window    => \$mw,
                                                                                    ini_ref        => \%ini,
                                                                                    auto_ini_ref   => \%auto_ini
                                                                                    );
                                                                }
                                                             }
                                                    )                                       ->grid (-row => 4, -column => 2, -sticky => 'w');

   $Pfam_button = $database_options->Label (-text     => 'Pfam db',
                                            -relief   => 'flat',
                                            -width    => '18',
                                           )                                                -> grid (-row => 6, -column => 0, -sticky => 'w');
   $database_options->Label(-justify     => 'left',
                           -relief       => 'groove',
                           -borderwidth  => 2,
                           -wraplength   => 170,
                           -textvariable => \$auto_ini{Pfam_db},
                           -width        => 30)                                             -> grid (-row => 6, -column => 1, -sticky => 'w');
   $Pfam_add = $database_options->Button(-text     => 'Add Pfam db',
                                         -width    => '18',
                                         -command => [sub {
                                                            &select_pfam_db(main_window    => \$mw,
                                                                            ini_ref        => \%ini,
                                                                            auto_ini_ref   => \%auto_ini
                                                                            );
                                                            #also test for selected Pfam db
                                                            if ($auto_ini{Pfam_selector} == 1) {
                                                               my $status = &test_pfam_db_selection(main_window    => \$mw,
                                                                                                    main_options   => \$main_options,
                                                                                                    ini_ref        => \%ini,
                                                                                                    auto_ini_ref   => \%auto_ini
                                                                                                     );
                                                               if ($status eq 'Pfam2error') {
                                                                  $mw->messageBox(-title   => 'Potential database conflict',
                                                                                  -message => 'It appears a Pfam3 database is selected for Hmmer2',
                                                                                  -type    => 'OK',
                                                                                  -icon    => 'info');
                                                               }
                                                               if ($status eq 'Pfam3error') {
                                                                  $mw->messageBox(-title   => 'Potential database conflict',
                                                                                  -message => 'It appears a non-Pfam3 database is selected for Hmmer3',
                                                                                  -type    => 'OK',
                                                                                  -icon    => 'info');
                                                               }
                                                            }
                                                           }])                              -> grid (-row => 6, -column => 2, -sticky => 'w');

   $Pfam_clear = $database_options->Button(-text    => 'Clear selection',
                                           -width   => '18',
                                           -command => [sub {
                                                               &clear_db_selection(main_window     => \$mw,
                                                                                   ini_ref         => \%ini,
                                                                                   auto_ini_ref    => \%auto_ini,
                                                                                   short_db        => 'Pfam_db',
                                                                                   full_db         => 'full_Pfam_db',
                                                                                  );
                                                            }])                            -> grid ( -row => 6, -column => 3, -sticky => 'w');

   $Rfam_button = $database_options->Label (-text     => 'Rfam db',
                                            -relief   => 'flat',
                                            -width    => '18',
                                            )                                              -> grid (-row => 7, -column => 0, -sticky => 'w');
   $database_options->Label(-justify      => 'left',
                            -relief       => 'groove',
                            -borderwidth  => 2,
                            -wraplength   => 170,
                            -textvariable => \$auto_ini{Rfam_db},
                            -width        => 30)                                           -> grid (-row => 7, -column => 1, -sticky => 'w');
   $Rfam_add = $database_options->Button(-text     => 'Add Rfam db',
                                         -width    => '18',
                                         -command  => [sub {
                                                            &select_rfam_db(main_window    => \$mw,
                                                                            ini_ref        => \%ini,
                                                                            auto_ini_ref   => \%auto_ini
                                                                            );
                                                           }])                             -> grid (-row => 7, -column => 2, -sticky => 'w');

   $Rfam_clear = $database_options->Button(-text     => 'Clear selection',
                                           -width    => '18',
                                           -command  => [sub {
                                                              &clear_db_selection(main_window     => \$mw,
                                                                                  ini_ref         => \%ini,
                                                                                  auto_ini_ref    => \%auto_ini,
                                                                                  short_db        => 'Rfam_db',
                                                                                  full_db         => 'full_Rfam_db',
                                                                                 );
                                                              }])                          -> grid ( -row => 7, -column => 3, -sticky => 'w');

   $TIGRfam_button = $database_options->Label (-text     => 'TIGRfam db',
                                               -relief   => 'flat',
                                               -width    => '18',
                                              )                                            -> grid (-row => 8, -column => 0, -sticky => 'w');
   $database_options->Label(-justify     => 'left',
                           -relief       => 'groove',
                           -borderwidth  => 2,
                           -wraplength   => 170,
                           -textvariable => \$auto_ini{TIGRfam_db},
                           -width        => 30)                                            -> grid (-row => 8, -column => 1, -sticky => 'w');

   $TIGRfam_add = $database_options->Button(-text    => 'Add TIGRfam db',
                                            -width   => '18',
                                            -command => [sub {
                                                              &select_tigrfam_db(main_window    => \$mw,
                                                                                 ini_ref        => \%ini,
                                                                                 auto_ini_ref   => \%auto_ini
                                                                                );
                                                             }])                           -> grid (-row => 8, -column => 2, -sticky => 'w');

   $TIGRfam_clear = $database_options->Button(-text    => 'Clear selection',
                                              -width   => '18',
                                              -command => [sub {
                                                                &clear_db_selection(main_window     => \$mw,
                                                                                    ini_ref         => \%ini,
                                                                                    auto_ini_ref    => \%auto_ini,
                                                                                    short_db        => 'TIGRfam_db',
                                                                                    full_db         => 'full_TIGRfam_db',
                                                                                   );
                                                               }])                         -> grid ( -row => 8, -column => 3, -sticky => 'w');
   #set initial state for buttons
   #Glimmer2/Glimmer3
   if (($auto_ini{use_glimmer2} == 1 || $auto_ini{use_glimmer3} == 1) && $auto_ini{internal_gm} == 1) {
      $glimmer_button       ->configure(-state => 'normal');
      $glimmer_add          ->configure(-state => 'normal');
      $glimmer_make_training->configure(-state => 'normal');
   } elsif ($auto_ini{use_glimmer2} == 1 && $auto_ini{use_glimmer3} == 1) {
      $glimmer_button       ->configure(-state => 'disabled');
      $glimmer_add          ->configure(-state => 'disabled');
      $glimmer_make_training->configure(-state => 'disabled');
   } elsif ($auto_ini{internal_gm} == 0) {
      $glimmer_button       ->configure(-state => 'disabled');
      $glimmer_add          ->configure(-state => 'disabled');
      $glimmer_make_training->configure(-state => 'disabled');
   }
   if ($auto_ini{make_training_file} == 1) {
      $glimmer_add          ->configure(-state => 'disabled');
   } elsif ($auto_ini{make_training_file} == 0) {
      $glimmer_add          ->configure(-state => 'normal');
   }
   #Critica
   if ($auto_ini{use_critica} == 1 && $auto_ini{internal_gm} == 1) {
      $critica_button->configure(-state => 'normal');
      $critica_add   ->configure(-state => 'normal');
      $critica_clear ->configure(-state => 'normal');
   } elsif ($auto_ini{use_critica} == 0 || $auto_ini{internal_gm} == 0) {
      $critica_button->configure(-state => 'disabled');
      $critica_add   ->configure(-state => 'disabled');
      $critica_clear ->configure(-state => 'disabled');
   }
   #Blast
   if ($auto_ini{blast_selector} == 1) {
      $blast_button    -> configure (-state => 'normal');
      $blast_sel       -> configure (-state => 'normal');
      $Blast_add       -> configure (-state => 'normal');
      $Blast_clear     -> configure (-state => 'normal');
      $gene_model_only -> configure (-state => 'disabled');
   } elsif ($auto_ini{blast_selector} == 0) {
      $blast_button    -> configure (-state => 'disabled');
      $blast_sel       -> configure (-state => 'disabled');
      $Blast_add       -> configure (-state => 'disabled');
      $Blast_clear       -> configure (-state => 'disabled');
      $gene_model_only -> configure (-state => 'normal');
   }
   #COG
   if ($auto_ini{COG_selector} == 1) {
      $COG_button   ->configure(-state => 'normal');
      $COG_sel      ->configure(-state => 'normal');
   } elsif ($auto_ini{COG_selector} == 0) {
      $COG_button   ->configure(-state => 'disabled');
      $COG_sel      ->configure(-state => 'disabled');
   }
   $mw->update;

   #PFam
   if ($auto_ini{Pfam_selector} == 1) {
      $Pfam_fam_sel -> configure (-state      => 'normal',
                                  -foreground => 'black');
      $Pfam_button  -> configure (-state      => 'normal');
      $Pfam_sel     -> configure (-state      => 'normal');
      $Pfam_add     -> configure (-state      => 'normal');
      $Pfam_clear   -> configure (-state      => 'normal');
      $Pfam_version2-> configure (-state      => 'normal');
      $Pfam_version3-> configure (-state      => 'normal');
      #if ($auto_ini{use_Pfam3} == 1) {
      $PFam3_thrshld-> configure (-state      => 'normal');
      #} elsif ($auto_ini{use_Pfam3} == 0) {
      #$PFam3_thrshld-> configure (-state      => 'disabled');
      #}
   } elsif ($auto_ini{Pfam_selector} == 0) {
      $Pfam_fam_sel -> configure (-state      => 'disabled',
                                  -foreground => 'gray');
      $Pfam_button  -> configure (-state      => 'disabled');
      $Pfam_sel     -> configure (-state      => 'disabled');
      $Pfam_add     -> configure (-state      => 'disabled');
      $Pfam_clear   -> configure (-state      => 'disabled');
      $Pfam_version2-> configure (-state      => 'disabled');
      $Pfam_version3-> configure (-state      => 'disabled');
      $PFam3_thrshld-> configure (-state      => 'disabled');
   }
   #RFam
   if ($auto_ini{ncrna_selector} == 1) {
      $Rfam_button->configure(-state => 'normal');
      $Rfam_add   ->configure(-state => 'normal');
      $Rfam_clear ->configure(-state => 'normal');
   } elsif ($auto_ini{ncrna_selector} == 0) {
      $Rfam_button->configure(-state => 'disabled');
      $Rfam_add   ->configure(-state => 'disabled');
      $Rfam_clear ->configure(-state => 'disabled');
   }
   #TIGRFam
   if ($auto_ini{TIGRfam_selector} == 1) {
      $TIGRfam_button-> configure (-state => 'normal');
      $TIGRfam_add   -> configure (-state => 'normal');
      $TIGRfam_clear -> configure (-state => 'normal');
      $tigr_anno     -> configure (-state => 'normal');
      if ($auto_ini{TIGR_annotation} == 1) {
         $TIGR_EConly    -> configure (-state => 'normal');
         $TIGR_ECCDS     -> configure (-state => 'normal');
         $TIGR_ECdiscard -> configure (-state => 'normal');
      } elsif ($auto_ini{TIGR_annotation} == 0) {
         $TIGR_EConly    -> configure (-state => 'disabled');
         $TIGR_ECCDS     -> configure (-state => 'disabled');
         $TIGR_ECdiscard -> configure (-state => 'disabled');
      }
   } elsif ($auto_ini{TIGRfam_selector} == 0) {
      $TIGRfam_button -> configure (-state => 'disabled');
      $TIGRfam_add    -> configure (-state => 'disabled');
      $TIGRfam_clear  -> configure (-state => 'disabled');
      $tigr_anno      -> configure (-state => 'disabled');
      $TIGR_EConly    -> configure (-state => 'disabled');
      $TIGR_ECCDS     -> configure (-state => 'disabled');
      $TIGR_ECdiscard -> configure (-state => 'disabled');
   }
}

#setup page4: structural options
{
   my $d_cutoff; #signalP trheshold
   $trna_source = $RNA->BrowseEntry(-choices  =>['Bacterial', 'Archaeal', 'Organellar', 'General_Model'],
                                                   -width    => 14,
                                                   -label    => 'Source Type',
                                                   -variable =>\$auto_ini{trna_org},
                                                   -command  => sub {
                                                                     if ($auto_ini{trna_org} eq 'Bacterial') {
                                                                        $auto_ini{trna_code} = '-B';
                                                                     } elsif ($auto_ini{trna_org} eq 'Archaeal') {
                                                                        $auto_ini{trna_code} = '-A';
                                                                     } elsif ($auto_ini{trna_org} eq 'Organellar') {
                                                                        $auto_ini{trna_code} = '-O';
                                                                     } elsif ($auto_ini{trna_org} eq 'General_Model') {
                                                                        $auto_ini{trna_code} = '-G';
                                                                     }
                                                                    }
                                                  )                        ->grid (-row => 0, -column => 1, -sticky => 'w');
   $trna_sensi  = $RNA->Checkbutton(-text     => 'max Sensitivity',
                                    -variable => \$auto_ini{trna_sensitivity},
                                    -command  => sub {
                                                      if ($auto_ini{trna_sensitivity} == 1) {
                                                         $auto_ini{trna_sensi} = '-C';
                                                      } elsif ($auto_ini{trna_sensitivity} == 0) {
                                                         $auto_ini{trna_sensi} = '';
                                                      }
                                                     }
                                   )                        ->grid (-row => 0, -column => 2, -sticky => 'w');
   $trna_save   = $RNA->Checkbutton(-text     => 'save secondary structure',
                                    -variable => \$auto_ini{trna_sec_struc},
                                   )                        ->grid (-row => 0, -column => 3, -sticky => 'w');
   $RNA->Checkbutton(-text     => 'tRNA search',
                     -variable => \$auto_ini{trna_selector},
                     -command  => sub {
                                       if ($auto_ini{trna_selector} == 0) {
                                          $trna_source->configure(-state => 'disabled');
                                          $trna_sensi ->configure(-state => 'disabled');
                                          $trna_save  ->configure(-state => 'disabled');
                                          $trna_source->configure(-labelForeground => 'grey');
                                       } elsif ($auto_ini{trna_selector} == 1) {
                                          $trna_source->configure(-state => 'normal');
                                          $trna_sensi ->configure(-state => 'normal');
                                          $trna_save  ->configure(-state => 'normal');
                                          $trna_source->configure(-labelForeground => 'black');
                                       }
                                      }
                    )                                          ->grid (-row => 0, -column => 0, -sticky => 'w');

   my $rrna_5S  = $RNA->Checkbutton(-text     => 'Search for 5S rRNAs',
                                    -variable => \$auto_ini{rrna_5S}
                                    )                                          ->grid (-row => 1, -column => 1, -sticky => 'w');
   my $rrna_16S = $RNA->Checkbutton(-text     => 'Search for 16S rRNAs',
                                    -variable => \$auto_ini{rrna_16S}
                                    )                                          ->grid (-row => 1, -column => 2, -sticky => 'w');
   my $rrna_23S = $RNA->Checkbutton(-text     => 'Search for 23S rRNAs',
                                    -variable => \$auto_ini{rrna_23S}
                                    )                                          ->grid (-row => 1, -column => 3, -sticky => 'w');
   $RNA->Checkbutton(-text     => 'Estimate position of ribosomal RNAs',
                     -variable => \$auto_ini{rrna_selector},
                     -command  => sub {
                                       if ($auto_ini{rrna_selector} == 0) {
                                          $rrna_5S  ->configure(-state => 'disabled');
                                          $rrna_16S ->configure(-state => 'disabled');
                                          $rrna_23S ->configure(-state => 'disabled');
                                       } elsif ($auto_ini{rrna_selector} == 1) {
                                          $rrna_5S  ->configure(-state => 'normal');
                                          $rrna_16S ->configure(-state => 'normal');
                                          $rrna_23S ->configure(-state => 'normal');
                                       }
                                      }
                    )                                                          ->grid (-row => 1, -column => 0, -sticky => 'w');

   my $ncrna_filter = $RNA->BrowseEntry(-label     => 'Set filter for pipeline    ',
                                         -choices  =>['--default',
                                                      '--max',
                                                      '--nohmm',
                                                      '--mid',
                                                      '--rfam',
                                                      '--hmmonly'
                                                     ],
                                         -width    => 10,
                                         -variable => \$auto_ini{ncrna_pipelinefilter}
                                       )                                       ->grid (-row => 2, -column => 1, -sticky => 'w');
   my $ncrna_threshold   = $RNA->BrowseEntry(-label    =>'Reporting threshold E-value',
                                             -choices  =>['10',
                                                          '1',
                                                          '1e-01',
                                                          '1e-02',
                                                          '1e-03',
                                                          '1e-04',
                                                          '1e-05',
                                                          '1e-06',
                                                          '1e-07',
                                                          '1e-08',
                                                          '1e-09',
                                                          '1e-10'
                                                          ],
                                             -width    => 10,
                                             -variable =>\$auto_ini{ncrna_threshold}
                                            )
                                                                               ->grid (-row => 3, -column => 1, -sticky => 'w', -columnspan => 2);
   my $ncrna_align       = $RNA->Checkbutton(-text       => 'Add hit alignments',
                                             -variable   => \$auto_ini{ncrna_align}
                                    )                                          ->grid (-row => 2, -column => 2, -sticky => 'w');
   my $ncrna_overlap     = $RNA->Checkbutton(-text     => "Allow overlapping RFam\'s",
                                             -variable => \$auto_ini{ncrna_overlap}
                                    )                                          ->grid (-row => 3, -column => 2, -sticky => 'w');
   my $ncrna_IGr_only    = $RNA->Checkbutton(-text     => "Keep only IG RFams",
                                             -variable => \$auto_ini{ncrna_IG_only}
                                    )                                          ->grid (-row => 3, -column => 3, -sticky => 'w');
   my $ncrna_ignore_rrna = $RNA->Checkbutton(-text     => 'Ignore rRNAs',
                                             -variable => \$auto_ini{ncrna_ignore_rrna}
                                    )                                          ->grid (-row => 2, -column => 3, -sticky => 'w');
   my $ncrna_ignore_trna = $RNA->Checkbutton(-text     => 'Ignore tRNAs',
                                             -variable => \$auto_ini{ncrna_ignore_trna}
                                    )                                          ->grid (-row => 2, -column => 4, -sticky => 'w');
   my $ncrna_algorithm   = $RNA->BrowseEntry(-label      =>'Alignment algorithm        ',
                                             -choices    =>['local',
                                                            'glocal'
                                                           ],
                                             -width      => 10,
                                             -variable   =>\$auto_ini{ncrna_local}
                                    )                                    ->grid (-row => 4, -column => 1, -sticky => 'w');
   my $ncrna_verbose = $RNA->Checkbutton(-text         => 'Verbose output',
                                         -variable     => \$auto_ini{ncrna_verbose}
                                        )                                      ->grid (-row => 4, -column => 2, -sticky => 'w');
   my $ncrna_blastlength = $RNA->BrowseEntry(-label    =>'Split sequences [knt]      ',
                                             -choices  =>[0..200],
                                             -width    => 10,
                                             -variable =>\$auto_ini{ncrna_splitseq}
                                             )
                                                                               ->grid (-row => 5, -column => 1, -sticky => 'w');
   $RNA->Checkbutton(-text     => 'non-coding RNA search',
                     -variable => \$auto_ini{ncrna_selector},
                     -command  => sub {
                                       if ($auto_ini{ncrna_selector} == 0) {
                                          $ncrna_filter         ->configure(-state           => 'disabled');
                                          $ncrna_filter         ->configure(-labelForeground => 'grey');
                                          $ncrna_algorithm      ->configure(-state           => 'disabled');
                                          $ncrna_algorithm      ->configure(-labelForeground => 'grey');
                                          $ncrna_align          ->configure(-state           => 'disabled');
                                          $ncrna_ignore_rrna    ->configure(-state           => 'disabled');
                                          $ncrna_ignore_trna    ->configure(-state           => 'disabled');
                                          $ncrna_overlap        ->configure(-state           => 'disabled');
                                          $ncrna_threshold      ->configure(-state           => 'disabled');
                                          $ncrna_threshold      ->configure(-labelForeground => 'grey');
                                          $ncrna_blastlength    ->configure(-state           => 'disabled');
                                          $ncrna_blastlength    ->configure(-labelForeground => 'grey');
                                          $Rfam_button          ->configure(-state           => 'disabled');
                                          $Rfam_add             ->configure(-state           => 'disabled');
                                          $Rfam_clear           ->configure(-state           => 'disabled');
                                          $ncrna_verbose        ->configure(-state           => 'disabled');
                                          $ncrna_IGr_only       ->configure(-state           => 'disabled');
                                       } elsif ($auto_ini{ncrna_selector} == 1) {
                                          $ncrna_filter         ->configure(-state           => 'normal');
                                          $ncrna_filter         ->configure(-labelForeground => 'black');
                                          $ncrna_algorithm      ->configure(-state           => 'normal');
                                          $ncrna_algorithm      ->configure(-labelForeground => 'black');
                                          $ncrna_align          ->configure(-state           => 'normal');
                                          $ncrna_ignore_rrna    ->configure(-state           => 'normal');
                                          $ncrna_ignore_trna    ->configure(-state           => 'normal');
                                          $ncrna_overlap        ->configure(-state           => 'normal');
                                          $ncrna_threshold      ->configure(-state           => 'normal');
                                          $ncrna_threshold      ->configure(-labelForeground => 'black');
                                          $ncrna_blastlength    ->configure(-state           => 'normal');
                                          $ncrna_blastlength    ->configure(-labelForeground => 'black');
                                          $Rfam_button          ->configure(-state           => 'normal');
                                          $Rfam_add             ->configure(-state           => 'normal');
                                          $Rfam_clear           ->configure(-state           => 'normal');
                                          $ncrna_verbose        ->configure(-state           => 'normal');
                                          $ncrna_IGr_only       ->configure(-state           => 'normal');
                                       }
                                      }
                                   )                                          ->grid (-row => 2, -column => 0, -sticky => 'w');
   $membrane->Checkbutton(-text     => 'Transmembrane prediction',
                          -variable => \$auto_ini{tmhmm_selector},
                          -command  => sub {
                                            if ($auto_ini{tmhmm_selector} == 0) {
                                               $tmm_domain->configure(-state => 'disabled');
                                            } elsif ($auto_ini{tmhmm_selector} == 1) {
                                               $tmm_domain->configure(-state => 'normal');
                                            }
                                           })                                  ->grid (-row => 4, -column => 0, -sticky => 'w');

   #initiate SignalP
   #detect SignalP version
   if (-e $ini{signalp_dir}.'/signalp') {
      my ($file_ref) = &slurp(main_window   => $mw,
                              auto_ini_ref  => \%auto_ini,
                              ini_ref       => \%ini,
                              filename      => 'signalp',
                              directory     => $ini{signalp_dir}
                             );
      if ($$file_ref =~ m/SignalP 4/si) {
         $signalp_version = 4;
      } elsif ($$file_ref =~ m/SignalP 3/si) {
         $signalp_version = 3;
      }
      undef $file_ref;
   }
   #assume latest version for now if nothing there yet
   unless (defined $signalp_version && $signalp_version =~ m/\d+/) {
      $signalp_version = 4;
   }

   $membrane->Checkbutton(-text     => 'SignalP prediction',
                          -variable => \$auto_ini{signalp_selector},
                          -command  => sub {
                                            if ($auto_ini{signalp_selector} == 0) {
                                               $signalp_source ->configure(-state => 'disabled');
                                               $signalp_trunc  ->configure(-state => 'disabled');
                                               $signalp_method ->configure(-state => 'disabled');
                                               if ($signalp_version == 4) {
                                                  $signalp_cutoff ->configure(-state => 'disabled');
                                                  $signalp_cutoff ->configure(-labelForeground => 'grey');
                                               }
                                               $signalp_source ->configure(-labelForeground => 'grey');
                                               $signalp_trunc  ->configure(-labelForeground => 'grey');
                                               $signalp_method ->configure(-labelForeground => 'grey');
                                            } elsif ($auto_ini{signalp_selector} == 1) {
                                               $signalp_source ->configure(-state => 'normal');
                                               $signalp_trunc  ->configure(-state => 'normal');
                                               $signalp_method ->configure(-state => 'normal');
                                               if ($signalp_version == 4) {
                                                  $signalp_cutoff ->configure(-state => 'normal');
                                                  $signalp_cutoff ->configure(-labelForeground => 'black');
                                               }
                                               $signalp_source ->configure(-labelForeground => 'black');
                                               $signalp_trunc  ->configure(-labelForeground => 'black');
                                               $signalp_method ->configure(-labelForeground => 'black');
                                            }
                                           }
                          )                                                    ->grid (-row => 5, -column => 0, -sticky => 'w');
   $DNAstructure->Checkbutton(-text     => 'Terminator structures',
                              -variable => \$auto_ini{terminator_selector}
                             )                                                 ->grid (-row => 6, -column => 0, -sticky => 'w');

   $DNAstructure->Checkbutton(-text     => 'CRISPR structures',
                              -variable => \$auto_ini{CRISPR_selector},
                              -command  => sub {
                                                if ($auto_ini{CRISPR_selector} == 0) {
                                                   $CRISPR_minNR   ->configure(-state => 'disabled');
                                                   $CRISPR_minRL   ->configure(-state => 'disabled');
                                                   $CRISPR_maxRL   ->configure(-state => 'disabled');
                                                   #$CRISPR_searchWL->configure(-state => 'disabled');
                                                   $CRISPR_minSL   ->configure(-state => 'disabled');
                                                   $CRISPR_maxSL   ->configure(-state => 'disabled');
                                                   $CRISPR_minNR   ->configure(-labelForeground => 'grey');
                                                   $CRISPR_minRL   ->configure(-labelForeground => 'grey');
                                                   $CRISPR_maxRL   ->configure(-labelForeground => 'grey');
                                                   #$CRISPR_searchWL->configure(-labelForeground => 'grey');
                                                   $CRISPR_minSL   ->configure(-labelForeground => 'grey');
                                                   $CRISPR_maxSL   ->configure(-labelForeground => 'grey');
                                                } elsif ($auto_ini{CRISPR_selector} == 1) {
                                                   $CRISPR_minNR   ->configure(-state => 'normal');
                                                   $CRISPR_minRL   ->configure(-state => 'normal');
                                                   $CRISPR_maxRL   ->configure(-state => 'normal');
                                                   #$CRISPR_searchWL->configure(-state => 'normal');
                                                   $CRISPR_minSL   ->configure(-state => 'normal');
                                                   $CRISPR_maxSL   ->configure(-state => 'normal');
                                                   $CRISPR_minNR   ->configure(-labelForeground => 'black');
                                                   $CRISPR_minRL   ->configure(-labelForeground => 'black');
                                                   $CRISPR_maxRL   ->configure(-labelForeground => 'black');
                                                   #$CRISPR_searchWL->configure(-labelForeground => 'black');
                                                   $CRISPR_minSL   ->configure(-labelForeground => 'black');
                                                   $CRISPR_maxSL   ->configure(-labelForeground => 'black');
                                                }
                                               }
                             )                                                 ->grid (-row => 7, -column => 0, -sticky => 'w');

   $tmm_domain  = $membrane->Checkbutton(-text     => 'Indicate individual membrane domains',
                                         -variable => \$auto_ini{tmhmm_domains}
                                        )
                                                                               ->grid (-row => 4, -column => 1, -sticky => 'w', -columnspan => 2);

   #define some original states here rather than below to set GUI
   if ($auto_ini{signalp_org} eq 'gram positive') {
      $auto_ini{signalp_code} = 'gram+';
   } elsif ($auto_ini{signalp_org} eq 'gram negative') {
      $auto_ini{signalp_code} = 'gram-';
   } elsif ($auto_ini{signalp_org} eq 'eukaryotic') {
      $auto_ini{signalp_code} = 'euk';
   } elsif ($auto_ini{signalp_org} eq 'gram pos+neg') {
      $auto_ini{signalp_code} = 'both';
   }
   if ($signalp_version == 3) {
      if ($auto_ini{signalp_method} eq 'neural networks') {
         $auto_ini{signalp_meth} = 'nn';
      } elsif ($auto_ini{signalp_method} eq 'hidden Markov') {
         $auto_ini{signalp_meth} = 'hmm';
      } elsif ($auto_ini{signalp_method} eq 'both') {
         $auto_ini{signalp_meth} = 'nn+hmm';
      }
   } else {
      if ($auto_ini{signalp4_method} eq 'SignalP-TM') {
         $auto_ini{signalp_meth} = 'best';
      } elsif ($auto_ini{signalp4_method} eq 'SignalP-noTM') {
         $auto_ini{signalp_meth} = 'notm';
      }
   }

   #adjust GUI to signalP version
   if ($signalp_version == 3) {
      $signalp_source = $membrane->BrowseEntry(-choices  =>['gram positive', 'gram negative', 'gram pos+neg', 'eukaryotic'],
                                               -width    => 14,
                                               -label    => 'Source Type',
                                               -variable =>\$auto_ini{signalp_org},
                                               -command  => sub {
                                                                 if ($auto_ini{signalp_org} eq 'gram positive') {
                                                                    $auto_ini{signalp_code} = 'gram+';
                                                                 } elsif ($auto_ini{signalp_org} eq 'gram negative') {
                                                                    $auto_ini{signalp_code} = 'gram-';
                                                                 } elsif ($auto_ini{signalp_org} eq 'eukaryotic') {
                                                                    $auto_ini{signalp_code} = 'euk';
                                                                 } elsif ($auto_ini{signalp_org} eq 'gram pos+neg') {
                                                                    $auto_ini{signalp_code} = 'both';
                                                                 }
                                                                }
                                              )                                   ->grid (-row => 5, -column => 1, -sticky => 'w');

      $signalp_trunc  = $membrane->BrowseEntry(-label   =>'Truncate to',
                                               -choices  =>[1..70],
                                               -width    => 5,
                                               -variable =>\$auto_ini{signalp_trunc}
                                              )                                   ->grid (-row => 5, -column => 2, -sticky => 'w');
      $signalp_method = $membrane->BrowseEntry(-choices  =>['neural networks', 'hidden Markov', 'both'],
                                               -width    => 14,
                                               -label    => 'Method',
                                               -variable =>\$auto_ini{signalp_method},
                                               -command  => sub {
                                                                 if ($auto_ini{signalp_method} eq 'neural networks') {
                                                                    $auto_ini{signalp_meth} = 'nn';
                                                                 } elsif ($auto_ini{signalp_method} eq 'hidden Markov') {
                                                                    $auto_ini{signalp_meth} = 'hmm';
                                                                 } elsif ($auto_ini{signalp_method} eq 'both') {
                                                                    $auto_ini{signalp_meth} = 'nn+hmm';
                                                                 }
                                                                }
                                              )                                   ->grid (-row => 5, -column => 3, -sticky => 'w');
   } elsif ($signalp_version == 4) {
      #define suggested D-cutoffs for TM (best) and noTM networks
      $d_cutoff->{'euk'}  ->{'notm'} = 0.45; $d_cutoff->{'euk'}  ->{'best'} = 0.50;
      $d_cutoff->{'gram+'}->{'notm'} = 0.57; $d_cutoff->{'gram+'}->{'best'} = 0.45;
      $d_cutoff->{'gram-'}->{'notm'} = 0.57; $d_cutoff->{'gram-'}->{'best'} = 0.51;

      #reset selection to valid option if required
      if ($auto_ini{signalp_org} eq 'gram pos+neg') {
         $auto_ini{signalp_org} = 'gram+';
      }

      $signalp_source = $membrane->BrowseEntry(-choices  =>['gram positive', 'gram negative', 'eukaryotic'],
                                               -width    => 14,
                                               -label    => 'Source Type',
                                               -variable =>\$auto_ini{signalp_org},
                                               -command  => sub {
                                                                 if ($auto_ini{signalp_org} eq 'gram positive') {
                                                                    $auto_ini{signalp_code} = 'gram+';
                                                                    $signalp_cutoff->configure(-variable=>\$d_cutoff->{$auto_ini{signalp_code}}->{$auto_ini{signalp_meth}});
                                                                    $auto_ini{signalp_cutoff} = $d_cutoff->{$auto_ini{signalp_code}}->{$auto_ini{signalp_meth}};
                                                                 } elsif ($auto_ini{signalp_org} eq 'gram negative') {
                                                                    $auto_ini{signalp_code} = 'gram-';
                                                                    $signalp_cutoff->configure(-variable=>\$d_cutoff->{$auto_ini{signalp_code}}->{$auto_ini{signalp_meth}});
                                                                    $auto_ini{signalp_cutoff} = $d_cutoff->{$auto_ini{signalp_code}}->{$auto_ini{signalp_meth}};
                                                                 } elsif ($auto_ini{signalp_org} eq 'eukaryotic') {
                                                                    $auto_ini{signalp_code} = 'euk';
                                                                    $signalp_cutoff->configure(-variable=>\$d_cutoff->{$auto_ini{signalp_code}}->{$auto_ini{signalp_meth}});
                                                                    $auto_ini{signalp_cutoff} = $d_cutoff->{$auto_ini{signalp_code}}->{$auto_ini{signalp_meth}};
                                                                 }
                                                                }
                                              )                                   ->grid (-row => 5, -column => 1, -sticky => 'w');

      $signalp_min  = $membrane->BrowseEntry(-label   =>'Minimum length',
                                               -choices  =>[1..70],
                                               -width    => 5,
                                               -variable =>\$auto_ini{signalp_min}
                                              )                                   ->grid (-row => 5, -column => 2, -sticky => 'w');
      $signalp_trunc  = $membrane->BrowseEntry(-label   =>'Truncate to   ',
                                               -choices  =>[1..70],
                                               -width    => 5,
                                               -variable =>\$auto_ini{signalp_trunc}
                                              )                                   ->grid (-row => 6, -column => 2, -sticky => 'w');
      $signalp_method = $membrane->BrowseEntry(-choices  =>['SignalP-TM', 'SignalP-noTM'],
                                               -width    => 14,
                                               -label    => 'Method',
                                               -variable =>\$auto_ini{signalp4_method},
                                               -command  => sub {
                                                                 if ($auto_ini{signalp4_method} eq 'SignalP-TM') {
                                                                    $auto_ini{signalp_meth} = 'best';
                                                                    $signalp_cutoff->configure(-variable=>\$d_cutoff->{$auto_ini{signalp_code}}->{$auto_ini{signalp_meth}});
                                                                    $auto_ini{signalp_cutoff} = $d_cutoff->{$auto_ini{signalp_code}}->{$auto_ini{signalp_meth}};
                                                                 } elsif ($auto_ini{signalp4_method} eq 'SignalP-noTM') {
                                                                    $auto_ini{signalp_meth} = 'notm';
                                                                    $signalp_cutoff->configure(-variable=>\$d_cutoff->{$auto_ini{signalp_code}}->{$auto_ini{signalp_meth}});
                                                                    $auto_ini{signalp_cutoff} = $d_cutoff->{$auto_ini{signalp_code}}->{$auto_ini{signalp_meth}};
                                                                 }
                                                                }
                                              )                                   ->grid (-row => 5, -column => 3, -sticky => 'w');

      $signalp_cutoff  = $membrane->BrowseEntry(-label    =>'Network cutoff ',
                                                -width    => 5,
                                                -variable =>\$d_cutoff->{$auto_ini{signalp_code}}->{$auto_ini{signalp_meth}},
                                                -command  => sub {
                                                                  $auto_ini{signalp_cutoff} = $d_cutoff->{$auto_ini{signalp_code}}->{$auto_ini{signalp_meth}};
                                                                 }
                                               )                                  ->grid (-row => 6, -column => 3, -sticky => 'w');
   }

   $CRISPR_minNR   = $DNAstructure->BrowseEntry(-label    =>'minimum number of repeats',
                                                -choices  =>[1..100],
                                                -width    => 5,
                                                -variable =>\$auto_ini{CRISPR_minNR}
                                              )
                                                                               ->grid (-row => 7, -column => 1, -sticky => 'w', -columnspan => 1);
   $CRISPR_minRL   = $DNAstructure->BrowseEntry(-label    =>'minimum length of repeats',
                                                -choices  =>[1..100],
                                                -width    => 5,
                                                -variable =>\$auto_ini{CRISPR_minRL}
                                              )
                                                                               ->grid (-row => 7, -column => 2, -sticky => 'w', -columnspan => 1);
   $CRISPR_maxRL   = $DNAstructure->BrowseEntry(-label    =>'maximum length of repeats',
                                                -choices  =>[1..100],
                                                -width    => 5,
                                                -variable =>\$auto_ini{CRISPR_maxRL}
                                              )
                                                                               ->grid (-row => 7, -column => 3, -sticky => 'w', -columnspan => 1);
   #$CRISPR_searchWL= $DNAstructure->BrowseEntry(-label    =>'length of search window___',
   #                                                   -choices  =>[6..9],
   #                                                   -width    => 5,
   #                                                   -variable =>\$auto_ini{CRISPR_maxWL}
   #                                                 )
   #                                                                      ->grid (-row => 8, -column => 1, -sticky => 'w', -columnspan => 1);
   $CRISPR_minSL   = $DNAstructure->BrowseEntry(-label    =>'minimum length of spacer_',
                                                -choices  =>[1..100],
                                                -width    => 5,
                                                -variable =>\$auto_ini{CRISPR_minSL}
                                              )
                                                                               ->grid (-row => 8, -column => 2, -sticky => 'w', -columnspan => 1);
   $CRISPR_maxSL   = $DNAstructure->BrowseEntry(-label    =>'maximum length of spacer_',
                                                -choices  =>[1..100],
                                                -width    => 5,
                                                -variable =>\$auto_ini{CRISPR_maxSL}
                                              )
                                                                               ->grid (-row => 8, -column => 3, -sticky => 'w', -columnspan => 1);
   my $vector_button = $others->Button(-text    => 'Select Vector db',
                                       -width   => '18',
                                       -command => [sub {
                                                         &select_blast_db(main_window    => \$mw,
                                                                          ini_ref        => \%ini,
                                                                          auto_ini_ref   => \%auto_ini,
                                                                          vector_db      => 1
                                                                         );
                                                        }])                    -> grid (-row => 7, -column => 1, -sticky => 'w');
   my $vector_path = $others->Label(-justify      => 'left',
                                    -relief       => 'groove',
                                    -borderwidth  => 2,
                                    -wraplength   => 170,
                                    -textvariable => \$auto_ini{selected_vector_db},
                                    -width        => 30)                       -> grid (-row => 7, -column => 2, -sticky => 'w');

   my $vector_hits = $others->BrowseEntry(-label    =>'Maximum hit number',
                                          -choices  =>[1..100],
                                          -width    => 5,
                                          -variable =>\$auto_ini{vector_hit_number}
                                         )                                     -> grid (-row => 7, -column => 3, -sticky => 'w');
   my $vector_penalty = $others->BrowseEntry(-label    =>'Penalty',
                                             -choices  =>[-100..100],
                                             -width    => 5,
                                             -variable =>\$auto_ini{vector_penalty}
                                           )
                                                                               -> grid (-row => 8, -column => 1, -sticky => 'w');
   my $vector_reward  = $others->BrowseEntry(-label    =>'Reward ',
                                             -choices  =>[-100..100],
                                             -width    => 5,
                                             -variable =>\$auto_ini{vector_reward}
                                           )
                                                                               -> grid (-row => 9, -column => 1, -sticky => 'w');
   my $vector_gapopen = $others->BrowseEntry(-label    =>'gapopen',
                                             -choices  =>[1..100],
                                             -width    => 5,
                                             -variable =>\$auto_ini{vector_gapopen}
                                           )
                                                                               -> grid (-row => 8, -column => 2, -sticky => 'w');
   my $vector_gapextend = $others->BrowseEntry(-label    =>'gapextend    ',
                                               -choices  =>[1..100],
                                               -width    => 10,
                                               -variable =>\$auto_ini{vector_gapextend}
                                              )
                                                                               -> grid (-row => 8, -column => 3, -sticky => 'w');
   my $vector_dust = $others->Checkbutton(-text     =>'dust ',
                                          -variable =>\$auto_ini{vector_dust}
                                         )
                                                                               -> grid (-row => 10, -column => 1, -sticky => 'w');
   my $vector_soft_masking = $others->Checkbutton(-text     =>'soft_masking',
                                                  -variable =>\$auto_ini{vector_soft_masking}
                                                  )
                                                                               -> grid (-row => 11, -column => 1, -sticky => 'w');
   my $vector_evalue = $others->BrowseEntry(-label     =>'evalue ',
                                             -choices  =>[1..1000],
                                             -width    => 5,
                                             -variable =>\$auto_ini{vector_evalue}
                                           )
                                                                               -> grid (-row => 9, -column => 2, -sticky => 'w');
   my $vector_xdrop = $others->BrowseEntry(-label    =>'xdrop_ungap  ',
                                           -choices  =>[1..100],
                                           -width    => 10,
                                           -variable =>\$auto_ini{vector_xdrop}
                                           )
                                                                               -> grid (-row => 9, -column => 3, -sticky => 'w');
   $others->Checkbutton(-text     => 'Screen for vector sequences',
                        -variable => \$auto_ini{vector_screen_selector},
                        -command  => sub {
                                          if ($auto_ini{vector_screen_selector} == 0) {
                                             $vector_button      ->configure(-state => 'disabled');
                                             $vector_path        ->configure(-state => 'disabled');
                                             $vector_penalty     ->configure(-state => 'disabled');
                                             $vector_gapopen     ->configure(-state => 'disabled');
                                             $vector_gapextend   ->configure(-state => 'disabled');
                                             $vector_dust        ->configure(-state => 'disabled');
                                             $vector_soft_masking->configure(-state => 'disabled');
                                             $vector_evalue      ->configure(-state => 'disabled');
                                             $vector_xdrop       ->configure(-state => 'disabled');
                                             $vector_hits        ->configure(-state => 'disabled');
                                             $vector_penalty     ->configure(-labelForeground => 'grey');
                                             $vector_gapopen     ->configure(-labelForeground => 'grey');
                                             $vector_gapextend   ->configure(-labelForeground => 'grey');
                                             $vector_evalue      ->configure(-labelForeground => 'grey');
                                             $vector_xdrop       ->configure(-labelForeground => 'grey');
                                          } elsif ($auto_ini{vector_screen_selector} == 1) {
                                             $vector_button      ->configure(-state => 'normal');
                                             $vector_path        ->configure(-state => 'normal');
                                             $vector_penalty     ->configure(-state => 'normal');
                                             $vector_gapopen     ->configure(-state => 'normal');
                                             $vector_gapextend   ->configure(-state => 'normal');
                                             $vector_dust        ->configure(-state => 'normal');
                                             $vector_soft_masking->configure(-state => 'normal');
                                             $vector_evalue      ->configure(-state => 'normal');
                                             $vector_xdrop       ->configure(-state => 'normal');
                                             $vector_hits        ->configure(-state => 'normal');
                                             $vector_penalty     ->configure(-labelForeground => 'black');
                                             $vector_gapopen     ->configure(-labelForeground => 'black');
                                             $vector_gapextend   ->configure(-labelForeground => 'black');
                                             $vector_evalue      ->configure(-labelForeground => 'black');
                                             $vector_xdrop       ->configure(-labelForeground => 'black');
                                          }
                                         }
                       )                                                       -> grid (-row => 7, -column => 0, -sticky => 'w');

   #set initial state for buttons
   #vector screen
   if ($auto_ini{vector_screen_selector} == 1) {
      $vector_button      ->configure(-state => 'normal');
      $vector_path        ->configure(-state => 'normal');
      $vector_penalty     ->configure(-state => 'normal');
      $vector_gapopen     ->configure(-state => 'normal');
      $vector_gapextend   ->configure(-state => 'normal');
      $vector_dust        ->configure(-state => 'normal');
      $vector_soft_masking->configure(-state => 'normal');
      $vector_evalue      ->configure(-state => 'normal');
      $vector_xdrop       ->configure(-state => 'normal');
      $vector_hits        ->configure(-state => 'normal');
      $vector_penalty     ->configure(-labelForeground => 'black');
      $vector_gapopen     ->configure(-labelForeground => 'black');
      $vector_gapextend   ->configure(-labelForeground => 'black');
      $vector_evalue      ->configure(-labelForeground => 'black');
      $vector_xdrop       ->configure(-labelForeground => 'black');
   } elsif ($auto_ini{vector_screen_selector} == 0) {
      $vector_button      ->configure(-state => 'disabled');
      $vector_path        ->configure(-state => 'disabled');
      $vector_penalty     ->configure(-state => 'disabled');
      $vector_gapopen     ->configure(-state => 'disabled');
      $vector_gapextend   ->configure(-state => 'disabled');
      $vector_dust        ->configure(-state => 'disabled');
      $vector_soft_masking->configure(-state => 'disabled');
      $vector_evalue      ->configure(-state => 'disabled');
      $vector_xdrop       ->configure(-state => 'disabled');
      $vector_hits        ->configure(-state => 'disabled');
      $vector_penalty     ->configure(-labelForeground => 'grey');
      $vector_gapopen     ->configure(-labelForeground => 'grey');
      $vector_gapextend   ->configure(-labelForeground => 'grey');
      $vector_evalue      ->configure(-labelForeground => 'grey');
      $vector_xdrop       ->configure(-labelForeground => 'grey');
   }
   #trna
   if ($auto_ini{trna_selector} == 1) {
      $trna_source->configure(-state => 'normal');
      $trna_sensi ->configure(-state => 'normal');
      $trna_save  ->configure(-state => 'normal');
      $trna_source->configure(-labelForeground => 'black');
   } elsif ($auto_ini{trna_selector} == 0) {
      $trna_source->configure(-state => 'disabled');
      $trna_sensi ->configure(-state => 'disabled');
      $trna_save  ->configure(-state => 'disabled');
      $trna_source->configure(-labelForeground => 'grey');
   }
   if ($auto_ini{trna_org} eq 'Bacterial') {
      $auto_ini{trna_code} = '-B';
   } elsif ($auto_ini{trna_org} eq 'Archaeal') {
      $auto_ini{trna_code} = '-A';
   } elsif ($auto_ini{trna_org} eq 'Organellar') {
      $auto_ini{trna_code} = '-O';
   } elsif ($auto_ini{trna_org} eq 'General_Model') {
      $auto_ini{trna_code} = '-G';
   }
   #rRNA
   if ($auto_ini{rrna_selector} == 0) {
      $rrna_5S  ->configure(-state => 'disabled');
      $rrna_16S ->configure(-state => 'disabled');
      $rrna_23S ->configure(-state => 'disabled');
   } elsif ($auto_ini{rrna_selector} == 1) {
      $rrna_5S  ->configure(-state => 'normal');
      $rrna_16S ->configure(-state => 'normal');
      $rrna_23S ->configure(-state => 'normal');
   }
   #TMM
   if ($auto_ini{tmhmm_selector} == 1) {
      $tmm_domain->configure(-state => 'normal');
   } elsif ($auto_ini{tmhmm_selector} == 0) {
      $tmm_domain->configure(-state => 'disabled');
   }
   #signalP
   if ($auto_ini{signalp_selector} == 1) {
      $signalp_source ->configure(-state => 'normal');
      $signalp_trunc  ->configure(-state => 'normal');
      $signalp_method ->configure(-state => 'normal');
      if ($signalp_version == 4) {
         $signalp_cutoff ->configure(-state => 'normal');
         $signalp_cutoff ->configure(-labelForeground => 'black');
         $auto_ini{signalp_cutoff} = $d_cutoff->{$auto_ini{signalp_code}}->{$auto_ini{signalp_meth}};
      }
      $signalp_source ->configure(-labelForeground => 'black');
      $signalp_trunc  ->configure(-labelForeground => 'black');
      $signalp_method ->configure(-labelForeground => 'black');
   } elsif ($auto_ini{signalp_selector} == 0) {
      $signalp_source ->configure(-state => 'disabled');
      $signalp_trunc  ->configure(-state => 'disabled');
      $signalp_method ->configure(-state => 'disabled');
      if ($signalp_version == 4) {
         $signalp_cutoff ->configure(-state => 'disabled');
         $signalp_cutoff ->configure(-labelForeground => 'grey');
      }
      $signalp_source ->configure(-labelForeground => 'grey');
      $signalp_trunc  ->configure(-labelForeground => 'grey');
      $signalp_method ->configure(-labelForeground => 'grey');
   }

   #ncRNA
   if ($auto_ini{ncrna_selector} == 0) {
      $ncrna_filter         ->configure(-state           => 'disabled');
      $ncrna_filter         ->configure(-labelForeground => 'grey');
      $ncrna_algorithm      ->configure(-state           => 'disabled');
      $ncrna_algorithm      ->configure(-labelForeground => 'grey');
      $ncrna_align          ->configure(-state           => 'disabled');
      $ncrna_ignore_rrna    ->configure(-state           => 'disabled');
      $ncrna_ignore_trna    ->configure(-state           => 'disabled');
      $ncrna_overlap        ->configure(-state           => 'disabled');
      $ncrna_threshold      ->configure(-state           => 'disabled');
      $ncrna_threshold      ->configure(-labelForeground => 'grey');
      $ncrna_blastlength    ->configure(-state           => 'disabled');
      $ncrna_blastlength    ->configure(-labelForeground => 'grey');
      $ncrna_verbose        ->configure(-state           => 'disabled');
      $ncrna_IGr_only       ->configure(-state           => 'disabled');
   } elsif ($auto_ini{ncrna_selector} == 1) {
      $ncrna_filter         ->configure(-state           => 'normal');
      $ncrna_filter         ->configure(-labelForeground => 'black');
      $ncrna_algorithm      ->configure(-state           => 'normal');
      $ncrna_algorithm      ->configure(-labelForeground => 'black');
      $ncrna_align          ->configure(-state           => 'normal');
      $ncrna_ignore_rrna    ->configure(-state           => 'normal');
      $ncrna_ignore_trna    ->configure(-state           => 'normal');
      $ncrna_overlap        ->configure(-state           => 'normal');
      $ncrna_threshold      ->configure(-state           => 'normal');
      $ncrna_threshold      ->configure(-labelForeground => 'black');
      $ncrna_blastlength    ->configure(-state           => 'normal');
      $ncrna_blastlength    ->configure(-labelForeground => 'black');
      $ncrna_verbose        ->configure(-state           => 'normal');
      $ncrna_IGr_only       ->configure(-state           => 'normal');
   }
   #CRISPR
   if ($auto_ini{CRISPR_selector} == 0) {
      $CRISPR_minNR   ->configure(-state => 'disabled');
      $CRISPR_minRL   ->configure(-state => 'disabled');
      $CRISPR_maxRL   ->configure(-state => 'disabled');
      #$CRISPR_searchWL->configure(-state => 'disabled');
      $CRISPR_minSL   ->configure(-state => 'disabled');
      $CRISPR_maxSL   ->configure(-state => 'disabled');
      $CRISPR_minNR   ->configure(-labelForeground => 'grey');
      $CRISPR_minRL   ->configure(-labelForeground => 'grey');
      $CRISPR_maxRL   ->configure(-labelForeground => 'grey');
      #$CRISPR_searchWL->configure(-labelForeground => 'grey');
      $CRISPR_minSL   ->configure(-labelForeground => 'grey');
      $CRISPR_maxSL   ->configure(-labelForeground => 'grey');
   } elsif ($auto_ini{CRISPR_selector} == 1) {
      $CRISPR_minNR   ->configure(-state => 'normal');
      $CRISPR_minRL   ->configure(-state => 'normal');
      $CRISPR_maxRL   ->configure(-state => 'normal');
      #$CRISPR_searchWL->configure(-state => 'normal');
      $CRISPR_minSL   ->configure(-state => 'normal');
      $CRISPR_maxSL   ->configure(-state => 'normal');
      $CRISPR_minNR   ->configure(-labelForeground => 'black');
      $CRISPR_minRL   ->configure(-labelForeground => 'black');
      $CRISPR_maxRL   ->configure(-labelForeground => 'black');
      #$CRISPR_searchWL->configure(-labelForeground => 'black');
      $CRISPR_minSL   ->configure(-labelForeground => 'black');
      $CRISPR_maxSL   ->configure(-labelForeground => 'black');
   }

   $structural_options->gridRowconfigure(0, -pad => 20);
   $structural_options->gridRowconfigure(1, -pad => 20);
   #$structural_options->gridRowconfigure(3, -pad => 20);
   $structural_options->gridRowconfigure(4, -pad => 20);
   $structural_options->gridRowconfigure(5, -pad => 20);
   $structural_options->gridRowconfigure(6, -pad => 20);
}

#setup page5: run
{
   #try as frame instead
   my $option_textbox = $run_options->Text(-width   => 65,
                                           -setgrid => 1,
                                           -relief  => 'flat')
                                  ->grid (-row => 2, -column => 0, -sticky => 'w');
   #insert buttons
   {

      $option_textbox->insert('end', "\n\n\n");
      $option_textbox->windowCreate('end', -window => $run_options->Checkbutton(-text     => "\tConcatenate input files.\n\tmsFASTA is automatically concatenated\n\tmsGenbank is concatenated into msFASTA",
                                                                                -variable => \$auto_ini{concatenate_input_files},
                                                                                -command  => sub {
                                                                                                  if ($auto_ini{concatenate_input_files} == 1) {
                                                                                                     $concatenate_files->configure(-state => 'normal');
                                                                                                     $conc_lb          ->configure(-state => 'normal');
                                                                                                  } elsif ($auto_ini{concatenate_input_files} == 0) {
                                                                                                     $concatenate_files->configure(-state => 'disabled');
                                                                                                     $conc_lb          ->configure(-state => 'disabled');
                                                                                                  }
                                                                                                 }
                                                                               )
                                    );
      $option_textbox->insert('end', "\n\n\n\n");
      $option_textbox->windowCreate('end', -window => $run_options->Checkbutton(-text     => "\tReplace internal N's with default spacer sequence",
                                                                                -variable => \$auto_ini{replace_N},
                                                                               )
                                    );
      $option_textbox->insert('end', "\n\n\n");

      $option_textbox->insert('end', "\n");
   }

   $run_options->Button(-text       => 'Run annotation',
                        -command    => sub {   if ($auto_ini{reuse_results} == 0) {
                                                  my $check = $mw->Dialog(-title          => 'Caution',
                                                                          -text           => 'All existing results will be deleted. Do you want to continue?',
                                                                          -default_button => 'Yes',
                                                                          -buttons        => ['Yes', 'Cancel'],
                                                                          -bitmap         => 'question'
                                                                         )->Show();
                                                  if ($check eq 'Cancel') {return};
                                               }
                                               #reset Error log
                                               unlink $Bin.'/Error.log';
                                               my $status = &run_analysis(main_window          => \$mw,
                                                                          progress_bar         => \$progress_bar,
                                                                          auto_ini_ref         => \%auto_ini,
                                                                          ini_ref              => \%ini,
                                                                          ini_comment          => \%ini_comment,
                                                                          concatenate_clusters => \%concatenate_clusters,
                                                                          conc_lb              => \$conc_lb
                                                                          );
                                               if ($status == 0) {
                                                  $mw->messageBox(-title => 'Error',
                                                                  -message => "Critical error encountered, aborting analysis",
                                                                  -type => 'OK',
                                                                  -icon => 'info'
                                                                  );
                                               }
                                           }
                        )                                                      ->grid (-row => 3, -column => 0, -sticky => 'w');

   $run_options                       ->Label(-text  => 'Current input files to be analysed')
                                                                               ->grid(-row => 0, -column => 2, -columnspan => 2, -sticky => 'nsew');
   #$select_files = $run_options       ->Button(-text => 'Add files',
   #                                            -command => sub {
   #                                                             &populate(main_window  => \$mw,
   #                                                                       auto_ini_ref => \%auto_ini,
   #                                                                       ini_ref      => \%ini,
   #                                                                       textbox      => $file_lb,
   #                                                                       frame        => \$run_options
   #                                                                      )
   #                                                            }
   #                                           )
   #                                                                            ->grid(-row => 1, -column => 2);
   $refresh_files = $run_options      ->Button(-text => 'Refresh content',
                                               -command => sub {
                                                                #clear cached files
                                                                %all_input_files = ();
                                                                @short_input_files = ();
                                                                $file_lb->delete(0, 'end');

                                                                #get new content
                                                                &grab_input_folder_content(ini_ref      => \%ini,
                                                                                           auto_ini_ref => \%auto_ini,
                                                                                          );

                                                                #repopulate
                                                                foreach my $entry (@short_input_files) {
                                                                   $file_lb->insert('end', $entry);
                                                                }
                                                              }
                                              )
                                                                               ->grid(-row => 1, -column => 2);
   $concatenate_files = $run_options  ->Button(-text    => 'Concatenate files',
                                               -command => sub {
                                                                #get new content
                                                                my ($grouping_ref) = &concatenate_files(main_window  => \$mw,
                                                                                                        ini_ref      => \%ini,
                                                                                                        auto_ini_ref => \%auto_ini,
                                                                                                        textbox      => \$file_lb,
                                                                                                        concbox      => \$conc_lb,
                                                                                                        frame        => \$run_options,
                                                                                                        width        => \$in_file_width
                                                                                                       );
                                                                if ($grouping_ref ne '0') {
                                                                   @concatenate_clusters{ keys %{$grouping_ref} } = values %{$grouping_ref};
                                                                }
                                                               }
                                              )
                                                                               ->grid(-row => 1, -column => 4);

   $file_lb         = $run_options    ->Scrolled("Listbox",
                                                 -width      => 5,
                                                 -selectmode => 'single',
                                                )
                                                                               ->grid(-row => 2, -column => 2, -columnspan => 2, -sticky => 'nsew');

   $conc_lb         = $run_options    ->Scrolled("Text",
                                                 -width      => 1,
                                                 )
                                                                               ->grid(-row => 2, -column => 4, -columnspan => 2, -sticky => 'nsew');
   #$run_options->gridColumnconfigure(0, -weight => 1);
   #$run_options->gridColumnconfigure(0, -minsize => 200);
   $run_options->gridColumnconfigure(2, -weight  => 2);
   #$run_options->gridColumnconfigure(2, -maxsize => 200);
   $run_options->gridColumnconfigure(4, -weight  => 1);

   $clear_selection = $run_options ->Button(-text    => 'Clear current Input folder',
                                            -command => sub {
                                                              #delete files
                                                              while (my ($key, $value) = (each %all_input_files)) {
                                                                 unlink $value;
                                                              }

                                                              %all_input_files = ();
                                                              @short_input_files = ();
                                                              $file_lb->delete(0, 'end');
                                                              $mw->update;
                                                             }
                                           )
                                                                               ->grid(-row => 3, -column => 2, -sticky => 'ew');
   $clear_clusters = $run_options ->Button(-text     => 'Clear current Cluster selection',
                                            -command => sub {
                                                              %concatenate_clusters            = ();
                                                              %concat_file_name                = ();
                                                              $auto_ini{concatenation_counter} = 0;
                                                              $conc_lb->selectAll;
                                                              $conc_lb->deleteSelected;
                                                              $mw->update;
                                                             }
                                           )
                                                                               ->grid(-row => 3, -column => 4, -sticky => 'ew');
   $read_clusters = $run_options ->Button(-text      => 'Read available Clusters',
                                            -command => sub {
                                                              &read_clusters(main_window       => \$mw,
                                                                             ini_ref           => \%ini,
                                                                             auto_ini_ref      => \%auto_ini,
                                                                             textbox           => \$file_lb,
                                                                             concbox           => \$conc_lb,
                                                                             frame             => \$run_options,
                                                                             width             => \$in_file_width,
                                                                             conc_clusters_ref => \%concatenate_clusters
                                                                            );
                                                              $mw->update;
                                                             }
                                           )
                                                                               ->grid(-row => 3, -column => 5, -sticky => 'ew');
   $run_options                       ->Label(-text => 'Double click entry to remove.')
                                                                               ->grid(-row => 3, -column => 3, -sticky => 'ew');
   $run_options->gridRowconfigure    (2, -weight => 1);
   #$run_options->gridColumnconfigure (0, -weight => 1);

   #define default states
   if ($auto_ini{concatenate_input_files} == 1) {
      $concatenate_files->configure(-state => 'normal');
   } elsif ($auto_ini{concatenate_input_files} == 0) {
      $concatenate_files->configure(-state => 'disabled');
   }

   #populate conclb if files are present and re-use is selected
   if ($auto_ini{reuse_results} == 1) {
      &read_clusters(main_window       => \$mw,
                     ini_ref           => \%ini,
                     auto_ini_ref      => \%auto_ini,
                     textbox           => \$file_lb,
                     concbox           => \$conc_lb,
                     frame             => \$run_options,
                     width             => \$in_file_width,
                     conc_clusters_ref => \%concatenate_clusters
                    );
   }

   #initially populate textbox with short files
   $file_lb->delete(0, 'end');
   foreach my $entry (@short_input_files) {
      $file_lb->insert('end', $entry);
   }


   #bind right mouse button to selection and removal of entry
   $file_lb->bind('<Double-ButtonRelease-1>' => sub {&remove_entry(main_window  => \$mw,
                                                                   auto_ini_ref => \%auto_ini,
                                                                   ini_ref      => \%ini,
                                                                   progress_bar => \$progress_bar,
                                                                   textbox      => \$file_lb,
                                                                   frame        => \$run_options
                                                                  )
                                                    });
}

#setup page6: miscallenous
{

   #$db_options->Button(-text    => 'Update databases via internet',
   #                      -command => sub {&update_db(main_window  => \$mw,
   #                                                  ini_ref      => \%ini,
   #                                                  auto_ini_ref => \%auto_ini
   #                                                 );
   #                                      }
   #                     )
   #   ->grid (-row => 0, -column => 0, -sticky => 'ew');
   $db_options->Button(-text    => 'Create custom Blast databases',
                         -padx    => 25,
                         -command => sub {&create_custom_db(main_window  => \$mw,
                                                            ini_ref      => \%ini,
                                                            auto_ini_ref => \%auto_ini);
                                         }
                         )
      ->grid (-row => 1, -column => 0, -sticky => 'ew');
   #$db_options->Button(-text    => 'Create custom COG databases',
   #                      -padx    => 25,
   #                      -command => sub {&create_custom_COG_db(main_window  => \$mw,
   #                                                             ini_ref      => \%ini,
   #                                                             auto_ini_ref => \%auto_ini);
   #                                      }
   #                      )
   #   ->grid (-row => 2, -column => 0, -sticky => 'ew');
   $db_options->gridRowconfigure   (0, -pad => 20);

   $gb_files_options->Button(-text    => 'Rotate single genbank files',
                         -padx    => 25,
                         -command => sub {&rotate_gb(main_window   => \$mw,
                                                      ini_ref      => \%ini,
                                                      auto_ini_ref => \%auto_ini,
                                                      progress_bar	=> \$progress_bar);
                                         }
                         )
      ->grid (-row => 0, -column => 0, -sticky => 'ew');

   $gb_files_options->Button(-text    => 'Prepare Genbank files for Sequin submission',
                         -padx    => 25,
                         -command => sub {&sequin(main_window      => \$mw,
                                                      ini_ref      => \%ini,
                                                      auto_ini_ref => \%auto_ini,
                                                      progress_bar	=> \$progress_bar);
                                         }
                         )
      ->grid (-row => 1, -column => 0, -sticky => 'ew');
   $gb_files_options->gridRowconfigure   (0, -pad => 20);

   #$annotation_options->Button(-text    => 'Transfer TIGRfams to annotation',
   #                      -padx    => 25,
   #                      -command => sub {&TIGRfam_annotation(main_window   => \$mw,
   #                                                           ini_ref       => \%ini,
   #                                                           auto_ini_ref  => \%auto_ini,
   #                                                           progress_bar	=> \$progress_bar);
   #                                      }
   #                      )
   #   ->grid (-row => 0, -column => 0, -sticky => 'ew');

   $annotation_options->Button(-text    => 'Transfer annotation between Genbank files',
                         -padx    => 25,
                         -command => sub {&transfer_annotation(main_window   => \$mw,
                                                               ini_ref       => \%ini,
                                                               auto_ini_ref  => \%auto_ini,
                                                               progress_bar	 => \$progress_bar);
                                         }
                         )
      ->grid (-row => 1, -column => 0, -sticky => 'ew');
   $annotation_options->gridRowconfigure   (0, -pad => 20);

   $metagenome_options->Button(-text    => 'Metagenome analysis',
                         -padx    => 25,
                         -command => sub {&metagenome_ui(main_window   => \$mw,
                                                         ini_ref       => \%ini,
                                                         auto_ini_ref  => \%auto_ini,
                                                         progress_bar  => \$progress_bar);
                                         }
                         )
      ->grid (-row => 0, -column => 0, -sticky => 'ew');
   $metagenome_options->gridRowconfigure   (0, -pad => 20);

   #resize proportionally
   #my ($misc_columns, $misc_rows) = $misc_options->gridSize();
   #for (my $misc_i = 0; $misc_i <= $misc_columns; $misc_i++) {
   #   $misc_options->gridColumnconfigure($misc_i, -weight => 1);
   #}
   #for (my $misc_i = 0; $misc_i <= $misc_rows; $misc_i++) {
   #   $misc_options->gridRowconfigure($misc_i, -weight => 1);
   #}
}

$mw->focus;

MainLoop;

exit;

sub restore_input_files {
   my %args = @_;
   my @input;

   opendir INPUTFILES, ${$args{ini_ref}}{input_files};
   @input = grep /_GAMOLAdna$/, readdir(INPUTFILES);
   closedir INPUTFILES;

   foreach my $file (@input) {
      $file =~ s/_GAMOLAdna$//;
      if (-e ${$args{ini_ref}}{move_msfasta}.'/'.$file) {
         unlink ${$args{ini_ref}}{input_files}.'/'.$file.'_GAMOLAdna';
         move (${$args{ini_ref}}{move_msfasta}.'/'.$file, ${$args{ini_ref}}{input_files}.'/'.$file);
      }
   }
   return (1);
}

sub grab_input_folder_content {
   my %args = @_;
   my @input;

   opendir INPUTFILES, ${$args{ini_ref}}{input_files};
   @input = grep !/^\./, readdir(INPUTFILES);
   closedir INPUTFILES;
   @input = sort(@input);

   foreach my $file (@input) {
      next if (-d ${$args{ini_ref}}{input_files}.'/'.$file); #skip directories

      #test if filename contains any non-acceptable characters, replace with '.', rename file
      if ($file =~ m/[\!\@\#\$\%\^\&\*\(\)\{\[\}\]\,\<\>\/\?\s]/) {
         my $valid_file = $file;
         $valid_file =~ s/[\!\@\#\$\%\^\&\*\(\)\{\[\}\]\,\<\>\/\?\s]/_/g;

         #exists already? try again
         if (-e ${$args{ini_ref}}{input_files}.'/'.$valid_file) {
            $valid_file = $file;
            $valid_file =~ s/[\!\@\#\$\%\^\&\*\(\)\{\[\}\]\,\<\>\/\?\s]/\./g;
         }

         #still exists? Complain and ask for manual intervention
         if (-e ${$args{ini_ref}}{input_files}.'/'.$valid_file) {
            return(0);
         }

         #rename file
         move(${$args{ini_ref}}{input_files}.'/'.$file, ${$args{ini_ref}}{input_files}.'/'.$valid_file);
         $file = $valid_file;
         undef $valid_file;
      }

      push (@short_input_files, $file);
      #configure new width for textbox
      if (length($file) > $in_file_width) {
         $in_file_width = length($file);
      }
      $all_input_files{$file} = ${$args{ini_ref}}{input_files}.'/'.$file;
   }
   @short_input_files = sort(@short_input_files);
   return (1);
}

sub populate {
   my %args = @_;

   my $types = [ ['Fasta files',   ['.fasta', '.msfasta', '.dna', '.seq']],
                 ['Genbank files', ['.gb', '.gbk', '.genbank']],
                 ['All Files',      '*'],
               ];

   (my $selected_files) = ${$args{main_window}}->getOpenFile(-initialdir => $current_in_dir,
                                                             -title      => 'Select files for custom database',
                                                             -multiple   => 1,
                                                             -filetypes => $types,
                                                             -defaultextension => '');

   if (defined(@{$selected_files})) {
      #generate short_list
      foreach my $file (@{$selected_files}) {
         'reset' =~ m/reset/;
         $file =~ m/(.+)\/([^\/]+)$/; #display only file_name
         $current_in_dir = $1;
         my $short = $2;
         push (@short_input_files, $short);
         #configure new width for textbox
         if (length($short) > $in_file_width) {
            $in_file_width = length($short);
         }

         #copy selected files into input folder if not already there
         if ($current_in_dir ne ${$args{ini_ref}}{input_files} && !-e ${$args{ini_ref}}{input_files}.'/'.$short) {
            copy ($file, ${$args{ini_ref}}{input_files}.'/'.$short);
         }

      }

      #remove duplicates
      &delete_duplicates(array => \@short_input_files);

      #sort input list
      @short_input_files = sort(@short_input_files);


      #populate textbox with short files
      $args{textbox}->delete(0, 'end');
      foreach my $entry (@short_input_files) {
         $args{textbox}->insert('end', $entry);
      }

      #merge selected files with overall selection array
      foreach (@{$selected_files}) {
         'reset' =~ m/reset/;
         m/\/([^\/]+)$/; #display only file_name
         my $short = $1;
         $all_input_files{$short} = ${$args{ini_ref}}{input_files}.'/'.$short;
      }

      #update display
      ${$args{frame}}->gridColumnconfigure(2, -minsize => ($in_file_width));
      $args{textbox} ->configure             (-width   => $in_file_width);
      ${$args{main_window}}->update;
   }
}

sub remove_entry {
   my %args = @_;
   return if ($#short_input_files < 0);

   my (@current_index) = ${$args{textbox}}->curselection();
   #remove from shortlist array
   splice(@short_input_files, $current_index[0], 1);

   #determine new width
   $in_file_width = 0;
   foreach my $entry (@short_input_files) {
      if (length($entry) > $in_file_width) {
         $in_file_width = length($entry);
      }
   }
   #remove from all files hash
   my $sel_key = ${$args{textbox}}->get($current_index[0]);
   foreach my $key (keys %all_input_files) {
      if ($key eq $sel_key) {
         delete $all_input_files{$key};
         ${$args{progress_bar}}->configure(-label=>"Moved entry $key to ${$args{auto_ini_ref}}{work_dir}");
         ${$args{main_window}}->update;
         #delete file
         #unlink ${$args{ini_ref}}{input_files}.'/'.$sel_key;
         #rather than delete, move to work directory
         move (${$args{ini_ref}}{input_files}.'/'.$sel_key , ${$args{auto_ini_ref}}{work_dir}.'/'.$sel_key);
         last;
      }
   }
   #remove from display
   ${$args{textbox}}->delete($current_index[0]);
   ${$args{main_window}}->update;

   #update display
   ${$args{frame}}      ->gridColumnconfigure(2, -minsize => ($in_file_width));
   ${$args{textbox}}    ->configure             (-width   => $in_file_width);
   ${$args{main_window}}->update;

}

sub delete_duplicates {
   my %args = @_;
   my %seen = ();
   @{$args{array}} = grep { ! $seen{$_} ++ } @{$args{array}};
}

sub concatenate_files {
   my %args = @_;
   my $wait;
   my @selection;
   my %grouping = ();
   #no files in listbox yet?
   return if ($#short_input_files < 0);

   #default color range for file groupings
   my %colour = (0=>'red', 1=>'orange', 2=>'dark green', 3=>'blue');


   #make new listbox in new pop up
   my $tl = ${$args{main_window}}->DialogBox(-title          => 'Select files to concatenate',
                                             -buttons        => [ 'OK', 'Cancel' ],
                                             -default_button => 'OK'
                                            );
   my $conc_lb = $tl->Scrolled("Listbox",
                               -setgrid    => 1,
                               -selectmode => 'extended',
                               -width      => ${$args{width}}
                               )->pack();
   #populate listbox
   $conc_lb->delete(0, 'end');
   foreach my $entry (@short_input_files) {
      $conc_lb->insert('end', $entry);
   }

   $wait = $tl->Show();

   #get selection-> return a list of indeces
   (@selection) = $conc_lb->curselection();

   if ($wait eq 'Cancel') {
      $tl->state('withdrawn');
      return(0);
   } elsif ($wait eq 'OK') {
      #something selected?
      if ($#selection >= 0) {
         #increase colour count
         $args{auto_ini_ref}{concatenation_counter}++;
         #define new max width
         unless (defined $conc_lb_max_width) {$conc_lb_max_width = 20};
         my $local_conc_width = 0;

         #iterate through selected files and colour accordingly; add to refernce hash for later
         foreach my $entry (@selection) {
            #set colour cycling
            my $color_number = $args{auto_ini_ref}{concatenation_counter} % 4;
            #set new max width
            $local_conc_width = length($short_input_files[$entry]);
            if ($local_conc_width > $conc_lb_max_width) {$conc_lb_max_width = $local_conc_width};

            #insert cluster into listbox
            ${$args{concbox}}->tagConfigure('tag'.$args{auto_ini_ref}{concatenation_counter}, -foreground => $colour{$color_number}, -font => "Courier 12 bold");
            ${$args{concbox}}->insert('end', $short_input_files[$entry]."\n", 'tag'.$args{auto_ini_ref}{concatenation_counter});

            #capture cluster into hash
            $grouping{$args{auto_ini_ref}{concatenation_counter}} .= $short_input_files[$entry]."\n";

            #set name hash to determine which file name to use. if all else fails, use appendices
            my @temp = split /\n/, $short_input_files[$entry];
            foreach (@temp) {
               unless (defined $concat_file_name{$_}) {
                  $concat_file_name{$_} = 0;
               }
            }
         }

         #determine proper file name for cluster. By default use first entry, then iterate through rest of files, then use appendices
         my $total_count = keys %concat_file_name;
         my $local_count = 0;
         my $chosen_name = '';
         #easy way, choose one of the selected files
         foreach my $entry (sort keys %concat_file_name) {
            $local_count++;
            $concat_file_name{$entry}++;
            if ($concat_file_name{$entry} == 1 && $chosen_name !~ /\w+/) {
               $chosen_name = $entry;
            }
         }
         #if already taken, increase local counter until unique name has been found.
         if ($local_count == $total_count) {
            $local_count = 0;
            while ($chosen_name !~ /\w+/) {
               foreach my $entry (sort keys %concat_file_name) {
                  $local_count++;
                  $concat_file_name{$entry.'_'.$local_count}++;
                  if ($concat_file_name{$entry.'_'.$local_count} == 1) {
                     $chosen_name = $entry.'_'.$local_count;
                     last;
                  }
               }
            }
         }

         $grouping{$args{auto_ini_ref}{concatenation_counter}} = "File_name:\t$chosen_name\n".$grouping{$args{auto_ini_ref}{concatenation_counter}};

         #set new width
         ${$args{concbox}}    ->configure          (   -width   => $conc_lb_max_width);
         ${$args{frame}}      ->gridColumnconfigure(4, -minsize => ($conc_lb_max_width * 10));
         ${$args{concbox}}    ->update;
         ${$args{main_window}}->update;
      }
      $tl->state('withdrawn');
      ${$args{main_window}}->update;
      return(\%grouping);
   }
}

sub read_clusters {
   my %args = @_;
   my (@input);

   #default color range for file groupings
   my %colour = (0=>'red', 1=>'orange', 2=>'dark green', 3=>'blue');

   #first, clear out existing selections
   ${$args{concbox}}->selectAll;
   ${$args{concbox}}->deleteSelected;

   #delete existing selection hash
   %{$args{conc_clusters_ref}} = ();

   #reset counter
   $args{auto_ini_ref}{concatenation_counter} = 0;

   #read available cluster files
   opendir INPUT, ${$args{auto_ini_ref}}{work_dir};
   @input = grep /\.cluster_files$/, readdir(INPUT);
   closedir INPUT;

   #increase colour count
   $args{auto_ini_ref}{concatenation_counter}++;
   #define new max width
   unless (defined $conc_lb_max_width) {$conc_lb_max_width = 20};
   my $local_conc_width = 0;

   #iterate through selected files and colour accordingly; add to refernce hash for later
   foreach my $entry (@input) {
      #capture cluster number
      'reset' =~ m/reset/;
      $entry =~ m/^(\d+)\./;
      my $c_number = $1;

      #read file content
      my ($file_ref) = &slurp(main_window   => $args{main_window},
                              auto_ini_ref  => $args{auto_ini_ref},
                              ini_ref       => $args{ini_ref},
                              filename      => $entry,
                              directory     => ${$args{auto_ini_ref}}{work_dir}
                             );

      if ($file_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Cannot read Contig_order file:\n$entry.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info'
                                                       );
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error concatenating input file ".
                        "\nCannot read Contig_order file:\n$entry in GAMOLA main module.\n\n";
         close ERRORLOG;
      }
      #set colour cycling
      my $color_number = $args{auto_ini_ref}{concatenation_counter} % 4;

      #set new max width
      my @temp = split/\n/, $$file_ref;
      foreach my $file (@temp) {
         next if ($file =~ m/^File_name:/);
         $local_conc_width = length($file);
         if ($local_conc_width > $conc_lb_max_width) {$conc_lb_max_width = $local_conc_width};
      }
      undef @temp;

      #insert cluster into listbox
      my $insert = $$file_ref;
      $insert =~ s/^File_name\:.+?\n//s;
      ${$args{concbox}}->tagConfigure('tag'.$args{auto_ini_ref}{concatenation_counter}, -foreground => $colour{$color_number}, -font => "Courier 12 bold");
      ${$args{concbox}}->insert('end', $insert."\n", 'tag'.$args{auto_ini_ref}{concatenation_counter});

      #capture cluster into hash
      ${$args{conc_clusters_ref}}{$c_number} = $$file_ref;

      $args{auto_ini_ref}{concatenation_counter}++;
      undef $file_ref;
      undef $insert;
   }
   #set new width
   ${$args{concbox}}    ->configure          (   -width   => $conc_lb_max_width);
   ${$args{frame}}      ->gridColumnconfigure(4, -minsize => ($conc_lb_max_width * 10));
   ${$args{concbox}}    ->update;
   ${$args{main_window}}->update;

   return;
}

sub welcome {
   my $hFont = '-*-helvetica-Medium-R-Normal--*-140-*-*-*-*-*-*';
   my $text = "GAMOLA2: Global Annotation of Multiplexed On-site bLasted DNA-sequences\n\n".
              "GAMOLA2 is a local annotation system for prokaryotic draft and complete genomes.\n\n".
              "Expert curated annotation remains one of the critical steps in achieving a\n".
              "reliable biological relevant annotation. \n".
              "GAMOLA2 is a user friendly and comprehensive software package to process, \n".
              "annotate and curate draft and complete bacterial and viral genomes.\n".
              "GAMOLA2 represents a wrapping tool to combine functional Blast, COG, Pfam, and \n".
              "TIGRfam analyses with structural predictions including detection of tRNAs, rRNA genes,\n".
              "non-coding RNAs, signal protein cleavage sites, transmembrane helices, CRISPR repeats and\n".
              "vector sequence contaminations.\n\n".
              "Results can best be viewed in the customised version of Sanger's Artemis,\n".
              "\(http:\/\/www.sanger.ac.uk\/science\/tools\/artemis\)".
              " that is provided with the GAMOLA2 software package.\n\n".
              "References for all external software packages can be viewed using the \"Reference\" button at the bottom of the main screen";

   my $welcome = $mw->Toplevel();
   #$welcome->minsize(400,400);
   #$welcome->maxsize(600,600);
   #$welcome->resizable(0,0);
   $welcome->title( 'Welcome to GAMOLA 2' );
   my $bitmap  = $welcome->Photo (-file => $Bin.'/lib/Archives/Gamola.gif');
   my $bottom_frame        = $welcome->Frame()-> pack (-side      => 'bottom',
                                                       -fill      => 'x',
                                                       -expand    => 'yes');
   my $welcome_left_frame  = $welcome->Frame()-> pack (-side      => 'left',
                                                       -anchor    => 'w');
   my $welcome_right_frame = $welcome->Frame()-> pack (-side      => 'right',
                                                       -expand    => 'yes',
                                                       -anchor    => 'e',
                                                       -fill      => 'both');

   $welcome_left_frame->Label(-image => $bitmap) -> pack (-side   => 'left',
                                                          -anchor => 'nw');
   my $message = $welcome_right_frame->Scrolled( 'ROText',
                                                -width  => 80,
                                                -height => 25,
                                                -font   => $hFont,
                                                -wrap   => 'word',
                                                -scrollbars  => 're',
                                              ) -> pack( -expand => 'yes', -fill => 'both' );
   $bottom_frame->Button(
                    -text    => 'OK',
                    -command => sub { $welcome->destroy() },
                   )->pack();
   $message->insert( 'end', join( '', $text ) );
   $welcome->focus();




}

sub view_references {
   #read reference list file
   (my $file_ref) = slurp(main_window => \$mw,
                          directory   => $auto_ini{work_dir}.'/lib/Archives',
                          filename    => 'References.ref'
                         );
   if ($file_ref eq '0') {
      return (0)
   };

   $$file_ref =~ s/[\r]+//gs;

   my $top_ref = $mw->Toplevel(-title=>'Cited References for external programs used in GAMOLA 2');
   my $text = $top_ref->Scrolled('ROText',
                                 -width  => 80,
                                 -height => 25,
                                 -wrap   => 'word',
                                 -font       => '-*-helvetica-Medium-R-Normal--*-140-*-*-*-*-*-*',
                                 -scrollbars => 'e',
                                )-> pack(-expand      => 1,
                                         -fill        => 'both'
                                         );
   $text->insert('end', $$file_ref);

   $top_ref -> Button (     -text => 'Close',
                            -command => sub {$top_ref -> withdraw}) -> pack();
   return;
}

