#!/opt/ActivePerl-5.8/bin/perl

#sequin: generate fasta and table format to submit via sequin
#input arguments: main_window, ini_ref, auto_ini_ref,

#optional agp_mapping file
#can have comment, strain line with '#'
#individual contigs with a scaffold must have same same as in Genbank file
#contigs must be in same order as in Genbank file
#scaffold boundaries are depicted by empty lines or lines containing '>'



package Miscellaneous::sequin;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&sequin);
use vars qw();

use Tk;
use Tk::Pane;
use initialise::read_me              qw(:DEFAULT);
use ProgrammeModules::sequence       qw(:DEFAULT);
use ProgrammeModules::genbank_parser qw(:DEFAULT);
use Basics::progress_bar             qw(:DEFAULT);
use Basics::Roman                    qw(:DEFAULT);
use Cwd;

#local variables
my (%args, $tl, @selected_files, %all_files, @short_files, $width, $current_dir,
    %seen_gene_annotation, %seen_CDS_annotation, $seen_gene_anno_count, $seen_CDS_anno_count);
our (@start, @stop, $start, $stop, $table_feature, @table_entry, $lineform);

#initialise counter in case of duplicated gene names
$seen_gene_anno_count = 2;
$seen_CDS_anno_count  = 2;

#setting up format
#format TABLE =
#^<<<<<<<<<~~ ____  ^<<<<<<<<~~ ____ ^<<<<<<<<<<<<<  ____  ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<~~
#shift(@start),   shift(@stop),  $table_feature,     shift (@table_entry)
#.

format FEATURES =
^<<<<<<<<<~~ ____ ^<<<<<<<<~~ ____ ^<<<<<<<<<<<<<<<<<<<<<<<
shift(@start),   shift(@stop),  $table_feature
.

format LINEFORM =
^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<~~
$lineform
.

#remember last directory for adding files
$current_dir = ${$args{auto_ini_ref}}{work_dir};

sub sequin {
   my %args = @_;
   my (@features, %selected_feature);
   my ($individual_files, $cluster_files, $cluster, $individual, $select_files,
       $file_lb, $locus_tag, $progress_bar, $frame_left, $frame_right, $frame_middle,
       $compile, $clear_selection, $pane, $reuse_locustag, $db_name, $project_name,
       $reset_locusid, $notes, $link_files, $file_combine, $separate_files,
       $file_individual, $force_product, $get_qualifiers, $qualifiers_ref,
       $genome_submission, $WGS_submission, $genome, $WGS, $agp_submission,
       $agp, $agp_setup, $agp_map_file, $gap, $no_gene_feature);


   $width = 40;

   #defined possible features for sequin
   @features = qw (attenuator C_region CAAT_signal CDS conflict D-loop D_segment
                   enhancer exon gap GC_signal gene iDNA intron J_segment LTR
                   mat_peptide misc_binding misc_difference misc_feature misc_recomb
                   misc_RNA misc_signal misc_structure modified_base mRNA ncRNA
                   N_region operon oriT polyA_signal polyA_site precursor_RNA
                   prim_transcript primer_bind promoter protein_bind RBS repeat_region
                   repeat_unit rep_origin rRNA S_region satellite sig_peptide
                   source stem_loop STS TATA_signal terminator tmRNA transit_peptide
                   tRNA unsure V_region V_segment variation 3'UTR 5'UTR -10_signal
                   -35_signal PFAM_match TIGR_match COG_match TMM signalP
                  );

   #set a few default values
   $selected_feature{'CDS'}  = 1;
   $selected_feature{'gene'} = 1;
   $cluster                  = 1;
   $genome_submission        = 1;
   $WGS_submission           = 0;
   $agp_submission           = 0;
   $gap                      = 100 - 32; #NCBI default unknown gap size - GAMOLA default contig spacer length
   $no_gene_feature          = 0;
   $reset_locusid            = '1';
   $file_individual          = '1';
   $force_product            = '0';
   $get_qualifiers           = '0';
   $db_name                  = 'ncbi';
   $project_name             = 'Genbank_submission';

   if (defined $tl && Tk::Exists($tl)) {
      $tl->state('normal');
      $progress_bar = $tl->Frame(-borderwidth => 2, -relief => 'groove');
      $frame_left   = $tl->Frame()                                      ;
      $frame_middle = $tl->Frame(-relief => 'raised')                   ;
      $frame_right  = $tl->Frame(-relief => 'groove')                   ;
   } else {
      $tl = ${$args{main_window}}->DialogBox(-title   => 'Prepare submission to NCBI/Genbank',
                                             -buttons => [ 'Exit' ],
                                            );
      $frame_left   = $tl->add("Frame", -borderwidth => 1, -relief => 'raised') -> grid(-row => 0, -column => 0,                   -sticky => 'nsw',  -ipady => 0);
      $frame_middle = $tl->add("Frame", -borderwidth => 1, -relief => 'raised') -> grid(-row => 0, -column => 1,                   -sticky => 'ns',   -ipady => 2);
      $frame_right  = $tl->add("Frame", -borderwidth => 1, -relief => 'raised') -> grid(-row => 0, -column => 2,                   -sticky => 'nse',  -ipady => 0);
      $progress_bar = $tl->add("Frame", -borderwidth => 2, -relief => 'groove') -> grid(-row => 1, -column => 0, -columnspan => 3, -sticky => 'swe',  -ipady => 0);
   }

   $progress_bar->configure(-label => " "); #create dummy for space
   $frame_right ->gridRowconfigure    (1, -weight => 1);
   $frame_right ->gridColumnconfigure (1, -weight => 1);

   #options into left frame
   {
      $frame_left                  -> Label (-text => "Options for Sequin Parser\n\n"
                                            )->grid(-row => 0, -column => 0, -columnspan => 2, -sticky => 'nsew');
      $agp       = $frame_left -> Checkbutton(-text     => '   Create AGP scaffold information',
                                              -variable => \$agp_submission,
                                             )->grid(-row => 3, -column => 0, -columnspan => 2, -sticky => 'w');
      $agp_setup = $frame_left -> Button     (-text     => 'Setup',
                                              -command  => sub {$agp_map_file = &agp_setup(main_window       => $args{main_window},
                                                                                           progress_bar      => \$progress_bar,
                                                                                           ini_ref           => $args{ini_ref},
                                                                                           auto_ini_ref      => $args{auto_ini_ref}
                                                                                          );
                                                                if (-e $agp_map_file) {
                                                                   $agp_setup -> configure(-foreground => 'green');
                                                                } else {
                                                                   $agp_setup -> configure(-foreground => 'red');
                                                                }

                                                               }
                                             )->grid(-row => 3, -column => 2, -sticky => 'w');

      #default state for agp
      $agp       -> configure(-state      => 'disabled');
      $agp_setup -> configure(-state      => 'disabled');
      $agp_setup -> configure(-foreground => 'red');

      $genome = $frame_left -> Checkbutton(-text     => 'Submit complete genome',
                                                  -variable => \$genome_submission,
                                                  -command  => sub {
                                                                     if ($genome_submission == 1) {
                                                                        $WGS_submission = 0;
                                                                        $agp       -> configure(-state => 'disabled');
                                                                        $agp_setup -> configure(-state => 'disabled');
                                                                     } elsif ($genome_submission == 0) {
                                                                        $WGS_submission = 1;
                                                                        $agp       -> configure(-state => 'normal');
                                                                        $agp_setup -> configure(-state => 'normal');
                                                                     }
                                                                   }
                                          )->grid(-row => 1, -column => 0, -columnspan => 2, -sticky => 'w');
      $WGS   = $frame_left -> Checkbutton(-text     => 'Submit WGS shotgun sequence',
                                                  -variable => \$WGS_submission,
                                                  -command  => sub {
                                                                     if ($WGS_submission == 1) {
                                                                        $genome_submission = 0;
                                                                        $agp       -> configure(-state => 'normal');
                                                                        $agp_setup -> configure(-state => 'normal');
                                                                     } elsif ($WGS_submission == 0) {
                                                                        $genome_submission = 1;
                                                                        $agp       -> configure(-state => 'disabled');
                                                                        $agp_setup -> configure(-state => 'disabled');
                                                                     }
                                                                   }
                                          )->grid(-row => 2, -column => 0, -columnspan => 2, -sticky => 'w');

      $cluster_files = $frame_left -> Checkbutton(-text     => 'Link selected files in Locus_tags',
                                                  -variable => \$cluster,
                                                  -command  => sub {
                                                                     if ($cluster == 1) {
                                                                        $individual = 0;
                                                                     } elsif ($cluster == 0) {
                                                                        $individual = 1;
                                                                     }
                                                                   }
                                          )->grid(-row => 4, -column => 0, -columnspan => 2, -sticky => 'w');

      $individual_files = $frame_left -> Checkbutton(-text     => 'Treat selected files individually in Locus_tags',
                                                     -variable => \$individual,
                                                     -command  => sub {
                                                                        if ($individual == 1) {
                                                                           $cluster = 0;
                                                                        } elsif ($individual == 0) {
                                                                           $cluster = 1;
                                                                        }
                                                                      }
                                             )->grid(-row => 5, -column => 0, -columnspan => 2, -sticky => 'w');

      $frame_left                 -> Checkbutton (-text     => 'Re-use existing Locus_tag numbers',
                                                  -variable => \$reuse_locustag
                                                  )->grid(-row => 6, -column => 0, -columnspan => 2, -sticky => 'w');
      $frame_left                 -> Checkbutton (-text     => 'Reset Locus_tag numbers for each file',
                                                  -variable => \$reset_locusid
                                                  )->grid(-row => 7, -column => 0, -columnspan => 2, -sticky => 'w');
      $frame_left                 -> Checkbutton (-text     => 'Force product qualifier in CDS',
                                                  -variable => \$force_product
                                                  )->grid(-row => 8, -column => 0, -columnspan => 2, -sticky => 'w');
      $frame_left                 -> Checkbutton (-text     => 'Allow additional CDS qualifiers',
                                                  -variable => \$get_qualifiers,
                                                  -command  => sub {
                                                                    if ($get_qualifiers == 1) {
                                                                       $qualifiers_ref = &determine_CDS_qualifiers(main_window       => $args{main_window},
                                                                                                                   ini_ref           => $args{ini_ref},
                                                                                                                   auto_ini_ref      => $args{auto_ini_ref},
                                                                                                                   );
                                                                    }
                                                                   }
                                                  )->grid(-row => 9, -column => 0, -columnspan => 2, -sticky => 'w');
      $frame_left                 -> Checkbutton (-text     => 'Allow notes in CDS features',
                                                  -variable => \$notes
                                                  )->grid(-row => 10, -column => 0, -columnspan => 2, -sticky => 'w');
      $frame_left                 -> Checkbutton (-text     => "Remove \'gene\' qualifier from \'gene\' feature",
                                                  -variable => \$no_gene_feature
                                                  )->grid(-row => 11, -column => 0, -columnspan => 2, -sticky => 'w');
      $frame_left                 ->LabEntry(-label        => "Enter Locus tag identifier  ",
                                             -labelPack    => [-side  => "left"],
                                             -width        => 20,
                                             -textvariable => \$locus_tag
                                            ) -> grid(-row => 12, -column => 0, -columnspan => 2, -sticky => 'w');
      $frame_left                 ->LabEntry(-label        => "Enter dbname string         ",
                                             -labelPack    => [-side  => "left"],
                                             -width        => 20,
                                             -textvariable => \$db_name
                                            ) -> grid(-row => 13, -column => 0, -columnspan => 2, -sticky => 'w');
      $frame_left                 ->LabEntry(-label        => "Enter Project name          ",
                                             -labelPack    => [-side  => "left"],
                                             -width        => 20,
                                             -textvariable => \$project_name
                                            ) -> grid(-row => 14, -column => 0, -columnspan => 2, -sticky => 'w');

      $compile   = $frame_left ->Button(-text => 'Convert to Sequin',
                                        -command => sub {
                                                          if ($genome_submission == 1) {
                                                             &compile(main_window       => $args{main_window},
                                                                      progress_bar      => \$progress_bar,
                                                                      ini_ref           => $args{ini_ref},
                                                                      auto_ini_ref      => $args{auto_ini_ref},
                                                                      locus_tag         => $locus_tag,
                                                                      cluster           => $cluster,
                                                                      reuse_locustag    => $reuse_locustag,
                                                                      reset_locusid     => $reset_locusid,
                                                                      force_product     => $force_product,
                                                                      qualifiers_ref    => $qualifiers_ref,
                                                                      notes             => $notes,
                                                                      db_name           => $db_name,
                                                                      project_name      => $project_name,
                                                                      file_individual   => $file_individual,
                                                                      selected_feature  => \%selected_feature,
                                                                      no_gene_feature   => $no_gene_feature
                                                                     );
                                                          } else {
                                                             &compile_WGS(main_window          => $args{main_window},
                                                                          progress_bar         => \$progress_bar,
                                                                          ini_ref              => $args{ini_ref},
                                                                          auto_ini_ref         => $args{auto_ini_ref},
                                                                          locus_tag            => $locus_tag,
                                                                          cluster              => $cluster,
                                                                          reuse_locustag       => $reuse_locustag,
                                                                          reset_locusid        => $reset_locusid,
                                                                          force_product        => $force_product,
                                                                          qualifiers_ref       => $qualifiers_ref,
                                                                          notes                => $notes,
                                                                          db_name              => $db_name,
                                                                          project_name         => $project_name,
                                                                          file_individual      => $file_individual,
                                                                          create_scaffold_info => $agp_submission,
                                                                          agp_map_file         => $agp_map_file,
                                                                          gap                  => $gap,
                                                                          selected_feature     => \%selected_feature,
                                                                          no_gene_feature   => $no_gene_feature
                                                                         );
                                                          }
                                                        }
                                          )->grid(-row => 15, -column => 0, -sticky => 'se');
      $frame_left->gridRowconfigure   (15, -weight => 1);
   }

   #selected fields for middle frame
   {
      $frame_middle         -> Label(-text => "Select features to include\n") -> pack (-side => 'top');
      $pane = $frame_middle -> Scrolled (qw/Pane -scrollbars osw/)            -> pack (-side => 'top', -anchor => 'w', -fill => 'both', -expand => 1);
      foreach my $field (sort @features) {
         $pane -> Checkbutton(-text     => $field,
                              -variable => \$selected_feature{$field}
                             ) -> pack (-side   => 'top',
                                        -anchor => 'w'
                                        );
      }
      $frame_middle         -> Button (-text => 'Select all',
                                       -command => sub {
                                                         foreach my $field (@features) {
                                                            $selected_feature{$field} = 1;
                                                         }
                                                         $pane->update;
                                                        }
                                      )                                       -> pack (-side => 'left', -expand => 1, -fill => 'x');
      $frame_middle         -> Button (-text => 'Deselect all',
                                       -command => sub {
                                                         foreach my $field (@features) {
                                                            $selected_feature{$field} = 0;
                                                         }
                                                         $pane->update;
                                                        }
                                      )                                       -> pack (-side => 'left', -expand => 1, -fill => 'x');
   }

   #Textbox and button in right frame
   {
      $select_files = $frame_right       ->Button(-text => 'Add files',
                                                  -command => sub {&populate(main_window  => \$tl,
                                                                             auto_ini_ref => $args{auto_ini_ref},
                                                                             textbox      => $file_lb,
                                                                             frame_right  => $frame_right
                                                                             )}
                                                 )                             ->grid(-row => 0, -column => 0, -sticky => 'e');
      $frame_right                       ->Label(-text => 'Double click entry to remove.')
                                                                               ->grid(-row => 0, -column => 1, -sticky => 'w');
      $clear_selection = $frame_right ->Button(-text => 'Clear current selection',
                                               -command => sub {
                                                                 %all_files = ();
                                                                 @short_files = ();
                                                                 $file_lb->delete(0, 'end');
                                                                 $tl->update;
                                                                }
                                              )                                ->grid(-row => 2, -column => 0, -columnspan => 2, -sticky => 'ew');
      $link_files     = $frame_right -> Checkbutton(-text     => 'Save as combined file',
                                                    -variable => \$file_combine,
                                                    -command  => sub {
                                                                       if ($file_combine == 1) {
                                                                          $file_individual = 0;
                                                                       } elsif ($file_combine == 0) {
                                                                          $file_individual = 1;
                                                                       }
                                                                     }
                                                   )                           ->grid(-row => 3, -column => 0, -sticky => 'w');

      $separate_files = $frame_right -> Checkbutton(-text     => 'Save individual files',
                                                    -variable => \$file_individual,
                                                    -command  => sub {
                                                                       if ($file_individual == 1) {
                                                                          $file_combine = 0;
                                                                       } elsif ($file_individual == 0) {
                                                                          $file_combine = 1;
                                                                       }
                                                                     }
                                                   )                           ->grid(-row => 3, -column => 1, -sticky => 'e');
      $file_lb         = $frame_right    ->Scrolled("Listbox",
                                                    -setgrid    => 0,
                                                    -selectmode => 'single',
                                                    -width      => $width
                                                    )                          ->grid(-row => 1, -column => 0, -columnspan => 2, -sticky => 'nsew');
   }
   #bind right mouse button to selection and removal of entry
   $file_lb->bind('<Double-ButtonRelease-1>' => sub {&remove_entry(textbox      => \$file_lb,
                                                                   main_window  => \$tl,
                                                                   progress_bar => \$progress_bar,
                                                                   ini_ref      => $args{ini_ref},
                                                                   auto_ini_ref => $args{auto_ini_ref}
                                                                  )
                                                    });
   my $wait = $tl->Show();
   if ($wait eq 'OK') {
      $file_lb->delete(0, 'end');
      $tl->update;
      undef @selected_files;
      undef %all_files;
      undef @short_files;
      undef $current_dir;
      undef %selected_feature;
      $tl->state('withdrawn');
   }
}

sub determine_CDS_qualifiers {
   my %args = @_;
   my ($pane, $db, @CDS_qualifiers, @selected_qualifier);

   @CDS_qualifiers = qw(prot_desc function EC_number
                        experiment inference go_component go_process
                        go_function db_xref pseudo exception transl_except
                        codon_start transl_table
                        );

   if (defined $db && Tk::Exists($db)) {
       $db->state('normal');
   } else {
      $db = ${$args{main_window}}->DialogBox(-title   => 'CDS qualifiers',
                                             -buttons => [ 'OK', 'Select all' ],
                                            );
   }

   $pane = $db->Scrolled(qw/Pane -scrollbars osw -height 280 -width 120/)->pack(-expand => 1, -fill => 'both');

   #set default value if not already defined
   foreach my $entry (sort @CDS_qualifiers) {
      unless (exists $args{selection}->{$entry} && $args{selection}->{$entry} == 1) {
         $args{selection}->{$entry} = 0;
      }
      $pane->Checkbutton(-text     => $entry,
                         -variable => \$args{selection}->{$entry}
                        )->pack(-side => 'top', -anchor => 'w');
   }

   my $wait = $db->Show();
   if ($wait eq 'OK') {
      while (my ($key, $value) = each %{$args{selection}}) {
         if ($value == 1) {
            push (@selected_qualifier, $key);
         }
      }
      $db->state('withdrawn');
   } elsif ($wait eq 'Select all') {
      foreach my $key (keys %{$args{selection}}) {
         $args{selection}->{$key} = 1;
      }
      &determine_CDS_qualifiers(main_window       => $args{main_window},
                                ini_ref           => $args{ini_ref},
                                auto_ini_ref      => $args{auto_ini_ref},
                                selection         => $args{selection}
                                );
      $db->state('withdrawn');
   }
   return (\@selected_qualifier);
}

sub compile {
   my %args = @_;
   my ($file_counter, $counter, $id_counter, @gb_order, $selected, $max_count, $local_path,
       $filename, @non_standard, %seen, @standard);

   #clear all variables
   $file_counter = '';
   $counter      = '';
   $id_counter   = '';
   $selected     = '';
   $max_count    = '';
   $local_path   = '';
   $filename     = '';
   @gb_order     = ();
   @non_standard = ();
   @standard     = ();
   %seen         = ();

   #create directory if not already there
   unless (-d ${$args{ini_ref}}{Sequin_submission}) {
      mkdir (${$args{ini_ref}}{Sequin_submission}, 0777) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create Sequin submission directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return;
      };
   }

   #create combined files
   if ($args{file_individual} == 0) {
      unlink ${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'.fasta';
      open FASTA, ">>".${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'.fasta' or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create Fasta Sequin file",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return;
      };

      unlink ${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'.tbl';
      open FEATURE, ">>".${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'.tbl' or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create Sequin Table file",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return;
      };

      #create agp file if selected
      if ($args{create_scaffold_info} == 1) {
         unlink ${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'.agp';
         open AGP, ">>".${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'.agp' or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot create Scaffold info file in directory $args{project_name}.",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return;
         };
      }
   }

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Processing Genbank files",
                   label        => 'Processing features'
                  );
   &show_pbar_2;

   #get defined features into search pattern
   $selected = '(';
   foreach my $key (keys (%{$args{selected_feature}})) {
      if ($args{selected_feature}{$key}) {
         $selected .= $key.'|';
         push (@standard, $key);
         #add _ to ~ replacement just in case
         if ($key =~ m/_/) {
            my $mod_key = $key;
            $mod_key =~ s/_/~/;
            $selected .= $mod_key.'|';
            push (@standard, $mod_key);
         }
      };
   }
   $selected =~ s/\|$//;
   $selected .= ')';

   #get standard descriptor of rest of selected features
   @non_standard = qw (gene mRNA CDS 5'UTR 3'UTR intron exon tRNA rRNA ncRNA misc_RNA misc~RNA);
   foreach my $entry (@standard, @non_standard) {$seen{$entry}++};
   @standard     = ();
   @non_standard = ();
   foreach my $entry (keys %seen) {
      if ($seen{$entry} == 1) {
         push (@standard, $entry);
      }
   }
   %seen = ();

   #split Genbank key order into array
   @gb_order = split /;/,${$args{auto_ini_ref}}{genbank_key_order};

   #iterate through each selected file
   $file_counter = 1;
   $id_counter   = 1;
   foreach my $file (@short_files) {
      my ($roman, %locus_counter, $fasta_header, $table_header, $max_number);
      #reset variables
      $roman         = '';
      $fasta_header  = '';
      $table_header  = '';
      $max_number    = '';
      %locus_counter = ();

      unless (exists $all_files{$file}) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "File $all_files{$file} cannot be accessed",
                                           -type    => 'OK',
                                           -icon    => 'error');
         &hide_pbar_2;
         return (0);
      }
      ${$args{progress_bar}} ->configure(-label => "Processing $file");
      ${$args{main_window}}  ->update;

      #pull filepath apart
      'reset' =~ m/reset/;
      $all_files{$file} =~ m/^(.+)\/([^\/]+)$/;
      $local_path = $1;
      $filename = $2;

      #create individual projects if selected
      if ($args{file_individual} == 1) {
         #create directory if not already there
         unless (-d ${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}) {
            mkdir (${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}, 0777) or do {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Cannot create individual Sequin project directory",
                                                 -icon    => 'error',
                                                 -type    => 'ok');
               &hide_pbar_2;
               return;
            };
         }

         #create individual files
         unlink ${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'/'.$filename.'.fasta';
         open FASTA, ">>".${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'/'.$filename.'.fasta' or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot create Fasta Sequin file",
                                              -icon    => 'error',
                                              -type    => 'ok');
            &hide_pbar_2;
            return;
         };

         unlink ${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'/'.$filename.'.tbl';
         open FEATURE, ">>".${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'/'.$filename.'.tbl' or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot create Sequin Table file",
                                              -icon    => 'error',
                                              -type    => 'ok');
            &hide_pbar_2;
            return;
         };

         #create agp file if selected
         if ($args{create_scaffold_info} == 1) {
            unlink ${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'/'.$filename.'.agp';
            open AGP, ">>".${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'/'.$filename.'.agp' or do {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Cannot create Scaffold info file in directory $args{project_name} for file $filename.",
                                                 -icon    => 'error',
                                                 -type    => 'ok');
               &hide_pbar_2;
               return;
            };
         }
      }

      #get current roman number, else empty
      if ($args{cluster} == 1 && $#short_files > 0) {
         $roman = roman($file_counter).'_';
      } else {
         $roman = '';
      }

      #determine format
      open READ, "<$all_files{$file}";
      $_ = <READ>;
      unless ($_ =~ /^locus/i) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Unrecognised file format for $all_files{$file}. Skipping file",
                                           -type    => 'OK',
                                           -icon    => 'error');
         close READ;
         next;
      }
      close READ;

      #parse gb file
      &update_pbar_2(title        => "Processing Genbank files",
                     label        => "Reading entries for: $filename",
                     progress     => 1,
                    );

      my ($header, $source, $core, $DNAseq, $start_ref,
          $feature_list_ref, $genbank_ref, $counter) = &gb_parser(main_window   => $args{main_window},
                                                                  progress_bar  => $args{progress_bar},
                                                                  ini_ref       => $args{ini_ref},
                                                                  auto_ini_ref  => $args{auto_ini_ref},
                                                                  gb_file       => $all_files{$file},
                                                                 );

      #define fasta header
      $fasta_header = &fasta(header => $header,
                             source => $source,
                             target => 'fasta'
                            );
      $fasta_header = '>'.$fasta_header;

      #print FASTA file
      print FASTA $fasta_header."\n".$DNAseq."\n";

      #define Table header
      $table_header = &fasta(header => $header,
                             source => $source,
                             target => 'table'
                            );
      print FEATURE '>Feature '.$table_header."\n";

      #sort feature_list_ref via start position and genbank_key_order
      $counter = 1;
      foreach my $key (@gb_order) {
         map { $_ =~ s/$key/$counter\-$key/ } @{$feature_list_ref};
         $counter++;
      }
      #add generic number to those entries not listed in key order
      foreach my $entry (@{$feature_list_ref}) {
         $entry =~ s/^(\d+_\d+_)(\D+\_.+)/$1$counter\-$2/;
      }
      my @sorted =
                  map  $_->[0] =>
                  sort { $a->[1] <=> $b->[1] ||
                         $a->[2] <=> $b->[2] }
                  map  [ $_, m/^(\d+)\_/, m/^\d+_\d+_(\d+)\-.*?_\d+/ ]
                  => @{$feature_list_ref};

      #get only features that are selected
      @sorted = grep{(/$selected\_/)} @sorted;

      #test for duplicated features: same key, same start-stop which will cause problems with Sequin
      %seen = ();
      foreach my $entry (@sorted) {
         'reset' =~ m/reset/;
         $entry =~ m/^(\d+_\d+_[^\_]+)_/;
         my $tag = $1;
         if (defined $seen{$tag}) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Likely duplicated feature detected, same bondaries with existing locus tag. Check Genbank source file.\nFeature: $entry\n",
                                              -type    => 'OK',
                                              -icon    => 'error');
            &hide_pbar_2;
            return;
         } else {
            $seen{$tag}++;
         }
      }
      %seen = ();

      #get maximum number of entries
      $max_number = length ($#sorted);

      #reset locus tag ID if selected
      if ($args{reset_locusid} == 1) {
         $id_counter = 1;
      }

      #create lookup hash for selected features; to be used in CDS annotation for PFAM, TIGRfam and COG hits
      my %sel_feat = ();
      foreach my $entry (@sorted) {$sel_feat{$entry}++};

      $max_count = $#sorted + 1;
      my $update_counter = 1;
      foreach my $entry (@sorted) {
         my ($annotation, $locus_tag, $gene, $CDS, $rRNA, $tRNA, $start, $stop, $orientation);

         #ignore PFAM, TIGRFAM and COG matches; to be included in CDS annotation if selected
         next if ($entry =~ /(COG|PFAM|TIGR)~match_/i);

         #create access key for maintaining previous annotation
         my $short = $entry;
         $short =~ s/_\d+$//;

         #remove genbank order key again
         $entry =~ s/_\d+\-(.+)$/\_$1/;

         #retrieve key
         'reset' =~ m/reset/;
         $entry =~ m/\d+_\d+_(.+?)_\d+/;
         $table_feature = $1;

         #retrieve boundaries
         ($start, $stop, $orientation) = &boundaries(feature => ${$genbank_ref}{$entry});

         #define locus tag for gene/CDS
         if ($table_feature =~ /(gene|CDS|5'UTR|3'UTR|intron|exon|RNA)$/) {
            unless (defined $locus_counter{$start.'_'.$stop}) {
               $locus_counter{$start.'_'.$stop} = $id_counter;
               $id_counter++;
            }

         }

         #determine locus tag
         #gene, mRNA, CDS, 5'UTR, 3'UTR, intron, exon, tRNA, rRNA and ncRNA require
         #their own locus_tags
         if ($table_feature =~ m/(gene|CDS|5'UTR|3'UTR|intron|exon|RNA)$/) {
            ($locus_tag   ) = &locus_tag (genbank_ref    => $genbank_ref,
                                          entry          => $entry,
                                          table_feature  => $table_feature,
                                          counter        => $id_counter,
                                          max_number     => $max_number,
                                          roman          => uc ($roman),
                                          db_name        => $args{db_name},
                                          locus_tag      => $args{locus_tag},
                                          reuse_locustag => $args{reuse_locustag},
                                          locus_counter  => $locus_counter{$start.'_'.$stop}
                                          );
         } else {
            $locus_tag = '';
         }

         #change feature back to original spelling
         $table_feature =~ s/~/_/;

         #update
         if (($update_counter % 100) == 0) {
            &update_pbar_2(title        => "Processing Genbank files",
                           label        => "Processing features",
                           progress     => ($update_counter / $max_count) * 100,
                          );
         }
         $update_counter++;

         #define gene if selected
         if ($table_feature eq 'gene') {
            my $annotation = &get_gene_annotation(genbank_ref    => $genbank_ref,
                                                  entry          => $entry
                                                  );
            &print_feature;
            print FEATURE "\t\t\t".$locus_tag;
            if ($annotation =~ /\w+/ && $args{no_gene_feature} == 0) { #only print gene qualifier if selected
               print FEATURE $annotation;
            }
            next;
         }

         #define CDS if selected
         if ($table_feature eq 'CDS') {
            my $annotation = &get_CDS_annotation(genbank_ref         => $genbank_ref,
                                                 entry               => $entry,
                                                 notes               => $args{notes},
                                                 force_product       => $args{force_product},
                                                 qualifiers_ref      => $args{qualifiers_ref},
                                                 sel_features_ref    => $args{selected_feature},
                                                 sorted_features_ref => \@sorted,
                                                 start               => $start,
                                                 stop                => $stop,
                                                 orientation         => $orientation
                                                );

            &print_feature;
            print FEATURE "\t\t\t$locus_tag";
            if ($annotation =~ /\w+/) {
               print FEATURE $annotation;
            }
            next;
         }

         #define rRNA if selected
         if ($table_feature eq 'rRNA') {
            #create additional gene feature
            $table_feature = 'gene';
            my $annotation    = &get_gene_annotation(genbank_ref    => $genbank_ref,
                                                     entry          => $entry
                                                    );
            #preserve boundaries for gene feature
            my @save_start = @start;
            my @save_stop = @stop;
            &print_feature;
            print FEATURE "\t\t\t".$locus_tag;
            if ($annotation =~ /\w+/) {
               print FEATURE $annotation;
            }

            #create additional rRNA feature
            $table_feature = 'rRNA';
            #restore boundaries
            @start = @save_start;
            @stop  = @save_stop;

            $annotation = &get_rRNA_annotation(genbank_ref    => $genbank_ref,
                                               entry          => $entry
                                              );
            &print_feature;
            if ($annotation =~ /\w+/) {
               print FEATURE $annotation;
            }
            next;
         }

         #define tRNA if selected
         if ($table_feature eq 'tRNA') {
            #create additional gene feature
            $table_feature = 'gene';

            my $annotation = &get_tRNA_annotation(genbank_ref    => $genbank_ref,
                                                  entry          => $entry
                                                 );
            #preserve boundaries for gene feature
            my @save_start = @start;
            my @save_stop = @stop;

            &print_feature;
            print FEATURE "\t\t\t$locus_tag";
            if ($annotation =~ /\w+/) {
               print FEATURE $annotation;
            }

            #create additional gene feature
            $table_feature = 'tRNA';
            #restore boundaries
            @start = @save_start;
            @stop  = @save_stop;

            #$annotation    = &get_gene_annotation(genbank_ref    => $genbank_ref,
            #                                      entry          => $entry
            #                                     );

            &print_feature;
            if ($annotation =~ /\w+/) {
               print FEATURE $annotation;
            }
            next;
         }

         #define ncRNA if selected
         #ncRNA is used as a feature for non-coding RNA (ncRNA)
         if ($table_feature eq 'ncRNA' || $table_feature eq 'misc~RNA') {
            #create additional gene feature
            $table_feature = 'gene';

            #my $annotation = &get_miscRNA_annotation(genbank_ref    => $genbank_ref,
            #                                         entry          => $entry
            #                                        );
            #preserve boundaries for gene feature
            my @save_start = @start;
            my @save_stop = @stop;

            &print_feature;
            print FEATURE "\t\t\t$locus_tag";
            #if ($annotation =~ /\w+/) {
            #   print FEATURE $annotation;
            #}

            #create additional gene feature
            $table_feature = 'ncRNA';
            #restore boundaries
            @start = @save_start;
            @stop  = @save_stop;
            $annotation    = &get_miscRNA_annotation(genbank_ref    => $genbank_ref,
                                                     entry          => $entry
                                                    );
            &print_feature;
            if ($annotation =~ /\w+/) {
               print FEATURE $annotation;
            }
            next;
         }

         #define remaining standard features
         {
            my $annotation = '';
            $annotation = &get_standard_annotation(genbank_ref    => $genbank_ref,
                                                   entry          => $entry,
                                                  );
            &print_feature;
            if ($annotation =~ /\w+/) {
               print FEATURE $annotation;
            }
         }

      }

      if ($args{file_individual} == 1) {
         close FASTA;
         close FEATURE;
      }

      undef %locus_counter;
      $file_counter++;
   }

   if ($args{file_individual} == 0) {
      close FASTA;
      close FEATURE;
   }
   &hide_pbar_2;

   ${$args{main_window}}->messageBox(-title   => 'Info',
                                     -message => "Successfully created Genbank submission project",
                                     -icon    => 'info',
                                     -type    => 'ok');

   return;
}

sub compile_WGS {
   my %args = @_;
   my ($file_counter, $counter, $id_counter, @gb_order, $selected, $max_count, $local_path,
       $filename, @non_standard, %seen, @standard, $WGS_contigs);
   my (%scaffold, %unique_gene_anno);
   my $status;

   #clear all variables
   $file_counter     = '';
   $counter          = '';
   $id_counter       = '';
   $selected         = '';
   $max_count        = '';
   $local_path       = '';
   $filename         = '';
   @gb_order         = ();
   @non_standard     = ();
   @standard         = ();
   %seen             = ();
   %unique_gene_anno = ();
   undef $WGS_contigs;

   #create directory if not already there
   unless (-d ${$args{ini_ref}}{Sequin_submission}) {
      mkdir (${$args{ini_ref}}{Sequin_submission}, 0777) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create Sequin submission directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return;
      };
   }

   #create combined files
   #this is for multiple input files in draft genome.
   #each contig will then be treated within the context of all files
   if ($args{file_individual} == 0) {
      unlink ${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'.fsa';
      open FASTA, ">>".${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'.fsa' or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create Fasta Sequin file",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return;
      };

      unlink ${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'.tbl';
      open FEATURE, ">>".${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'.tbl' or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create Sequin Table file",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return;
      };
   }

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Processing Genbank files",
                   label        => 'Processing features'
                  );
   &show_pbar_2;

   #get defined features into search pattern
   $selected = '(';
   foreach my $key (keys (%{$args{selected_feature}})) {
      if ($args{selected_feature}{$key}) {
         $selected .= $key.'|';
         push (@standard, $key);
         #add _ to ~ replacement just in case
         if ($key =~ m/_/) {
            my $mod_key = $key;
            $mod_key =~ s/_/~/;
            $selected .= $mod_key.'|';
            push (@standard, $mod_key);
         }
      };
   }
   $selected =~ s/\|$//;
   $selected .= ')';

   #get standard descriptor of rest of selected features
   @non_standard = qw (gene mRNA CDS 5'UTR 3'UTR intron exon tRNA rRNA ncRNA misc_RNA misc~RNA);
   foreach my $entry (@standard, @non_standard) {$seen{$entry}++};
   @standard     = ();
   @non_standard = ();
   foreach my $entry (keys %seen) {
      if ($seen{$entry} == 1) {
         push (@standard, $entry);
      }
   }
   %seen = ();

   #split Genbank key order into array
   @gb_order = split /;/,${$args{auto_ini_ref}}{genbank_key_order};

   #iterate through each selected file
   $file_counter = 1;
   $id_counter   = 1;
   foreach my $file (@short_files) {
      my ($roman, %locus_counter, $fasta_header, $table_header, $max_number);
      #reset variables
      $roman         = '';
      $fasta_header  = '';
      $table_header  = '';
      $max_number    = '';
      %locus_counter = ();

      unless (exists $all_files{$file}) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "File $all_files{$file} cannot be accessed",
                                           -type    => 'OK',
                                           -icon    => 'error');
         &hide_pbar_2;
         return (0);
      }
      ${$args{progress_bar}} ->configure(-label => "Processing $file");
      ${$args{main_window}}  ->update;


      #pull filepath apart
      'reset' =~ m/reset/;
      $all_files{$file} =~ m/^(.+)\/([^\/]+)$/;
      $local_path = $1;
      $filename = $2;

      #create individual projects if selected
      if ($args{file_individual} == 1) {
         #create directory if not already there
         unless (-d ${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}) {
            mkdir (${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}, 0777) or do {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Cannot create individual Sequin project directory",
                                                 -icon    => 'error',
                                                 -type    => 'ok');
               &hide_pbar_2;
               return;
            };
         }

         #create combined files
         unlink ${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'/'.$filename.'.fsa';
         open FASTA, ">>".${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'/'.$filename.'.fsa' or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot create Fasta Sequin file",
                                              -icon    => 'error',
                                              -type    => 'ok');
            &hide_pbar_2;
            return;
         };

         unlink ${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'/'.$filename.'.tbl';
         open FEATURE, ">>".${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'/'.$filename.'.tbl' or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot create Sequin Table file",
                                              -icon    => 'error',
                                              -type    => 'ok');
            &hide_pbar_2;
            return;
         };

         #create AGP file if selected
         if ($args{create_scaffold_info} == 1) {
            unlink ${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'/'.$filename.'.agp';
            open AGP, ">>".${$args{ini_ref}}{Sequin_submission}.'/'.$args{project_name}.'/'.$filename.'.agp' or do {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Cannot create Scaffold info file in directory $args{project_name} for file $filename.",
                                                 -icon    => 'error',
                                                 -type    => 'ok');
               &hide_pbar_2;
               return;
            };
         }
      }

      #get current roman number, else empty
      if ($args{cluster} == 1 && $#short_files > 0) {
         $roman = roman($file_counter);
      } else {
         $roman = '';
      }

      #determine format
      open READ, "<$all_files{$file}";
      $_ = <READ>;
      unless ($_ =~ /^locus/i) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Unrecognised file format for $all_files{$file}. Skipping file",
                                           -type    => 'OK',
                                           -icon    => 'error');
         close READ;
         next;
      }
      close READ;

      #parse gb file
      &update_pbar_2(title        => "Processing Genbank and AGP mapping files",
                     label        => "Reading entries for: $filename",
                     progress     => 1,
                    );

      #read selected AGP mapping file is selected
      if (defined $args{agp_map_file} && -e $args{agp_map_file}) {
         my $scf_counter = 1;
         my $position    = 1;
         open READ, "<$args{agp_map_file}";
         while (<READ>) {
            next if ($_ =~ m/^\#/);
            if ($_ =~ m/^(>|\s)/) {
               $scf_counter++;
               $position = 1;
               next;
            }
            chomp $_;
            $_ =~ s/\s+//g;
            $scaffold{$_}   = {(scf_counter => $scf_counter,
                                start       => '0',
                                position    => $position
                              )};
            $position++;
         }
         close READ;
      } else {
         %scaffold = ();
      }

      my ($header, $source, $core, $DNAseq, $start_ref,
          $feature_list_ref, $genbank_ref, $counter) = &gb_parser(main_window   => $args{main_window},
                                                                  progress_bar  => $args{progress_bar},
                                                                  ini_ref       => $args{ini_ref},
                                                                  auto_ini_ref  => $args{auto_ini_ref},
                                                                  gb_file       => $all_files{$file},
                                                                 );

      #get fasta header
      $fasta_header = &fasta(header => $header,
                             source => $source,
                             target => 'fasta'
                            );

      #define fasta entries for each contig in each file
      $status = &fasta_WGS(main_window  => $args{main_window},
                           header       => $header,
                           source       => $source,
                           core         => $core,
                           sequence     => $DNAseq,
                           fasta_header => $fasta_header,
                           target       => 'fasta',
                           wgs_hash_ref => \$WGS_contigs,
                           scaffold_ref => \%scaffold
                          );
      if ($status == 0) {
         return (0);
      }

      #print FASTA file for each entry; no features yet, therefore feature_counter = 1
      foreach my $counter_contig (sort {$a<=>$b} keys %{ $WGS_contigs->{$fasta_header} }) {
         print FASTA '>'.$WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'ctg_name'}.
                     ' '.$WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'fasta_header'}.' [gcode=11] '."\n".
                         $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'ctg_sequence'}."\n";
      }

      #define Table header
      $table_header = &fasta(header => $header,
                             source => $source,
                             target => 'table'
                            );

      #sort feature_list_ref via start position and genbank_key_order
      $counter = 1;
      foreach my $key (@gb_order) {
         map { $_ =~ s/$key/$counter\-$key/ } @{$feature_list_ref};
         $counter++;
      }
      #add generic number to those entries not listed in key order
      foreach my $entry (@{$feature_list_ref}) {
         $entry =~ s/^(\d+_\d+_)(\D+\_.+)/$1$counter\-$2/;
      }
      my @sorted =
                  map  $_->[0] =>
                  sort { $a->[1] <=> $b->[1] ||
                         $a->[2] <=> $b->[2] }
                  map  [ $_, m/^(\d+)\_/, m/^\d+_\d+_(\d+)\-.*?_\d+/ ]
                  => @{$feature_list_ref};

      #get only features that are selected
      @sorted = grep{(/$selected\_/)} @sorted;

      #test for duplicated features: same key, same start-stop which will cause problems with Sequin
      %seen = ();
      foreach my $entry (@sorted) {
         'reset' =~ m/reset/;
         $entry =~ m/^(\d+_\d+_[^\_]+)_/;
         my $tag = $1;
         if (defined $seen{$tag}) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Likely duplicated feature detected, same bondaries with existing locus tag. Check Genbank source file $all_files{$file}.\nFeature: $entry\n",
                                              -type    => 'OK',
                                              -icon    => 'error');
            return;
         } else {
            $seen{$tag}++;
         }
      }
      %seen = ();

      #get maximum number of entries
      $max_number = length ($#sorted);

      #add each feature to WGS hash
      my %feature_counter = ();
      my $check           = 0;
      foreach my $entry (@sorted) {
         #get start and stop position
         $entry =~ m/^(\d+)_(\d+)_/;
         my ($feature_start, $feature_stop) = ($1, $2);

         #now test which contig it belongs to and add feature
         foreach my $counter_contig (sort {$a<=>$b} keys %{ $WGS_contigs->{$fasta_header} }) {
            if ($feature_start >= $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'ctg_left_bd'} - 22 &&
                $feature_stop  <= $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'ctg_right_bd'} + 22 ) {
               $feature_counter{$counter_contig}++;
               $WGS_contigs->{$fasta_header}->{$counter_contig}->{$feature_counter{$counter_contig}}->{'feature'}        = $entry;
               $WGS_contigs->{$fasta_header}->{$counter_contig}->{$feature_counter{$counter_contig}}->{'ctg_name'}       = $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'ctg_name'};
               $WGS_contigs->{$fasta_header}->{$counter_contig}->{$feature_counter{$counter_contig}}->{'ctg_left_bd'}    = $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'ctg_left_bd'};
               $WGS_contigs->{$fasta_header}->{$counter_contig}->{$feature_counter{$counter_contig}}->{'ctg_right_bd'}   = $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'ctg_right_bd'};
               $WGS_contigs->{$fasta_header}->{$counter_contig}->{$feature_counter{$counter_contig}}->{'scf_counter'}    = $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'scf_counter'};
               $WGS_contigs->{$fasta_header}->{$counter_contig}->{$feature_counter{$counter_contig}}->{'scf_start'}      = $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'scf_start'};
               $WGS_contigs->{$fasta_header}->{$counter_contig}->{$feature_counter{$counter_contig}}->{'ctg_pos_in_scf'} = $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'ctg_pos_in_scf'};
               $check = 1;
               last;
            }
         }
         if ($check == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Info',
                                              -message => "Error assigning feature $entry to any contig in $fasta_header",
                                              -icon    => 'info',
                                              -type    => 'ok');
            &hide_pbar_2;
            return (0);
         }
      }

      #reset locus tag ID if selected
      if ($args{reset_locusid} == 1) {
         $id_counter = 1;
      }

      #create lookup hash for selected features; to be used in CDS annotation for PFAM, TIGRfam and COG hits
      my %sel_feat = ();
      foreach my $entry (@sorted) {$sel_feat{$entry}++};

      $max_count           = $#sorted + 1;
      my $update_counter   = 1;
      #this corrects for the difference between real GAMOLA spacer sequence length and default unknown gap size of 100 nt by NCBI
      # will be multiplied by contig number - 1 for each consecutive contig
      #if AGP mapping file is selected, then this will also be taken into account, resetting gene model coordinate for each new scaffold

      foreach my $counter_contig (sort {$a<=>$b} keys %{ $WGS_contigs->{$fasta_header} }) {
         print FEATURE '>Features '.$WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'ctg_name'}."\n";
         foreach my $counter_feature (sort {$a<=>$b} keys %{ $WGS_contigs->{$fasta_header}->{$counter_contig} }) {
            my ($annotation, $locus_tag, $gene, $CDS, $rRNA, $tRNA, $start, $stop, $orientation);

            #ignore PFAM, TIGRFAM and COG matches; to be included in CDS annotation if selected
            next if ($WGS_contigs->{$fasta_header}->{$counter_contig}->{$counter_feature}->{'feature'} =~ /(COG|PFAM|TIGR)~match_/i);

            #create access key for maintaining previous annotation
            my $short = $WGS_contigs->{$fasta_header}->{$counter_contig}->{$counter_feature}->{'feature'};
            $short =~ s/_\d+$//;

            #remove genbank order key again
            $WGS_contigs->{$fasta_header}->{$counter_contig}->{$counter_feature}->{'feature'} =~ s/_\d+\-(.+)$/\_$1/;

            #retrieve key
            'reset' =~ m/reset/;
            $WGS_contigs->{$fasta_header}->{$counter_contig}->{$counter_feature}->{'feature'} =~ m/\d+_\d+_(.+?)_\d+/;
            $table_feature = $1;

            #retrieve boundaries

            ($start, $stop, $orientation) = &boundaries(feature => ${$genbank_ref}{$WGS_contigs->{$fasta_header}->{$counter_contig}->{$counter_feature}->{'feature'}});
            #define locus tag for gene/CDS
            if ($table_feature =~ /(gene|CDS|5'UTR|3'UTR|intron|exon|RNA)$/) {
               unless (defined $locus_counter{$start.'_'.$stop}) {
                  $locus_counter{$start.'_'.$stop} = $id_counter;
                  $id_counter++;
               }
            }

            #determine locus tag
            #gene, mRNA, CDS, 5'UTR, 3'UTR, intron, exon, tRNA, rRNA and ncRNA require
            #their own locus_tags
            if ($table_feature =~ m/(gene|CDS|5'UTR|3'UTR|intron|exon|RNA)$/) {
               ($locus_tag   ) = &locus_tag (genbank_ref    => $genbank_ref,
                                             entry          => $WGS_contigs->{$fasta_header}->{$counter_contig}->{$counter_feature}->{'feature'},
                                             table_feature  => $table_feature,
                                             counter        => $id_counter,
                                             max_number     => $max_number,
                                             roman          => uc ($roman),
                                             db_name        => $args{db_name},
                                             locus_tag      => $args{locus_tag},
                                             reuse_locustag => $args{reuse_locustag},
                                             locus_counter  => $locus_counter{$start.'_'.$stop}
                                             );
            } else {
               $locus_tag = '';
            }

            #change feature back to original spelling
            $table_feature =~ s/~/_/;

            #update
            if (($update_counter % 100) == 0) {
               &update_pbar_2(title        => "Processing Genbank files",
                              label        => "Processing features",
                              progress     => ($update_counter / $max_count) * 100,
                             );
            }
            $update_counter++;

            #define gene if selected
            if ($table_feature eq 'gene') {
               my $annotation = &get_gene_annotation(genbank_ref          => $genbank_ref,
                                                     unique_gene_anno_ref => \%unique_gene_anno,
                                                     entry                => $WGS_contigs->{$fasta_header}->{$counter_contig}->{$counter_feature}->{'feature'}
                                                     );
               &print_feature_WGS (ctg_counter      => $counter_contig,
                                   gap              => $args{gap},
                                   fasta_header     => $fasta_header,
                                   scaffold_ref     => \%scaffold,
                                   counter_feature  => $counter_feature,
                                   scf_counter      => $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'scf_counter'},
                                   wgs_hash_ref     => \$WGS_contigs,
                                  );
               print FEATURE "\t\t\t".$locus_tag;
               if ($annotation =~ /\w+/ && $args{no_gene_feature} == 0) { #print only if selected
                  print FEATURE $annotation;
               }
               next;
            }

            #define CDS if selected
            if ($table_feature eq 'CDS') {
               my $annotation = &get_CDS_annotation(genbank_ref         => $genbank_ref,
                                                    entry               => $WGS_contigs->{$fasta_header}->{$counter_contig}->{$counter_feature}->{'feature'},
                                                    notes               => $args{notes},
                                                    force_product       => $args{force_product},
                                                    qualifiers_ref      => $args{qualifiers_ref},
                                                    sel_features_ref    => $args{selected_feature},
                                                    sorted_features_ref => \@sorted,
                                                    start               => $start,
                                                    stop                => $stop,
                                                    orientation         => $orientation
                                                   );

               &print_feature_WGS (ctg_counter     => $counter_contig,
                                   gap             => $args{gap},
                                   fasta_header    => $fasta_header,
                                   scaffold_ref    => \%scaffold,
                                   counter_feature => $counter_feature,
                                   scf_counter     => $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'scf_counter'},
                                   wgs_hash_ref    => \$WGS_contigs,
                                  );
               print FEATURE "\t\t\t$locus_tag";
               if ($annotation =~ /\w+/) {
                  print FEATURE $annotation;
               }
               next;
            }

            #define rRNA if selected
            if ($table_feature eq 'rRNA') {
               #create additional gene feature
               $table_feature = 'gene';
               my $annotation    = &get_rRNA_annotation(genbank_ref    => $genbank_ref,
                                                        entry          => $WGS_contigs->{$fasta_header}->{$counter_contig}->{$counter_feature}->{'feature'}
                                                       );
               #preserve boundaries for gene feature
               my @save_start = @start;
               my @save_stop = @stop;
               &print_feature_WGS (ctg_counter     => $counter_contig,
                                   gap             => $args{gap},
                                   fasta_header    => $fasta_header,
                                   scaffold_ref    => \%scaffold,
                                   counter_feature => $counter_feature,
                                   scf_counter     => $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'scf_counter'},
                                   wgs_hash_ref    => \$WGS_contigs,
                                  );
               print FEATURE "\t\t\t".$locus_tag;
               if ($annotation =~ /\w+/) {
                  print FEATURE $annotation;
               }

               #create additional gene feature
               $table_feature = 'rRNA';
               #restore boundaries
               @start = @save_start;
               @stop  = @save_stop;

               $annotation = &get_gene_annotation(genbank_ref          => $genbank_ref,
                                                  unique_gene_anno_ref => \%unique_gene_anno,
                                                  entry                => $WGS_contigs->{$fasta_header}->{$counter_contig}->{$counter_feature}->{'feature'}
                                                 );
               &print_feature_WGS (ctg_counter     => $counter_contig,
                                   gap             => $args{gap},
                                   fasta_header    => $fasta_header,
                                   scaffold_ref    => \%scaffold,
                                   counter_feature => $counter_feature,
                                   scf_counter     => $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'scf_counter'},
                                   wgs_hash_ref    => \$WGS_contigs,
                                  );
               if ($annotation =~ /\w+/) {
                  print FEATURE $annotation;
               }
               next;
            }

            #define tRNA if selected
            if ($table_feature eq 'tRNA') {
               #create additional gene feature
               $table_feature = 'gene';

               my $annotation = &get_tRNA_annotation(genbank_ref    => $genbank_ref,
                                                     entry          => $WGS_contigs->{$fasta_header}->{$counter_contig}->{$counter_feature}->{'feature'}
                                                    );
               #preserve boundaries for gene feature
               my @save_start = @start;
               my @save_stop = @stop;

               &print_feature_WGS (ctg_counter     => $counter_contig,
                                   gap             => $args{gap},
                                   fasta_header    => $fasta_header,
                                   scaffold_ref    => \%scaffold,
                                   counter_feature => $counter_feature,
                                   scf_counter     => $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'scf_counter'},
                                   wgs_hash_ref    => \$WGS_contigs,
                                  );
               print FEATURE "\t\t\t$locus_tag";
               if ($annotation =~ /\w+/) {
                  print FEATURE $annotation;
               }

               #create additional gene feature
               $table_feature = 'tRNA';
               #restore boundaries
               @start = @save_start;
               @stop  = @save_stop;

               $annotation    = &get_gene_annotation(genbank_ref    => $genbank_ref,
                                                     entry          => $WGS_contigs->{$fasta_header}->{$counter_contig}->{$counter_feature}->{'feature'}
                                                    );
               &print_feature_WGS (ctg_counter     => $counter_contig,
                                   gap             => $args{gap},
                                   fasta_header    => $fasta_header,
                                   scaffold_ref    => \%scaffold,
                                   counter_feature => $counter_feature,
                                   scf_counter     => $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'scf_counter'},
                                   wgs_hash_ref    => \$WGS_contigs,
                                  );
               if ($annotation =~ /\w+/) {
                  print FEATURE $annotation;
               }
               next;
            }

            #define ncRNA if selected
            #ncRNA is used as a feature for non-coding RNA (ncRNA)
            if ($table_feature eq 'ncRNA' || $table_feature eq 'misc~RNA') {
               #create additional gene feature
               $table_feature = 'gene';

               #my $annotation = &get_miscRNA_annotation(genbank_ref    => $genbank_ref,
               #                                         entry          => $entry
               #                                        );
               #preserve boundaries for gene feature
               my @save_start = @start;
               my @save_stop = @stop;

               &print_feature_WGS (ctg_counter     => $counter_contig,
                                   gap             => $args{gap},
                                   fasta_header    => $fasta_header,
                                   scaffold_ref    => \%scaffold,
                                   counter_feature => $counter_feature,
                                   scf_counter     => $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'scf_counter'},
                                   wgs_hash_ref    => \$WGS_contigs,
                                  );
               print FEATURE "\t\t\t$locus_tag";
               #if ($annotation =~ /\w+/) {
               #   print FEATURE $annotation;
               #}

               #create additional gene feature
               $table_feature = 'ncRNA';
               #restore boundaries
               @start = @save_start;
               @stop  = @save_stop;
               $annotation    = &get_miscRNA_annotation(genbank_ref    => $genbank_ref,
                                                        entry          => $WGS_contigs->{$fasta_header}->{$counter_contig}->{$counter_feature}->{'feature'}
                                                       );
               &print_feature_WGS (ctg_counter     => $counter_contig,
                                   gap             => $args{gap},
                                   fasta_header    => $fasta_header,
                                   scaffold_ref    => \%scaffold,
                                   counter_feature => $counter_feature,
                                   scf_counter     => $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'scf_counter'},
                                   wgs_hash_ref    => \$WGS_contigs,
                                  );
               if ($annotation =~ /\w+/) {
                  print FEATURE $annotation;
               }
               next;
            }

            #define remaining standard features
            {
               my $annotation = '';
               $annotation = &get_standard_annotation(genbank_ref    => $genbank_ref,
                                                      entry          => $WGS_contigs->{$fasta_header}->{$counter_contig}->{$counter_feature}->{'feature'},
                                                     );
               &print_feature_WGS (ctg_counter     => $counter_contig,
                                   gap             => $args{gap},
                                   fasta_header    => $fasta_header,
                                   scaffold_ref    => \%scaffold,
                                   counter_feature => $counter_feature,
                                   scf_counter     => $WGS_contigs->{$fasta_header}->{$counter_contig}->{'1'}->{'scf_counter'},
                                   wgs_hash_ref    => \$WGS_contigs,
                                  );
               if ($annotation =~ /\w+/) {
                  print FEATURE $annotation;
               }
            }
         }
      }

      #create scaffold_info file if selected
      if ($args{create_scaffold_info} == 1) {
         &create_scaffold_info(main_window     => $args{main_window},
                               header          => $header,
                               source          => $source,
                               core            => $core,
                               sequence        => $DNAseq,
                               fasta_header    => $fasta_header,
                               locus_tag       => $args{locus_tag},
                               target          => 'fasta',
                               filename        => $filename,
                               project_name    => $args{project_name},
                               file_individual => $args{file_individual},
                               wgs_hash_ref    => \$WGS_contigs
                              );
      }

      if ($args{file_individual} == 1) {
         close FASTA;
         close FEATURE;
      }

      undef %locus_counter;
      $file_counter++;
   }

   if ($args{file_individual} == 0) {
      close FASTA;
      close FEATURE;
   }
   &hide_pbar_2;

   ${$args{main_window}}->messageBox(-title   => 'Info',
                                     -message => "Successfully created Genbank submission project",
                                     -icon    => 'info',
                                     -type    => 'ok');

   return;
}

sub fasta {
   my %args = (target => 'fasta',
               @_
              );
   my ($accession, $organism, $strain, @mods);
   my @modifiers = qw (acronym 	country 	lat-lon 	specific-host
                       anamorph 	cultivar 	map 	specimen-voucher
                       authority 	dev-stage 	metagenomic 	strain
                       biotype 	ecotype 	note 	sub-species
                       biovar 	endogenous-virus-name 	organism 	subclone
                       breed 	forma 	pathovar 	subgroup
                       cell-line 	forma-specialis 	plasmid-name 	substrain
                       cell-type 	fwd-PCR-primer-seq 	plastid-name 	subtype
                       chemovar 	genotype 	pop-variant 	synonym
                       chromosome 	group 	rev-PCR-primer-seq 	teleomorph
                       clone 	haplotype 	segment 	tissue-lib
                       clone-lib 	identified-by 	serogroup 	tissue-type
                       collected-by 	isolate 	serotype 	type
                       collection-date 	isolation-source 	serovar 	variety
                       common 	lab-host 	sex
                      );
   'reset' =~ m/reset/;
   $args{header} =~ m/^locus\s+(\S+)\s+/i;
   $accession = $1;
   unless (defined $accession) {
      $args{header} =~ m/\n\r?accession\s+(\S+)\s*\n\r?/si;
      $accession = $1;
   }
   unless (defined $accession) {
      $args{source} =~ m/\/strain\=\"([^\"]*?)\"/si;
      $accession = $1;
   }
   unless (defined $accession) {
      $accession = 'anonymous';
   }
   $accession =~ s/\s+/ /gs;
   $accession .= ' ';

   #return accession for Sequin table
   if ($args{target} eq 'table') {return $accession};

   foreach my $modifier (@modifiers) {
      'reset' =~ m/reset/;
      $args{source} =~ m/\/$modifier\=\"([^\"]*?)\"/si;
      my $hit = $1;
      if (defined $hit) {
         $hit = '['.$modifier.'='.$hit.']';
      } else {
         $hit = '';
      }
      push (@mods, $hit);
   }

   #merge all hits
   foreach my $hit (@mods) {
      $accession .= $hit.' ';
   }
   $accession =~ s/\s+/ /g;
   return $accession;
}

sub fasta_WGS {
   my %args = (target => 'fasta',
               @_
              );
   my ($accession, $organism, $strain, $counter_contig, $counter_feature);
   my @contigs = ();

   #remove trailing info from fasta header
   my $fsa_head = $args{fasta_header};
   $fsa_head =~ s/^[^\[]+//;

   #get contig entries across the genbank file
   @contigs = ($args{core} =~ m/\s+(fasta_record\s+\d+\.\.\d+.+?)\/colour/gsi);

   #iterate through contigs and get boundaries, name and sequence
   $counter_contig  = 1;
   $counter_feature = 1;
   foreach my $entry (@contigs) {
      my ($orientation, $left_bd, $right_bd, $ctg_name, $seq);
      'reset' =~ m/reset/;
      $entry =~ m/fasta_record\s+(complement)?(\d+)\.\.(\d+).+?\/product\=\"Contig\:\s+([^\,]+?)\,/s;
      ($orientation, $left_bd, $right_bd, $ctg_name) = ($1, $2, $3, $4);
      #different format? try capture all until first whitespace
      unless (defined $ctg_name) {
         'reset' =~ m/reset/;
         $entry =~ m/fasta_record\s+(complement)?(\d+)\.\.(\d+).+?\/product\=\"([^\s\"]+)[\s\"]/s;
         ($orientation, $left_bd, $right_bd, $ctg_name) = ($1, $2, $3, $4);
      }
      unless (defined $ctg_name) {
         print "\nError parsing entry $entry in fasta_WGS entry $args{header}.\n";
         return(0);
      };
      if (defined $orientation) {
         $orientation = '-';
      } else {
         $orientation = '+';
      }
      $seq = substr($args{sequence}, ($left_bd - 1), ($right_bd - $left_bd + 1));

      #test if contig is > 199 nt
      if (length ($seq) < 200 ) {
         ${$args{main_window}}->messageBox(-title   => 'Info',
                                           -message => "Error: contig $ctg_name is smaller than 200nt in $args{fasta_header}.\nRemove this contig and rebuild assembly for NCBI WGS submission",
                                           -icon    => 'info',
                                           -type    => 'ok');
         return (0);
      }

      if (${$args{scaffold_ref}}{$ctg_name}->{'position'} == 1) {
         ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$counter_contig}->{'1'} = {( ctg_name        => $ctg_name,
                                                                                      ctg_orientation => $orientation,
                                                                                      ctg_left_bd     => $left_bd,
                                                                                      ctg_right_bd    => $right_bd,
                                                                                      ctg_sequence    => $seq,
                                                                                      fasta_header    => $fsa_head,
                                                                                      scf_counter     => ${$args{scaffold_ref}}{$ctg_name}->{'scf_counter'},
                                                                                      scf_start       => $left_bd,
                                                                                      ctg_pos_in_scf  => ${$args{scaffold_ref}}{$ctg_name}->{'position'},
                                                                                      feature         => ''
                                                                                   )};
         ${$args{scaffold_ref}}{$ctg_name}->{'start'} = $left_bd;
      } else {
         ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$counter_contig}->{'1'} = {( ctg_name        => $ctg_name,
                                                                                      ctg_orientation => $orientation,
                                                                                      ctg_left_bd     => $left_bd,
                                                                                      ctg_right_bd    => $right_bd,
                                                                                      ctg_sequence    => $seq,
                                                                                      fasta_header    => $fsa_head,
                                                                                      scf_counter     => ${$args{scaffold_ref}}{$ctg_name}->{'scf_counter'},
                                                                                      scf_start       => '',
                                                                                      ctg_pos_in_scf  => ${$args{scaffold_ref}}{$ctg_name}->{'position'},
                                                                                      feature         => ''
                                                                                   )}; #scaffold start needs to be back referenced to first contig in scaffold
      }
      $counter_contig++;
   }
   return (1) ;
}

sub create_scaffold_info {
   my %args          = @_;
   my $part_number   = 1;
   my $left_bd       = 1;
   my $right_bd      = 0;


   #generate default AGP header
   print AGP '##compatible to agp-version 2.0'."\n".
             "\# ORGANISM: ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{'1'}->{'1'}->{'fasta_header'}\n";
   #iterate through contigs
   #this currently assumes that contigs are ordered, and separated by default GAMOLA non-bleeding separator
   #gaps between contigs are set to default of 100 nt as per NCBI guideline
   #no linkage bewteen contigs is assumed
   #no actual scaffold information is used, as this information is not captured by the Genbank annotation
   my $scaffold_tracker = 0;
   my $scf              = '';
   my $gap              = '';
   foreach my $counter_contig (sort {$a <=> $b} keys %{ ${$args{wgs_hash_ref}}->{$args{fasta_header}} }) {
      #contig
      $part_number++;
      #print "\n... ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$counter_contig}->{'1'}->{'scf_counter'} ...";
      if (${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$counter_contig}->{'1'}->{'scf_counter'} != $scaffold_tracker) {
         $left_bd = 1;
         $scaffold_tracker = ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$counter_contig}->{'1'}->{'scf_counter'};
         $gap = '';
         $part_number = 1;
      } else {
         print AGP $gap;
         $gap = '';
      }
      $right_bd =  $left_bd + ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$counter_contig}->{'1'}->{'ctg_right_bd'} - ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$counter_contig}->{'1'}->{'ctg_left_bd'};
      print AGP $args{locus_tag} . '_scaffold'. ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$counter_contig}->{'1'}->{'scf_counter'} ."\t";
      print AGP $left_bd  ."\t";
      print AGP $right_bd ."\t";
      print AGP $part_number."\t";
      print AGP 'W'."\t";
      print AGP ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$counter_contig}->{'1'}->{'ctg_name'}."\t";
      print AGP '1'."\t";
      print AGP ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$counter_contig}->{'1'}->{'ctg_right_bd'} - ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$counter_contig}->{'1'}->{'ctg_left_bd'} + 1 ."\t";
      print AGP ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$counter_contig}->{'1'}->{'ctg_orientation'}."\n";
      #gap
      $part_number++;
      $left_bd  = $right_bd + 1;
      $right_bd += 100;
      $gap .= $args{locus_tag} . '_scaffold'. ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$counter_contig}->{'1'}->{'scf_counter'} ."\t";
      $gap .= $left_bd  . "\t";
      $gap .= $right_bd . "\t";
      $gap .= $part_number."\t";
      $gap .= 'U'."\t";
      $gap .= '100'."\t";
      $gap .= 'contig'."\t";
      $gap .= 'no'."\t";
      $gap .= 'na'."\n";
      $left_bd = $right_bd + 1;
   }
   close AGP;
   return (1);

}


sub get_tRNA_annotation {
   my %args = @_;
   my ($annotation, $tRNA, $anticodon);

   'reset' =~ m/reset/;
   ${$args{genbank_ref}}{$args{entry}} =~ m/\/product\=\"([^\"]*?)\"/s;
   $annotation = $1;
   if (!defined $annotation) {
      ${$args{genbank_ref}}{$args{entry}} =~ m/\/gene\=\"([^\"]*?)\"/s;
      $annotation = $1;
   }
   if (!defined $annotation || $annotation !~ /\w+/) {
      $annotation = 'Pseudo';
   }
   $annotation =~ s/\s+/ /gs;
   $annotation =~ s/[\[\]]//g;

   #get tRNA and anticodon
   #does it come from tRNScan?
   if ($annotation =~ m/Anticodon/) {
      'reset' =~ m/reset/;
      $annotation =~ m/tRNA:\s+(\w+)\;\s+(Anticodon\:\s+\w+)/;
      ($tRNA, $anticodon) = ($1, $2);
   }
   #assume ncRNA scan results
   else {
      'reset' =~ m/reset/;
      $annotation =~ m/RF(\d+)\,\s+(Function\:\s+\w+)/;
      ($tRNA, $anticodon) = ($1, $2);
      $tRNA = 'RF'.$tRNA;
   }

   #$annotation = "product\ttRNA\-$annotation\n";
   $annotation = "product\ttRNA\-$tRNA\n\t\t\tnote\t$anticodon\n";

   return ("\t\t\t".$annotation);
}

sub get_rRNA_annotation {
   my %args = @_;
   my ($annotation);

   'reset' =~ m/reset/;
   ${$args{genbank_ref}}{$args{entry}} =~ m/\/product\=\"([^\"]*?)\"/s;
   $annotation = $1;
   if (!defined $annotation) {
      ${$args{genbank_ref}}{$args{entry}} =~ m/\/gene\=\"([^\"]*?)\"/s;
      $annotation = $1;
   }
   if (!defined $annotation || $annotation !~ /\w+/) {
      $annotation = 'ribosomal RNA';
   }
   $annotation =~ s/\s+/ /gs;
   $annotation =~ s/[\[\]]//g;

   $annotation = "\t\t\tproduct\t$annotation\n";
   return ($annotation);
}

sub get_miscRNA_annotation {
   my %args = @_;
   my ($annotation);

   'reset' =~ m/reset/;
   ${$args{genbank_ref}}{$args{entry}} =~ m/\/product\=\"([^\"]*?)\"/s;
   $annotation = $1;
   if (!defined $annotation) {
      ${$args{genbank_ref}}{$args{entry}} =~ m/\/gene\=\"([^\"]*?)\"/s;
      $annotation = $1;
   }
   if (!defined $annotation || $annotation !~ /\w+/) {
      $annotation = '';
   }
   $annotation =~ s/\s+/ /gs;
   $annotation =~ s/[\[\]]//g;

   $annotation = "\t\t\tproduct\tncRNA: $annotation\n";
   return ($annotation);
}

sub get_CDS_annotation {
   my %args = @_;
   my (@features, @mods, $annotation);
   $annotation = '';

   #includes notes if seleted
   if ($args{notes} == 1) {
      push (@{$args{qualifiers_ref}}, 'note');
   }

   #removing duplicates from start and feature lists
   my %seen = ();
   @{$args{qualifiers_ref}} = grep { ! $seen{$_} ++ } @{$args{qualifiers_ref}};

   #get product if forced
   if ($args{force_product} == 1) {
      'reset' =~ m/reset/;
      ${$args{genbank_ref}}{$args{entry}} =~ m/\/product\=\"([^\"]*?)\"/s;
      $annotation = $1;
      if (!defined $annotation) {
         ${$args{genbank_ref}}{$args{entry}} =~ m/\/gene\=\"([^\"]*?)\"/s;
         $annotation = $1;
      }
   } else {
      'reset' =~ m/reset/;
      ${$args{genbank_ref}}{$args{entry}} =~ m/\/gene\=\"([^\"]*?)\"/s;
      $annotation = $1;
      if (!defined $annotation) {
         ${$args{genbank_ref}}{$args{entry}} =~ m/\/product\=\"([^\"]*?)\"/s;
         $annotation = $1;
      }
   }
   if (!defined $annotation) {
      $annotation = 'hypothetical protein';
   }
   $annotation =~ s/\s+/ /gs;
   $annotation =~ s/_\s*\d+\s*$//gs;
   $annotation =~ s/[\[\]]//g;

   #remove Length, Score and evalue; not allowed in SequIn
   $annotation =~ s/Length\s*\=\s*\d+.+$//;

   #rephrase annotation for hypotheticals
   'reset' =~ m/reset/;
   if ($annotation =~ m/(cons\.?\s*hypo|hyp\.?\s*prot\.?|conserved hypothetical|unknown|hypothetical protein|predicted protein|putative protein|conserved domain protein|putative uncharacterized protein)/i) {
      $annotation = 'hypothetical protein';
   }
   $annotation =~ s/prot\./protein/gi;
   $annotation =~ s/hyp\./hypothetical/gi;
   $annotation =~ s/ dna / DNA /gi;
   $annotation =~ s/ rna / RNA /gi;
   $annotation =~ s/[\.\,]\s*(putative)?\s*$//g;

   #eliminate some other words
   $annotation =~ s/(putative|put\.?\s+|uncharacterised|Chloroplast|\_)/ /i;
   $annotation =~ s/\,\s*partial\s*/ /gs;
   $annotation =~ s/\s+/ /gs;
   $annotation =~ s/^\s+//;
   $annotation =~ s/\(Role\:\d+\)\s*\:?\s*//g;
   $annotation =~ s/(domain)\s*$/$1 protein/;
   $annotation =~ s/(repeat)\s*$/$1 protein/;
   $annotation =~ s/\s*\,?\s*\.\s*\,?\s*/ /g;
   $annotation =~ s/genes/gene/g;
   $annotation =~ s/;/,/g;
   $annotation =~ s/[\.\,]\s*$//;
   $annotation =~ s/[\<\>]//g;
   $annotation =~ s/(possible|probable|uncharacterized)/putative/g;
   unless ($annotation =~ /hypothetical protein/i) {
      $annotation =~ s/hypothetical//i;
      $annotation =~ s/\s+/ /gs;
      $annotation =~ s/^\s+//;
   }

   #check for string length
   if (length ($annotation) > 100 ) {
      $annotation =~ s/\,[^\,]+\,//;
   }
   if (length ($annotation) > 100 ) {
      $annotation =~ s/\,.+$//;
   }
   while (length ($annotation) > 100 ) {
      $annotation =~ s/\s+\S+\s*$//;
   }

   #format annotation
   #$lineform = $annotation;
   #$annotation = &print_line;
   #$annotation =~ s/\n/\n\t/gs;

   #set normal annotation
   $annotation = "\t\t\tproduct\t$annotation\n";

   #generate special case for predicted secreted and membrane proteins
   if ($annotation =~ m/^product\t\s*(secreted protein|membrane protein|extracellular protein|transmembrane protein)\s*$/) {
      $annotation =~ s/^\t\t\tproduct/\t\t\tnote/;
      $annotation = "\t\t\tproduct\thypothetical protein\n".$annotation;
   }

   #other qualifiers present?
   foreach my $modifier (@{$args{qualifiers_ref}}) {
      'reset' =~ m/reset/;
      ${$args{genbank_ref}}{$args{entry}} =~ m/\/$modifier\=\"([^\"]*?)\"/s; #catch regular entry with ""
      my $hit = $1;
      unless (defined $hit) {
         'reset' =~ m/reset/;
         ${$args{genbank_ref}}{$args{entry}} =~ m/\/$modifier\=(.*?)\s\//s; #catch entries without ""
         $hit = $1;
         if (defined $hit) {
            $hit =~ s/\s+/ /gs;
            $hit =~ s/\///s;
            $hit =~ s/\s+$//g;
         }
      }
      unless (defined $hit) {
         'reset' =~ m/reset/;
         ${$args{genbank_ref}}{$args{entry}} =~ m/\/$modifier\=(.+)/s; #catch entries at the end of a feature
         $hit = $1;
         if (defined $hit) {
            $hit =~ s/\s+/ /gs;
            $hit =~ s/\///s;
            $hit =~ s/\s+$//g;
         }
      }

      if (defined $hit && $hit =~ /\w+/) {
         $hit =~ s/\s+/ /gs;
         $hit =~ s/\///s;
         $hit =~ s/\s+$//g;
         $hit = "\t\t\t$modifier\t$hit\n";
      } else {
         $hit = '';
      }
      push (@mods, $hit);
   }

   #merge all hits
   foreach my $hit (@mods) {
      $annotation .= $hit;
   }

   #add notes for PFam, TIGRfam and COG matches if selected
   if (defined ${$args{sel_features_ref}}{'COG_match'}  ||
       defined ${$args{sel_features_ref}}{'PFAM_match'} ||
       defined ${$args{sel_features_ref}}{'TIGR_match'} ) {
      my ($CDS_left_bd, $CDS_right_bd);
      #get proper left and right boundaries from CDS feature
      if ($args{orientation} eq 'sense') {
         $CDS_left_bd  = $args{start};
         $CDS_right_bd = $args{stop};
      } else {
         $CDS_left_bd  = $args{stop};
         $CDS_right_bd = $args{start};
      }
      #get amximum range in case of joint CDS
      $CDS_left_bd =~ s/^(\d+).*/$1/;
      $CDS_right_bd =~ s/.*_(\d+)$/$1/;
      my @functional_entries = ();
      #grep functional hits
      if (defined ${$args{sel_features_ref}}{'COG_match'}) {
         my @temp = grep { /COG~match/ } @{$args{'sorted_features_ref'}};
         @functional_entries = (@functional_entries, @temp);
         undef @temp;
      }
      if (defined ${$args{sel_features_ref}}{'PFAM_match'}) {
         my @temp = grep { /PFAM~match/ } @{$args{'sorted_features_ref'}};
         @functional_entries = (@functional_entries, @temp);
         undef @temp;
      }
      if (defined ${$args{sel_features_ref}}{'TIGR_match'}) {
         my @temp = grep { /TIGR~match/ } @{$args{'sorted_features_ref'}};
         @functional_entries = (@functional_entries, @temp);
         undef @temp;
      }

      #iterate through entries
      my %min_evalue       = ();
      my %local_annotation = ();
      $min_evalue{'COG~match'}  = 1000;
      $min_evalue{'PFAM~match'} = 1000;
      $min_evalue{'TIGR~match'} = 1000;
      foreach my $entry (@functional_entries) {
         my ($local_left_bd, $local_right_bd, $local_feature, $local_label, $local_evalue);
         'reset' =~ m/reset/;
         $entry =~ m/(\d+)_(\d+)_(\d+\-)?([^_]+)_/;
         ($local_left_bd, $local_right_bd, $local_feature) = ($1, $2, $4);

         #remove counter in entry
         $entry =~ s/_\d+\-/_/;

         next if (($local_left_bd < $CDS_left_bd  && $local_right_bd < $CDS_left_bd) ||
                  ($local_left_bd > $CDS_right_bd && $local_right_bd > $CDS_right_bd));
         next unless (defined ${$args{genbank_ref}}{$entry});

         #parse entry
         'reset' =~ m/reset/;
         ${$args{genbank_ref}}{$entry} =~ m/\/label\=\s*([^\/]+)\//s;
         $local_label = $1;
         $local_label =~ s/\s{2,}/ /gs;
         $local_label =~ s/\"//g;

         #get evalue
         'reset' =~ m/reset/;
         ${$args{genbank_ref}}{$entry} =~ m/\/product\=.*Expect\s*\=\s*([^\"]+)\"/s;
         $local_evalue = $1;
         unless (defined $local_evalue) {$local_evalue = 1000};
         $local_evalue =~ s/\s+//;
         if ($local_evalue < $min_evalue{$local_feature}) {
            $min_evalue{$local_feature} = $local_evalue;
            my $label = $local_feature;
            $label =~ s/~match/ label/;
            $local_label =~ s/;\s+GO\:.+//;
            $local_annotation{$local_feature} = "\t\t\tnote\t; $label\: $local_label\n";

            #test if TIGRfam has GO annotation; add only if CDS down not have GO qualifiers on its own
            if ($entry =~ m/TIGR/i &&
                ${$args{genbank_ref}}{$entry} =~ /\/go_from_interpro/ &&
                ${$args{genbank_ref}}{$entry} !~ /\/(go_component|go_process|go_function)/i) {
               my ($GO, $go_component, $go_process, $go_function);
               'reset' =~ m/reset/;
               ${$args{genbank_ref}}{$entry} =~ m/\/go_from_interpro\=\"([^\"]+)\"/;
               $GO = $1;
               unless (defined ($GO)) {$GO = ''};

               #get component
               'reset' =~ m/reset/;
               $GO =~ m/(GO\:\d+)\,/;
               $go_component = $1;
               if (defined $go_component && $go_component =~ /\w+/) {
                  $local_annotation{$local_feature} .= "\t\t\tgo_component\t$go_component\n";
               }

               #get process
               'reset' =~ m/reset/;
               $GO =~ m/GO\:\d+\,([^\;]+)\;/;
               $go_process = $1;
               if (defined $go_process && $go_process =~ /\w+/) {
                  $go_process =~ s/\n//gs;
                  $go_process =~ s/\s{2,}/ /g;
                  $local_annotation{$local_feature} .= "\t\t\tgo_process\t$go_process\n";
               }

               #get function
               'reset' =~ m/reset/;
               $GO =~ m/GO\:\d+\,[^\;]+\;(.+)/;
               $go_function = $1;
               if (defined $go_function && $go_function =~ /\w+/) {
                  $go_function =~ s/\n//gs;
                  $go_function =~ s/\s{2,}/ /g;
                  $local_annotation{$local_feature} .= "\t\t\tgo_function\t$go_function\n";
               }
            }
         }
      }
      foreach my $key (sort keys %local_annotation) {
         $annotation .= $local_annotation{$key};
      }
      undef %local_annotation;
      undef %min_evalue;
   }

   return ($annotation);
}

sub get_gene_annotation {
   my %args = @_;
   my ($annotation);

   'reset' =~ m/reset/;
   ${$args{genbank_ref}}{$args{entry}} =~ m/\/gene\=\"([^\"]*?)\"/s;
   $annotation = $1;
   if (!defined $annotation) {
      ${$args{genbank_ref}}{$args{entry}} =~ m/\/product\=\"([^\"]*?)\"/s;
      $annotation = $1;
   }
   if (!defined $annotation) {
      ${$args{genbank_ref}}{$args{entry}} =~ m/\/note\=\"([^\"]*?)\"/s;
      $annotation = $1;
   }
   unless (defined $annotation) {return ('')};
   $annotation =~ s/\s+/ /gs;
   $annotation =~ s/_\s*\d+\s*$//gs;
   $annotation =~ s/[\[\]]//g;
   $annotation =~ s/[\<\>]//g;
   $annotation =~ s/[\.\,]\s*(putative)?\s*$//g;
   if (length $annotation > 5) {return ('')};
   if ($annotation !~ /\w+/) {return ('')};

   #change first letter to lower case
   $annotation = lcfirst($annotation);

   #add some extra parsing expection
   if ($annotation =~ /^(hyp|hypo|put|prot|con|cons)$/) {return ('')}; #no annotation if hypothetical
   ${$args{unique_gene_anno_ref}}{$annotation}++;
   if (${$args{unique_gene_anno_ref}}{$annotation} == 1) {
      $annotation = "\t\t\tgene\t$annotation\n";
   } else {
      $annotation = "\t\t\tgene\t$annotation"."${$args{unique_gene_anno_ref}}{$annotation}\n";
   }

   #return ("\t".$annotation);
   return ($annotation);
}

sub get_standard_annotation {
   my %args = @_;
   my (@qualifiers, $annotation);
   $annotation = '';

   #get all qualifiers
   'reset' =~ m/reset/;
   @qualifiers = (${$args{genbank_ref}}{$args{entry}} =~ m/  \/([^\=]*?)\="/gs);

   #removing duplicates from start and feature lists
   my %seen = ();
   @qualifiers = grep { ! $seen{$_} ++ } @qualifiers;
   undef %seen;

   foreach my $entry (@qualifiers) {
      my (@local_qualifiers, $label);
      'reset' =~ m/reset/;
      @local_qualifiers = (${$args{genbank_ref}}{$args{entry}} =~ m/\/$entry\=\"([^\"]+)\"/sg);
      foreach my $label (@local_qualifiers) {
         #$label = $1;
         #unless (defined $label) {return ('')};
         $label =~ s/\s+/ /gs;
         $label =~ s/_\s*\d+\s*$//gs;
         $label =~ s/[\[\]]//g;
         $label =~ s/\<\>//g;
         if ($label !~ /\w+/) {return ('')};

         $annotation .= "\t\t\t$entry\t$label\n";
      }
   }

   return ($annotation);
}

sub locus_tag {
   my %args = @_;
   my ($locus_tag, $locus_id);

   #re-use existing locus tag in entry?
   if ($args{reuse_locustag} == 1 && ${$args{genbank_ref}}{$args{entry}} =~ /\/locus_tag/) {
      'reset' =~ m/reset/;
      ${$args{genbank_ref}}{$args{entry}} =~ m/\/locus_tag\=\".*?(\d+)\s*\"/s;
      $locus_id = $1;
   }

   #new locus tag for gene and CDS
   if ($args{reuse_locustag} == 0) { # && $args{table_feature} =~ m/(gene|CDS)/) {
      $locus_id = $args{locus_counter};
   }
   #new locus tag for all entries without a defined locus tag
   unless ($locus_id =~ m/\d+/) {
      $locus_id = $args{locus_counter};
   }

   #remove any preceding '0's
   $locus_id =~ s/^0//g;

   #add sufficient leading '0's for proper format
   $locus_id = '0' x ($args{max_number} - length($locus_id)).$locus_id;

   #create locus tag for CDS
   if ($args{table_feature} =~ /CDS/) {
      $locus_tag = 'protein_id'."\t".'gnl|'.$args{db_name}.'|'.$args{locus_tag}.'_'.$args{roman}.$locus_id."\n";
   }
   #create locus tags for everything else
   else {
      $locus_tag = 'locus_tag'."\t".$args{locus_tag}.'_'.$args{roman}.$locus_id."\n";
   }
   return ($locus_tag);
}

sub boundaries {
   my %args = @_;
   my ($boundary, $bd_values, @temp, $orientation);

   #reset boundaries
   $start = '';
   $stop  = '';
   @start = ();
   @stop  = ();

   #get boundary line
   'reset' =~ m/reset/;
   $args{feature} =~ m/^\s*\S+\s+([^\/]*?)\//s;
   $boundary = $1;
   #if not defined try again in case there are no other qualifiers present
   unless (defined $boundary) {
      $args{feature} =~ m/^\s*\S+\s+(\S.+)$/s;
      $boundary = $1;
   }
   $boundary =~ s/\s+/ /gs;

   #get boundary numbers
   'reset' =~ m/reset/;
   $boundary =~ m/([<>]?\d+.+)/;
   $bd_values = $1;
   $bd_values =~ s/\s//gs;
   $bd_values =~ s/\)//;
   @temp = split /\,/,$bd_values;

   #orientation
   if ($boundary =~ /complement/) {
      $orientation = 'antisense';
   } else {
      $orientation = 'sense';
   }

   #iterate through individual boundaries
   foreach my $entry (@temp) {
      my ($left_bd, $right_bd);
      'reset' =~ m/reset/;
      $entry =~ m/([><]?\d+)\.\.([><]?\d+)/;
      $left_bd = $1;
      $right_bd = $2;
      if ($orientation eq 'sense') {
         push (@start, $left_bd);
         push (@stop,  $right_bd);
         $start .= $left_bd.'_';
         $stop  .= $right_bd.'_';
      } else {
         push (@start, $right_bd);
         push (@stop,  $left_bd);
         $start .= $right_bd.'_';
         $stop  .= $left_bd.'_';
      }
   }
   $start =~ s/_$//s;
   $stop  =~ s/_$//s;
   return ($start, $stop, $orientation);
}


sub populate {
   my %args = @_;

   (my $selected_files) = ${$args{main_window}}->getOpenFile(-initialdir => $current_dir,
                                                             -title      => 'Select files for custom database',
                                                             -multiple   => 1);
   if (defined(@{$selected_files})) {
      #generate short_list
      foreach my $file (@{$selected_files}) {
         'reset' =~ m/reset/;
         $file =~ m/(.+)\/([^\/]+)$/; #display only file_name
         $current_dir = $1;
         my $short = $2;
         push (@short_files, $short);
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
         'reset' =~ m/reset/;
         m/\/([^\/]+)$/; #display only file_name
         my $short = $1;
         $all_files{$short} = $_;
      }

      #update display
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

sub delete_duplicates {
   my %args = @_;
   my %seen = ();
   @{$args{array}} = grep { ! $seen{$_} ++ } @{$args{array}};
}

sub print_feature {
   my ($format);
   open FEATURES, '>', \$format;
   write FEATURES;
   close FEATURES;

   #$format =~ s/\t?____\t?/\t/gs;
   $format =~ s/[ \t]?____[ \t]?/\t/gs;
   $format =~ s/\t\n/\n/gs;
   $format =~ s/ {2,}//gs;
   print FEATURE $format;
}

sub print_feature_WGS {
   my %args = @_;
   my ($format, $correction, $scf_start);

   #find correct scaffold start
   #my $ctg_carryover;
   #foreach my $ctg_name (sort {$a<=>$b} keys %{ $args{scaffold_ref} }) {
   #   if (${$args{scaffold_ref}}{$ctg_name}->{'scf_counter'} == $args{scf_counter} &&
   #       ${$args{scaffold_ref}}{$ctg_name}->{'position'}    == 1) {
   #         $scf_start = ${$args{scaffold_ref}}{$ctg_name}->{'start'};
   #         $ctg_carryover = $ctg_name;
   #         last;
   #       }
   #}

   open FEATURES, '>', \$format;
   write FEATURES;
   close FEATURES;
   #print "\n...$format ...\n"; <STDIN>;

   #$format =~ s/\t?____\t?/\t/gs;
   $format =~ s/[ \t]?____[ \t]?/\t/gs;
   $format =~ s/\t\n/\n/gs;
   $format =~ s/ {2,}//gs;
   $format =~ m/(\d+)\D+(\d+)/;
   my ($left_bd, $right_bd) = ($1, $2);
   my ($lefty, $righty) = ($1, $2);
   $left_bd  = $left_bd  - ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$args{ctg_counter}}->{$args{counter_feature}}->{'ctg_left_bd'} + 1;
   $right_bd = $right_bd - ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$args{ctg_counter}}->{$args{counter_feature}}->{'ctg_left_bd'} + 1;
   #$left_bd  = $left_bd  - ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$args{ctg_counter}}->{$args{counter_feature}}->{'scf_start'} + 1 + ($args{gap} * (${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$args{ctg_counter}}->{$args{counter_feature}}->{'ctg_pos_in_scf'} - 1));
   #$right_bd = $right_bd - ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$args{ctg_counter}}->{$args{counter_feature}}->{'scf_start'} + 1 + ($args{gap} * (${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$args{ctg_counter}}->{$args{counter_feature}}->{'ctg_pos_in_scf'} - 1));
   #$left_bd  += ($args{gap} * ($args{ctg_counter} - 1 ));
   #$right_bd += ($args{gap} * ($args{ctg_counter} - 1 ));
   #if ORFs are bleeding into spacer region
   if ($left_bd < 1) {
      $left_bd = '<'.1;
   }
   if ($right_bd < 1) {
      $right_bd = '>'.1;
   }
   if ($left_bd > (${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$args{ctg_counter}}->{$args{counter_feature}}->{'ctg_right_bd'} - ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$args{ctg_counter}}->{$args{counter_feature}}->{'ctg_left_bd'} + 1) ) {
      $left_bd = '<'.(${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$args{ctg_counter}}->{$args{counter_feature}}->{'ctg_right_bd'} - ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$args{ctg_counter}}->{$args{counter_feature}}->{'ctg_left_bd'} + 1);
   }
   if ($right_bd > (${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$args{ctg_counter}}->{$args{counter_feature}}->{'ctg_right_bd'} - ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$args{ctg_counter}}->{$args{counter_feature}}->{'ctg_left_bd'} + 1) ) {
      $right_bd = '>'.(${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$args{ctg_counter}}->{$args{counter_feature}}->{'ctg_right_bd'} - ${$args{wgs_hash_ref}}->{$args{fasta_header}}->{$args{ctg_counter}}->{$args{counter_feature}}->{'ctg_left_bd'} + 1);
   }

   $format =~ s/\d+(\D+)\d+(.+)/$left_bd$1$right_bd$2/;
   print FEATURE $format;
}

sub print_line {
   my ($format);
   open LINEFORM, '>', \$format;
   write LINEFORM;
   close LINEFORM;
   return $format;
}

sub agp_setup {
   my %args = @_;
   my ($agp_setup, $pane, $agp_mapping_file);

   #if (defined $agp_setup && Tk::Exists($agp_setup)) {
   #    $agp_setup->state('normal');
   #} else {
   #   $agp_setup = ${$args{main_window}}->DialogBox(-title   => 'AGP setup',
   #                                                 -buttons => [ 'OK' ],
   #                                                 );
   #}

   #$pane = $agp_setup->Scrolled(qw/Pane -scrollbars osw -height 400 -width 150/)->pack(-expand => 1, -fill => 'both');
   my $types = [ ['txt', ['.txt']],
                 ['All', ['.*'  ]]
               ];
   ($agp_mapping_file) = ${$args{main_window}}->getOpenFile(-initialdir       => $current_dir,
                                                            -title            => 'Select scaffold mapping file',
                                                            -filetypes        => $types,
                                                            -multiple         => 0);
   unless (-e $agp_mapping_file) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Error finding AGP mapping file $agp_mapping_file.\nPlease try again",
                                                    -bitmap  => 'error',
                                                    -buttons => ['ok']);
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      return (0);
   }
   return ($agp_mapping_file);
}

1;
