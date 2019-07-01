#!/opt/ActivePerl-5.8/bin/perl

#verify GUI selection and all variables

package initialise::verify;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&verify_all);
use vars qw();

use Basics::MesgBox;
use Basics::StatusBar;
use GeneModel::blast_contigs;
use GeneModel::make_blastpairs;
use Tk::LabEntry;
use initialise::read_me       qw(:DEFAULT);
use Cwd;
use ProgrammeModules::blast   qw(:DEFAULT);
use Proc::Background;

#local variables
my (%args, $key, $value, $msg, $manual, $message, $pid);

sub verify_all {
   my %args     = @_;
   my $status   = 0;
   my $program  = '';
   my $database = '';
   my $message  = '';

   #create dummy files
   &create_dummy(main_window  => $args{main_window},
                 progress_bar => $args{progress_bar},
                 auto_ini_ref => $args{auto_ini_ref},
                 ini_ref      => $args{ini_ref}
                );

   #check if all directories exist
   while ( ($key, $value) = each %{$args{ini_ref}}) {
      if (${$args{ini_comment}}{$key} =~ /^\#?[A-Z\s]+\:.*?(directory|location)/) {
         $status = &verify_directory(current_dir          => ${$args{ini_ref}}{$key},
                                     main_window          => $args{main_window},
                                     progress_bar         => $args{progress_bar},
                                     directory_descriptor => ${$args{ini_comment}}{$key});

         if ($status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Invalid directory selection',
                                              -message => 'Correct directory pathways in setup',
                                              -icon    => 'error',
                                              -type    => 'ok',
                                              -parent  => ${$args{main_window}});
            &delete_dummy(main_window  => $args{main_window},
                          progress_bar => $args{progress_bar},
                          auto_ini_ref => $args{auto_ini_ref},
                          ini_ref      => $args{ini_ref}
                         );
            return;
         }
         $status = 0;
      }
   }

   #check if selected options do have their required executables
   {
      ($status, $program) = &verify_selected_executables(main_window  => $args{main_window},
                                                         progress_bar => $args{progress_bar},
                                                         auto_ini_ref => $args{auto_ini_ref},
                                                         ini_ref      => $args{ini_ref}
                                                        );
      if ($status == 0) {
         ${$args{main_window}}->messageBox(-title   => 'Missing executable',
                                           -message => "$program cannot be found in specified location",
                                           -icon    => 'error',
                                           -type    => 'ok',
                                           -parent  => ${$args{main_window}});
         &delete_dummy(main_window  => $args{main_window},
                       progress_bar => $args{progress_bar},
                       auto_ini_ref => $args{auto_ini_ref},
                       ini_ref      => $args{ini_ref}
                      );
         return;
      }
      $status = 0;
      $program = "";
   }

   #check if all databases and parsing files are present for selected options
   {
      ($status, $message) = &verify_database_dependency(main_window  => $args{main_window},
                                                        progress_bar => $args{progress_bar},
                                                        auto_ini_ref => $args{auto_ini_ref},
                                                        ini_ref      => $args{ini_ref},
                                                        directory    => $args{directory}
                                                       );

      if ($status == 0) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "$message",
                                           -icon    => 'error',
                                           -type    => 'ok',
                                           -parent  => ${$args{main_window}});
         &delete_dummy(main_window  => $args{main_window},
                       progress_bar => $args{progress_bar},
                       auto_ini_ref => $args{auto_ini_ref},
                       ini_ref      => $args{ini_ref}
                      );
         return;
      }
      $status = 0;
      $database = "";
   }

   #check for structural analyses
   {
      ($status, $program) = &verify_structural_analyses(main_window  => $args{main_window},
                                                        progress_bar => $args{progress_bar},
                                                        auto_ini_ref => $args{auto_ini_ref},
                                                        ini_ref      => $args{ini_ref},
                                                        directory    => $args{directory}
                                                       );

      if ($status == 0) {
         ${$args{main_window}}->messageBox(-title   => 'Missing structural executable',
                                           -message => "$program cannot be found in specified location",
                                           -icon    => 'error',
                                           -type    => 'ok',
                                           -parent  => ${$args{main_window}});
         &delete_dummy(main_window  => $args{main_window},
                       progress_bar => $args{progress_bar},
                       auto_ini_ref => $args{auto_ini_ref},
                       ini_ref      => $args{ini_ref}
                      );
         return;
      }
      $status = 0;
      $program = "";
   }

   #delete dummy files
   &delete_dummy(main_window  => $args{main_window},
                 progress_bar => $args{progress_bar},
                 auto_ini_ref => $args{auto_ini_ref},
                 ini_ref      => $args{ini_ref}
                );

}

sub verify_directory {
   my %args = @_;
   my $status = "";
   my $char_length = length($args{current_dir}) + 20;

   #verify directory exists
   if (-e $args{current_dir}) {
      ${$args{progress_bar}}->configure(-label=>"Verifying existence of directory $args{current_dir}...OK");
      ${$args{progress_bar}}->after(50);
      ${$args{main_window}}->update;
      return (1);
   } else {
      ${$args{progress_bar}}->configure(-label=>"Verifying existence of directory $args{current_dir}...not found");
      my $msg = ${$args{main_window}}->MesgBox(
                                            -title     => 'Directory not found',
                                            -text      => "Enter path manually?",
                                            -icon      => 'QUESTION',
                                            -buttons   => ['Yes', 'No'],
                                            -defbutton => 'Yes',
                                           );
      if ($msg->Show eq 'Yes') {
         $manual = ${$args{main_window}}->DialogBox(-title          => "Pathway for $args{directory_descriptor}",
                                                    -buttons        => ['OK', 'Cancel'],
                                                    -default_button => 'OK');
         $manual -> add('LabEntry', -textvariable => \$args{current_dir},
                                    -width        => $char_length,
                                    -label        => 'Pathway',
                                    -labelPack    => [-side => 'left']) -> pack();
         #if cancel, then no new pathway selected, generate error message and restart verification
         if ($manual->Show eq 'Cancel') {
            ${$args{main_window}}->MesgBox(
                                         -title     => 'Invalid selection',
                                         -text      => "No pathway selected",
                                         -icon      => 'ERROR',
                                         -buttons   => ['OK'],
                                         -defbutton => 'OK',
                                        );
         }
         #verify new selection or give chance to re-enter in case of cancel
         &verify_directory(%args);
      } else {
         #return bad status
         return (0);
      }
   }
}

sub verify_selected_executables {
   my %args = @_;

   #Glimmer2
   if (${$args{auto_ini_ref}}{internal_gm} == 1 && ${$args{auto_ini_ref}}{use_glimmer2} == 1) {
      unless (-e ${$args{ini_ref}}{glprog2}.'/glimmer2') {
         return (0, "Glimmer2");
      }
   }

   #Glimmer3
   if (${$args{auto_ini_ref}}{internal_gm} == 1 && ${$args{auto_ini_ref}}{use_glimmer3} == 1) {
      unless (-e ${$args{ini_ref}}{glprog3}.'/glimmer3') {
         return (0, "Glimmer3");
      }
   }

   #Critica
   if (${$args{auto_ini_ref}}{internal_gm} == 1 && ${$args{auto_ini_ref}}{use_critica} == 1) {
      unless (-e ${$args{auto_ini_ref}}{critica_bin}.'/critica' && -e ${$args{auto_ini_ref}}{critica_scripts}.'/blast-contigs') {
         return (0, "Critica");
      }
   }

   #Blast
   if (${$args{auto_ini_ref}}{blast_selector} == 1 && ${$args{auto_ini_ref}}{blast_type} =~ /(BlastP|gappedBlastP|BlastN|tBlastX)/) {
      unless (-e ${$args{ini_ref}}{blast_executables}.'/blastall') {
         return (0, "Blast: blastall");
      }
   }
   if (${$args{auto_ini_ref}}{blast_selector} == 1 && ${$args{auto_ini_ref}}{blast_type} =~ /PSI-Blast/) {
      unless (-e ${$args{ini_ref}}{blast_executables}.'/blastall' && -e ${$args{ini_ref}}{blast_executables}.'/blastpgp') {
         return (0, "PSI-Blast\: blastall or blastpgp");
      }
   }

   #COG
   if (${$args{auto_ini_ref}}{COG_selector} == 1) {
      unless (-e ${$args{ini_ref}}{blast_executables}.'/blastall') {
         return (0, "COG\: blastall");
      }
   }

   #PFam
   if (${$args{auto_ini_ref}}{Pfam_selector} == 1) {
      unless (-e ${$args{ini_ref}}{pfam_executable}.'/hmmpfam') {
         return (0, "Pfam");
      }
   }

   #TIGRFam
   if (${$args{auto_ini_ref}}{TIGRfam_selector} == 1) {
      unless (-e ${$args{ini_ref}}{pfam_executable}.'/hmmpfam') {
         return (0, "TIGRfam");
      }
   }

   #tRNAscan
   if (${$args{auto_ini_ref}}{trna} == 1) {
      unless (-e ${$args{ini_ref}}{trnascan_dir}.'/tRNAscan-SE' || -e ${$args{ini_ref}}{trnascan_dir}.'/bin/tRNAscan-SE') {
         return (0, "tRNAScan-SE");
      }
   }

   #Transmembrane predictions
   if (${$args{auto_ini_ref}}{tmhmm} == 1) {
      unless (-e ${$args{ini_ref}}{tmhmm_dir}.'/bin/decodeanhmm') {
         return (0, "TMHMM");
      }
   }

   #SignalP predictions
   if (${$args{auto_ini_ref}}{signalp} == 1) {
      unless (-e ${$args{ini_ref}}{signalP_dir}.'/signalp') {
         return (0, "SignalP");
      }
   }

   #Terminator predictions
   if (${$args{auto_ini_ref}}{terminator_selector} == 1) {
      unless (-e ${$args{ini_ref}}{transterm_dir}.'/transterm') {
         return (0, "Transterm");
      }
   }

   #CIRSPR prediction
   if (${$args{auto_ini_ref}}{CRISPR_selector} == 1) {
      unless (-e ${$args{ini_ref}}{CRISPR_dir}.'/CRT_java.jar') {
         return (0, "CRISPR");
      }
   }

   return(1);
}

sub verify_database_dependency {
   my %args = (status => 0, @_);
   my ($current_dir, @results, @infofile, $status);

   $status = 0;

   #Glimmer2 model file, if not training on the fly
   if (${$args{auto_ini_ref}}{internal_gm}        == 1 &&
       ${$args{auto_ini_ref}}{use_glimmer2}       == 1 &&
       ${$args{auto_ini_ref}}{make_training_file} == 0) {
      ${$args{progress_bar}}->configure(-label=>"Verifying Glimmer2 setup");
      ${$args{main_window}}->update;

      #Model file still there?
      unless (-e ${$args{auto_ini_ref}}{selected_gl_model}) {
         return (0, "Selected Glimmer2 model file ${$args{auto_ini_ref}}{selected_gl_model} is missing.");
      }

      #model valid for glimmer?
      `${$args{ini_ref}}{glprog2}/glimmer2 ${$args{ini_ref}}{tempfile}  ${$args{auto_ini_ref}}{selected_gl_model} -l +f > ${$args{ini_ref}}{gm_output}/dummy.glcoord 2>&1`;

      local( $/, *ENTRY ) ;
      open( ENTRY, ${$args{ini_ref}}{genemodel_output}.'/dummy.glcoord' ) or do
      {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not open file dummy.glcoord.",
                                           -icon    => 'error',
                                           -type    => 'OK');
         return (0, "Selected Glimmer2 model file is missing.");
      };
      my $value = <ENTRY>;
      close ENTRY;

      unless ($value =~ /\w+/) {
         return (0, "Selected Glimmer2 model file ${$args{auto_ini_ref}}{selected_gl_model} is corrupt.");
      }
      ${$args{progress_bar}}->configure(-label=>"Verifying Glimmer2 setup...OK");
      ${$args{main_window}}->update;
   }

   #Glimmer3 model file
   if (${$args{auto_ini_ref}}{internal_gm}        == 1 &&
       ${$args{auto_ini_ref}}{use_glimmer3}       == 1 &&
       ${$args{auto_ini_ref}}{make_training_file} == 0) {
      ${$args{progress_bar}}->configure(-label=>"Verifying Glimmer3 setup");
      ${$args{main_window}}->update;

      #Model file still there?
      unless (-e ${$args{auto_ini_ref}}{selected_gl_model}) {
         return (0, "Selected Glimmer3 model file is missing.");
      }

      #model valid for glimmer?
      `${$args{ini_ref}}{glprog3}/glimmer3 -z 11 -g 100 -l ${$args{ini_ref}}{tempfile}  ${$args{auto_ini_ref}}{selected_gl_model} ${$args{ini_ref}}{gm_output}/dummy 2>&1`;

      local( $/, *ENTRY ) ;
      open( ENTRY, ${$args{ini_ref}}{gm_output}.'/dummy.detail' ) or do
      {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                          -message => "Could not open file dummy.glcoord.",
                                          -icon    => 'error',
                                          -type    => 'OK');
         return (0, "Selected Glimmer3 model file is missing.");
      };
      my $value = <ENTRY>;
      close ENTRY;

      unless ($value =~ /\w+/) {
         return (0, "Selected Glimmer3 model file ${$args{auto_ini_ref}}{selected_gl_model} is corrupt.");
      }
      ${$args{progress_bar}}->configure(-label=>"Verifying Glimmer3 setup...OK");
      ${$args{main_window}}->update;
   }

   #Critica and Critica database
   if (${$args{auto_ini_ref}}{internal_gm} == 1 && ${$args{auto_ini_ref}}{use_critica} == 1) {
      ${$args{progress_bar}}->configure(-label=>"Verifying Critica setup");
      ${$args{main_window}}->update;

      #test for file presence
      {
         unless (-e ${$args{ini_ref}}{critica}.'/scripts/blast-contigs') {
            return (0, "Missing Critica component: blast-contigs.");
         }
         unless (-e ${$args{ini_ref}}{critica}.'/scripts/make-blastpairs') {
            return (0, "Missing Critica component: make-blastpairs.");
         }
         unless (-e ${$args{ini_ref}}{critica}.'/scanblastpairs') {
            return (0, "Missing Critica component: scanblastpairs.");
         }
         unless (-e ${$args{ini_ref}}{critica}.'/scripts/iterate-critica') {
            return (0, "Missing Critica component: iterate-critica.");
         }
         unless (-e ${$args{ini_ref}}{critica}.'/critica') {
            return (0, "Missing Critica component: critica.");
         }
         unless (-e ${$args{ini_ref}}{critica}.'/dicodontable') {
            return (0, "Missing Critica component: dicodontable.");
         }
         unless (-e ${$args{ini_ref}}{critica}.'/critica') {
            return (0, "Missing Critica component: critica.");
         }
         unless (-e ${$args{ini_ref}}{critica}.'/addlongorfs') {
            return (0, "Missing Critica component: addlongorfs.");
         }
         unless (-e ${$args{ini_ref}}{critica}.'/removeoverlaps') {
            return (0, "Missing Critica component: removeoverlaps.");
         }
      }

      #test blast-contigs
      {
         #temporarily set dummy critica blastn db
         my $temp_critica_db = ${$args{auto_ini_ref}}{full_critica_db};
         ${$args{auto_ini_ref}}{full_critica_db} = ${$args{auto_ini_ref}}{work_dir}.'/dummy.nt';

         &blastcontigs(main_window   => $args{main_window},
                       progress_bar  => $args{progress_bar},
                       auto_ini_ref  => $args{auto_ini_ref},
                       ini_ref       => $args{ini_ref},
                       filename      => 'temp',
                       directory     => ${$args{auto_ini_ref}}{work_dir},
                      );

         ${$args{auto_ini_ref}}{full_critica_db} = $temp_critica_db;

         local( $/, *ENTRY ) ;
         open( ENTRY, ${$args{auto_ini_ref}}{work_dir}.'/temp.blast' ) or do
         {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                             -message => "Could not open critica blast-contigs results.",
                                             -icon    => 'error',
                                             -type    => 'OK');
            return (0, "Error in Critica execution. Something is wrong with blast-contigs results.");
         };
         my $value = <ENTRY>;
         close ENTRY;

         unless ($value =~ /gb\|\d\|\s+DUMMY/) {
            return (0, "Error in Critica execution. Something is wrong with blast-contigs results.");
         }
      }

      #test make blast pairs
      {
         &blastpair_critica(main_window   => $args{main_window},
                            progress_bar  => $args{progress_bar},
                            auto_ini_ref  => $args{auto_ini_ref},
                            ini_ref       => $args{ini_ref},
                            filename      => 'temp',
                            directory     => ${$args{auto_ini_ref}}{work_dir}
                           );

         local( $/, *ENTRY ) ;
         open( ENTRY, ${$args{auto_ini_ref}}{work_dir}.'/temp.blast.pairs' ) or do
         {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                             -message => "Could not open critica make-blastpairs results.",
                                             -icon    => 'error',
                                             -type    => 'OK');
            return (0, "Error in Critica execution. Something is wrong with make-blastpairs results.");
         };
         my $value = <ENTRY>;
         close ENTRY;

         unless ($value =~ /\>/) {
            return (0, "Error in Critica execution. Something is wrong with make-blastpairs results.");
         }
      }

      #make test for scanblastpairs
      {

         #change to critica directory
         $current_dir = cwd(); #get current working directory
         chdir ${$args{ini_ref}}{critica}.'/scripts'; #change to critica directory for working environment

         `${$args{ini_ref}}{critica}/scanblastpairs  ${$args{auto_ini_ref}}{work_dir}/temp.blast ${$args{auto_ini_ref}}{work_dir}/temp.blast.pairs ${$args{auto_ini_ref}}{work_dir}/temp.triplets`;

         local( $/, *ENTRY ) ;
         open( ENTRY, ${$args{auto_ini_ref}}{work_dir}.'/temp.triplets' ) or do
         {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                             -message => "Could not open critica scanblastpairs results.",
                                             -icon    => 'error',
                                             -type    => 'OK');
            #delete tempfiles if error

            return (0, "Error in Critica execution. Something is wrong with scanblastpairs results.");
         };
         my $value = <ENTRY>;
         close ENTRY;

         unless ($value =~ /DUMMY/) {
            return (0, "Error in Critica execution. Something is wrong with scanblastpairs results.");
         }
      }

      #change back to current directory
      chdir $current_dir;

      ${$args{progress_bar}}->configure(-label=>"Verifying Critica setup...OK");
      ${$args{main_window}}->update;
   }

   #Blast files and databases
   if (${$args{auto_ini_ref}}{blast_selector} == 1) {
      ${$args{progress_bar}}->configure(-label=>"Verifying Blast setup");
      ${$args{main_window}}->update;

      #make dummy gene model
      my @dummy_gm = ();
      push (@dummy_gm, 'temp___1___1___300___sense');
      push (@dummy_gm, 'temp___2___303___900___sense');



      my ($blast_ref) = &blast_file(main_window      => $args{main_window},
                                    progress_bar     => $args{progress_bar},
                                    auto_ini_ref     => $args{auto_ini_ref},
                                    ini_ref          => $args{ini_ref},
                                    directory        => ${$args{auto_ini_ref}}{work_dir},
                                    input_file       => 'temp',
                                    combined_orf_ref => \@dummy_gm
                                   );

      #are all blast executables present?
      unless (-e ${$args{ini_ref}}{blast_executables}.'/blastall') {
         return (0, "Blast executable 'blastall' is missing.");
      }
      unless (-e ${$args{ini_ref}}{blast_executables}.'/blastpgp') {
         return (0, "Blast executable 'blastpgp' is missing.");
      }

      #test if selected database produces valid results
      if (${$args{auto_ini_ref}}{blast_type} =~ /(blastp|gappedblastp|psi-blast)/i) {
         {
            local ($/, *ENTRY);
            open (ENTRY, "<${$args{ini_ref}}{tempfile}_aa");
            $args{sequence} = <ENTRY>;
            close ENTRY;
         }
      } elsif (${$args{auto_ini_ref}}{blast_type} =~ /(blastn|tblastx)/i) {
         {
            local ($/, *ENTRY);
            open (ENTRY, "<${$args{ini_ref}}{tempfile}");
            $args{sequence} = <ENTRY>;
            close ENTRY;
         }
      }

      #start background process for independent progress bar
      my $status = Proc::Background->new('perl', $args{directory}."\/lib\/Basics\/progress_bar.pl", "Blast status");

      my @childs = ();
      my $ref = "";
      #start forking
      my $pid = fork();
      if ($pid) {
         # parent
         push(@childs, $pid);
      } elsif ($pid == 0) {
         # child
         if (${$args{auto_ini_ref}}{blast_type}  =~ /psi-blast/) {
            $ref = &blast_seq(main_window   => $args{main_window},
                              progress_bar  => $args{progress_bar},
                              auto_ini_ref  => $args{auto_ini_ref},
                              ini_ref       => $args{ini_ref},
                              sequence      => $args{sequence},
                              aa_input_file => "${$args{ini_ref}}{tempfile}_aa"
                             );
         } else {
            $ref = &blast_seq(main_window   => $args{main_window},
                              progress_bar  => $args{progress_bar},
                              auto_ini_ref  => $args{auto_ini_ref},
                              ini_ref       => $args{ini_ref},
                              sequence      => $args{sequence}
                             );
         }

         #write temp results file
         open WRITE, "+>${$args{ini_ref}}{tempfile}_dummy";
         print WRITE $ref;
         close WRITE;

         #exit safely
         CORE::exit();
      } else {
         die "Couldn\'t fork\: $!\n";
      }

      #wait
      foreach (@childs) {
         waitpid($_, 0);
      }

      #kill status bar
      $status->die;

      #check status of test blast
      {
         local ($/, *ENTRY);
         open (ENTRY, "<${$args{ini_ref}}{tempfile}_dummy");
         $status = <ENTRY>;
         close ENTRY;
      }

      #error in Blast?
      if ($status =~ /Error in Blast analysis using database/) {
         unlink "${$args{ini_ref}}{tempfile}_dummy";
         return (0, "Selected Blast database  ${$args{auto_ini_ref}}{blast_db} is corrupt.");
      }

      ${$args{progress_bar}}->configure(-label=>"Verifying Blast setup...OK");
      ${$args{main_window}}->update;
   }

   #PFam files and databases
   if (${$args{auto_ini_ref}}{Pfam_selector} == 1) {
      ${$args{progress_bar}}->configure(-label=>"Verifying PFam setup");
      ${$args{main_window}}->update;

      #are all Pfam executables present?
      unless (-e ${$args{ini_ref}}{pfam_executable}.'/hmmpfam') {
         return (0, "Pfam executable 'hmmpfam' is missing.");
      }
      #are local Pfam databases present?
      #generate local array for selected PFam databases
      my @local_db = split/ ; /,${$args{auto_ini_ref}}{full_Pfam_db};
      foreach my $pfam_db (@local_db) {
         unless (-e $pfam_db) {
            return (0, "Selected Pfam database $pfam_db is missing.");
         }
      }

      #are databases for PFAM verbose present?
      if (${$args{auto_ini_ref}}{pfam_verbose} eq 'verbose') {
         unless (-e ${$args{ini_ref}}{pfam_dir}.'/pfamA.txt') {
            return (0, "Pfam code database 'pfamA.txt' is missing.");
         }
         unless (-e ${$args{ini_ref}}{pfam_dir}.'/interpro.txt') {
            return (0, "InterPro reference database 'interpro.txt' is missing.");
         }
      }

      my @childs = ();
      my $ref = "";
      #start forking
      my $pid = fork();
      if ($pid) {
         # parent
         push(@childs, $pid);
      } elsif ($pid == 0) {
         # child

         #iterate over local db's
         foreach my $pfam_db (@local_db) {
            #start background process for independent progress bar
            my $status = Proc::Background->new('perl', $args{directory}."\/lib\/Basics\/progress_bar.pl", "PFam status", "Checking Pfam database $pfam_db");
            my $local_PFAM;
            my $write_result = $pfam_db.'_dummy';
            system "${$args{ini_ref}}{pfam_executable}/hmmpfam $pfam_db ${$args{ini_ref}}{tempfile}_aa > $write_result";
            #kill status bar
            $status->die;
         }
         #exit safely
         CORE::exit();
      } else {
         die "Couldn\’t fork\: $!\n";
      }

      #wait
      foreach (@childs) {
         waitpid($_, 0);
      }

      #iterate over local db's
      foreach my $pfam_db (@local_db) {
         my $local_PFAM = "";
         ${$args{progress_bar}}->configure(-label=>"Verifying PFam setup...verifying results for $pfam_db");
         ${$args{main_window}}->update;
         #check status of test blast
         {
            local( $/, *WRITE_TEST ) ;
            open( WRITE_TEST, "<$pfam_db\_dummy" ) or do {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Could not open Pfam test file $pfam_db\_dummy",
                                                 -icon    => 'error',
                                                 -type    => 'ok',
                                                 -parent  => ${$args{main_window}});
            };
            $local_PFAM = <WRITE_TEST>;
            close WRITE_TEST;
         }

         #correct outcome?
         unless ($local_PFAM =~ /Scores for sequence family classification([^\n]*?\n){4,}\n?Parsed for domains/si) {
            return (0, "Error in checking $pfam_db.\nCheck Pfam setup");
         }
      }
      ${$args{progress_bar}}->configure(-label=>"Verifying Pfam setup...OK");
      ${$args{main_window}}->update;
   }

   #TIGRFam files and databases
   if (${$args{auto_ini_ref}}{TIGRfam_selector} == 1) {
      ${$args{progress_bar}}->configure(-label=>"Verifying TIGRFam setup");
      ${$args{main_window}}->update;

      #are all TIGRfam executables present?
      unless (-e ${$args{ini_ref}}{pfam_executable}.'/hmmpfam') {
         return (0, "TIGRfam executable 'hmmpfam' is missing.");
      }

      #are local TIGRfam databases present?
      #generate local array for selected TIGRFam databases
      my @local_db = split/ ; /,${$args{auto_ini_ref}}{full_TIGRfam_db};
      foreach my $tigrfam_db (@local_db) {
         unless (-e $tigrfam_db) {
            return (0, "Selected TIGRfam database $tigrfam_db is missing.");
         }
      }

      #are databases for TIGRFAM verbose present?
      if (${$args{auto_ini_ref}}{tigrfam_verbose} eq 'verbose + GO') {
         unless (-e ${$args{ini_ref}}{TIGRfam_dir}.'/TIGRFAMS_GO_LINK') {
            return (0, "TIGRfam GO links are missing.");
         }
         unless (-e ${$args{ini_ref}}{TIGRfam_dir}.'/TIGRFAMS_ROLE_LINK') {
            return (0, "TIGRfam role links are missing.");
         }
         unless (-e ${$args{ini_ref}}{TIGRfam_dir}.'/TIGR_ROLE_NAMES') {
            return (0, "TIGRfam role names are missing.");
         }
         #TIGRfam_info present?
         opendir SEQINPUT, ${$args{ini_ref}}{TIGRfam_info} or do {
            return (0, "TIGRfam info directory is missing at ${$args{ini_ref}}{TIGRfam_info}.");
         };
         (@infofile) = grep !/^\./, readdir(SEQINPUT);
         if ($infofile[0] eq "" || $infofile[0] !~ /\.info/i) {
            return (0, "TIGRfam info directory ${$args{ini_ref}}{TIGRfam_info} appears to be empty.");
         }
      }

      my @childs = ();
      my $ref = "";
      #start forking
      my $pid = fork();
      if ($pid) {
         # parent
         push(@childs, $pid);
      } elsif ($pid == 0) {
         # child

         #iterate over local db's
         foreach my $tigrfam_db (@local_db) {
            #start background process for independent progress bar
            my $status = Proc::Background->new('perl', $args{directory}."\/lib\/Basics\/progress_bar.pl", "TIGRFam status", "Checking TIGRfam database $tigrfam_db");
            my $local_TIGRFAM;
            my $write_result = $tigrfam_db.'_dummy';
            system "${$args{ini_ref}}{pfam_executable}/hmmpfam $tigrfam_db ${$args{ini_ref}}{tempfile}_aa > $write_result";
            #kill status bar
            $status->die;
         }
         #exit safely
         CORE::exit();
      } else {
         die "Couldn\’t fork\: $!\n";
      }

      #wait
      foreach (@childs) {
         waitpid($_, 0);
      }
      #iterate over local db's
      foreach my $tigrfam_db (@local_db) {
         my $local_TIGRFAM = "";
         ${$args{progress_bar}}->configure(-label=>"Verifying TIGRFam setup...verifying results for $tigrfam_db");
         ${$args{main_window}}->update;
         #check status of test blast
         {
            local( $/, *WRITE_TEST ) ;
            open( WRITE_TEST, "<$tigrfam_db\_dummy" ) or do {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Could not open TIGRfam test file $tigrfam_db\_dummy",
                                                 -icon    => 'error',
                                                 -type    => 'ok',
                                                 -parent  => ${$args{main_window}});
            };
            $local_TIGRFAM = <WRITE_TEST>;
            close WRITE_TEST;
         }

         #correct outcome?
         unless ($local_TIGRFAM =~ /Scores for sequence family classification([^\n]*?\n){4,}\n?Parsed for domains/si) {
            return (0, "Error in checking $tigrfam_db.\nCheck Pfam setup");
         }
      }
      ${$args{progress_bar}}->configure(-label=>"Verifying TIGRfam setup...OK");
      ${$args{main_window}}->update;
   }

   return (1);
}

sub verify_structural_analyses {
   my %args = @_;

   #tRNA
   if (${$args{auto_ini_ref}}{trna} == 1 ) {
      unless (-e ${$args{ini_ref}}{trnascan_dir}.'/bin/tRNAscan-SE') {
         return (0, "tRNAscan-SE");
      }
   }

   #tmhmm2
   if (${$args{auto_ini_ref}}{tmhmm} == 1 ) {
      unless (-e ${$args{ini_ref}}{tmhmm_dir}.'/bin/decodeanhmm') {
         return (0, "decodeanhmm");
      }
      unless (-e ${$args{ini_ref}}{tmhmm_dir}.'/bin/tmhmmformat.pm') {
         return (0, "tmhmmformat");
      }
      unless (-e ${$args{ini_ref}}{tmhmm_dir}.'/bin/tmhmm.pm') {
         return (0, "tmhmm");
      }
      unless (-e ${$args{ini_ref}}{tmhmm_dir}.'/lib/TMHMM1.0.model') {
         return (0, "TMHMM1.0.model");
      }
      unless (-e ${$args{ini_ref}}{tmhmm_dir}.'/lib/TMHMM2.0.model') {
         return (0, "TMHMM2.0.model");
      }
   }

   #signalp
   if (${$args{auto_ini_ref}}{signalp} == 1 ) {
      unless (-e ${$args{ini_ref}}{signalP_dir}.'/signalp') {
         return (0, "signalp");
      }
   }

   #terminator
   if (${$args{auto_ini_ref}}{terminator} == 1 ) {
      unless (-e ${$args{ini_ref}}{transterm_dir}.'/transterm') {
         return (0, "transterm");
      }
      unless (-e ${$args{ini_ref}}{transterm_dir}.'/expterm.dat') {
         return (0, "expterm.dat");
      }
   }
   return (1);
}


sub create_dummy {
   my %args = @_;
   #generate dummy DNA file for nt database
   open  WRITE, "+>".${$args{ini_ref}}{tempfile} or return (0, "Couldn't write to tempfile.");
   print WRITE "\>gi\|1\|gb\|1\| DUMMY\n".
                  "TTTCTATTATCTTTTTTCCACAAACTTGTGAATTTGTTAGTTATTTTTGTTCCTTGTGGAAAACTTTAAA\n".
                  "AATTCATGCTAAACTTATTTTGACTTAATCTGTGGATAATACAGTACTTGGAGGAATTTTTTTGTTTGAT\n".
                  "TTAGATAAATTTTGGCAATTTTTCAATGCTGAGATGAAAAAAAGCTACAGCACCGTTGCCTATAACGCTT\n".
                  "GGTTTAAAAATACTAAACCTATTTCCTTTAATAAAAAGACAAAAGAAATGATAATTGCTGTTGAATCTCC\n".
                  "AGTTGCAAAAGGGTATTGGGAAAAGAACTTGGCTTCTCAACTGATTCAAGAGGCGTATGCTTATGCAGAC\n".
                  "ATGGAAATTCAACCAAAATTCGAAGTGGCAGGTAAAGAAGGACCTGAACGTCTAGTTACGCCAAAACCAA\n".
                  "TAAGTACACTTTTGATACCTTTGTTCAGGGTGAAGGAAATAAGCTAGCAGCCGGAGCCGCTTTAGCAGTT\n".
                  "GCCGATAACCCAGGAAGTTTTTATAATCCTTTATTTATCTTTGGTGGTGTGGGTTTAGGTAAGACGCACT\n".
                  "TAATGCAAGCTATCGGCCATCAAATGCTAGCCGAAAAACCCCATGCAAAAGTGGTCTACATTCAAAGTGA\n".
                  "AATTGTGATCTGCTATTAGTAGATGATATTCAATTCTTCTCTAAAAAAGAGGGGATTCAAGAGGAGTTTT\n".
                  "TCCATACTTTTGAAACTCTTTATAATGATCAAAAACAGATTGTTATGACTAGTGATCGTCTGCCAACTGA\n".
                  "AATTCCTGAATTATCAGAACGATTGGTTTCTCGATTCGCCTGGGGCTTACAAGTTGAAATTACGCCACCT\n".
                  "GATTTAGAAACTAGAATTGCAATTTTACGTAAAAAAGCTGAAACAGATGGTTTGGCAATTGATGACAGTA\n".
                  "CCCTTGACTATATTGCTTCGCAAGTTGATACCAATATTCGAGAACTCGAAGGCGCTTTAGTTAAAGTTCA\n".
                  "AGCCCATGCAACAATTGAACGTGAAGATATTAACGTTGATTTAGCAAAAGAAGCTTTAGCAGACCTAAAA\n".
                  "CTTGTTCAAAAAAATAGAGGACTTCAAATCTCAAAAATCCAAGAAGTTGTCGCAAATTACTTTCAAACTT\n".
                  "CAACTACAGAACTAAAGGGAAAGAAACGGGTAAAACAAATTGTTGTCCCGCGGCAAATTGCGATGTATTT\n".
                  "GTCACGTGAGTTAACAGATTCCAGCTTACCTAAAATTGGACAAGAGTTTGGTGGCAAAGATCATACCACG\n".
                  "GTTATGCATGCCTGCGATAAAATTTCTCGTGCCTTAAAAACCGATGCCGAAATTAAAGCAGCTGTGTATG\n".
                  "ATTTAAAAGCAATGCTTGAACACTAAAAGTTAATAACATGTGGATAAGTTCTTCTTTTGTCCACATCCAT\n".
                  "TAAACAGACAAAAGTCCGTTTTATTCAGGTGTTTGACTACTTTCCCACAGATTCCACTGGCCCTACTACT\n".
                  "ACTTCTAAGTATTTAAATATATAATTAAATAAATAAGACAAAGTACATATAACTTCACAGGAGGTAAACG\n".
                  "ATGCAGTTTACAATTAATCGTAATTTATTCCTCGAAAACCTAAATAATGTAATGCGTGCAATTTCTTCAC\n".
                  "GTGCAACTATTCCAATTTTAAGTGGTATAAAACTTAACCTTACTGATGAGATGCTAACTTTAACCGGTAG\n".
                  "TGACACTGACATTTCGATCGAAATTCAGATTCCTGTAAACGATGATCTAGTTGTTCAATCTACAGGTTCG\n".
                  "ATTGTTTTACCTGCACGCTTTTTTAGTGAGATCGTTAAAAAATTACCTGGTAAAGACTTTTCATTTGAAG\n".
                  "TAAGAGAAAGCTTCCAGACTAAAATTGTTTCTGAAAATACTGAATTTATGATCAATGGCTTGGATGCAAA\n".
                  "TAACTATCCTCATCTACCAGAAATTTCTACTGATGCCTCATTTAAGATTTCAGGTAAAACTTTTAGAGAA\n".
                  "ATTATTAATGAAACTGTTTTTGCTGTTGCTACTCAAGAAAGTCGTCCTACTTTAACAGGTGTTAATTTCA\n".
                  "TTTTTAATAACTCATCAATTAAGGCAGTTGCTACCGACAGTCACCGTTTATCACAGCGTCAAATTTCTTT\n".
                  "AGAAAACGGCCCACAAACTAGTACTGACTTAATAATTCCAGGAAAGAGTTTAGTAGAATTAGCTAGAATT\n".
                  "ATTGGTGAAAACGATCCTGAAATTACAGTAAATCCAGGTGAGAACCAAGTTTTATTTGAAGTCGGAAATA\n".
                  "TTGCATTCTATTCCCGTTTGCTTGATGGTCAATATCCAGATACTGATCGTTTAATTCCAACTGAATCTAC\n".
                  "TACTTCTGTTGAATTTGAACTGCCAGTTTTAGCTAGGTCTCTAGAACGTGCAAGTCTTCTTACGCATGAA\n".
                  "AGTCGAAATAACGTAGTAAAAATGACTCTTGATGTTCAAAATCAATTAGTTAAGCTTCAAGGTGATTCTC\n".
                  "CAGAGATAGGTAATGTGGAAGAAGAAATTAGTTTTAAGAAGCTTGAAGGCGAAGGTTTGACAATTTCATT\n".
                  "TAACCCTGACTATTTAAGAGAGGCTTTGCGCGCATCAATTACTGATTCCATTATTATGAACTTTACGCAG\n".
                  "CCATTGAGACCATTTACAGTTGTACCTGCAAAGCAGGATGTCAACTTTACGCAACTGATTACTCCAGTGA\n".
                  "GAACATTTTAAATGCCTGAAAAAGAGCTACTTTGGTAGCTCTTTTTTTCTTGTTGTACTTAAAAGCTAAA\n".
                  "TTTAGCCTGATTAAACGAATTTTTCTCTAGTTAAGTATAGTCATATTAAATTAGGAGAAAATGGCTTAAA\n".
                  "AAAGGCTATAATGCCTCATAAAATTGAATTTCACATAAAAAAACGGTATAATAACCAGGTATGACAAGGT\n".
                  "AAATGGGGTGATATTATCAAAAAATTCACAATTAAGGGAGAATATATTACGCTTACTCAATTTTTAAAAG\n".
                  "AAGAAAGCGTGATTTCTTCTGGTGGTCAAGCTAAATGGTATCTAAAAGAAAATCCTGTTAAGTTAAATGG\n".
                  "CGAACTTGAAGATCGTCGCGGAAAGAAAATACATGCAGGCGATGTATTGAACTTAGCTGGAGAAGAGTAT\n".
                  "GAATTTGTTTCAGAGTAGCCAATGTATCTCGCAAACTTTGAATTGAAAGACTTTAGGAATTTTAAAGAAT\n".
                  "TAAAAACAGATTTTGATCCTCATGTAAATATTTTTATTGGTCCGAATGCCCAAGGAAAAACTAATTTACT\n".
                  "TGAGGCAATTTATTTTTTAGCCCTAACACGATCCCATCGAACTAATAGCGATAAAGAATTAATTCGCTTT\n".
                  "GGCAGTAAATTTGCCGGTTTACAAGGTCGAGTTCATAAAAGCCAACTGCAAGTAGAGCTTAAATTACGTT\n".
                  "TAACAGCTAATGGAAAAAAGGCTTGGGTAAATCGCTTGGAGCAAAAGAAGCTTTCTGCTTATGTTGGGCA\n".
                  "AATGAATGCAATTTTGTTTTCTCCAGAAGATTTAGCCCTAGTTAAGGGAGCACCATCTGTTAGACGGAGA\n".
                  "TTTATGGATTTAGAGTTTGGGCAGATTAATTCAGAATATCTATATTTTTCAAGTCAGTATCGTCAGGTAC\n".
                  "TGCAACAGAGAAACAATTATTTAAAGCAACTAAGCATCAAAAAAGCAAATGATCAAGTTTTTTTAGATGT\n".
                  "CTTATCTGATCAGCTGGCAGGAATTGCAGCAGAGATAATATCGCGAAGAATAAAATACATCAAAAAACTC\n".
                  "AATTCCTACGCTAAAGCTGCTCATAGCGAGATTAGTGGTCAAGCTGAAAAATTGCAGATTTTTTACCGTC\n".
                  "CATCAGTTAAGGAAATTATCCCAGAGGATAATGTAGAGACCATTTATCGAAAAGTTATTACTAGTTATAA\n".
                  "GAAAAATCGTCCTAATGAAATTCGAAAGGGCACTACTCTAAGTGGCCCGCATCGAGATGACTTGGAGTTT\n".
                  "TTAATTAATGAAAAAAATGCCCATGATTTTGCATCTCAAGGACAGCAACGGACGATTTCTCTAAGCGTCA\n".
                  "AGTTAGCCGAGATTCAATTAGTACATGAATTGACGCAAGAATATCCAATTCTGTTATTAGATGATGTGAT\n".
                  "GAGTGAGTTAGATCATCGAAGACAGAGCCGTTTATTGAACTATATTCATGGCAAAACGCAAACATTTATA\n".
               "\>gi\|2\|gb\|2\| DUMMY_A\n".
                  "ATGACCGACAGCTCATTACCTAAAATTGGCCAAGAGTTTGGCGGTAAAGATCATACGACCGTCATGCATG\n".
                  "CTCATGAACGGATTAGTCAAAGTTTGACAACCGATCAAAACCTAAAAGATGCCATTTTGGATTTAAAGAA\n".
                  "GTCTTGCACAAAGAATACTGGGAAAATAACTTGGCGACTAAAGTCGTTGAAGGGGCCTATGAGTTTGCGG\n".
                  "CACAATGAAGTCTTAAAAATGGGTTGTGGATAAAGGTAAAAGTTATCCTAAGGATTTCTCACATTTTATC\n".
                  "CACAGGTGGAAAACTATGAGATTCTTGGCTTTCAATACTTATCCACAGGTTCCACAGGCCCTACTACTAT\n".
                  "TACTTTAAAATCTTTTAATATATATATAATATTCACCCTCGGAGGTAGCACCTATGAAATTTTCAATTAA\n".
                  "CCGCGCCCTATTTATCAAAACGCTCAATGATGTGCAACGGGCCATTTCAGCTAAAACAACGATTCCGATC\n".
                  "CAAGTATCTCTATAATATTTCGCAATACAAATCAATTCTGAAACAACGTAATCAATATCTTAAGCAACTC\n".
                  "CTAACAGGGCTTAAACTTGTTCTATCGCAAGAAGGCTTGACGCTAACCGGTAGTGATGCGGATATCTCAA\n".
                  "TTGAAGCAACATTGTCAGCCTCAGATGAAAAAAACCAATTAGTCGTCGAAGAAGTTGGCGAAATTGTGTT\n".
                  "ACCAGCACGTTTCTTTAGTGAAATCGTTAAAAAATTACCAGAAGATCAGATGCGAATTGCAGTTGGAAAT\n".
                  "AATTTCCAAACAACAATTCGCTCCGGTAAATCAGAATTTACAATCAACGGGTTAGACGCTAATAGTTATC\n";

   close WRITE;

   #generate dummy nt db
   if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
      chdir ${$args{ini_ref}}{blast_executables};
      `.formatdb -i ${$args{ini_ref}}{tempfile} -p F -o T -n ${$args{auto_ini_ref}}{work_dir}/dummy.nt`;
      chdir ${$args{auto_ini_ref}}{work_dir};
   } elsif (${$args{auto_ini_ref}}{blast_plus} == 1) {
      chdir ${$args{ini_ref}}{blast_plus_executables};
      `makeblastdb -dbtype nucl -in ${$args{ini_ref}}{tempfile} -parse_seqids -hash_index -out ${$args{auto_ini_ref}}{work_dir}/dummy.nt`;
      chdir ${$args{auto_ini_ref}}{work_dir};
   }

   #generate dummy DNA file for critica
   unlink ${$args{ini_ref}}{tempfile};
   open  WRITE, "+>".${$args{ini_ref}}{tempfile} or return (0, "Couldn't write to tempfile.");
   print WRITE "\>DUMMY_2\n".
                  "TGTATTTAGTCTTTCTTTATTGAGTATTTTTCAATAAAGTGGTGATTTTATGTGAGGAGGACACGCTACA\n".
                  "TGCCTGATCTTCAGGAACTTTGGAATTATCTGCAAGAAAAATTCCAAACCGATCTAACCACGGTTGGCTT\n".
                  "TAACGCTTGGATTAAAACCGCTAAGCCGCTTGCTTTTAGAGCTAACGAACTACTAATTGAAGTGCCTTCT\n".
                  "GTCTTGCACAAAGAATACTGGGAAAATAACTTGGCGACTAAAGTCGTTGAAGGGGCCTATGAGTTTGCGG\n".
                  "AAATCGAGCTGACCCCTATTTTTGTGTTACCTGGTGAATCTGACAACTTAACACCGCTTGAACCAGAAGA\n".
                  "AGAACACGTCCTAACGAAGGCTGAAACACCAACCTTTTTACGTGAAACGCATTTAAATAGTAAATATACA\n".
                  "TTTGATACGTTTGTCACTGGGAAAGGCAATCAAATGGCCCACGCTGCGGCGTTGGTGGTTTCCGAAGAAC\n".
                  "CCGGTGTCTTATACAATCCCCTCTTCTTATATGGTGGTGTTGGTTTAGGTAAAACCCATTTGATGCAAGC\n".
                  "AATTGGTCACCAATTACTCGAATCAAAACCCGAAACGAACGTTAAGTACGTGACGAGTGAGGCTTTTGCC\n".
                  "AACGATTTTATCAATTCGATTCAAACGAAGAATCAAGAGAAATTTAGACAAGAATACCGGAACGTCGACC\n".
                  "TTTTATTAGTGGATGATATTCAATTTTTCGCTGATAAAGAAGGAACGCAAGAAGAATTCTTCCATACATT\n".
                  "TAATGATCTCTATAATGATAAGAAACAAATTGTGTTGACCTCTGATCGCCTCCCAAATGAGATTCCCAAG\n".
                  "CTGCAAGAACGCTTGGTTTCTCGGTTTAAATGGGGGTTACCAGTCGATATCACGCCACCTGATCTTGAGA\n".
                  "CGCGAATTGCGATTTTAAGAAATAAAGCCGATGCGGAACATCTGGAAATTCCAGAAGATACGCTCAGTTA\n".
                  "TATTGCTGGGCAAATCGATTCCAACGTTCGGGAATTAGAAGGCGCACTCGTGCGGGTCCAAGCCTATGCA\n".
                  "ACGATGCAAAACGCTGAGATTACAACCAGTCTAGCGGCCGACGCGCTCAAAGGCTTGAAACTCAACGGTA\n".
                  "AATCTAGTCAGCTCTCGATTGCCAAGATTCAATCGGTGGTTGCTAAATATTACAGCTTAACCGTTGCTGA\n".
                  "TCTTAAAGGTCGCAAGCGGGTCAAAGAGATTGTGTTGCCGCGCCAAATCGCAATGTACCTTGCTCGCGAA\n".
                  "ATGACCGACAGCTCATTACCTAAAATTGGCCAAGAGTTTGGCGGTAAAGATCATACGACCGTCATGCATG\n".
                  "CTCATGAACGGATTAGTCAAAGTTTGACAACCGATCAAAACCTAAAAGATGCCATTTTGGATTTAAAGAA\n".
                  "CACAATGAAGTCTTAAAAATGGGTTGTGGATAAAGGTAAAAGTTATCCTAAGGATTTCTCACATTTTATC\n".
                  "CACAGGTGGAAAACTATGAGATTCTTGGCTTTCAATACTTATCCACAGGTTCCACAGGCCCTACTACTAT\n".
                  "TACTTTAAAATCTTTTAATATATATATAATATTCACCCTCGGAGGTAGCACCTATGAAATTTTCAATTAA\n".
                  "CCGCGCCCTATTTATCAAAACGCTCAATGATGTGCAACGGGCCATTTCAGCTAAAACAACGATTCCGATC\n".
                  "CTAACAGGGCTTAAACTTGTTCTATCGCAAGAAGGCTTGACGCTAACCGGTAGTGATGCGGATATCTCAA\n".
                  "TTGAAGCAACATTGTCAGCCTCAGATGAAAAAAACCAATTAGTCGTCGAAGAAGTTGGCGAAATTGTGTT\n".
                  "ACCAGCACGTTTCTTTAGTGAAATCGTTAAAAAATTACCAGAAGATCAGATGCGAATTGCAGTTGGAAAT\n".
                  "AATTTCCAAACAACAATTCGCTCCGGTAAATCAGAATTTACAATCAACGGGTTAGACGCTAATAGTTATC\n".
                  "CACACCTTCCTGAAATTGCAGGGCAAAATACATTTAGCATCAGCACTGATGTTTTGAAACAACTCATTCA\n".
                  "CCAAACCGTCATCGCTGTTTCTAATCAAGAAAGTCGCCCAATCTTAACTGGGGTCCACTTGGTCCTAGCC\n".
                  "AATGGCGAAATGTTAGCCGTTGCGACTGATAGCCATCGTTTAGCACAACGTAAATTACCAGTTGCGCTTC\n".
                  "CAAGTGATGCCAATTATGACGTCATCATTCCTGGTAAGAGTTTGGTAGAACTTTCAAGAACGCTTGCGGA\n".
                  "TGATGCCGAAGACGTTGAAATTCGGATTGCTGAAAACCAAGTCCTCTTTAAAGCTGGCAACTTAGCCTTT\n".
                  "TATTCACGCTTATTAGAAGGTAACTATCCTGATACAGCACGTTTAATTCCAACTTCATCAAGTACCCAAA\n".
                  "TCGAATTTAACGCGCCAGTCTTATTGGCCGCTATCGAACGGGCGTCATTGTTATCCCATGAAAGTCGGAA\n".
                  "TAACGTTGTCCGTTTAACACTCGACACAGAGGCTAACACAGCGATTATTTATGGTAATTCACCAGATGTG\n".
                  "GGGAATGTTGAAGAAGCCCTTCAATTTGAAAACCTCACGGGTGAATCCCTTGAAATTTCATTTAACCCTG\n".
                  "ACTACATGAAGGATGCCTTGCGTTCATTTGGCCAAACAACAATCGTGATTAACTTTACGGCTGCTTTACG\n".
                  "ACCATTTACGTTGGTTCCAAGTGAAGATCAAGAACACTATATTCAATTGATCACACCAGTCCGGACATTC\n".
                  "TAAAATTAAATAACAACTCAAAAACACCTCGATGAGAGGCGTTTTTTTGTGAACTAGAATTAAAATAAAT\n".
                  "TCGATTACAGCGATTTAAGATATTTCTAGTAGGAGAGGGGTAATTTTGTCAAAAACGCTGAGAAAGCATT\n".
                  "TTAACCATGTTTATTTGTTTTTAGTGCTAAAAAAGTCTATAATTACTAGGTATAACGGATCGTTATTTTT\n".
                  "CCAGTAATTAATATCGTTAAGAATGGTGGGTTAATATGACACAAGAAATTAAACTTGAGGCTGAATTCGT\n".
                  "CACATTAGGACAACTCCTAAAAGAAGCCGGCATTATTGAAACGGGTGGCAAAGCCAAATGGTTCCTAAGA\n".
                  "GAAAATACGGTTTTAGTCAACGGCAAACCTGATGATCGACGTGGCCGGAAGTTATATCCAGAGGATGTCA\n".
                  "TTGAAGTCCCAGATAACGGCCAATTCATCGTTAAACAACAAGAAAGATTGTAGGCAGCTATGTATCTAAG\n".
                  "CGAACTACAACTGAATCATTATCGTAATTATGAGTCAGTTGATGTGCACTTTTCACCGGATACGAACGTC\n".
                  "TTGATTGGTGAAAATGCGCAGGGGAAAACGAACTTACTAGAAGCAATCTACGTTTTAGCCTTGGCCCGCT\n".
                  "CACATCGCACGAATACCGATCGCGAATTAATTCAATGGCATGAAGATTTCGCTAAAATTACCGGCTTGGT\n".
                  "CCAACGTTCTGCAGGGAAAACCCCGTTAGAACTCGTTTTAAGCCAAAAAGGTAAAAAGGCGAAGGTCAAT\n".
                  "CATTTAGAACAGGCCAAGTTATCACAATATATTGGTCAGCTAAACGTCGTGTTATTTGCGCCGGAGGATT\n".
                  "TAAACATCGTTAAGGGTTCGCCCGCTGTGAGACGTCATTTCATTGATATGGAATTCGGGCAGATGAGTAG\n".
                  "CAAGTATCTCTATAATATTTCGCAATACAAATCAATTCTGAAACAACGTAATCAATATCTTAAGCAACTC\n".
                  "CAACGGCGCCAAGCCAAGGATCTTGTTTATTTAGGTGTTTTATCGGATCAATTAGCGGCTTACGGGGCGG\n".
                  "AAGTGACCGTGGCCCGGCGCCAGTTTCTCCAACAAATGGAAAAATGGGCCCAAAAATTGCATCAGGAGAT\n".
                  "TACGAAGGATCGGGAAGTCTTGACGTTTAAATATCAAAGTCAGATTCCAGAAGAACAATTGGATCAAAGT\n".
                  "GTCGAAGAACTCTATCAACAATTTCAAACCTTGTATGAGAAACAACAAATCCGTGAGGTTGAACAAGGCA\n";

   close WRITE;

   #generate dummy aa file
   open  WRITE, "+>".${$args{ini_ref}}{tempfile}.'_aa' or return (0, "Couldn't write to aa tempfile.");
   print WRITE "\>aa_dummy\n".
               "MRRTARFTAAELAKINHKQLEAVYNQFSTFARYRLFNSPIPDISLIDGNAIIPRIRRVGNLLNQITHTVN\n".
               "LTGTVSRKQVGAVKELVGELSKVLKKHFLERVLALCKTNKSDE\n";
   close WRITE;

}

sub delete_dummy {
   my %args = @_;
   unlink ${$args{ini_ref}}{tempfile};
   unlink ${$args{ini_ref}}{gm_output}.'/dummy.glcoord';

   unlink ${$args{ini_ref}}{gm_output}.'/dummy.detail';
   unlink ${$args{ini_ref}}{gm_output}.'/dummy.predict';

   unlink ${$args{ini_ref}}{tempfile};
   unlink ${$args{auto_ini_ref}}{work_dir}.'/dummy.nt.nhr';
   unlink ${$args{auto_ini_ref}}{work_dir}.'/dummy.nt.nin';
   unlink ${$args{auto_ini_ref}}{work_dir}.'/dummy.nt.nnd';
   unlink ${$args{auto_ini_ref}}{work_dir}.'/dummy.nt.nni';
   unlink ${$args{auto_ini_ref}}{work_dir}.'/dummy.nt.nsd';
   unlink ${$args{auto_ini_ref}}{work_dir}.'/dummy.nt.nsi';
   unlink ${$args{auto_ini_ref}}{work_dir}.'/dummy.nt.nsq';
   unlink ${$args{auto_ini_ref}}{work_dir}.'/dummy.nt.nog';
   unlink ${$args{ini_ref}}{critica}.'/EC.blast';
   unlink ${$args{ini_ref}}{critica}.'/EC.blast.pairs';
   unlink ${$args{ini_ref}}{critica}.'/EC.triplets';

   unlink ${$args{ini_ref}}{tempfile}.'_aa';

   #generate local array for selected PFam databases
   my @local_db = split/ ; /,${$args{auto_ini_ref}}{full_Pfam_db};
   foreach my $pfam_db (@local_db) {
      unlink $pfam_db.'_dummy';
   }

}



1;