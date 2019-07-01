#!/opt/ActivePerl-5.8/bin/perl

#recompile programs if run for first time or if relocated
#input arguments: main_window, progress_bar, auto_ini_ref, ini_ref,

package initialise::recompile;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&recompile);
use vars qw();

use initialise::read_me        qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);
use Basics::Recursive          qw(dircopy);
use Cwd;
use File::Find;
use File::Copy;
use File::Path;

#local variables
my (%args, $newdir);

sub recompile {
   my %args = (force_compile => 0,
               @_
              );

   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'Recompiling programs',
                   label        => 'Recompiling'
                  );

   #just to make sure, change access right to all archives:
   `chmod -R 777 ${$args{auto_ini_ref}}{work_dir}/lib/Archives/*`;

   #always test for zip
   {
      unless (-e ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives/zip') {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling programs',
                        label        => 'Recompiling Zip executables',
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_zip(main_window  => $args{main_window},
                                                 progress_bar => $args{progress_bar},
                                                 auto_ini_ref => $args{auto_ini_ref},
                                                 ini_ref      => $args{ini_ref},
                                                 directory    => ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives',
                                                );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not compile ZIP program.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }
   }

   #test for recompiling legacy Blast; compile anyway, essential programme
   if ((${$args{auto_ini_ref}}{legacy_blast}       == 1   ||
        ${$args{auto_ini_ref}}{blast_plus}         == 1)  &&
       ($args{force_compile}                       == 1   ||
        ${$args{auto_ini_ref}}{blast_selector}     == 1   ||
        ${$args{auto_ini_ref}}{COG_selector}       == 1   ||
        ${$args{auto_ini_ref}}{rrna_selector}      == 1   ||
        ${$args{auto_ini_ref}}{runintergenicblast} == 1)
       ) {
      my $blast_dir = ${$args{ini_ref}}{blast_executables};
      $blast_dir =~ s/\/bin$//;
      my ($status) = &test_recompile(main_window  => $args{main_window},
                                     progress_bar => $args{progress_bar},
                                     auto_ini_ref => $args{auto_ini_ref},
                                     ini_ref      => $args{ini_ref},
                                     directory    => $blast_dir,
                                     executable   => 'blastall'
                                    );
      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling programs',
                        label        => 'Recompiling Blast executables',
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_blast(main_window  => $args{main_window},
                                                   progress_bar => $args{progress_bar},
                                                   auto_ini_ref => $args{auto_ini_ref},
                                                   ini_ref      => $args{ini_ref},
                                                   directory    => $blast_dir,
                                                   );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not compile Blast programs.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }
   }

   #test for recompiling Blast plus; compile anyway, essential programme
   if ((${$args{auto_ini_ref}}{blast_plus} == 1     ||
        ${$args{auto_ini_ref}}{legacy_blast} == 1)  &&
       ($args{force_compile} == 1                   ||
        ${$args{auto_ini_ref}}{blast_selector} == 1 ||
        ${$args{auto_ini_ref}}{COG_selector} == 1   ||
        ${$args{auto_ini_ref}}{rrna_selector} == 1  ||
        ${$args{auto_ini_ref}}{runintergenicblast} == 1)
       ) {
      my $blast_dir = ${$args{ini_ref}}{blast_plus_executables};
      $blast_dir =~ s/\/bin$//;
      my ($status) = &test_recompile(main_window  => $args{main_window},
                                     progress_bar => $args{progress_bar},
                                     auto_ini_ref => $args{auto_ini_ref},
                                     ini_ref      => $args{ini_ref},
                                     directory    => $blast_dir,
                                     executable   => 'blastp'
                                    );
      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling programs',
                        label        => 'Recompiling Blast executables',
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_blast_plus (main_window  => $args{main_window},
                                                         progress_bar => $args{progress_bar},
                                                         auto_ini_ref => $args{auto_ini_ref},
                                                         ini_ref      => $args{ini_ref},
                                                         directory    => $blast_dir,
                                                         );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not compile Blast programs.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }
   }

   #test for recompiling fixed COG databases
   if ($args{force_compile} == 1 || ${$args{auto_ini_ref}}{COG_selector} == 1) {
      # set parent COG directory, archive already contains folder
      my $parent_COG_dir = ${$args{ini_ref}}{COG_db_path};
      #abbreviate COG_db_path if required
      if ($parent_COG_dir =~ m/(arCOG|arCOG2014|COG2003|COG2008|COG2014|POG2013)/) {
         $parent_COG_dir =~ s/\/\w+$//;
      }
      #standard COG
      my ($status) = &test_recompile(main_window  => $args{main_window},
                                     progress_bar => $args{progress_bar},
                                     auto_ini_ref => $args{auto_ini_ref},
                                     ini_ref      => $args{ini_ref},
                                     directory    => $parent_COG_dir.'/COG2003',
                                     type         => 'COG2003',
                                     executable   => 'COG.phr'
                                    );
      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling databases',
                        label        => 'Extracting COG database files',
                        progress     => 1,
                       );

         my ($recompile_status) = &recompile_COG(main_window  => $args{main_window},
                                                 progress_bar => $args{progress_bar},
                                                 auto_ini_ref => $args{auto_ini_ref},
                                                 ini_ref      => $args{ini_ref},
                                                 directory    => ${$args{ini_ref}}{COG_db_path},
                                                 type         => 'COG2003'
                                                );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not extract standard COG database.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }

      #COG2008
      ($status) = &test_recompile(main_window  => $args{main_window},
                                  progress_bar => $args{progress_bar},
                                  auto_ini_ref => $args{auto_ini_ref},
                                  ini_ref      => $args{ini_ref},
                                  directory    => $parent_COG_dir.'/COG2008',
                                  type         => 'COG2008',
                                  executable   => 'COG2008.phr'
                                 );
      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling databases',
                        label        => 'Extracting COG2008 database files',
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_COG(main_window  => $args{main_window},
                                                 progress_bar => $args{progress_bar},
                                                 auto_ini_ref => $args{auto_ini_ref},
                                                 ini_ref      => $args{ini_ref},
                                                 directory    => ${$args{ini_ref}}{COG_db_path},
                                                 type         => 'COG2008'
                                                );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not extract COG2008 database.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }

      #COG2014
      ($status) = &test_recompile(main_window  => $args{main_window},
                                  progress_bar => $args{progress_bar},
                                  auto_ini_ref => $args{auto_ini_ref},
                                  ini_ref      => $args{ini_ref},
                                  directory    => $parent_COG_dir.'/COG2014',
                                  type         => 'COG2014',
                                  executable   => 'COG2014.phr'
                                 );
      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling databases',
                        label        => 'Extracting COG2014 database files',
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_COG(main_window  => $args{main_window},
                                                 progress_bar => $args{progress_bar},
                                                 auto_ini_ref => $args{auto_ini_ref},
                                                 ini_ref      => $args{ini_ref},
                                                 directory    => ${$args{ini_ref}}{COG_db_path},
                                                 type         => 'COG2014'
                                                );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not extract COG2014 database.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }

      #arCOG
      ($status) = &test_recompile(main_window  => $args{main_window},
                                  progress_bar => $args{progress_bar},
                                  auto_ini_ref => $args{auto_ini_ref},
                                  ini_ref      => $args{ini_ref},
                                  directory    => $parent_COG_dir.'/arCOG',
                                  type         => 'arCOG',
                                  executable   => 'arCOG.phr'
                                 );
      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling databases',
                        label        => 'Extracting arCOG database files',
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_COG(main_window  => $args{main_window},
                                                 progress_bar => $args{progress_bar},
                                                 auto_ini_ref => $args{auto_ini_ref},
                                                 ini_ref      => $args{ini_ref},
                                                 directory    => ${$args{ini_ref}}{COG_db_path},
                                                 type         => 'arCOG'
                                                );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not extract arCOG database.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }

      #arCOG2014
      ($status) = &test_recompile(main_window  => $args{main_window},
                                  progress_bar => $args{progress_bar},
                                  auto_ini_ref => $args{auto_ini_ref},
                                  ini_ref      => $args{ini_ref},
                                  directory    => $parent_COG_dir.'/arCOG2014',
                                  type         => 'arCOG2014',
                                  executable   => 'arCOG2014.phr'
                                 );
      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling databases',
                        label        => 'Extracting arCOG2014 database files',
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_COG(main_window  => $args{main_window},
                                                 progress_bar => $args{progress_bar},
                                                 auto_ini_ref => $args{auto_ini_ref},
                                                 ini_ref      => $args{ini_ref},
                                                 directory    => ${$args{ini_ref}}{COG_db_path},
                                                 type         => 'arCOG2014'
                                                );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not extract arCOG2014 database.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }

      #POG2013
      ($status) = &test_recompile(main_window  => $args{main_window},
                                  progress_bar => $args{progress_bar},
                                  auto_ini_ref => $args{auto_ini_ref},
                                  ini_ref      => $args{ini_ref},
                                  directory    => $parent_COG_dir.'/POG2013',
                                  type         => 'POG2013',
                                  executable   => 'POGseqs_HighVQ.phr'
                                 );
      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling databases',
                        label        => 'Extracting POG2013 database files',
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_COG(main_window  => $args{main_window},
                                                 progress_bar => $args{progress_bar},
                                                 auto_ini_ref => $args{auto_ini_ref},
                                                 ini_ref      => $args{ini_ref},
                                                 directory    => ${$args{ini_ref}}{COG_db_path},
                                                 type         => 'POG2013'
                                                );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not extract POG2013 database.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }
   }

   #test for recompiling hmmer
   if ($args{force_compile} == 1 || ${$args{auto_ini_ref}}{Pfam_selector} == 1 || ${$args{auto_ini_ref}}{TIGRfam_selector} == 1) {
      #hmmer2
      my ($status) = &test_recompile(main_window  => $args{main_window},
                                     progress_bar => $args{progress_bar},
                                     auto_ini_ref => $args{auto_ini_ref},
                                     ini_ref      => $args{ini_ref},
                                     directory    => ${$args{ini_ref}}{pfam_executable},
                                     executable   => 'hmmpfam'
                                    );

      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling programs',
                        label        => 'Recompiling hmmer executables',
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_hmmer(main_window  => $args{main_window},
                                                   progress_bar => $args{progress_bar},
                                                   auto_ini_ref => $args{auto_ini_ref},
                                                   ini_ref      => $args{ini_ref},
                                                   directory    => ${$args{ini_ref}}{pfam_executable},
                                                   );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not compile HMMER.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }

      #hmmer3
      ($status) = &test_recompile(main_window  => $args{main_window},
                                  progress_bar => $args{progress_bar},
                                  auto_ini_ref => $args{auto_ini_ref},
                                  ini_ref      => $args{ini_ref},
                                  directory    => ${$args{ini_ref}}{pfam3_executable},
                                  executable   => 'hmmscan'
                                 );

      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling programs',
                        label        => 'Recompiling hmmer3 executables',
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_hmmer3(main_window  => $args{main_window},
                                                    progress_bar => $args{progress_bar},
                                                    auto_ini_ref => $args{auto_ini_ref},
                                                    ini_ref      => $args{ini_ref},
                                                    directory    => ${$args{ini_ref}}{pfam3_executable},
                                                    );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not compile HMMER3.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }

      #test for presence of relnote.txt
      #unless (-e ${$args{ini_ref}}{pfam_db_path}.'/relnotes.txt') {
      #   ${$args{main_window}}->messageBox(-title   => 'Error',
      #                                     -message => "Could not find release notes file \'relnotes.txt\' in directory ${$args{ini_ref}}{pfam_db_path} \.\nAborting",
      #                                     -type    => 'OK',
      #                                     -icon    => 'info');
      #   &hide_pbar_2;
      #   return (0);
      #}

   }

   #test for recompiling glimmer2 + 3
   if ($args{force_compile} == 1 ||
       ${$args{auto_ini_ref}}{internal_gm} == 1 &&
       (${$args{auto_ini_ref}}{use_glimmer2} == 1 || ${$args{auto_ini_ref}}{use_glimmer3} == 1)) {
      my (%glimmer);
      $glimmer{'glimmer2'} = ${$args{ini_ref}}{glprog2};
      $glimmer{'glimmer3'} = ${$args{ini_ref}}{glprog3};
      while (my ($glimmer_exe, $glimmer_dir) = (each %glimmer)) {
         my ($status) = &test_recompile(main_window  => $args{main_window},
                                        progress_bar => $args{progress_bar},
                                        auto_ini_ref => $args{auto_ini_ref},
                                        ini_ref      => $args{ini_ref},
                                        directory    => $glimmer_dir,
                                        executable   => $glimmer_exe
                                       );

         if ($status == 0) {
            &show_pbar_2;
            &update_pbar_2(title        => 'Recompiling programs',
                           label        => "Recompiling $glimmer_exe executables",
                           progress     => 1,
                          );
            my ($recompile_status) = &recompile_glimmer(main_window  => $args{main_window},
                                                        progress_bar => $args{progress_bar},
                                                        auto_ini_ref => $args{auto_ini_ref},
                                                        ini_ref      => $args{ini_ref},
                                                        directory    => $glimmer_dir,
                                                        version      => $glimmer_exe
                                                        );
            if ($recompile_status == 0) {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Could not compile $glimmer_exe.\nAborting",
                                                 -type    => 'OK',
                                                 -icon    => 'info');
               &hide_pbar_2;
               return (0);
            }
            &hide_pbar_2;
         }
      }
   }

   #test for copying prodigal
   if ($args{force_compile} == 1 ||
       ${$args{auto_ini_ref}}{internal_gm} == 1 && ${$args{auto_ini_ref}}{use_prodigal} == 1 ) {
      my $recompile_prodigal = 0;
      my $executable;
      #test if prodigal directory exists
      if (!-d  ${$args{ini_ref}}{prodigal}) {$recompile_prodigal = 1};
      #test if there is a prodigal exe
      if ($recompile_prodigal == 0) {
         find (sub {$executable = $File::Find::name if (-f && $_ =~ m/^prodigal/i) }, ${$args{ini_ref}}{prodigal});
         unless (defined $executable) {$recompile_prodigal = 1};
      }
      if ($recompile_prodigal == 1) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Copying programs',
                        label        => "Copying Prodigal executable",
                        progress     => 1,
                       );

         if (!-d  ${$args{ini_ref}}{prodigal}) {
            mkdir ${$args{ini_ref}}{prodigal};
         }

         #find prodigal executable in archive
         find (sub {$executable = $File::Find::name if (-f && $_ =~ m/^prodigal/i) }, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives');

         #copy file
         copy ($executable, ${$args{ini_ref}}{prodigal});

         #change access rights just to be sure
         find (sub {$executable = $File::Find::name if (-f && $_ =~ m/^prodigal/i) }, ${$args{ini_ref}}{prodigal});
         chmod 0777, $executable;

         #set compiler flag
         open WRITE, '+>'.${$args{ini_ref}}{prodigal}.'/recompiled' or do {return (0)};
         print WRITE ${$args{ini_ref}}{prodigal};
         close WRITE;

         &hide_pbar_2;
      }
   }

   #test for recompiling critica
   if ($args{force_compile} == 1 || ${$args{auto_ini_ref}}{internal_gm} == 1 && ${$args{auto_ini_ref}}{use_critica} == 1) {
      my ($status) = &test_recompile(main_window  => $args{main_window},
                                     progress_bar => $args{progress_bar},
                                     auto_ini_ref => $args{auto_ini_ref},
                                     ini_ref      => $args{ini_ref},
                                     directory    => ${$args{ini_ref}}{critica},
                                     executable   => 'critica'
                                    );
      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling programs',
                        label        => "Recompiling Critica executables",
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_critica(main_window  => $args{main_window},
                                                     progress_bar => $args{progress_bar},
                                                     auto_ini_ref => $args{auto_ini_ref},
                                                     ini_ref      => $args{ini_ref},
                                                     directory    => ${$args{ini_ref}}{critica},
                                                     );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not compile Critica.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         #add critica paths to environment
         if ($ENV{'PATH'} !~ /${$args{auto_ini_ref}}{critica_bin}/) {
            $ENV{'PATH'} = ${$args{auto_ini_ref}}{critica_bin}.':'.${$args{auto_ini_ref}}{critica_scripts}.':'.$ENV{'PATH'};
         }
         &hide_pbar_2;
      }
   }

   #test for recompiling SignalP3 or 4, but test is only if v3 or v4 is present
   if ($args{force_compile} == 1 || ${$args{auto_ini_ref}}{signalp_selector} == 1) {
      my ($status) = &test_recompile(main_window  => $args{main_window},
                                     progress_bar => $args{progress_bar},
                                     auto_ini_ref => $args{auto_ini_ref},
                                     ini_ref      => $args{ini_ref},
                                     directory    => ${$args{ini_ref}}{signalp_dir},
                                     executable   => 'signalp'
                                    );

      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling programs',
                        label        => "Recompiling SignalP executables",
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_signalp(main_window  => $args{main_window},
                                                     progress_bar => $args{progress_bar},
                                                     auto_ini_ref => $args{auto_ini_ref},
                                                     ini_ref      => $args{ini_ref},
                                                     directory    => ${$args{ini_ref}}{signalp_dir},
                                                     );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not compile SignalP.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }
   }

   #test for recompiling TMHMM
   if ($args{force_compile} == 1 || ${$args{auto_ini_ref}}{tmhmm_selector} == 1) {
      my ($status) = &test_recompile(main_window  => $args{main_window},
                                     progress_bar => $args{progress_bar},
                                     auto_ini_ref => $args{auto_ini_ref},
                                     ini_ref      => $args{ini_ref},
                                     directory    => ${$args{ini_ref}}{tmhmm_dir},
                                     executable   => 'decodeanhmm'
                                    );
      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling programs',
                        label        => "Recompiling TMHMM2 executables",
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_tmhmm(main_window  => $args{main_window},
                                                   progress_bar => $args{progress_bar},
                                                   auto_ini_ref => $args{auto_ini_ref},
                                                   ini_ref      => $args{ini_ref},
                                                   directory    => ${$args{ini_ref}}{tmhmm_dir},
                                                  );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not compile TMHMM2.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }
   }

   #test for recompiling tRNAscan
   if ($args{force_compile} == 1 || ${$args{auto_ini_ref}}{trna_selector} == 1) {
      my ($status) = &test_recompile(main_window  => $args{main_window},
                                     progress_bar => $args{progress_bar},
                                     auto_ini_ref => $args{auto_ini_ref},
                                     ini_ref      => $args{ini_ref},
                                     directory    => ${$args{ini_ref}}{trnascan_dir},
                                     executable   => 'tRNAscan-SE'
                                    );
      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling programs',
                        label        => "Recompiling tRNAscan executables",
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_trnascan(main_window  => $args{main_window},
                                                      progress_bar => $args{progress_bar},
                                                      auto_ini_ref => $args{auto_ini_ref},
                                                      ini_ref      => $args{ini_ref},
                                                      directory    => ${$args{ini_ref}}{trnascan_dir},
                                                     );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not compile tRNAscan.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         #add tRNAscan paths to environment
         if ($ENV{'PATH'} !~ /${$args{ini_ref}}{trnascan_dir}/) {
            $ENV{'PATH'} = ${$args{ini_ref}}{trnascan_dir}.':'.$ENV{'PATH'};
         }
         &hide_pbar_2;
      }
   }

   #test for recompiling Transterm
   if ($args{force_compile} == 1 || ${$args{auto_ini_ref}}{terminator_selector} == 1) {
      my ($status) = &test_recompile(main_window  => $args{main_window},
                                     progress_bar => $args{progress_bar},
                                     auto_ini_ref => $args{auto_ini_ref},
                                     ini_ref      => $args{ini_ref},
                                     directory    => ${$args{ini_ref}}{transterm_dir},
                                     executable   => 'transterm'
                                    );
      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling programs',
                        label        => "Recompiling Transterm executables",
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_transterm(main_window  => $args{main_window},
                                                       progress_bar => $args{progress_bar},
                                                       auto_ini_ref => $args{auto_ini_ref},
                                                       ini_ref      => $args{ini_ref},
                                                       directory    => ${$args{ini_ref}}{transterm_dir},
                                                      );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not compile Transterm.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }
   }

   #test for recompiling rFAM/infernal
   if ($args{force_compile} == 1 || ${$args{auto_ini_ref}}{ncrna_selector} == 1) {
      my ($status) = &test_recompile(main_window  => $args{main_window},
                                     progress_bar => $args{progress_bar},
                                     auto_ini_ref => $args{auto_ini_ref},
                                     ini_ref      => $args{ini_ref},
                                     directory    => ${$args{ini_ref}}{infernal_dir},
                                     executable   => 'cmsearch'
                                    );
      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling programs',
                        label        => "Recompiling Infernal executables",
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_infernal(main_window  => $args{main_window},
                                                      progress_bar => $args{progress_bar},
                                                      auto_ini_ref => $args{auto_ini_ref},
                                                      ini_ref      => $args{ini_ref},
                                                      directory    => ${$args{ini_ref}}{infernal_dir},
                                                     );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not compile Infernal.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }
   }

   #test for recompiling rRNA databases
   if ($args{force_compile} == 1 || ${$args{auto_ini_ref}}{rrna_selector} == 1 ) {
      my $status;
      my @db = qw(5S_rRNA.nin 16S_rRNA.nin 23S_rRNA.nin);

      #test for all 3 databases
      foreach my $set (@db) {
         $status = &test_recompile(main_window  => $args{main_window},
                                   progress_bar => $args{progress_bar},
                                   auto_ini_ref => $args{auto_ini_ref},
                                   ini_ref      => $args{ini_ref},
                                   directory    => ${$args{ini_ref}}{rrna_db_path},
                                   executable   => $set
                                  );
         last unless $status;
      }

      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling rRNA database',
                        label        => 'Recompiling rRNA databases',
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_rRNA(main_window  => $args{main_window},
                                                  progress_bar => $args{progress_bar},
                                                  auto_ini_ref => $args{auto_ini_ref},
                                                  ini_ref      => $args{ini_ref},
                                                  directory    => ${$args{ini_ref}}{rrna_db_path},
                                                  );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not creat rRNA databases.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }
   }

   #test for recompiling CRISPR CRT
   if ($args{force_compile} == 1 || ${$args{auto_ini_ref}}{CRISPR_selector} == 1) {
      my ($status) = &test_recompile(main_window  => $args{main_window},
                                     progress_bar => $args{progress_bar},
                                     auto_ini_ref => $args{auto_ini_ref},
                                     ini_ref      => $args{ini_ref},
                                     directory    => ${$args{ini_ref}}{CRISPR_dir},
                                     executable   => 'CRT'
                                    );
      if ($status == 0) {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling programs',
                        label        => "Recompiling CRT executable",
                        progress     => 1,
                       );
         my ($recompile_status) = &recompile_CRT(main_window  => $args{main_window},
                                                 progress_bar => $args{progress_bar},
                                                 auto_ini_ref => $args{auto_ini_ref},
                                                 ini_ref      => $args{ini_ref},
                                                 directory    => ${$args{ini_ref}}{CRISPR_dir},
                                                );
         if ($recompile_status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not compile CRT.\nAborting",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0);
         }
         &hide_pbar_2;
      }
   }

   #test for recompiling COGnitor
   # if ($args{force_compile} == 1 || $args{test_cognitor} == 1) {
      # my %cognitor = ();
      # $cognitor{COGcognitor}  = ${$args{ini_ref}}{COGcognitor_path};
      # $cognitor{COGlse}       = ${$args{ini_ref}}{COGcognitor_path};
      # $cognitor{COGmakehash}  = ${$args{ini_ref}}{COGcognitor_path};
      # $cognitor{COGreadblast} = ${$args{ini_ref}}{COGcognitor_path};
      # $cognitor{COGtriangles} = ${$args{ini_ref}}{COGcognitor_path};

      # while (my ($cog_exe, $cog_dir) = (each %cognitor)) {
         # my ($status) = &test_recompile(main_window  => $args{main_window},
                                        # progress_bar => $args{progress_bar},
                                        # auto_ini_ref => $args{auto_ini_ref},
                                        # ini_ref      => $args{ini_ref},
                                        # directory    => $cog_dir,
                                        # executable   => $cog_exe
                                       # );
         # if ($status == 0) {
            # &show_pbar_2;
            # &update_pbar_2(title        => 'Recompiling programs',
                           # label        => "Recompiling COGnitor executables",
                           # progress     => 1,
                          # );
            # my ($recompile_status) = &recompile_COGnitor(main_window  => $args{main_window},
                                                         # progress_bar => $args{progress_bar},
                                                         # auto_ini_ref => $args{auto_ini_ref},
                                                         # ini_ref      => $args{ini_ref},
                                                         # directory    => $cog_dir,
                                                        # );
            # if ($recompile_status == 0) {
               # ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 # -message => "Could not compile COGnitor.\nAborting",
                                                 # -type    => 'OK',
                                                 # -icon    => 'info');
               # &hide_pbar_2;
               # return (0);
            # }
            # &hide_pbar_2;
         # }
      # }
   # }

   #test if VecScreen Database is present in default Blast_db dir
   if ($args{force_compile} == 1 || ${$args{auto_ini_ref}}{vector_screen_selector} == 1) {
      chdir ${$args{ini_ref}}{blast_plus_executables};

      #test for UniVec
      unless (-e ${$args{ini_ref}}{vector_db_path}.'/UniVec.nhr') {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling UniVec',
                        label        => "Recompiling UniVec database",
                        progress     => 1,
                       );
         copy (${$args{auto_ini_ref}}{work_dir}.'/lib/Archives/UniVec', ${$args{ini_ref}}{vector_db_path}.'/UniVec');
         #format db
         `./makeblastdb -dbtype nucl -in ${$args{ini_ref}}{vector_db_path}/UniVec -parse_seqids -hash_index -out ${$args{ini_ref}}{vector_db_path}/UniVec`;
         unlink ${$args{ini_ref}}{vector_db_path}.'/UniVec';
         &hide_pbar_2;
      }
      #test for UniVec_core
      unless (-e ${$args{ini_ref}}{vector_db_path}.'/UniVec_Core.nhr') {
         &show_pbar_2;
         &update_pbar_2(title        => 'Recompiling UniVec_Core',
                        label        => "Recompiling UniVec_Core database",
                        progress     => 1,
                       );
         copy (${$args{auto_ini_ref}}{work_dir}.'/lib/Archives/UniVec_Core', ${$args{ini_ref}}{vector_db_path}.'/UniVec_Core');
         #format db
         `./makeblastdb -dbtype nucl -in ${$args{ini_ref}}{vector_db_path}/UniVec_Core -parse_seqids -hash_index -out ${$args{ini_ref}}{vector_db_path}/UniVec_Core`;
         unlink ${$args{ini_ref}}{vector_db_path}.'/UniVec_Core';
         &hide_pbar_2;
      }
      chdir ${$args{auto_ini_ref}}{work_dir};
      unless (-e ${$args{ini_ref}}{vector_db_path}.'/UniVec_Core.nhr' &&
              -e ${$args{ini_ref}}{vector_db_path}.'/UniVec.nhr' ) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Failed to compile UniVec databases in directory ${$args{ini_ref}}{vector_db_path}",
                                           -icon    => 'error',
                                           -type    => 'ok');
      }

   }

   return (1);
}

sub test_recompile {
   my %args = (type => '',
               @_);
   my ($test_dir, $executable);

   #recompile on first GAMOLA use
   unless (-e $args{directory}.'/'.$args{type}.'recompiled') {
      return (0); #recompile
   }

   #test if relocated
   open READ, $args{directory}.'/'.$args{type}.'recompiled';
   while (<READ>) {$test_dir = $_};
   close READ;

   $test_dir =~ s/[\n\r]//gs;
   #if ($args{type} eq '') {
      if ($test_dir ne $args{directory}) {
         return (0); #recompile
      }
   #} else {
   #   if ($test_dir ne $args{directory}.'/'.$args{type} ) {
   #      return (0); #recompile
   #   }
   #}
   #test if executable exists
   if ($args{type} eq '') {
      unless (-e $args{directory}.'/'.$args{executable}) {
         #try finding executable if not in standard path
         find (sub {$executable = $File::Find::name if (-f && $_ =~ m/^$args{executable}/) }, $args{directory});
         #if exe is found then assume everything is OK
         if (defined $executable && $executable =~ /\w+/) {
            return (1); #everything's OK
         }
         return (0); #recompile
      }
   } else {
      unless (-e $args{directory}.'/'.$args{type}.'/'.$args{executable}) {
         #try finding executable if not in standard path
         find (sub {$executable = $File::Find::name if (-f && $_ =~ m/^$args{executable}/) }, $args{directory});
         #if exe is found then assume everything is OK
         if (defined $executable && $executable =~ /\w+/) {
            return (1); #everything's OK
         }
         return (0); #recompile
      }
   }
   return (1); #everything's OK
}


sub recompile_zip {
   my %args = @_;
   my ($zip_archive, $zip_dir, $short_name, $version);

   #find current zip file in Archive directory
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Finding Zip archive',
                  progress     => 25,
                 );

   #does directory and Makefile already exists?
   unless (-e $args{directory}.'/zipsrc/unix/Makefile') {
      find (sub {$zip_archive = $File::Find::name if (/^zip.*?\.(zip|Z|tar|gz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );

      unless (defined $zip_archive && $zip_archive =~ /\w+/) {
         return (0);
      }

      $short_name = $zip_archive;
      $short_name =~ s/.+\/(.+)$/$1/;

      unless (-d $args{directory}.'/zipsrc') {
         mkdir ($args{directory}.'/zipsrc', 0777)  or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot create Zip source directory",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return (0);
         };
      }
      chdir $args{directory}.'/zipsrc';

      &update_pbar_2(title        => 'Recompiling programs',
                     label        => 'Unzipping Zip archive',
                     progress     => 55,
                    );

      `unzip $zip_archive`;
   }

   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Compiling Zip',
                  progress     => 75,
                 );
   chdir $args{directory}.'/zipsrc';
   `make -f unix/Makefile generic`;

   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Moving Zip executables',
                  progress     => 95,
                 );
   move ($args{directory}.'/zipsrc/zip'     , $args{directory}.'/zip');
   #move ($args{directory}.'/zipsrc/zipcloak', $args{directory}.'/zipcloak');
   #move ($args{directory}.'/zipsrc/zipnote' , $args{directory}.'/zipnote');
   #move ($args{directory}.'/zipsrc/zipsplit', $args{directory}.'/zipsplit');

   #make sure everything is accessible and executable
   chdir $args{directory};
   `chmod -R 777 *`;

   #cleanup
   chdir ${$args{auto_ini_ref}}{work_dir};
   cleanup ($args{directory}.'/zipsrc');

   return (1);
}

sub recompile_blast {
   my %args = @_;
   my ($blast_archive, $blast_dir, $short_name, $version);
   my @blast_dirs = qw(bin data doc);

   #determine bit version of Linux
   my $li_version = `uname -a`;
   if ($li_version =~ / x86_64 /) {
      $li_version = 'x86_64';
   } else {
      $li_version = 'ia32';
   }

   #find current Blast zip file in Archive directory
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Finding Blast archive',
                  progress     => 25,
                 );

   #find archive for legacy Blast (32bit only for Blast 2.2.18)
   find (sub {$blast_archive = $File::Find::name if (/^blast\-.*?\.(Z|tar|gz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );
   unless (defined $blast_archive && $blast_archive =~ /\w+/) {
      return (0);
   }
   $short_name = $blast_archive;
   $short_name =~ s/.+\/(.+)$/$1/;
   $version    = $short_name;
   $version =~ s/^(blast\-[^\-]*?)\-.+/$1/i;

   #check if directory exists
   unless (-d $args{directory}) {
      mkdir $args{directory};
   }

   #copy Blast archive to Blast exe directory
   copy($blast_archive , $args{directory}.'/'.$short_name);

   #uncompress file
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Uncompressing Blast archive',
                  progress     => 50,
                 );
   chdir $args{directory};
   if ($short_name =~ m/\.tar.Z$/i) {
      `gzip -d $short_name`;
      $short_name =~ s/\.Z$//;
      `tar -xf $short_name`;
   } elsif ($short_name =~ m/\.tar.(gz|zip)$/i) {
      `tar -xzf $short_name`;
   } elsif ($short_name =~ m/\.tar$/i) {
      `tar -xf $short_name`;
   }

   #copy directories to proper Blast path
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Moving Blast software and cleanup',
                  progress     => 95,
                 );
   find (sub {$blast_dir = $File::Find::dir if ($_ =~ m/version/i) }, $args{directory});
   foreach my $dir (@blast_dirs) {
      dircopy ($blast_dir.'/'.$dir , $args{directory}.'/'.$dir);
   }

   #set path to blast executable:
   ${$args{ini_ref}}{blast_executables} = $args{directory}.'/bin';

   #cleanup
   cleanup ($blast_dir);
   unlink $args{directory}.'/'.$short_name;


   #make sure everything is accessible and executable
   `chmod -R 777 *`;

   #set compiler flag
   open WRITE, '+>'.$args{directory}.'/recompiled';
   print WRITE $args{directory};
   close WRITE;
   chdir ${$args{auto_ini_ref}}{work_dir};
   return (1);
}

sub recompile_blast_plus {
   my %args = @_;
   my ($blast_archive, $blast_dir, $short_name, $version);
   my @blast_dirs = qw(bin doc);

   #determine bit version of Linux
   my $li_version = `uname -a`;
   if ($li_version =~ / x86_64 /) {
      $li_version = 'x86_64';
   } else {
      $li_version = 'ia32';
   }

   #find current Blast zip file in Archive directory
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Finding Blast Plus archive',
                  progress     => 25,
                 );

   #find archive for legacy Blast (32bit only for Blast 2.2.18)
   #or find Blast Plus archive for 32 or 64 bit systems
   if ($li_version eq 'ia32') {
      find (sub {$blast_archive = $File::Find::name if (/ncbi\-blast\-.*?ia32.*?\.(Z|tar|gz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );
   } elsif ($li_version eq 'x86_64') {
      find (sub {$blast_archive = $File::Find::name if (/ncbi\-blast\-.*?x64.*?\.(Z|tar|gz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );
   } else {
      return (0);
   }
   unless (defined $blast_archive && $blast_archive =~ /\w+/) {
      return (0);
   }
   $short_name = $blast_archive;
   $short_name =~ s/.+\/(.+)$/$1/;
   $version    = $short_name;
   $version    =~ s/^(ncbi\-blast\-[^\-\+]*?)\-.+/$1/i;

   #copy Blast archive to Blast exe directory
   unless (-d $args{directory}) { mkdir $args{directory} };
   copy($blast_archive , $args{directory}.'/'.$short_name);

   #uncompress file
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Uncompressing Blast archive',
                  progress     => 50,
                 );
   chdir $args{directory};
   if ($short_name =~ m/\.tar.Z$/i) {
      `gzip -d $short_name`;
      $short_name =~ s/\.Z$//;
      `tar -xf $short_name`;
   } elsif ($short_name =~ m/\.tar.(gz|zip)$/i) {
      `tar -xzf $short_name`;
   } elsif ($short_name =~ m/\.tar$/i) {
      `tar -xf $short_name`;
   }

   #copy directories to proper Blast path
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Moving Blast software and cleanup',
                  progress     => 95,
                 );
   find (sub {$blast_dir = $File::Find::dir if ($_ =~ m/ncbi_package_info/i) }, $args{directory});
   foreach my $dir (@blast_dirs) {
      dircopy ($blast_dir.'/'.$dir , $args{directory}.'/'.$dir);
   }

   #set path to blast executable:
   ${$args{ini_ref}}{blast_plus_executables} = $args{directory}.'/bin';

   #cleanup
   cleanup ($blast_dir);
   unlink $args{directory}.'/'.$short_name;


   #make sure everything is accessible and executable
   `chmod -R 777 *`;

   #set compiler flag
   open WRITE, '+>'.$args{directory}.'/recompiled';
   print WRITE $args{directory};
   close WRITE;

   return (1);
}

sub recompile_COG {
   my %args = @_;
   my ($cog_archive, $cog_type_dir);

   #find current COG zip file in Archive directory
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Finding $args{type} COG archive",
                  progress     => 15,
                 );

   $cog_archive = ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives/'. $args{type} .'.tar.gz';
   unless (-e $cog_archive) {
      return (0);
   }

   #test if directory exists
   $cog_type_dir = $args{directory}.'/'.$args{type};
   unless (-d $cog_type_dir) {
      mkdir ($cog_type_dir, 0777) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create $cog_type_dir directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
   }

   #copy COG archive to correct cog directory
   copy($cog_archive , $args{directory});

   #uncompress file
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Uncompressing $args{type} COG archive",
                  progress     => 50,
                 );
   chdir $args{directory};
   `tar -xzf $args{directory}/$args{type}.tar.gz`;

   #make sure everything is accessible and executable
   `chmod -R 777 $args{directory}/$args{type}/*`;

   #cleanup source distro
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Cleanup',
                  progress     => 95,
                 );
   chdir ${$args{auto_ini_ref}}{work_dir};
   unlink $args{directory}.'/'.$args{type}.'.tar.gz';

   #set compiler flag
   open WRITE, '+>'.$cog_type_dir.'/'.$args{type}.'recompiled' or do {return (0)};
   print WRITE $cog_type_dir;
   close WRITE;

   return (1);
}

sub recompile_hmmer {
   my %args = @_;
   my ($hmmer_archive, @hmmer_archive, $hmmer_dir, $short_name, $version);

   #find current HMMER zip file in Archive directory
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Finding HMMER archive',
                  progress     => 15,
                 );

   find (sub {$hmmer_archive = $File::Find::name if (/hmmer\-2.*?\.(Z|tar|gz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );
   unless (defined $hmmer_archive && $hmmer_archive =~ /\w+/) {
      return (0);
   }

   $short_name = $hmmer_archive;
   $short_name =~ s/.+\/(.+)$/$1/;

   $version    = $short_name;
   $version    =~ s/^(hmmer\-.*?)\.(tar|gz|Z).*/$1/i;

   #test if directory exists
   unless (-d $args{directory}) {
      mkdir ($args{directory}, 0777) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create $args{directory} directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
   }

   #copy hmmer archive to hmmer exe directory
   copy($hmmer_archive , $args{directory}.'/'.$short_name);

   #uncompress file
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Uncompressing HMMER archive version $version",
                  progress     => 30,
                 );
   chdir $args{directory};
   if ($short_name =~ m/\.tar.Z$/i) {
      `gzip -d $short_name`;
      $short_name =~ s/\.Z$//;
      `tar -xf $short_name`;
   } elsif ($short_name =~ m/\.tar.(gz|zip)$/i) {
      `tar -xzf $short_name`;
   } elsif ($short_name =~ m/\.tar$/i) {
      `tar -xf $short_name`;
   }

   #change to proper directory to compile
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Configuring HMMER archive version $version",
                  progress     => 45,
                 );
   find (sub {$hmmer_dir = $File::Find::dir if ($_ =~ m/license/i) }, $args{directory});
   chdir $hmmer_dir or do {return (0)};

   #configure to proper path
   `./configure --mandir=${$args{ini_ref}}{pfam_executable} --bindir=${$args{ini_ref}}{pfam_executable}`;

   #run Makefile
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Compiling HMMER archive version $version",
                  progress     => 65,
                 );
   `make`;

   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Moving HMMER distribution",
                  progress     => 75,
                 );
   `make install`;

   #make sure everything is accessible and executable
   `chmod -R 777 ${$args{ini_ref}}{pfam_executable}/*`;

   #cleanup source distro
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Cleanup',
                  progress     => 95,
                 );
   chdir ${$args{auto_ini_ref}}{work_dir};
   cleanup ($hmmer_dir);
   unlink $args{directory}.'/'.$short_name;

   #set compiler flag
   open WRITE, '+>'.$args{directory}.'/recompiled' or do {return (0)};
   print WRITE $args{directory};
   close WRITE;

   return (1);
}

sub recompile_hmmer3 {
   my %args = @_;
   my ($hmmer_archive, @hmmer_archive, $hmmer_dir, $short_name, $version);

   #determine bit version of Linux
   my $li_version = `uname -a`;
   if ($li_version =~ / x86_64 /) {
      $li_version = 'x86_64';
   } else {
      $li_version = 'ia32';
   }

   #find current HMMER zip file in Archive directory
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Finding HMMER3 archive',
                  progress     => 15,
                 );

   find (sub {push (@hmmer_archive, $File::Find::name) if (/hmmer\-3.*?\.(Z|tar|gz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );
   unless ($#hmmer_archive >= 0) {
      return (0);
   }

   #get hmmer3 tarball
   foreach my $entry (@hmmer_archive){
      $hmmer_archive = $entry;
      last if ($hmmer_archive =~ /$li_version/);
   }

   unless (defined $hmmer_archive && $hmmer_archive =~ /$li_version/) {
      return (0);
   }

   $short_name = $hmmer_archive;
   $short_name =~ s/.+\/(.+)$/$1/;

   $version    = $short_name;
   $version    =~ s/^(hmmer\-.*?)\.(tar|gz|Z).*/$1/i;

   #does the directory still exist?
   unless (-d $args{directory}) {
      mkdir $args{directory};
   }
   #copy hmmer archive to hmmer exe directory
   copy($hmmer_archive , $args{directory});#.'/'.$short_name);

   #uncompress file
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Uncompressing HMMER3 archive version $version",
                  progress     => 30,
                 );
   chdir $args{directory};
   if ($short_name =~ m/\.tar.Z$/i) {
      `gzip -d $short_name`;
      $short_name =~ s/\.Z$//;
      `tar -xf $short_name`;
   } elsif ($short_name =~ m/\.tar.(gz|zip)$/i) {
      `tar -xzf $short_name`;
   } elsif ($short_name =~ m/\.tar$/i) {
      `tar -xf $short_name`;
   }

   #change to proper directory to compile
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Configuring HMMER3 archive version $version",
                  progress     => 45,
                 );
   find (sub {$hmmer_dir = $File::Find::dir if ($_ =~ m/RELEASE-NOTES/i) }, $args{directory});
   chdir $hmmer_dir or do {return (0)};

   #configure to proper path
   `./configure --prefix=$args{directory}`;# CC=icc LDFLAGS=\"-static\"`;

   #run Makefile
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Compiling HMMER3 archive version $version",
                  progress     => 65,
                 );
   `make`;

   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Moving HMMER3 distribution",
                  progress     => 75,
                 );
   `make install`;

   #test if copying worked
   my @files = qw(hmmalign hmmbuild hmmconvert hmmemit hmmfetch hmmpress hmmscan hmmsearch hmmsim hmmstat jackhmmer phmmer);
   if (-e $args{directory}.'/bin/hmmsearch') {
      foreach my $entry (@files) {
         copy($args{directory}."/bin/".$entry , $args{directory});
      }
   } else {
      foreach my $entry (@files) {
         copy($hmmer_dir."/src/".$entry , $args{directory});
      }
   }

   #set permissions
   chdir $args{directory};
   system "chmod 777 *";

   #cleanup source distro
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Cleanup',
                  progress     => 95,
                 );
   chdir ${$args{auto_ini_ref}}{work_dir};

   #remove problem child manually
   unlink $hmmer_dir.'/src/impl';
   cleanup ($hmmer_dir);

   if (-e $args{directory}.'/bin') {
      cleanup ($args{directory}.'/bin');
   }
   if (-e $args{directory}.'/lib') {
      cleanup ($args{directory}.'/lib');
   }
   if (-e $args{directory}.'/include') {
      cleanup ($args{directory}.'/include');
   }
   unlink $args{directory}.'/'.$short_name;

   #set compiler flag
   open WRITE, '+>'.$args{directory}.'/recompiled' or do {return (0)};
   print WRITE $args{directory};
   close WRITE;

   return (1);
}

sub recompile_glimmer {
   my %args = @_;
   my ($glimmer_archive, $glimmer_dir, $short_name, $version);
   my @glimmer2 = qw(adjust anomaly build-icm check codon-usage compare-lists
                     extract generate get-len get-putative glimmer2 long-orfs
                    );


   #find current glimmer zip file in Archive directory
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Finding Glimmer archive',
                  progress     => 15,
                 );

   if ($args{version} eq 'glimmer2') {
      find (sub {$glimmer_archive = $File::Find::name if (/glimmer2.*?\.(Z|tar|gz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );
   } elsif ($args{version} eq 'glimmer3') {
      find (sub {$glimmer_archive = $File::Find::name if (/glimmer3.*?\.(Z|tar|gz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );
   }
   unless (defined $glimmer_archive && $glimmer_archive =~ /\w+/) {
      return (0);
   }
   $short_name = $glimmer_archive;
   $short_name =~ s/.+\/(.+)$/$1/;

   $version    = $short_name;
   $version    =~ s/^(glimmer.*?)\.(tar|gz|Z).*/$1/i;
   $version    =~ s/^(glimmer\d)(\d+)/$1\.$2/i;

   #clear out directory in case of preciously corrupted installation
   rmtree ($args{directory}, {keep_root => 1} );

   #test if directory exists
   unless (-d $args{directory}) {
      mkdir ($args{directory}, 0777) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create $args{directory} directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
   }

   #copy glimmer archive to glimmer exe directory
   copy($glimmer_archive , $args{directory}.'/'.$short_name);

   #uncompress filegz
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Uncompressing Glimmer version $version",
                  progress     => 30,
                 );
   chdir $args{directory};
   if ($short_name =~ m/\.tar.Z$/i) {
      `gzip -d $short_name`;
      $short_name =~ s/\.Z$//;
      `tar -xf $short_name`;
   } elsif ($short_name =~ m/\.tar.(gz|zip)$/i) {
      `tar -xzf $short_name`;
   } elsif ($short_name =~ m/\.tar$/i) {
      `tar -xf $short_name`;
   }
   `chmod -R 777 *`;

   #change to proper directory to compile
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Configuring Glimmer version $version",
                  progress     => 45,
                 );

   if ($args{version} eq 'glimmer2') {
      find (sub {$glimmer_dir = $File::Find::dir if ($_ =~ m/Makefile/i) }, $args{directory});
      chdir $glimmer_dir or do {return (0)};
      #run Makefile
      &update_pbar_2(title        => 'Recompiling programs',
                     label        => "Compiling Glimmer version $version",
                     progress     => 65,
                    );
      `make`;

      #check if compiled. Tere may be a problem in delcher.h that we can correct
      unless (-e $glimmer_dir.'/glimmer2') {
         print "\nError encountered while compiling Glimmer2. Trying to patch delcher.h\n";
         my $delcher = slurp(directory => $glimmer_dir,
                             filename  => 'delcher.h');

         $$delcher =~ s/#include  <iostream.h>/#include  <iostream>/;
         $$delcher =~ s/#include  <iomanip.h>/#include  <iomanip>/;
         $$delcher =~ s/#include  <fstream.h>/#include  <fstream>/;
         open WRITE, "+>$glimmer_dir\/delcher.h";
         print WRITE $$delcher;
         close WRITE;
         undef $delcher;

         $delcher = slurp(directory => $glimmer_dir,
                          filename  => 'adjust.cc');

         $$delcher =~ s/cin /std::cin /gs;
         $$delcher =~ s/cout /std::cout /gs;
         $$delcher =~ s/setw /std::setw /gs;
         $$delcher =~ s/setiosflags/std::setiosflags/gs;
         $$delcher =~ s/ios /std::ios /gs;
         $$delcher =~ s/endl/std::endl/gs;
         open WRITE, "+>$glimmer_dir\/adjust.cc";
         print WRITE $$delcher;
         close WRITE;
         undef $delcher;

         $delcher = slurp(directory => $glimmer_dir,
                          filename  => 'check.cc');

         $$delcher =~ s/cin /std::cin /gs;
         $$delcher =~ s/cout /std::cout /gs;
         $$delcher =~ s/setw /std::setw /gs;
         $$delcher =~ s/setiosflags/std::setiosflags/gs;
         $$delcher =~ s/ios /std::ios /gs;
         $$delcher =~ s/endl/std::endl/gs;
         open WRITE, "+>$glimmer_dir\/check.cc";
         print WRITE $$delcher;
         close WRITE;
         undef $delcher;

         $delcher = slurp(directory => $glimmer_dir,
                          filename  => 'rnabin.h');

         $$delcher =~ s/cin /std::cin /gs;
         $$delcher =~ s/cout /std::cout /gs;
         $$delcher =~ s/setw /std::setw /gs;
         $$delcher =~ s/setiosflags/std::setiosflags/gs;
         $$delcher =~ s/ios /std::ios /gs;
         $$delcher =~ s/endl/std::endl/gs;
         open WRITE, "+>$glimmer_dir\/rnabin.h";
         print WRITE $$delcher;
         close WRITE;
         undef $delcher;

         $delcher = slurp(directory => $glimmer_dir,
                          filename  => 'rnabin.cc');

         $$delcher =~ s/cin /std::cin /gs;
         $$delcher =~ s/cout /std::cout /gs;
         $$delcher =~ s/setw /std::setw /gs;
         $$delcher =~ s/setiosflags/std::setiosflags/gs;
         $$delcher =~ s/ios /std::ios /gs;
         $$delcher =~ s/endl/std::endl/gs;
         open WRITE, "+>$glimmer_dir\/rnabin.cc";
         print WRITE $$delcher;
         close WRITE;
         undef $delcher;

         `make`;
      }

      #are files there? yes, move; no, use static Glimmer instead
      if (-e $glimmer_dir.'/glimmer2') {
         #move files
         foreach my $file (@glimmer2) {
            move ($glimmer_dir.'/'.$file , $args{directory}.'/'.$file);
         }
      } else {
         chdir $args{directory};
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not compile Glimmer2, using static archive instead.",
                                           -icon    => 'error',
                                           -type    => 'ok');
         copy (${$args{auto_ini_ref}}{work_dir}.'/lib/Archives/glimmer213_static.zip', $args{directory}.'/glimmer213_static.zip');
         `unzip glimmer213_static.zip`;
         unlink 'glimmer213_static.zip';
      }
   } elsif ($args{version} eq 'glimmer3') {
      find (sub {$glimmer_dir = $File::Find::dir if ($_ =~ m/license/i) }, $args{directory});

      #change to source directory
      chdir $glimmer_dir.'/src' or do {return (0)};

      #run Makefile
      &update_pbar_2(title        => 'Recompiling programs',
                     label        => "Compiling Glimmer version $version",
                     progress     => 65,
                    );
      `make`;

      #did compiling go wrong? start patching
      #First try patching the code
      unless (-e $glimmer_dir.'/bin/glimmer3') {
         my $delcher = slurp(directory => $glimmer_dir.'/src/Common',
                             filename  => 'delcher.hh');
         $$delcher =~ s/#include  <string>/#include  <cstring>\n#include  <climits>\n#include  <memory>\n/s;

         open WRITE, "+>$glimmer_dir\/src\/Common\/delcher.hh";
         print WRITE $$delcher;
         close WRITE;
         undef $delcher;

         $delcher = slurp(directory => $glimmer_dir.'/src/Common',
                          filename  => 'gene.cc');
         $$delcher =~ s/char  \* tag \= \"acgt\"/char  \* tag \= \(char\*\)\"acgt\"/s;
         $$delcher =~ s/p \= strchr \(CONVERSION_STRING\, tolower \(ch\)\)/p = \(char\*\)strchr \(CONVERSION_STRING, tolower \(ch\)\)/s;
         open WRITE, "+>$glimmer_dir\/src\/Common\/gene.cc";
         print WRITE $$delcher;
         close WRITE;
         undef $delcher;

         $delcher = slurp(directory => $glimmer_dir.'/src/ICM',
                          filename  => 'icm.cc');
         $$delcher =~ s/p \= strchr \(ALPHA_STRING\, tolower \(Filter \(ch\)\)\)/p \= \(char\*\) strchr \(ALPHA_STRING\, tolower \(Filter \(ch\)\)\)/s;
         open WRITE, "+>$glimmer_dir\/src\/ICM\/icm.cc";
         print WRITE $$delcher;
         close WRITE;
         undef $delcher;

      }

      #remove some features
      unless (-e $glimmer_dir.'/bin/glimmer3') {
         print "\nError encountered while compiling Glimmer3. Trying to patch delcher.hh\n";
         my $delcher = slurp(directory => $glimmer_dir.'/src/Common',
                             filename  => 'delcher.hh');

         $$delcher =~ s/#define  ALLOW_LONG_OPTIONS  1/#define  ALLOW_LONG_OPTIONS  1/;
         open WRITE, "+>$glimmer_dir\/src\/Common\/delcher.hh";
         print WRITE $$delcher;
         close WRITE;
         undef $delcher;
         `make`;
      }

      unless (-e $glimmer_dir.'/bin/glimmer3') {
         print "\nError encountered while compiling Glimmer3. Trying to patch delcher.hh more.\n";
         my $delcher = slurp(directory => $glimmer_dir.'/src/Common',
                             filename  => 'delcher.hh');

         $$delcher =~ s/#include  <getopt.h>//;
         open WRITE, "+>$glimmer_dir\/src\/Common\/delcher.hh";
         print WRITE $$delcher;
         close WRITE;
         undef $delcher;
         `make`;
      }

      unless (-e $glimmer_dir.'/bin/glimmer3') {
         print "\nError encountered while compiling Glimmer3. One last try.\n";
         chdir $glimmer_dir.'/SimpleMake';
         `make`;
      }

      #are files there? yes, move; no, use static Glimmer instead
      if (-e $glimmer_dir.'/bin/glimmer3') {
         #copy exes to proper dir
         chdir $args{directory};
         dircopy($glimmer_dir.'/bin', $args{directory});
      } else {
         chdir $args{directory};
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not compile Glimmer3, using static archive instead.",
                                           -icon    => 'error',
                                           -type    => 'ok');
         copy (${$args{auto_ini_ref}}{work_dir}.'/lib/Archives/glimmer302_static.zip', $args{directory}.'/glimmer302_static.zip');
         `unzip glimmer302_static.zip`;
         unlink 'glimmer302_static.zip';
      }
   }

   #cleanup source distro
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Cleanup',
                  progress     => 95,
                 );

   #make sure everything is accessible and executable
   chdir $args{directory};
   `chmod -R 777 *`;

   chdir ${$args{auto_ini_ref}}{work_dir};
   cleanup ($glimmer_dir);
   unlink $args{directory}.'/'.$short_name;

   #set compiler flag
   open WRITE, '+>'.$args{directory}.'/recompiled' or do {return (0)};
   print WRITE $args{directory};
   close WRITE;

   return (1);
}

sub recompile_critica {
   my %args = @_;
   my ($critica_archive, $critica_dir, $short_name, $version);
   my @critica = qw(lookat motiffind translate extractcoding intergenic
                    empirical-matrix dicodontable addlongorfs removeoverlaps
                    scanblastpairs critica reportscore
                    );


   #find current critica zip file in Archive directory
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Finding Critica archive',
                  progress     => 15,
                 );

   find (sub {$critica_archive = $File::Find::name if (/critica.*?\.(Z|tar|gz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );

   unless (defined $critica_archive && $critica_archive =~ /\w+/) {
      return (0);
   }
   $short_name = $critica_archive;
   $short_name =~ s/.+\/(.+)$/$1/;

   $version    = $short_name;
   $version    =~ s/^(critica\d+)\.(tar|gz|Z).*/$1/i;

   #test if directory exists
   unless (-d $args{directory}) {
      mkdir ($args{directory}, 0777) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create $args{directory} directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
   }

   #copy critica archive to critica exe directory
   copy($critica_archive , $args{directory}.'/'.$short_name);

   #uncompress file
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Uncompressing Critica version $version",
                  progress     => 30,
                 );
   chdir $args{directory};
   if ($short_name =~ m/\.tar.Z$/i) {
      `gzip -d $short_name`;
      $short_name =~ s/\.Z$//;
      `tar -xf $short_name`;
   } elsif ($short_name =~ m/\.tar.(gz|zip)$/i) {
      `tar -xzf $short_name`;
   } elsif ($short_name =~ m/\.tar$/i) {
      `tar -xf $short_name`;
   }

   #change to proper directory to compile
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Configuring Critica version $version",
                  progress     => 45,
                 );
   find (sub {$critica_dir = $File::Find::dir if ($_ =~ m/Makefile/i) }, $args{directory});
   chdir $critica_dir or do {return (0)};

   #run Makefile
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Compiling Critica version $version",
                  progress     => 75,
                 );
   `make`;

   #move files
   foreach my $file (@critica) {
      move ($critica_dir.'/'.$file , $args{directory}.'/'.$file);
   }
   #move script directory
   dircopy ($critica_dir.'/scripts' , $args{directory}.'/scripts');

   #make sure everything is accessible and executable
   chdir $args{directory};
   `chmod -R 777 *`;

   #cleanup source distro
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Cleanup',
                  progress     => 95,
                 );
   chdir ${$args{auto_ini_ref}}{work_dir};
   cleanup ($critica_dir);
   unlink $args{directory}.'/'.$short_name;

   #set compiler flag
   open WRITE, '+>'.$args{directory}.'/recompiled' or do {return (0)};
   print WRITE $args{directory};
   close WRITE;

   return (1);
}

sub recompile_signalp {
   my %args = @_;
   my (@signalp_archive, $signalp_dir, $short_name, $version, $latest);

   #find current signalp zip file in Archive directory
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Finding SignalP archive',
                  progress     => 15,
                 );

   find (sub {push (@signalp_archive, $File::Find::name) if (/signalp.*?\.(Z|tar|gz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );

   unless ($#signalp_archive >= 0) {
      return (0);
   }

   #find latest version
   foreach my $version (@signalp_archive) {
      my ($number, $max_number, );
      'reset' =~ m/reset/;
      $version =~ m/\/signalp\-(.+)\.linux/i;
      $number = $1;
      $number =~ s/\D//g;
      if ($number > $max_number) {$max_number = $number; $latest = $version};
   }

   $short_name = $latest;
   $short_name =~ s/.+\/(.+)$/$1/;

   $version    = $short_name;
   $version    =~ s/^(signalp.+)\.linux\.(tar|gz|Z).*/$1/i;

   #test if directory exists
   unless (-d $args{directory}) {
      mkdir ($args{directory}, 0777) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create $args{directory} directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
   }

   #copy signalp archive to signalp exe directory
   copy($latest , $args{directory}.'/'.$short_name);

   #uncompress file
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Uncompressing SignalP version $version",
                  progress     => 30,
                 );
   chdir $args{directory};
   if ($short_name =~ m/\.tar.Z$/i) {
      `gzip -d $short_name`;
      $short_name =~ s/\.Z$//;
      `tar -xf $short_name`;
   } elsif ($short_name =~ m/\.tar.(gz|zip)$/i) {
      `tar -xzf $short_name`;
   } elsif ($short_name =~ m/\.tar$/i) {
      `tar -xf $short_name`;
   }

   #just to make sure, change access right to all archives:
   `chmod -R 777 $args{directory}/*`;

   #change to proper directory to configuration
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Configuring SignalP version $version",
                  progress     => 45,
                 );
   find (sub {$signalp_dir = $File::Find::dir if ($_ =~ m/readme/i) }, $args{directory});
   my ($file_ref) = &slurp(main_window   => $args{main_window},
                           auto_ini_ref  => $args{auto_ini_ref},
                           ini_ref       => $args{ini_ref},
                           filename      => 'signalp',
                           directory     => $signalp_dir
                          );
   if ($file_ref eq '0') {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => 'Could not read SignalP script in folder '.$signalp_dir,
                                        -type    => 'OK',
                                        -icon    => 'info');
      return (0); #return to main
   }

   #modify filepath
   if ($short_name =~ m/signalp\-3/i) {
      $$file_ref =~ s/SIGNALP=\/usr\/opt\/signalp-3\.0/SIGNALP=$args{directory}/s;
   } elsif ($short_name =~ m/signalp\-4/i) {
      $$file_ref =~ s/\$ENV\{SIGNALP\} = \'.+?\'/\$ENV\{SIGNALP\} = \'$args{directory}\'/i;
      $$file_ref =~ s/\$outputDir = \".+?\"/\$outputDir = \"${$args{ini_ref}}{signalp_results}\"/i;
   } else {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => 'Could not identify SignalP version based on the main script in folder '.$signalp_dir,
                                        -type    => 'OK',
                                        -icon    => 'info');
   }

   #update file
   open WRITE, '+>'.$signalp_dir.'/signalp' or do {return (0)};
   print WRITE $$file_ref;
   close WRITE;
   undef $file_ref;

   #move files
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Moving SignalP version $version",
                  progress     => 75,
                 );
   dircopy ($signalp_dir , $args{directory});

   #make sure everything is accessible and executable
   chdir $args{directory};
   `chmod -R 777 *`;

   #cleanup source distro
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Cleanup',
                  progress     => 95,
                 );
   chdir ${$args{auto_ini_ref}}{work_dir};
   cleanup ($signalp_dir);
   unlink $args{directory}.'/'.$short_name;

   #set compiler flag
   open WRITE, '+>'.$args{directory}.'/recompiled' or do {return (0)};
   print WRITE $args{directory};
   close WRITE;

   return (1);
}

sub recompile_tmhmm {
   my %args = @_;
   my ($tmhmm_archive, $tmhmm_dir, $short_name, $version, $latest);


   #find current tmhmm zip file in Archive directory
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Finding TMHMM archive',
                  progress     => 15,
                 );

   find (sub {$tmhmm_archive = $File::Find::name if (/tmhmm.*?\.(Z|tar|gz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );

   unless (defined $tmhmm_archive && $tmhmm_archive =~ /\w+/) {
      return (0);
   }

   $short_name = $tmhmm_archive;
   $short_name =~ s/.+\/(.+)$/$1/;

   $version    = $short_name;
   $version    =~ s/^(tmhmm.*)\.(tar|gz|Z).*/$1/i;

   #test if directory exists
   unless (-d $args{directory}) {
      mkdir ($args{directory}, 0777) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create $args{directory} directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
   }

   #copy tmhmm archive to tmhmm exe directory
   copy($tmhmm_archive , $args{directory}.'/'.$short_name);
   unless (-e $args{directory}.'/'.$short_name) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Failed to copy the THMHMM archive to directory $args{directory}",
                                        -icon    => 'error',
                                        -type    => 'ok');
   }

   #uncompress file
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Uncompressing TMHMM",
                  progress     => 30,
                 );
   chdir $args{directory};
   if ($short_name =~ m/\.tar.Z$/i) {
      `gzip -d $short_name`;
      $short_name =~ s/\.Z$//;
      `tar -xf $short_name`;
   } elsif ($short_name =~ m/\.tar.(gz|zip)$/i) {
      `tar -xzf $short_name`;
   } elsif ($short_name =~ m/\.tar$/i) {
      `tar -xf $short_name`;
   }

   #change to proper directory to configuration
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Configuring TMHMM",
                  progress     => 45,
                 );
   find (sub {$tmhmm_dir = $File::Find::dir if ($_ =~ m/^README/) }, $args{directory});

   #move files
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Moving TMHMM",
                  progress     => 75,
                 );
   unless (-d $args{directory}.'/bin') {
      mkdir ($args{directory}.'/bin', 0777) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create $args{directory}\/bin directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         hide_pbar_2;
         return (0);
      };
   }
   copy ($tmhmm_dir.'/bin/decodeanhmm' , $args{directory}.'/bin/decodeanhmm');
   dircopy ($tmhmm_dir.'/lib' , $args{directory}.'/lib');

   #make sure everything is accessible and executable
   chdir $args{directory};
   `chmod -R 777 *`;

   #cleanup source distro
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Cleanup',
                  progress     => 95,
                 );
   chdir ${$args{auto_ini_ref}}{work_dir};
   cleanup ($tmhmm_dir);
   unlink $args{directory}.'/'.$short_name;

   #set compiler flag
   open WRITE, '+>'.$args{directory}.'/recompiled' or do {return (0)};
   print WRITE $args{directory};
   close WRITE;

   return (1);
}

sub recompile_transterm {
   my %args = @_;
   my ($transterm_archive, $transterm_dir, $short_name, $version);

   #find current transterm zip file in Archive directory
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Finding Transterm archive',
                  progress     => 15,
                 );

   find (sub {$transterm_archive = $File::Find::name if (/transterm.*?\d\.(zip|Z|tar|gz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );
   unless (defined $transterm_archive && $transterm_archive =~ /\w+/) {
      return (0);
   }
   $short_name = $transterm_archive;
   $short_name =~ s/.+\/(.+)$/$1/;

   $version    = $short_name;
   $version    =~ s/^(transterm.*?)\.(zip|tar|gz|Z).*/$1/i;

   #test if directory exists
   unless (-d $args{directory}) {
      mkdir ($args{directory}, 0777) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create $args{directory} directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
   }

   #copy transterm archive to transterm exe directory
   copy($transterm_archive , $args{directory}.'/'.$short_name);
   unless (-e $args{directory}.'/'.$short_name) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Failed to copy the Transterm archive to directory $args{directory}",
                                        -icon    => 'error',
                                        -type    => 'ok');
   }

   #uncompress file
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Uncompressing Transterm archive version $version",
                  progress     => 30,
                 );
   chdir $args{directory};
   if ($short_name =~ m/\.tar.Z$/i) {
      `gzip -d $short_name`;
      $short_name =~ s/\.Z$//;
      `tar -xf $short_name`;
   } elsif ($short_name =~ m/\.tar.(gz|zip)$/i) {
      `tar -xzf $short_name`;
   } elsif ($short_name =~ m/\.tar$/i) {
      `tar -xf $short_name`;
   } elsif ($short_name =~ m/\.zip$/i) {
      `unzip $transterm_archive`;
   }

   #change to proper directory to compile
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Compiling Transterm archive version $version",
                  progress     => 65,
                 );

   find (sub {$transterm_dir = $File::Find::dir if ($_ =~ m/Makefile/) }, $args{directory});
   chdir $transterm_dir;

   `make clean transterm`;
   #are files there? yes, move; no, use static Transterm instead
   if (-e $transterm_dir.'/transterm') {
      copy ($transterm_dir.'/transterm', $args{directory}.'/transterm');
      copy ($transterm_dir.'/expterm.dat', $args{directory}.'/expterm.dat');
   } else {
      chdir $args{directory};
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Could not compile TransTerm, using static archive instead.",
                                        -icon    => 'error',
                                        -type    => 'ok');
      copy (${$args{auto_ini_ref}}{work_dir}.'/lib/Archives/transterm_static.zip', $args{directory}.'/transterm_static.zip');
      `unzip transterm_static.zip`;
      unlink 'transterm_static.zip';
   }

   #make sure everything is accessible and executable
   chdir $args{directory};
   `chmod -R 777 *`;

   #cleanup source distro
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Cleanup',
                  progress     => 95,
                 );
   chdir ${$args{auto_ini_ref}}{work_dir};
   cleanup ($transterm_dir);
   unlink $args{directory}.'/'.$short_name;

   #set compiler flag
   open WRITE, '+>'.$args{directory}.'/recompiled' or do {return (0)};
   print WRITE $args{directory};
   close WRITE;

   return (1);
}

sub recompile_trnascan {
   my %args = (change_sqio => 0,
               @_);
   my ($trnascan_archive, $trnascan_dir, $short_name, $version);

   #find current trnascan zip file in Archive directory
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Finding tRNAscan archive',
                  progress     => 15,
                 );

   #find tar or zip file
   find (sub {$trnascan_archive = $File::Find::name if (/trnascan.*?\.(Z|tar|gz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );

   unless (defined $trnascan_archive && $trnascan_archive =~ /\w+/) {
      return (0);
   }
   $short_name = $trnascan_archive;
   $short_name =~ s/.+\/(.+)$/$1/;

   $version    = $short_name;
   $version    =~ s/^(trnascan.*?)\.(tar|gz|Z).*/$1/i;

   #test if directory exists
   unless (-d $args{directory}) {
      mkdir ($args{directory}, 0777) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create $args{directory} directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
   }

   #copy trnascan archive to trnascan exe directory
   copy($trnascan_archive , $args{directory}.'/'.$short_name);
   unless (-e $args{directory}.'/'.$short_name) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Failed to copy the tRNAscan-SE archive to directory $args{directory}",
                                        -icon    => 'error',
                                        -type    => 'ok');
   }

   #uncompressing
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Uncompressing tRNAscan archive version $version",
                  progress     => 30,
                 );
   chdir $args{directory};
   if ($short_name =~ m/\.tar.Z$/i) {
      `gzip -dqf $short_name`;
      $short_name =~ s/\.Z$//;
      `tar -xf $short_name`;
      unlink $args{directory}.'/'.$short_name.'.Z';
      unlink $args{directory}.'/'.$short_name;
   } elsif ($short_name =~ m/\.tar.(gz|zip)$/i) {
      `tar -xzf $short_name`;
   } elsif ($short_name =~ m/\.tar$/i) {
      `tar -xf $short_name`;
   }

   #reconfigure Makefile
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Configuring tRNAscan version $version",
                  progress     => 45,
                 );
   find (sub {$trnascan_dir = $File::Find::dir if ($_ =~ m/Makefile/) }, $args{directory});
   my $file_ref = &slurp(main_window => $args{main_window},
                         directory   => $trnascan_dir,
                         filename    => 'Makefile'
                        );

   if ($file_ref eq '0') {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => 'Could not read Makefile in folder '.$args{directory}.'/'.$short_name,
                                        -type    => 'OK',
                                        -icon    => 'info');
      return (0); #return to main
   }

   $$file_ref =~ s/\n\r?BINDIR  = .*?\n\r?/\nBINDIR  = ${$args{ini_ref}}{trnascan_dir}\n/s;
   $$file_ref =~ s/\n\r?LIBDIR  = .*?\n\r?/\nLIBDIR  = ${$args{ini_ref}}{trnascan_dir}\n/s;
   $$file_ref =~ s/\n\r?MANDIR  = .*?\n\r?/\nMANDIR  = ${$args{ini_ref}}{trnascan_dir}\n/s;
   $$file_ref =~ s/\n\r?TEMPDIR = .*?\n\r?/\nTEMPDIR = ${$args{ini_ref}}{trnascan_results}\n/s;

   #delete testrun
   $$file_ref =~ s/testrun:.*?\n\r?\#/\n\n\#/s;

   #write new Makefile
   open WRITE, "+>".$trnascan_dir.'/Makefile' or do {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => 'Error writing tRNAscan-SE Makefile in directory '.$trnascan_dir,
                                        -type    => 'OK',
                                        -icon    => 'info');
   };
   print WRITE $$file_ref;
   close WRITE;
   undef $$file_ref;

   #reconfigure sqio.c on 1st compilation failure
   if ($args{change_sqio} == 1) {
      find (sub {$trnascan_dir = $File::Find::dir if ($_ =~ m/sqio\.c/) }, $args{directory});
      $file_ref = &slurp(main_window => $args{main_window},
                         directory   => $trnascan_dir,
                         filename    => 'sqio.c'
                        );

      if ($file_ref eq '0') {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => 'Could not read sqio.c in folder '.$trnascan_dir,
                                           -type    => 'OK',
                                           -icon    => 'info');
         return (0); #return to main
      }

      $$file_ref =~ s/getline/getlineRSV/gs;

      #write new sqio.c
      open WRITE, "+>".$trnascan_dir.'/sqio.c' or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => 'Error writing tRNAscan-SE sqio.c in directory '.$trnascan_dir,
                                           -type    => 'OK',
                                           -icon    => 'info');
      };
      print WRITE $$file_ref;
      close WRITE;
      undef $$file_ref;
   }

   #run Makefile
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Compiling tRNAscan version $version",
                  progress     => 65,
                 );
   chdir $trnascan_dir;
   `make`;

   #move files
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Moving tRNAscan version $version",
                  progress     => 75,
                 );
   `make install`;

   #are files there? if not, modify makefile further and try again
   unless (-e $trnascan_dir.'/tRNAscan-SE' && $args{change_sqio} == 1) {
      cleanup ($trnascan_dir);
      chdir ${$args{auto_ini_ref}}{work_dir};
      #restart recompile, but change 'sqio.c' along the way;
      &recompile_trnascan(main_window  => $args{main_window},
                          progress_bar => $args{progress_bar},
                          auto_ini_ref => $args{auto_ini_ref},
                          ini_ref      => $args{ini_ref},
                          directory    => $args{directory},
                          change_sqio  => 1
                         );
   }

   #are files there? yes, move; no, use static Glimmer instead
   unless (-e $trnascan_dir.'/tRNAscan-SE') {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Could not compile tRNAscan-SE. Please check manually for problems; the software MUST be compiled using pre-set file paths.\nThe rest of the installation will continue.",
                                        -icon    => 'error',
                                        -type    => 'ok');
      return (1);
   }
   #make sure everything is accessible and executable
   #chdir $args{directory};
   `chmod -R 777 $args{directory}`;

   #cleanup source distro
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Cleanup',
                  progress     => 95,
                 );
   chdir ${$args{auto_ini_ref}}{work_dir};
   #cleanup ($trnascan_dir);
   #unlink $args{directory}.'/'.$short_name;

   #set compiler flag
   open WRITE, '+>'.$args{directory}.'/recompiled' or do {return (0)};
   print WRITE $args{directory};
   close WRITE;

   return (1);
}

sub recompile_infernal {
   my %args = @_;
   my ($infernal_archive, $infernal_dir, $short_name, $version);
   my @infernal = qw(cmalign cmbuild cmemit cmscore cmsearch cm2hmm cm2hmmsearch);

   #find current infernal zip file in Archive directory
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Finding Infernal archive',
                  progress     => 15,
                 );

   find (sub {$infernal_archive = $File::Find::name if (/infernal.*?\.(Z|tar|gz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );
   unless (defined $infernal_archive && $infernal_archive =~ /\w+/) {
      return (0);
   }
   $short_name = $infernal_archive;
   $short_name =~ s/.+\/(.+)$/$1/;

   $version    = $short_name;
   $version    =~ s/^(infernal.*?)\.(tar|gz|Z).*/$1/i;

   #test if directory exists
   unless (-d $args{directory}) {
      mkdir ($args{directory}, 0777) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create $args{directory} directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
   }

   #copy infernal archive to infernal exe directory
   copy($infernal_archive , $args{directory}.'/'.$short_name);

   #uncompress file
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Uncompressing Infernal archive",
                  progress     => 30,
                 );
   chdir $args{directory};
   if ($short_name =~ m/\.tar.Z$/i) {
      `gzip -d $short_name`;
      $short_name =~ s/\.Z$//;
      `tar -xf $short_name`;
   } elsif ($short_name =~ m/\.tar.(gz|zip)$/i) {
      `tar -xzf $short_name`;
   } elsif ($short_name =~ m/\.tar$/i) {
      `tar -xf $short_name`;
   }

   #change to proper directory to compile
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Configuring Infernal archive",
                  progress     => 45,
                 );
   find (sub {$infernal_dir = $File::Find::dir if ($_ =~ m/Userguide.pdf/i) }, $args{directory});

   chdir $infernal_dir;
   `./configure`;

   #compile infernal
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => "Compiling Infernal archive",
                  progress     => 65,
                 );
   `make`;

   #move files
   foreach my $file (@infernal) {
      move ($infernal_dir.'/src/'.$file , $args{directory}.'/'.$file);
   }

   #make sure everything is accessible and executable
   chdir $args{directory};
   `chmod -R 777 *`;

   #cleanup source distro
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Cleanup',
                  progress     => 95,
                 );
   chdir ${$args{auto_ini_ref}}{work_dir};
   cleanup ($infernal_dir);
   unlink $args{directory}.'/'.$short_name;

   #set compiler flag
   open WRITE, '+>'.$args{directory}.'/recompiled' or do {return (0)};
   print WRITE $args{directory};
   close WRITE;

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

sub recompile_rRNA {
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

   #test if directory exists
   unless (-d $args{directory}) {
      mkdir ($args{directory}, 0777) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create $args{directory} directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
   }

   #copy rRNA archive to rRNA db directory
   copy($rRNA_archive , $args{directory}.'/'.$short_name);
   unless (-e $args{directory}.'/'.$short_name) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Failed to copy the rRNA databases archive to directory $args{directory}",
                                        -icon    => 'error',
                                        -type    => 'ok');

   }

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
      if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
         chdir ${$args{ini_ref}}{blast_executables};
         `./formatdb -i $args{directory}/$db -p F -o T -n $args{directory}/$db`;
      } else {
         `${$args{ini_ref}}{blast_plus_executables}/makeblastdb -dbtype nucl -in $args{directory}/$db -parse_seqids -hash_index -out $args{directory}/$db`;
      }
      #remove source fasta file
      unlink $args{directory}.'/'.$db;
   }
   chdir $args{directory};

   #cleanup
   unlink $args{directory}.'/'.$short_name;
   unlink $args{directory}.'/formatdb.log';

   unless (-e $args{directory}.'/16S_rRNA.nhr' &&
           -e $args{directory}.'/23S_rRNA.nhr' &&
           -e $args{directory}.'/5S_rRNA.nhr') {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Failed to compile rRNA databases in directory $args{directory}",
                                        -icon    => 'error',
                                        -type    => 'ok');
   }

   #make sure everything is accessible and executable
   `chmod -R 777 *`;

   #set compiler flag
   open WRITE, '+>'.$args{directory}.'/recompiled';
   print WRITE $args{directory};
   close WRITE;

   return (1);
}

sub recompile_CRT {
   my %args = @_;
   my ($CRT_archive, $CRT_dir, $short_name);

   #find current CRT file in Archive directory
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Finding CRT archive',
                  progress     => 25,
                 );
   find (sub {$CRT_archive = $File::Find::name if (/^CRT.*?\.(zip|Z|tar|gz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );
   $CRT_archive =~ m/.+\/(CRT.+)$/;
   $short_name = $1;

   #does directory and Makefile already exists?
   unless (-d $args{directory}) {
      mkdir ($args{directory}, 0777)  or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create CRISPR directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
   }

   #copy CRT archive to proper directory
   copy ($CRT_archive, $args{directory});
   chdir $args{directory};

   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Unzipping CRT archive',
                  progress     => 55,
                 );

   `unzip $short_name`;
   unlink $short_name;

   #rename to default name
   find (sub {$CRT_archive = $File::Find::name if (/^CRT.*/i)}, $args{directory} );
   rename $CRT_archive, $args{directory}.'/CRT_java.jar';

   #make sure everything is accessible and executable
   `chmod -R 777 *`;

   #test if Results directory exists
   unless (-d ${$args{ini_ref}}{CRISPR_results}) {
      mkdir (${$args{ini_ref}}{CRISPR_results}, 0777)  or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create CRISPR results directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
   }

   #cleanup
   chdir ${$args{auto_ini_ref}}{work_dir};

   #test if Java is installed
   my $java_test = `java -version 2>&1` ;
   unless ($java_test =~ /version/i) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "It seems Java is not installed.\n".
                                                    "Installation continues normal, but CRT will not be avalable unless Java is accessible",
                                        -icon    => 'error',
                                        -type    => 'ok');
   };
   #set compiler flag
   open WRITE, '+>'.$args{directory}.'/recompiled';
   print WRITE $args{directory};
   close WRITE;

   return (1);
}

sub recompile_COGnitor {

   ###################
   return (1);
   ###################

   my %args = @_;
   my ($COGnitor_archive, $short_name, $extension, %COG_subdir);

   #define COGsubdirectories for compiling
   $COG_subdir{COGcognitor}  = 1;
   $COG_subdir{COGlse}       = 1;
   $COG_subdir{COGmakehash}  = 1;
   $COG_subdir{COGreadblast} = 1;
   $COG_subdir{COGtriangles} = 1;

   #find current CRT file in Archive directory
   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Finding COGnitor archive',
                  progress     => 25,
                 );
   find (sub {$COGnitor_archive = $File::Find::name if (/^COGsoft.*?\.(zip|Z|tar|gz|tgz)$/i)}, ${$args{auto_ini_ref}}{work_dir}.'/lib/Archives' );
   $COGnitor_archive =~ m/.+\/(COGsoft.+)$/;
   $short_name = $1;

   #does directory and Makefile already exists?
   unless (-d $args{directory}) {
      mkdir ($args{directory}, 0777)  or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create COGnitor directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      };
   }
   #copy COGnitor archive to proper directory
   copy ($COGnitor_archive, $args{directory});
   chdir $args{directory};

   &update_pbar_2(title        => 'Recompiling programs',
                  label        => 'Uncompressing COGnitor archive',
                  progress     => 55,
                 );
   `tar zxvf $short_name`;

   #remove extension
   $short_name =~ s/(\.[^\.]+?)$//;
   $extension = $1;

   unless (-d $args{directory}.'/'.$short_name) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Error uncompressing COGnitor archive",
                                        -icon    => 'error',
                                        -type    => 'ok');
      return (0);
   }
   unlink $short_name.$extension;

   #compile individual COG progs
   foreach my $COG_prog (keys %COG_subdir) {
      chdir $args{directory}.'/'.$short_name.'/'.$COG_prog;
      `make`;
      unless (-e $args{directory}.'/'.$short_name.'/'.$COG_prog.'/'.$COG_prog) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Error compiling $COG_prog",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return (0);
      }
      move ($args{directory}.'/'.$short_name.'/'.$COG_prog.'/'.$COG_prog, $args{directory}.'/'.$COG_prog);
   }
   #move wrapper script
   move ($args{directory}.'/'.$short_name.'/COGtriangles/COGtriangles.reformat.pl', $args{directory}.'/COGtriangles.reformat.pl');

   #make sure everything is accessible and executable
   chdir $args{directory};
   `chmod -R 777 *`;

   #cleanup
   cleanup ($args{directory}.'/'.$short_name);
   chdir ${$args{auto_ini_ref}}{work_dir};

   #set compiler flag
   open WRITE, '+>'.$args{directory}.'/recompiled';
   print WRITE $args{directory};
   close WRITE;

   return (1);
}

1;