#!/opt/ActivePerl-5.8/bin/perl

#non-coding RNA uses Infernal program calls
#input arguments: main_window, progress_bar, auto_ini_ref, ini_ref,
#                 directory, input_array_ref, output_folder, input_type
#takes the 'rfam.txt.gz' descritor from 'ftp://ftp.sanger.ac.uk/pub/databases/Rfam/11.0/database_files'
#and Rfam.cm.1_1.gz from 'ftp://ftp.sanger.ac.uk/pub/databases/Rfam/xx.0'

package ProgrammeModules::ncrna;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&ncrna &recompile_Rfam_fasta);
use vars qw();

use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);
use Cwd;
use File::Copy;
use File::Find;
use File::Path;

#local variables
my (%args, $newdir);

sub ncrna {
   my %args = @_;
   my ($Rfam_cm, $status);

   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'Non-coding RNA analysis',
                   label        => 'Non-coding RNA analysis'
                  );
   &show_pbar_2;

   &update_pbar_2(title        => "Non-coding RNA analysis",
                  label        => "Testing databases",
                  progress     => 1,
                 );
   ($status, $Rfam_cm) = &test_databases(main_window  => $args{main_window},
                                         progress_bar => $args{progress_bar},
                                         auto_ini_ref => $args{auto_ini_ref},
                                         ini_ref      => $args{ini_ref},
                                        );
   if ($status == 0) {
      &hide_pbar_2;
      return(0);
   }

   RFam_10(main_window      => $args{main_window},
           progress_bar     => $args{progress_bar},
           auto_ini_ref     => $args{auto_ini_ref},
           ini_ref          => $args{ini_ref},
           directory        => ${$args{ini_ref}}{input_files},
           input_array_ref  => $args{input_array_ref},
           combined_orf_ref => $args{gene_model_ref},
           Rfam_cm          => $Rfam_cm,
           input_type       => 'fasta'
          );

   &hide_pbar_2;
   return (1);
}




sub RFam_10 {
   my %args = @_;
   my ($chunksize, $Rfam_jobs_ref, $delete_chunks_ref, $Rfam_jobs_count);
   my (%cm_models);

   #read threshold file
   #&update_pbar_2(title        => "Non-coding RNA analysis",
   #               label        => "Reading Rfam CM models",
   #               progress     => 1,
   #              );

   #save memory space and write individual CMs into their respective CM directory
   #&write_CMs(main_window      => $args{main_window},
   #           auto_ini_ref     => $args{auto_ini_ref},
   #           ini_ref          => $args{ini_ref},
   #           Rfam_cm          => $args{Rfam_cm}
   #          );

   #set correct chunksize
   $chunksize = ${$args{auto_ini_ref}}{ncrna_splitseq} * 1000;

   #test if input sequences need to be split up into chunks
   ($Rfam_jobs_ref, $delete_chunks_ref,
    $Rfam_jobs_count) = &split_into_chunks(main_window     => $args{main_window},
                                           auto_ini_ref    => $args{auto_ini_ref},
                                           ini_ref         => $args{ini_ref},
                                           input_array_ref => $args{input_array_ref},
                                           chunksize       => $chunksize
                                          );

   #max number of fasta entries
   my $max_count = $Rfam_jobs_count + 1;
   if ($max_count < 1) {$max_count = 1};

   #start multithreading Blast part of ncRNA
   &multithread_Rfam(main_window        => $args{main_window},
                      auto_ini_ref      => $args{auto_ini_ref},
                      ini_ref           => $args{ini_ref},
                      Rfam_jobs_count   => $Rfam_jobs_count,
                      max_count         => $max_count,
                      Rfam_jobs_ref     => $Rfam_jobs_ref,
                      delete_chunks_ref => $delete_chunks_ref,
                      Rfam_cm           => $args{Rfam_cm}
                     );
}

sub multithread_Rfam {
   my %args = (format => '',
               @_
              );
   my $align_method;

   #define alignment algorithm
   if (${$args{auto_ini_ref}}{ncrna_local} eq 'local') {
      $align_method = '';
   } else {
      $align_method = '-g';
   }

   #start multithreading Rfam
   for (my $i=0; $i<=$args{Rfam_jobs_count} + 1; $i = $i+${$args{auto_ini_ref}}{CPU}) {
      my $count = $i + ${$args{auto_ini_ref}}{CPU};
      if ($count > $args{Rfam_jobs_count} + 1) {$count = $args{Rfam_jobs_count} + 1};

      &update_pbar_2(title        => "Non-coding RNA analysis",
                     label        => "Non-coding RNA analysis: $count of $args{max_count}",
                     progress     => ($count / $args{max_count}) * 100,
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
            my (@list, $seq_ref, $Rfam_result, $Rfam_hits_ref);

            #test if files exists and is complete;
            #skip if re-use
            if (-e ${$args{ini_ref}}{rfam_results}.'/'.${$args{Rfam_jobs_ref}}[$j].'.Rfam') {
               my $last_line = `tail -1 ${$args{ini_ref}}{rfam_results}\/${$args{Rfam_jobs_ref}}[$j]\.Rfam`;
               if ($last_line =~ m/\[ok\]/i && -e ${$args{ini_ref}}{input_files}.'/'.${$args{delete_chunks_ref}}[$j]) {
                  unlink ${$args{ini_ref}}{input_files}.'/'.${$args{delete_chunks_ref}}[$j];
                  #exit safely
                  CORE::exit();
               } else {
                  #uncomplete file, delete
                  unlink ${$args{ini_ref}}{rfam_results}.'/'.${$args{Rfam_jobs_ref}}[$j].'.Rfam';
               }
            }

            chdir ${$args{ini_ref}}{infernal_dir};
            system "./cmsearch $align_method ${$args{auto_ini_ref}}{ncrna_pipelinefilter} --notextw -o ${$args{ini_ref}}{rfam_results}\/${$args{Rfam_jobs_ref}}[$j]\.Rfam $args{Rfam_cm} ${$args{ini_ref}}{input_files}\/${$args{Rfam_jobs_ref}}[$j] ";
            if (-e ${$args{ini_ref}}{input_files}.'/'.${$args{delete_chunks_ref}}[$j] && ${$args{delete_chunks_ref}}[$j] =~ m/_chunk_\d+$/) {
               unlink ${$args{ini_ref}}{input_files}.'/'.${$args{delete_chunks_ref}}[$j];
            }
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
}

sub split_into_chunks {
   my %args          = @_;
   my @blast_jobs    = ();
   my @delete_chunks = ();
   #test if input sequences need to be split up into chunks

   foreach my $entry (@{$args{input_array_ref}}) {
      my ($seq_ref, $header, $sequence);

      #read fasta input file and match header to filename
      ($seq_ref) = &slurp_cmd(auto_ini_ref => $args{auto_ini_ref},
                              ini_ref      => $args{ini_ref},
                              directory    => ${$args{ini_ref}}{input_files},
                              filename     => $entry
                             );
      if ($seq_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not read sequence file $entry\nSkipping entry",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error in ncRNA:".
                        "\nCould not read input file $entry in directory ${$args{ini_ref}}{input_files}\n\n";
         close ERRORLOG;
         next;
      }
      #grab header and sequence
      $$seq_ref =~ m/^>([^\n]+?)\n(.+)/;
      ($header, $sequence) = ($1, $2);
      unless (defined $sequence) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse sequence file $entry\nSkipping entry",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error in ncRNA:".
                        "\nCould not grab content of input file $entry in directory ${$args{ini_ref}}{input_files}\n\n";
         close ERRORLOG;
         next;
      }

      #cleanup sequence
      $sequence =~ s/\s//gs;

      #test if sequence is larger than defined chunklength and split into chunks if necessary
      #if (length($sequence) > $args{chunksize}) {
         my $offset = 0;
         #iterate through long sequence
         while ($offset <= length($sequence)) {
            my $chunk = substr ($sequence, $offset, $args{chunksize});
            open WRITE, "+>${$args{ini_ref}}{input_files}\/$entry\_chunk_$offset";
            print WRITE '>'.$header.' offset: '.$offset."\n".$chunk;
            close WRITE;
            push (@blast_jobs,    $entry.'_chunk_'.$offset);
            push (@delete_chunks, $entry.'_chunk_'.$offset); #add separate delete chunks array for cleanup later
            $offset += $args{chunksize};
         }
      #}
      #if sequence is smaller then chunklength, then just add to joblist
      #else {
      #   push (@blast_jobs, $entry);
      #}
   }
   return (\@blast_jobs, \@delete_chunks, $#blast_jobs);
}


sub write_CMs {
   my %args = @_;

   #save memory space and write individual CMs into their respective CM directory
   unless (-d ${$args{ini_ref}}{rfam_cm_path}) {
      mkdir ( ${$args{ini_ref}}{rfam_cm_path}, 0777);
   }
   unless (-e ${$args{ini_ref}}{rfam_cm_path}.'/RF00001.cm') { #assume that if one is there all are there
      open READ, ">$args{Rfam_cm}";
      while (<READ>) {
         my ($acc, $content);
         if (/^\/\//) {
            open WRITE, "+>${$args{ini_ref}}{rfam_cm_path}\/$acc\.cm";
            print WRITE $content."\n\/\/";
            close WRITE;
         }
         $content .= $_;
         if (m/^accession\s+(\S+)/i) {
            $acc = $1;
         }
      }
      close READ;
   }
}

sub test_databases {
   my %args = @_;
   my ($RFam_db);
   #test if database exists
   my $recompile = 0;
   find (sub {$RFam_db = $File::Find::name if (/Rfam.*\.cm$/i)}, ${$args{ini_ref}}{rfam_db_path} );

   unless (defined $RFam_db && -e $RFam_db) {$recompile = 1};

   if ($recompile) {
      &update_pbar_2(title        => "Non-coding RNA analysis",
                     label        => "Recompiling Rfam database",
                     progress     => 1,
                    );
      my ($status, $RFam_db) = &recompile_Rfam_fasta(main_window  => $args{main_window},
                                                     progress_bar => $args{progress_bar},
                                                     auto_ini_ref => $args{auto_ini_ref},
                                                     ini_ref      => $args{ini_ref}
                                                    );
      if ($status == 0) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not compile Rfam.cm nt database in directory ${$args{ini_ref}}{rfam_db_path}",
                                           -type    => 'OK',
                                           -icon    => 'info');
         &hide_pbar_2;
         return (0); #exit to main
      }
   }
   return (1, $RFam_db);
}

sub recompile_Rfam_fasta {
   my %args = (Rfam_file => '0',
               @_);
   my ($RFam_cm_archive, $RFam_cm);

   #find Rfam CM archive in GAMOLA working directory
   find (sub {$RFam_cm_archive = $File::Find::name if (/Rfam.*\.cm.*?\.gz$/i)}, ${$args{auto_ini_ref}}{work_dir} );
   unless (defined $RFam_cm_archive && $RFam_cm_archive =~ /\w+/) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => 'Could not find Rfam cm database archive',
                                        -type    => 'OK',
                                        -icon    => 'info');
      my $types = [ ['gz',  ['.gz']],
                    ['All', ['.*' ]]
                  ];
      $RFam_cm_archive = ${$args{main_window}}->getOpenFile(-initialdir => ${$args{ini_ref}}{work_dir},
                                                            -title      => 'Select RFam cm archive',
                                                            -filetypes  => $types
                                                           );
      unless (defined $RFam_cm_archive) {
         &hide_pbar_2;
         return (0,0)
      };
   }

   #copy archive into RFam_db directory
   copy($RFam_cm_archive, ${$args{ini_ref}}{rfam_db_path}.'/Rfam.cm.gz');
   #uncompress
   `gunzip -f -q -d ${$args{ini_ref}}{rfam_db_path}/Rfam.cm.gz`;

   #input files present?
   find (sub {$RFam_cm = $File::Find::name if (/Rfam.*\.cm$/i)}, ${$args{ini_ref}}{rfam_db_path} );
   unless (defined $RFam_cm && $RFam_cm =~ /\w+/) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Input file\n\>$RFam_cm_archive\<\n was not found for compiling Rfam",
                                        -type    => 'OK',
                                        -icon    => 'info');
      &hide_pbar_2;
      return (0,0);
   }

   chdir ${$args{auto_ini_ref}}{work_dir};
   hide_pbar_2;
   return (1, $RFam_cm);
}


1;