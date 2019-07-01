#!/opt/ActivePerl-5.8/bin/perl

#tRNAscan program calls and reformatting of results
#input arguments: main_window, progress_bar, auto_ini_ref, ini_ref,
#                 directory, input_array_ref, output_folder, input_type

package ProgrammeModules::trnascan;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&trna);
use vars qw();

use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);
use Cwd;
use File::Find;

#local variables
my (%args, $newdir, $DNAseq);

#define format
format DNA_SEQ =
^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<~~
$DNAseq
.

sub trna {
   my %args = @_;
   my (@list, $covariance_model, $executable);

   #max number of fasta entries
   my $max_count = $#{$args{input_array_ref}} + 1;
   if ($max_count < 1) {$max_count = 1};

   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "tRNAscan analysis",
                   label        => 'tRNAscan analysis'
                  );
   &show_pbar_2;

   #remove existing entries if re-use results
   if (${$args{auto_ini_ref}}{reuse_results} == 1) {
      my (@trna_results, %exist_trna_file);
      &update_pbar_2(title        => 'tRNAscan analysis',
                     label        => "Testing for existing results",
                     progress     => 1,
                    );

      #read existing results
      opendir SEQINPUT, ${$args{ini_ref}}{trnascan_results};
      @trna_results = grep /_tRNAscan$/, readdir(SEQINPUT);
      closedir SEQINPUT;
      #create hash
      foreach (@trna_results) {
         $exist_trna_file{$_} = '1';
      }

      #iterate over file model
      foreach my $entry (@{ $args{input_array_ref} }) {
         #skip if result file exists and re-use is activated
         next if (-e ${$args{ini_ref}}{trnascan_results}.'/'.$entry.'_tRNAscan' && -s ${$args{ini_ref}}{trnascan_results}.'/'.$entry.'_tRNAscan' > 0);
         #add if required
         push (@list, $entry);
         #delete empty file
         if (-e ${$args{ini_ref}}{trnascan_results}.'/'.$entry.'_tRNAscan' && -s ${$args{ini_ref}}{trnascan_results}.'/'.$entry.'_tRNAscan' < 1) {
            unlink ${$args{ini_ref}}{trnascan_results}.'/'.$entry.'_tRNAscan';
         }
      }
      undef @trna_results;
      undef %exist_trna_file;
   } else {
      #prepare full job list
      foreach my $entry (@{ $args{input_array_ref} }) {
         push (@list, $entry);
      }
   }

   #start multithreading tRNAscan
   for (my $i=0; $i<=$#list + 1; $i = $i+${$args{auto_ini_ref}}{CPU}) {
      my $count = $i + ${$args{auto_ini_ref}}{CPU};
      if ($count > $#list + 1) {$count = $#list + 1};

      &update_pbar_2(title        => "tRNAscan analysis",
                     label        => "tRNAscan analysis $count of $max_count",
                     progress     => ($count / $max_count) * 100,
                    );

      my @childs = ();
      for (my $j=$i; $j<$count; $j++) {
         #start forking
         my $pid = fork();
         if ($pid) {
            # parent
            push(@childs, $pid);
         } elsif ($pid == 0) {
            # child
            my ($max_sensitivity);
            #define input/output filename
            my $output_file        = ${$args{ini_ref}}{trnascan_results}.'/'.$list[$j].'_tRNAscan';
            my $sec_structure_file = ${$args{ini_ref}}{trnascan_results}.'/'.$list[$j].'_tRNAscan.structure';
            my $input              = ${$args{ini_ref}}{trnascan_dir}.    '/'.$list[$j];

            #reformat and check input file
            my $file_ref = &slurp_cmd(directory => ${$args{ini_ref}}{input_files},
                                      filename  => $list[$j]
                                     );
            my ($header, $dnasequence);
            'reset'    =~ m/reset/;
            $$file_ref =~ m/^(>[^\n\r]*?)\n\r?(.+)/;
            $header    = $1;
            $DNAseq    = $2;

            unless (defined $DNAseq) {
               print "\nError parsing input file $list[$j]\n";
               CORE::exit();
            }

            $header    =~ s/^(>\s*\S*?)/$1/;
            $DNAseq    =~ s/\s+//gs;
            $DNAseq    =~ s/[^acgtnACGTN]/N/g;

            open DNA_SEQ, '>', \$dnasequence;
            write DNA_SEQ;
            close DNA_SEQ;

            open WRITE, "+>$input" or do {
               print "\nError writing temp input file $list[$j]\n";
               CORE::exit();
            };
            print WRITE $header."\n".$dnasequence;
            close WRITE;

            #properly define max_sensitivity parameter
            if (${$args{auto_ini_ref}}{sensitivity} == 1) {
               $max_sensitivity = '-C';
            } else {
               $max_sensitivity = '';
            }

            #run tRNAscan for selected sequence
            chdir ${$args{ini_ref}}{trnascan_dir};
            if (${$args{auto_ini_ref}}{trna_sec_struc} == 0) {
               `tRNAscan-SE -q ${$args{auto_ini_ref}}{trna_code} $max_sensitivity -o $output_file $input`;
            } elsif (${$args{auto_ini_ref}}{trna_sec_struc} == 1) {
               `tRNAscan-SE -q ${$args{auto_ini_ref}}{trna_code} $max_sensitivity -o $output_file -f $sec_structure_file $input`;
            }
            chdir ${$args{auto_ini_ref}}{work_dir};

            #delete temp file
            unlink $input;

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

   &update_pbar_2(title        => 'tRNAscan analysis',
                  label        => "Finished tRNAscan",
                  progress     => 1,
                 );
   &hide_pbar_2;
   return;
}

sub test_recompile {
   my %args = @_;
   my ($test_dir);
   #recompile tRNAscan on first GAMOLA use
   unless (-e ${$args{ini_ref}}{trnascan_dir}.'/recompiled') {
      my ($status) = &recompile(main_window  => $args{main_window},
                                progress_bar => $args{progress_bar},
                                auto_ini_ref => $args{auto_ini_ref},
                                ini_ref      => $args{ini_ref}
                               );
      if ($status == 0) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not compile tRNAscan.\nAborting",
                                           -type    => 'OK',
                                           -icon    => 'info');
         return (0);
      }
   }
   #test if relocated
   open READ, ${$args{ini_ref}}{trnascan_dir}.'/recompiled';
   while (<READ>) {$test_dir = $_};
   close READ;
   $test_dir =~ s/\n//gs;
   if ($test_dir ne ${$args{ini_ref}}{trnascan_dir}) {
      my ($status) = &recompile(main_window  => $args{main_window},
                                progress_bar => $args{progress_bar},
                                auto_ini_ref => $args{auto_ini_ref},
                                ini_ref      => $args{ini_ref}
                               );
      if ($status == 0) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not compile tRNAscan.\nAborting",
                                           -type    => 'OK',
                                           -icon    => 'info');
         return (0);
      }
   }
   return (1);
}


sub recompile {
   my %args = @_;
   my ($compressed_file, $compressed_dir);

   ${$args{progress_bar}}->configure(-label=>"Recompiling tRNAscan before first use");
   ${$args{main_window}}->update;

   #find tar or zip file
   find (sub {$compressed_file = $File::Find::name if (/tRNAscan.*?\.(Z|tar)$/)}, ${$args{ini_ref}}{trnascan_dir} );
   find (sub {$compressed_dir  = $File::Find::dir  if (/tRNAscan.*?\.(Z|tar)$/)}, ${$args{ini_ref}}{trnascan_dir} );

   if (defined $compressed_file && -e $compressed_file) {
      chdir $compressed_dir;
      if ($compressed_file =~ m/\.Z$/i) {
         `uncompress $compressed_file`;
         `tar -xf $compressed_file`;
      } elsif ($compressed_file =~ m/\.tar$/i) {
         `tar -xf $compressed_file`;
      }
      chdir ${$args{auto_ini_ref}}{work_dir};
   }


   #find the proper directory to compile
   find (sub {$newdir =  $File::Find::dir if (/Makefile/) }, ${$args{ini_ref}}{trnascan_dir} );

   #reconfigure Makefile
   my $makefile_ref = &slurp(main_window => $args{main_window},
                             directory   => $newdir,
                             filename    => 'Makefile'
                            );

   if ($makefile_ref eq '0') {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Could not READ Makfile in folder $newdir",
                                        -type    => 'OK',
                                        -icon    => 'info');
      return (0); #return to main
   }

   ${$makefile_ref} =~ s/\n\r?BINDIR  = .*?\n\r?/\nBINDIR  = ${$args{ini_ref}}{trnascan_dir}\n/s;
   ${$makefile_ref} =~ s/\n\r?LIBDIR  = .*?\n\r?/\nLIBDIR  = ${$args{ini_ref}}{trnascan_dir}\n/s;
   ${$makefile_ref} =~ s/\n\r?MANDIR  = .*?\n\r?/\nMANDIR  = ${$args{ini_ref}}{trnascan_dir}\n/s;
   ${$makefile_ref} =~ s/\n\r?TEMPDIR = .*?\n\r?/\nTEMPDIR = ${$args{ini_ref}}{trnascan_dir}\n/s;

   #delete testrun
   ${$makefile_ref} =~ s/testrun:.*?\n\r?\#/\n\n\#/s;

   #write new Makefile
   open WRITE, "+>".$newdir.'/Makefile' or warn "\nError writing Makefile\n";;
   print WRITE ${$makefile_ref};
   close WRITE;
   undef ${$makefile_ref};

   #run Makefile
   chdir $newdir;
   `make`;
   `make install`;
   `make clean`;
   chdir ${$args{auto_ini_ref}}{work_dir};

   #set compiler flag
   open WRITE, '+>'.${$args{ini_ref}}{trnascan_dir}.'/recompiled';
   print WRITE ${$args{ini_ref}}{trnascan_dir};
   close WRITE;
   return (1);
}



1;