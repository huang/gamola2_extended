#!/opt/ActivePerl-5.8/bin/perl
#Central station: main hub for directing all analyses
#input arguments: main_window, ini_ref, auto_ini_ref,


package ProgrammeModules::CentralStation;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&run_analysis);
use vars qw();

use initialise::read_me               qw(:DEFAULT);
use ProgrammeModules::sequence        qw(:DEFAULT);
use ProgrammeModules::genbank_parser  qw(:DEFAULT);
use ProgrammeModules::blast           qw(:DEFAULT);
use ProgrammeModules::COG_blast       qw(:DEFAULT);
use ProgrammeModules::pfam            qw(:DEFAULT);
use ProgrammeModules::tigrfam         qw(:DEFAULT);
use ProgrammeModules::trnascan        qw(:DEFAULT);
use ProgrammeModules::transterm       qw(:DEFAULT);
use ProgrammeModules::tmhmm_run       qw(:DEFAULT);
use ProgrammeModules::signalp         qw(:DEFAULT);
use ProgrammeModules::ncrna           qw(:DEFAULT);
use ProgrammeModules::rrna            qw(:DEFAULT);
use ProgrammeModules::vector_screen   qw(:DEFAULT);
use ProgrammeModules::CRISPR          qw(:DEFAULT);
use CompilerModules::blast_parser     qw(:DEFAULT);
use CompilerModules::cog_parser       qw(:DEFAULT);
use CompilerModules::pfam_parser      qw(:DEFAULT);
use CompilerModules::tigrfam_parser   qw(:DEFAULT);
use CompilerModules::signalp_parser   qw(:DEFAULT);
use CompilerModules::tmhmm_parser     qw(:DEFAULT);
use CompilerModules::trnascan_parser  qw(:DEFAULT);
use CompilerModules::transterm_parser qw(:DEFAULT);
use CompilerModules::ncrna_parser     qw(:DEFAULT);
use CompilerModules::rrna_parser      qw(:DEFAULT);
use CompilerModules::CRISPR_parser    qw(:DEFAULT);
use CompilerModules::vector_parser    qw(:DEFAULT);
use CompilerModules::rbs_parser       qw(:DEFAULT);
use GeneModel::run_and_process        qw(:DEFAULT);
use GeneModel::read_genemodel         qw(:DEFAULT);
use initialise::gb_header             qw(:DEFAULT);
use initialise::recompile             qw(:DEFAULT);
use Basics::progress_bar              qw(:DEFAULT);
use Basics::Rename;
use Basics::Recursive                 qw(dircopy);
use Cwd;
use File::Copy;
use File::Find;

#local variables
my (%args, $status, $input_list_ref, $orf_list_ref, @fasta_input, @genbank_input,
    @combined_gm, $file_ref, $progress, @rename_files, $rename_filename, $rename_id,
    $COGcode_to_number, $COGcode_to_header, $COGcode_to_letter, $COGletter_to_family,$COGcode_to_phylo,
    $arCOGcode_to_COG, $arCOGacc_to_org, $COG2008genomes, $COG2008_to_def, $COG2008_to_acc,
    $COG2014genomes, $COG2014_to_def, $COG2014_to_acc, $arCOGcode2014_to_COG, $arCOGacc2014_to_org,
    $COG2014_to_refseq,
    $POG2013_to_gi,
    $genbank_header, $genbank_source, $conc_files_hash);


#make concatenate_file-hash global to catch all single file names
my %local_name_hash = ();

#central hub
sub run_analysis {
   my %args = @_;
   my ($input, $count, $gm_hash_ref, $fasta_gene_model_ref, $fasta_input_list_ref,
       $gb_gene_model_ref, $gb_input_list_ref, $gene_model_ref, $input_list_ref,
       @msgenbank_move);

   #resetting system
   @{$input_list_ref} = ();
   @fasta_input       = ();
   @genbank_input     = ();
   @msgenbank_move    = ();

   #clear all output folders if so selected
   &clear_results(progress_bar  => $args{progress_bar},
                  main_window   => $args{main_window},
                  auto_ini_ref  => $args{auto_ini_ref},
                  ini_ref       => $args{ini_ref},
                 );

   #test for recompiling
   &recompile(progress_bar  => $args{progress_bar},
              main_window   => $args{main_window},
              auto_ini_ref  => $args{auto_ini_ref},
              ini_ref       => $args{ini_ref},
             );

   #test for concatenated file clusters
   if (${$args{auto_ini_ref}}{concatenate_input_files} == 1) {
      #re-use results?
      if (${$args{auto_ini_ref}}{reuse_results} == 1) {
         #if reuse, then cluster files have been saved already and concatenate_clusters hash has been build
         #check if all files relevant to each cluster are present
         foreach my $entry (sort {$a<=>$b} keys %{ $args{concatenate_clusters} }) {
            my ($conc_file_name);
            #test for concatenated file
            'reset' =~ m/reset/;
            $args{concatenate_clusters}->{$entry} =~ m/^File_name:\t([^\n]+)\n/;
            $conc_file_name = $1;

            #file exists and contig order exists and has size? if not, treat as new
            unless (-e ${$args{ini_ref}}{input_files}.'/'.$conc_file_name.'.cb'                         &&
                    -s ${$args{ini_ref}}{input_files}.'/'.$conc_file_name.'.cb' > 0                     &&
                    -e ${$args{ini_ref}}{input_files}.'/contig_order/'.$conc_file_name.'.contig_order'    &&
                    -s ${$args{ini_ref}}{input_files}.'/contig_order/'.$conc_file_name.'.contig_order' > 0) {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Concatenation system for cluster $entry could not be established.\nReturning files to original state",
                                                             -buttons => ['OK'],
                                                             -bitmap  => 'info'
                                                            );
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error re-using existing concatenation selections".
                              "\nConcatenation system for cluster $entry could not be established.\nReturning files to original state\n\n";
               close ERRORLOG;

               #move original files to input directory
               my @conc_files = split/\n/, $args{concatenate_clusters}->{$entry};
               foreach my $file (@conc_files) {
                  next if ($file =~ m/^File_name\:/);
                  if (-e ${$args{ini_ref}}{move_msfasta}.'/'.$file) {
                     move (${$args{ini_ref}}{move_msfasta}.'/'.$file , ${$args{ini_ref}}{input_files}.'/'.$file);
                  }
               }
               #delete conc file, contig order and cluster file
               unlink ${$args{ini_ref}}{input_files}.'/'.$conc_file_name.'.cb' ;
               unlink ${$args{ini_ref}}{input_files}.'/contig_order/'.$conc_file_name.'.contig_order';
            }
         }
      } elsif (${$args{auto_ini_ref}}{reuse_results} == 0) {
         # all that is required is to delete existing contig_order files
         # current concatenate_clusters hash contains the selected clusters and does not need to be modified
         my @input = ();
         opendir INPUT, ${$args{auto_ini_ref}}{work_dir};
         @input = grep /\.cluster_files$/, readdir(INPUT);
         closedir INPUT;
         foreach my $entry (@input) {
            unlink ${$args{auto_ini_ref}}{work_dir}.'/'.$entry;
         }
      }

      &concatenate(progress_bar         => $args{progress_bar},
                   main_window          => $args{main_window},
                   auto_ini_ref         => $args{auto_ini_ref},
                   ini_ref              => $args{ini_ref},
                   concatenate_clusters => $args{concatenate_clusters}
                  ); #concatenate_clusters hash contains concat_file name plus individual files

      #move source files to msfasta directory
      foreach my $entry (keys %{ $args{concatenate_clusters} }) {
         my @conc_files = split/\n/, $args{concatenate_clusters}->{$entry};
         undef $file_ref;
         foreach my $file (@conc_files) {
            next if ($file =~ m/^File_name\:/);
            if (-e ${$args{ini_ref}}{input_files}.'/'.$file) {
               my ($status) = move (${$args{ini_ref}}{input_files}.'/'.$file , ${$args{ini_ref}}{move_msfasta}.'/'.$file);
               if ($status eq '0') {
                  my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                                -text    => "Error moving original file $file for concatenate cluster $entry.\nOriginal file will remain in file queue.",
                                                                -buttons => ['OK'],
                                                                -bitmap  => 'info'
                                                               );
                  $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
                  $error_msg-> Show();
                  open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
                  print ERRORLOG "Error moving concatenation selections".
                                 "\nError moving original file $file for concatenate cluster $entry.\nOriginal file will remain in file queue.\n\n";
                  close ERRORLOG;
               }
            }
         }
      }
   }

   #grab all input files from input folder
   #this automatically updates the files that were concatenated
   ${$args{progress_bar}}->configure(-label=>"Reading input data");
   ${$args{main_window}}->update;
   ($input_list_ref, $status) = &read_dir(main_window   => $args{main_window},
                                          auto_ini_ref  => $args{auto_ini_ref},
                                          ini_ref       => $args{ini_ref},
                                         );
   if ($status eq '0') {
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error processing input files".
                     "\nCould not read from directory ${$args{ini_ref}}{input_files}\n\n";
      close ERRORLOG;
   #   ${$args{main_window}}->messageBox(-title   => 'Error',
   #                                     -message => "Could not read from directory ${$args{ini_ref}}{input_files}",
   #                                     -type    => 'OK',
   #                                     -icon    => 'info');
      return (0); #exit to main
   }

   #no input files?
   if ($#{$input_list_ref} < 0) {
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error processing input files".
                     "\nNo input files in directory ${$args{ini_ref}}{input_files}\n\n";
      close ERRORLOG;
   #   ${$args{main_window}}->messageBox(-title   => 'Error',
   #                                     -message => "No input files in directory ${$args{ini_ref}}{input_files}",
   #                                     -type    => 'OK',
   #                                     -icon    => 'info');
      return (0); #return to main
   };

   #check format of input files: Fasta and GB OK, everyting else gets error message and will be skipped
   foreach my $input (@{$input_list_ref}) {
      open READ, '<'.${$args{ini_ref}}{input_files}.'/'.$input or do {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not open file $input in directory ${$args{ini_ref}}{input_files}",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nCould not open file $input in directory ${$args{ini_ref}}{input_files}\n\n";
         close ERRORLOG;
         return (0); #return to main
      };
      while (<READ>) {
         next if ($_ =~ m/^\s+$/);
         if ($_ =~ /^\>/) {
            push (@fasta_input, $input); #fasta format
         } elsif ($_ =~ /^locus/i) {
            push (@genbank_input, $input); #GenBank format
         } else {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Format for input file $input is not recognised. Skipping for analysis.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info'
                                                          );
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error processing input files".
                           "\nFormat for input file $input is not recognised. Skipping for analysis.\n\n";
            close ERRORLOG;
         }
         last;
      }
      close READ;
   }

   #check for msfasta files and modify as necessary
   &msfasta(main_window    => $args{main_window},
            progress_bar   => $args{progress_bar},
            auto_ini_ref   => $args{auto_ini_ref},
            ini_ref        => $args{ini_ref},
            input_list_ref => \@fasta_input,
            directory      => ${$args{ini_ref}}{input_files}
           );

   #check for msGenbank files and modify as necessary
   #will be converted to fasta and concatenated, require gene model rebuild!
   &msGenbank(main_window        => $args{main_window},
              progress_bar       => $args{progress_bar},
              auto_ini_ref       => $args{auto_ini_ref},
              ini_ref            => $args{ini_ref},
              input_list_ref     => \@genbank_input,
              fasta_list_ref     => \@fasta_input,
              directory          => ${$args{ini_ref}}{input_files},
              msgenbank_move_ref => \@msgenbank_move
             );

   #create Genbank headers for fasta files
   foreach my $fasta_file (@fasta_input) {
      my ($header, $seq);
      #read file into memory
      ($file_ref) = &slurp(main_window   => $args{main_window},
                           auto_ini_ref  => $args{auto_ini_ref},
                           ini_ref       => $args{ini_ref},
                           filename      => $fasta_file,
                           directory     => ${$args{ini_ref}}{input_files}
                          );
      if ($file_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not find fasta input file $fasta_file. Terminating analysis",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info'
                                                      );
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error reading input files".
                        "\nCould not find fasta input file $fasta_file. Terminating analysis.\n\n";
         close ERRORLOG;
         return (0);
      }
      'reset'         =~ /reset/;
      $$file_ref      =~ m/^\>([^\n\r]*?)\n\r?(.+)/s;
      ($header, $seq) = ($1, $2);
      $header         =~ s/\s+Length\s*\=.*//;
      $header         =~ s/\W+/_/g;
      $seq            =~ s/\s+//gs;

      #submit to Genbank header for individual genbank headers
      if (${$args{auto_ini_ref}}{fasta_header_individual} == 1) {
         ($genbank_header->{$fasta_file}) = &gb_header(main_window       => $args{main_window},
                                                       auto_ini_ref      => $args{auto_ini_ref},
                                                       ini_ref           => $args{ini_ref},
                                                       filename          => $fasta_file,
                                                       seq_length        => length($seq),
                                                       definition        => $header
                                                      );
         #use default header if no new header defined
         if ($genbank_header->{$fasta_file} !~ /\w+/) {
            ($genbank_header->{$fasta_file}) = &default_header(main_window       => $args{main_window},
                                                            auto_ini_ref      => $args{auto_ini_ref},
                                                            ini_ref           => $args{ini_ref},
                                                            filename          => $fasta_file,
                                                            seq_length        => length($seq),
                                                            definition        => $header
                                                           );
         }
      }
      #use default Genbank header if selected
      elsif (${$args{auto_ini_ref}}{fasta_header_individual} == 0) {
         ($genbank_header->{$fasta_file}) = &default_header(main_window       => $args{main_window},
                                                            auto_ini_ref      => $args{auto_ini_ref},
                                                            ini_ref           => $args{ini_ref},
                                                            filename          => $fasta_file,
                                                            seq_length        => length($seq),
                                                            definition        => $header
                                                           );
      }

      $genbank_source->{$fasta_file} = '     source          1..'.length($seq)."\n".
                                       '                     /organism="'.$header.'"';
   }

   #pre-define Genbank headers first if selected
   foreach my $gb_file (@genbank_input) {
      #get header and source
      my ($header, $source) = &gb_parser (main_window   => $args{main_window},
                                          progress_bar  => $args{progress_bar},
                                          auto_ini_ref  => $args{auto_ini_ref},
                                          ini_ref       => $args{ini_ref},
                                          gb_file       => ${$args{ini_ref}}{input_files}.'/'.$gb_file,
                                          process_gb    => 'sequence'
                                         );

      #submit to individual Genbank header
      if (${$args{auto_ini_ref}}{gb_header_new} == 1) {
         ($genbank_header->{$gb_file}) = &gb_header(main_window       => $args{main_window},
                                                    auto_ini_ref      => $args{auto_ini_ref},
                                                    ini_ref           => $args{ini_ref},
                                                    predefined_header => $header,
                                                    filename          => $gb_file
                                                   );
         #re-use old header if no new header defined
         if ($genbank_header->{$gb_file} !~ /\w+/) {
            $genbank_header->{$gb_file} = $header;
         }
      }
      #re-use existing Genbank header
      elsif (${$args{auto_ini_ref}}{gb_header_new} == 0) {
         ($genbank_header->{$gb_file}) = $header;
      }

      $genbank_source->{$gb_file} = $source;

   }

   #analyse Genbank files first and retrieve gene model and generated fasta file lists
   if ($#genbank_input >= 0) {
      ($gb_gene_model_ref, $gb_input_list_ref) = &from_gb(main_window    => $args{main_window},
                                                          progress_bar   => $args{progress_bar},
                                                          auto_ini_ref   => $args{auto_ini_ref},
                                                          ini_ref        => $args{ini_ref},
                                                          input_list_ref => \@genbank_input,
                                                         );

      if ($gb_gene_model_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while analysing Genbank files.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nError while analysing Genbank files.\n\n";
         close ERRORLOG;
         return (0); #return to main
      }
   }
   #define empty arrays
   else {
      @{$gb_gene_model_ref} = ();
      @{$gb_input_list_ref} = ();
   }

   #then analyse remaining fasta files
   if ($#fasta_input >= 0) {
      ($fasta_gene_model_ref, $fasta_input_list_ref) = &from_fasta(main_window    => $args{main_window},
                                                                   progress_bar   => $args{progress_bar},
                                                                   auto_ini_ref   => $args{auto_ini_ref},
                                                                   ini_ref        => $args{ini_ref},
                                                                   input_list_ref => \@fasta_input,
                                                                   source         => 'fasta'
                                                                  );
      if ($fasta_gene_model_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while analysing Fasta files.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nError while analysing Fasta files.\n\n";
         close ERRORLOG;
         return (0); #return to main
      }
   }
   #define empty arrays
   else {
      @{$fasta_gene_model_ref} = ();
      @{$fasta_input_list_ref} = ();
   }

   #merge fasta and gb lists for compiler
   @{$gene_model_ref} = (@{$gb_gene_model_ref}, @{$fasta_gene_model_ref});
   @{$input_list_ref} = (@{$gb_input_list_ref}, @{$fasta_input_list_ref});

   #run program modules for each file, then return
   #reset status frame
   ${$args{progress_bar}}->configure(-label=>" Running analyses ");
   ${$args{main_window}}->update;
   &run_programmes (main_window     => $args{main_window},
                    progress_bar    => $args{progress_bar},
                    auto_ini_ref    => $args{auto_ini_ref},
                    ini_ref         => $args{ini_ref},
                    gene_model_ref  => $gene_model_ref,
                    input_list_ref  => $input_list_ref
                   );

   #run result compiler modules on all input files
   ${$args{progress_bar}}->configure(-label=>" Parsing results ");
   ${$args{main_window}}->update;
   &run_compiler (main_window      => $args{main_window},
                  progress_bar     => $args{progress_bar},
                  auto_ini_ref     => $args{auto_ini_ref},
                  ini_ref          => $args{ini_ref},
                  gene_model_ref   => $gene_model_ref,
                  input_list_ref   => $input_list_ref,
                  genbank_source   => $genbank_source,
                  genbank_header   => $genbank_header,
                 );

   #enter contig boundaries if selected
   if (${$args{auto_ini_ref}}{concatenate_input_files} == 1) {
      #enter contig boundaries
      &contig_boundaries(progress_bar         => $args{progress_bar},
                         main_window          => $args{main_window},
                         auto_ini_ref         => $args{auto_ini_ref},
                         ini_ref              => $args{ini_ref},
                        );
      #cleanup and move all files
      &move_cb_files(progress_bar         => $args{progress_bar},
                     main_window          => $args{main_window},
                     auto_ini_ref         => $args{auto_ini_ref},
                     ini_ref              => $args{ini_ref},
                     concatenate_clusters => $args{concatenate_clusters},
                     msgenbank_move_ref   => \@msgenbank_move
                    );
   }

   #remove GB based fasta files and move original Genbank files back
   foreach my $corename (@{$gb_input_list_ref}) {
      $corename =~ s/_GAMOLAdna$//;
      unlink ${$args{ini_ref}}{input_files}.'/'.$corename.'_GAMOLAdna';
      move (${$args{ini_ref}}{move_gb}.'/'.$corename , ${$args{ini_ref}}{input_files}.'/'.$corename);
   }

   #compress results if selected
   if (${$args{auto_ini_ref}}{compress_results} == 1) {
      &compress_results(main_window    => $args{main_window},
                        progress_bar   => $args{progress_bar},
                        auto_ini_ref   => $args{auto_ini_ref},
                        ini_ref        => $args{ini_ref},
                        input_list_ref => $input_list_ref
                       );
   }


   ${$args{progress_bar}}->configure(-label=>"Finished");
   ${$args{main_window}}->update;
   if (${$args{auto_ini_ref}}{sort_results} == 1) {
      ${$args{main_window}}->messageBox(-title   => 'Finished',
                                        -message => "Finished annotation Pipeline.\nAnnotations are stored in $args{auto_ini_ref}{work_dir}\/Consolidated_results",
                                        -type    => 'OK',
                                        -icon    => 'info');
   } else {
      ${$args{main_window}}->messageBox(-title   => 'Finished',
                                        -message => "Finished annotation Pipeline.\nAnnotations are stored in $args{ini_ref}{results}",
                                        -type    => 'OK',
                                        -icon    => 'info');
   }
   if (-s ${$args{auto_ini_ref}}{work_dir}.'/Error.log' > 0) {
      ${$args{main_window}}->messageBox(-title   => 'Errors generated',
                                        -message => "This annotation run created runtime errors.\nCheck Error.log file for more details\n",
                                        -type    => 'OK',
                                        -icon    => 'info');
   }

   #clear variables
   undef @fasta_input;
   undef @genbank_input;
   undef @combined_gm;
   undef $file_ref;
   undef $progress;
   undef @rename_files;
   undef $rename_filename;
   undef $rename_id;
   undef $genbank_header;
   undef $genbank_source;
   undef $conc_files_hash;
   %{$args{concatenate_clusters}} = ();
   undef $args{concatenate_clusters};
   ${$args{auto_ini_ref}}{concatenation_counter} = 0;
   undef %local_name_hash;
   undef $input;
   undef $count;
   undef $gm_hash_ref;
   undef $fasta_gene_model_ref;
   undef $fasta_input_list_ref;
   undef $gb_gene_model_ref;
   undef $gb_input_list_ref;
   undef $gene_model_ref;
   undef $input_list_ref;



   &hide_pbar_1;
   &hide_pbar_2;
   &hide_pbar_3;
   return (1);

}

sub read_dir {
   my %args = @_;
   my @inputfile = ();
   my @temp_dir_array = ();

   opendir SEQINPUT, ${$args{ini_ref}}{input_files} or do {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Could not open directory ${$args{ini_ref}}{input_files}",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error processing input files".
                     "\nCould not open directory ${$args{ini_ref}}{input_files}\n\n";
      close ERRORLOG;
      return (undef, 0);
   };
   @inputfile = grep !/^\./, readdir(SEQINPUT);
   if ($inputfile[0] eq "") {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No files present in ${$args{ini_ref}}{input_files}",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error processing input files".
                     "\nNo files present in ${$args{ini_ref}}{input_files}\n\n";
      close ERRORLOG;
      return (undef, 0);
   }

   foreach my $dir_check (@inputfile) {
      my $jjj = ${$args{ini_ref}}{input_files}.'/'.$dir_check;
      if (-d $jjj) {next};
      push (@temp_dir_array, $dir_check);
   }
   @inputfile = @temp_dir_array;
   @temp_dir_array = ();
   return (\@inputfile, 1);
}

sub clear_results {
   my %args = @_;

   if (${$args{auto_ini_ref}}{reuse_results} == 0) {
      my @inputfile = ();
      my $curdir = getcwd();
      my @clear_folder = qw[results blast_results COG_results pfam_results TIGRfam_results
                            annotation_transfer_dir session_comparison
                            trnascan_results transterm_results tmhmm_results CRISPR_results
                            rfam_results rrna_results genemodel_output vector_results
                           ];

      foreach my $folder (@clear_folder) {
         ${$args{progress_bar}}->configure(-label=>"Cleaning up $folder folder");
         ${$args{main_window}}->update;
         chdir ${$args{ini_ref}}{$folder} or do {&error_message(dir         => ${$args{ini_ref}}{$folder},
                                                                main_window => $args{main_window});
                                                 return(0);
                                                };
         unlink <*>;
      }

      #cleaning SignalP and TMHMM folder
      @clear_folder = qw[tmhmm_results signalp_results];
      foreach my $folder (@clear_folder) {
         ${$args{progress_bar}}->configure(-label=>"Cleaning up $folder folder");
         ${$args{main_window}}->update;
         opendir SEQINPUT, ${$args{ini_ref}}{$folder} or do {&error_message(dir         => ${$args{ini_ref}}{$folder},
                                                                            main_window => $args{main_window});
                                                             return(0);
                                                            };
         @inputfile = grep !/^\./, readdir(SEQINPUT);
         foreach my $dir (@inputfile) {
            if (-d ${$args{ini_ref}}{$folder}.'/'.$dir) {
               chdir ${$args{ini_ref}}{$folder}.'/'.$dir or do {&error_message(dir         => ${$args{ini_ref}}{$folder},
                                                                               main_window => $args{main_window});
                                                                return(0);
                                                               };
               unlink <*>;
               chdir ${$args{ini_ref}}{$folder} or do {&error_message(dir         => ${$args{ini_ref}}{$folder},
                                                                      main_window => $args{main_window});
                                                       return(0);
                                                      };
               rmdir ${$args{ini_ref}}{$folder}.'/'.$dir;
            }
         }
         chdir ${$args{ini_ref}}{$folder} or do {&error_message(dir         => ${$args{ini_ref}}{$folder},
                                                                main_window => $args{main_window});
                                                 return(0);
                                                };
         unlink <*>;
         chdir $curdir;
      }

   }
   #clear IG blast results if selected
   #requires both flags to be set. IG can take a long time to run and this acts as an additional safety layer
   if (${$args{auto_ini_ref}}{reuseIGresults} == 0 && ${$args{auto_ini_ref}}{reuse_results} == 0) {
      my @inputfile = ();
      my $curdir = getcwd();
      ${$args{progress_bar}}->configure(-label=>"Cleaning up IG Blast result folder");
      ${$args{main_window}}->update;
      chdir ${$args{ini_ref}}{ig_results} or do {&error_message(dir         => ${$args{ini_ref}}{ig_results},
                                                                main_window => $args{main_window});
                                                 return(0);
                                                };
      unlink <*>;
      chdir $curdir;
   } else {
      ${$args{progress_bar}}->configure(-label=>"Keeping IG Blast results");
      ${$args{main_window}}->update;
   }
}

sub compress_results {
   my %args = @_;
   my ($count, $max_count);
   my @result_folders = (${$args{ini_ref}}{results},
                         ${$args{ini_ref}}{ig_results},
                         ${$args{ini_ref}}{genemodel_output},
                         ${$args{ini_ref}}{blast_results},
                         ${$args{ini_ref}}{COG_results},
                         ${$args{ini_ref}}{pfam_results},
                         ${$args{ini_ref}}{TIGRfam_results},
                         ${$args{ini_ref}}{trnascan_results},
                         ${$args{ini_ref}}{tmhmm_results},
                         ${$args{ini_ref}}{signalp_results},
                         ${$args{ini_ref}}{transterm_results},
                         ${$args{ini_ref}}{CRISPR_results},
                         ${$args{ini_ref}}{rfam_results},
                         ${$args{ini_ref}}{rrna_results},
                         ${$args{ini_ref}}{vector_results},
                        );
   $max_count = $#result_folders + 2;

   #set up progress bar
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'Compressing results',
                   label        => 'Compressing results'
                  );
   &show_pbar_2;

   #first compress object results
   &update_pbar_1(title        => 'Compressing results',
                  label        => "Compressing object related results in folder ${$args{ini_ref}}{results}"
                 );
   ${$args{progress_bar}}->configure(-label=>"Compressing result data...object related data");
   ${$args{main_window}}->update;

   chdir ${$args{auto_ini_ref}}{work_dir};
   #set name if empty
   if (${$args{auto_ini_ref}}{project_name} !~ /\w+/) {${$args{auto_ini_ref}}{project_name} = 'genome'};

   #remove existing zip file
   unlink ${$args{auto_ini_ref}}{project_name}.'.object_results.zip';

   $count = 1;
   foreach my $directory (@result_folders) {
      my $short_name;
      $directory =~ m/.*?\/(.+)$/;
      $short_name = $1;
      &update_pbar_2(title        => "Compressing results",
                     label        => "Compressing object related folder $short_name",
                     progress     => ($count / $max_count) * 100,
                    );
      $count++;
      `zip -q -u -r -9 ${$args{auto_ini_ref}}{project_name}\.object_results\.zip $directory`;
   }

   #compress consolidated results if selected
   if (${$args{auto_ini_ref}}{sort_results} == 1) {
      &update_pbar_1(title        => 'Compressing results',
                     label        => "Compressing consolidated results in folder ${$args{auto_ini_ref}}{work_dir}\/Consolidated_results"
                    );
      ${$args{progress_bar}}->configure(-label=>"Compressing result data...consolidated data");
      ${$args{main_window}}->update;


      chdir ${$args{auto_ini_ref}}{work_dir};
      #set name if empty
      if (${$args{auto_ini_ref}}{project_name} !~ /\w+/) {${$args{auto_ini_ref}}{project_name} = 'genome'};

      #remove existing zip file
      unlink ${$args{auto_ini_ref}}{project_name}.'.consolidated_results.zip';

      $count = 1;
      $max_count = $#{$args{input_list_ref}} + 2;
      foreach my $input_file (@{$args{input_list_ref}}) {
         #clean up name if required
         $input_file =~ s/_GAMOLAdna$//;
         $input_file =~ s/\.(gb|gbk)$//i;
         &update_pbar_2(title        => "Compressing results",
                        label        => "Compressing consolidated results for $input_file",
                        progress     => ($count / $max_count) * 100,
                       );
         $count++;
         `zip -q -u -r -9 ${$args{auto_ini_ref}}{project_name}\.consolidated_results\.zip ${$args{auto_ini_ref}}{work_dir}\/Consolidated_results\/$input_file`;
      }
   }
   &hide_pbar_2;
   ${$args{progress_bar}}->configure(-label=>" ");
   ${$args{main_window}}->update;
   return (1);
}


sub error_message {
   my %args = @_;
   my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                 -text    => "Error while changing to directory $args{dir} for clearing results.\n Aborting.",
                                                 -buttons => ['OK'],
                                                 -bitmap  => 'info');
   $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
   $error_msg-> Show();
   open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
   print ERRORLOG "Error processing input files".
                  "\nError while changing to directory $args{dir} for clearing results.\n Aborting.\n\n";
   close ERRORLOG;
}

sub get_genemodel_fasta {
   my %args = (genbank_new => '0',
               @_);

   if (${$args{auto_ini_ref}}{internal_gm} == 1 || $args{genbank_new} == 1) {
      #generate gene models for all fasta entries
      ${$args{progress_bar}}->configure(-label=>"Generating gene models for fasta files");
      ${$args{main_window}}->update;
      my $status = &make_genemodel(main_window   => $args{main_window},
                                   progress_bar  => $args{progress_bar},
                                   auto_ini_ref  => $args{auto_ini_ref},
                                   ini_ref       => $args{ini_ref},
                                   input_file    => $args{input_file}
                                  );
      if ($status == 0) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while generating gene model for $args{input_file}.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nError while generating gene model for $args{input_file}.\n\n";
         close ERRORLOG;
         return (0); #return to main
      }
      #read gene model to memory model will be in hash reference
      ($orf_list_ref) = &read_internal_gene_model(main_window   => $args{main_window},
                                                  progress_bar  => $args{progress_bar},
                                                  auto_ini_ref  => $args{auto_ini_ref},
                                                  ini_ref       => $args{ini_ref},
                                                  input_file    => $args{input_file}
                                                 );
      if ($orf_list_ref eq 0) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while reading internal gene models.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nError while reading internal gene models for $args{input_file}.\n\n";
         close ERRORLOG;
         return (0); #return to main
      }
      ${$args{main_window}}->update;

      return ($orf_list_ref);

   } elsif (${$args{auto_ini_ref}}{external_gm} == 1) {
      #read external gene models to memory
      ($orf_list_ref) = &read_external_gene_model(main_window   => $args{main_window},
                                                  progress_bar  => $args{progress_bar},
                                                  auto_ini_ref  => $args{auto_ini_ref},
                                                  ini_ref       => $args{ini_ref},
                                                  input_file    => $args{input_file}
                                                 );
      if ($orf_list_ref eq 0) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while reading external gene model for $args{input_file}.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nError while reading external gene model for $args{input_file}.\n\n";
         close ERRORLOG;
         return (0); #return to main
      }
      ${$args{main_window}}->update;
      return ($orf_list_ref);

   }
}

sub from_fasta {
   my %args = (input_files => '', @_);
   my $orf_list_ref;
   my @combined_gm = ();

   #initialise progress bar
   &progress_bar_1(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'Processing Fasta files',
                   label        => 'FASTA'
                  );
   &show_pbar_1;


   #processing each fasta input file for gene model
   #max number of fasta entries
   my $max_count = $#fasta_input + 1;
   my $label = "Determining Gene models";
   $progress = 0;
   &update_pbar_1(title        => 'Gene model',
                  label        => $label
                 );

   foreach my $input_file (@{$args{input_list_ref}}) {
      #update progress bar
      $progress++;
      &update_pbar_1(title       => 'Gene model',
                     label       => "Determining gene model for $input_file",
                     progress    => ($progress / $max_count) * 100
                    );


      #generate gene models for all FASTA entries or read external gene model files
      my ($orf_list_ref) = &get_genemodel_fasta(main_window   => $args{main_window},
                                                progress_bar  => $args{progress_bar},
                                                auto_ini_ref  => $args{auto_ini_ref},
                                                ini_ref       => $args{ini_ref},
                                                input_file    => $input_file
                                                );
      if ($orf_list_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while generating and reading gene model for $input_file.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nError while generating and reading gene model for $input_file.\n\n";
         close ERRORLOG;
         &hide_pbar_1;
         return (0); #return to main
      }

      #merge gene models into one array
      @combined_gm = (@combined_gm, @{$orf_list_ref});

      my %seen = ();
      @combined_gm = grep { ! $seen{$_} ++ } @combined_gm;

      #undef $orf_list_ref;
   }
   &hide_pbar_1;
   ${$args{main_window}}->update;

   return (\@combined_gm, $args{input_list_ref});
}

sub run_programmes {
   my %args = @_;
   #run Blast analysis for each file if selected
   if (${$args{auto_ini_ref}}{blast_selector} == 1) {
      foreach my $input_file (@{$args{input_list_ref}}) {
         #reset status frame
         ${$args{progress_bar}}->configure(-label=>" Blasting ");
         ${$args{main_window}}->update;

         my ($blast_ref) = &blast_file(main_window      => $args{main_window},
                                       progress_bar     => $args{progress_bar},
                                       auto_ini_ref     => $args{auto_ini_ref},
                                       ini_ref          => $args{ini_ref},
                                       directory        => ${$args{ini_ref}}{input_files},
                                       input_file       => $input_file,
                                       combined_orf_ref => $args{gene_model_ref}
                                      );
         if ($blast_ref eq '0') {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Error while Blasting $input_file.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error processing input files".
                           "\nError while Blasting $input_file.\n\n";
            close ERRORLOG;
            #return (0); #return to main
            #keep analysis going - solve problem at re-run
         }
      }
   }

   #run COG analysis if selected
   if (${$args{auto_ini_ref}}{COG_selector} == 1) {
      undef $COGcode_to_number;
      undef $COGcode_to_header;
      undef $COGcode_to_letter;
      undef $COGcode_to_phylo;
      undef $COGletter_to_family;
      undef $COG2008genomes;
      undef $COG2008_to_def;
      undef $COG2008_to_acc;
      undef $arCOGcode_to_COG;
      undef $arCOGacc_to_org;
      undef $COG2014genomes;
      undef $COG2014_to_def;
      undef $COG2014_to_acc;
      undef $COG2014_to_refseq;
      undef $arCOGcode2014_to_COG;
      undef $arCOGacc2014_to_org;
      undef $POG2013_to_gi;
      undef $COGletter_to_family;

      #reading translation tables for original COG
      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
         ${$args{progress_bar}}->configure(-label=>" Parsing COG2003 reference data ");
         ${$args{main_window}}->update;
         &translate_COG2003(main_window         => $args{main_window},
                            progress_bar        => $args{progress_bar},
                            auto_ini_ref        => $args{auto_ini_ref},
                            ini_ref             => $args{ini_ref},
                            COGcode_to_number   => \$COGcode_to_number,
                            COGcode_to_header   => \$COGcode_to_header,
                            COGcode_to_letter   => \$COGcode_to_letter,
                            COGcode_to_phylo    => \$COGcode_to_phylo,
                            COGletter_to_family => \$COGletter_to_family
                           );

      }

      #get additional codes for arCOG if selected
      if (${$args{auto_ini_ref}}{COG_type} eq 'arCOG') {
         ${$args{progress_bar}}->configure(-label=>" Parsing arCOG reference data ");
         ${$args{main_window}}->update;
         &translate_arCOG(main_window         => $args{main_window},
                          progress_bar        => $args{progress_bar},
                          auto_ini_ref        => $args{auto_ini_ref},
                          ini_ref             => $args{ini_ref},
                          arCOGcode_to_COG    => \$arCOGcode_to_COG,
                          arCOGacc_to_org     => \$arCOGacc_to_org,
                         );
      }

      #get additional codes for arCOG2014 if selected
      if (${$args{auto_ini_ref}}{COG_type} eq 'arCOG2014') {
         ${$args{progress_bar}}->configure(-label=>" Parsing arCOG2014 reference data ");
         ${$args{main_window}}->update;
         &translate_arCOG2014(main_window          => $args{main_window},
                              progress_bar         => $args{progress_bar},
                              auto_ini_ref         => $args{auto_ini_ref},
                              ini_ref              => $args{ini_ref},
                              arCOGcode2014_to_COG => \$arCOGcode2014_to_COG,
                              arCOGacc2014_to_org  => \$arCOGacc2014_to_org,
                              COGletter_to_family  => \$COGletter_to_family,
                             );
      }

      #get additional codes for COG2008 is selected
      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2008') {
         ${$args{progress_bar}}->configure(-label=>" Parsing COG2008 reference data ");
         ${$args{main_window}}->update;
         &translate_COG2008(main_window         => $args{main_window},
                            progress_bar        => $args{progress_bar},
                            auto_ini_ref        => $args{auto_ini_ref},
                            ini_ref             => $args{ini_ref},
                            COG2008genomes      => \$COG2008genomes,
                            COG2008_to_def      => \$COG2008_to_def,
                            COG2008_to_acc      => \$COG2008_to_acc,
                           );
      }

      #get additional codes for COG2014 is selected
      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2014') {
         ${$args{progress_bar}}->configure(-label=>" Parsing COG2014 reference data ");
         ${$args{main_window}}->update;
         &translate_COG2014(main_window         => $args{main_window},
                            progress_bar        => $args{progress_bar},
                            auto_ini_ref        => $args{auto_ini_ref},
                            ini_ref             => $args{ini_ref},
                            COG2014genomes      => \$COG2014genomes,
                            COG2014_to_def      => \$COG2014_to_def,
                            COG2014_to_acc      => \$COG2014_to_acc,
                            COGletter_to_family => \$COGletter_to_family,
                            COG2014_to_refseq   => \$COG2014_to_refseq,
                           );
      }

      #get additional codes for POG2013 is selected
      if (${$args{auto_ini_ref}}{COG_type} eq 'POG2013') {
         ${$args{progress_bar}}->configure(-label=>" Parsing POG2013 reference data ");
         ${$args{main_window}}->update;
         &translate_POG2013(main_window         => $args{main_window},
                            progress_bar        => $args{progress_bar},
                            auto_ini_ref        => $args{auto_ini_ref},
                            ini_ref             => $args{ini_ref},
                            POG2013_to_gi       => \$POG2013_to_gi,
                           );
      }

      foreach my $input_file (@{$args{input_list_ref}}) {
         #reset status frame
         ${$args{progress_bar}}->configure(-label=>" COG Blast ");
         ${$args{main_window}}->update;

         my ($COG_ref) = &COG_file(main_window      => $args{main_window},
                                   progress_bar     => $args{progress_bar},
                                   auto_ini_ref     => $args{auto_ini_ref},
                                   ini_ref          => $args{ini_ref},
                                   directory        => ${$args{ini_ref}}{input_files},
                                   input_file       => $input_file,
                                   combined_orf_ref => $args{gene_model_ref},
                                   COGcode_to_number    => \$COGcode_to_number,
                                   COGcode_to_header    => \$COGcode_to_header,
                                   COGcode_to_letter    => \$COGcode_to_letter,
                                   COGcode_to_phylo     => \$COGcode_to_phylo,
                                   COGletter_to_family  => \$COGletter_to_family,
                                   COG2008genomes       => \$COG2008genomes,
                                   COG2008_to_def       => \$COG2008_to_def,
                                   COG2008_to_acc       => \$COG2008_to_acc,
                                   arCOGcode_to_COG     => \$arCOGcode_to_COG,
                                   arCOGacc_to_org      => \$arCOGacc_to_org,
                                   COG2014genomes       => \$COG2014genomes,
                                   COG2014_to_def       => \$COG2014_to_def,
                                   COG2014_to_acc       => \$COG2014_to_acc,
                                   COG2014_to_refseq    => \$COG2014_to_refseq,
                                   arCOGcode2014_to_COG => \$arCOGcode2014_to_COG,
                                   arCOGacc2014_to_org  => \$arCOGacc2014_to_org,
                                   POG2013_to_gi        => \$POG2013_to_gi,
                                   COGletter_to_family  => \$COGletter_to_family,
                                  );
         if ($COG_ref eq '0') {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Error while Blasting $input_file with COG database.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error processing input files".
                           "\nError while Blasting $input_file with COG database.\n\n";
            close ERRORLOG;
            #return (0); #return to main
            #keep analysis going - solve problem at re-run
         }
      }
   }

   #run PFam analysis if selected
   if (${$args{auto_ini_ref}}{Pfam_selector} == 1) {
      #test if PFam3 databases are in binary format
      if (${$args{auto_ini_ref}}{use_Pfam3} == 1) {
         my @local_db = split/ ; /,${$args{auto_ini_ref}}{full_Pfam_db};
         foreach my $db (@local_db) {
            unless (-e $db.'.h3m') {
               ${$args{progress_bar}}->configure(-label=>" Pressing Pfam database into binary format ");
               ${$args{main_window}}->update;
               system "${$args{ini_ref}}{pfam3_executable}/hmmpress $db";
            }
         }
      }

      foreach my $input_file (@{$args{input_list_ref}}) {
         #reset status frame
         ${$args{progress_bar}}->configure(-label=>" PFam ");
         ${$args{main_window}}->update;

         my ($pfam_ref) = &Pfam_file(main_window      => $args{main_window},
                                     progress_bar     => $args{progress_bar},
                                     auto_ini_ref     => $args{auto_ini_ref},
                                     ini_ref          => $args{ini_ref},
                                     directory        => ${$args{ini_ref}}{input_files},
                                     input_file       => $input_file,
                                     combined_orf_ref => $args{gene_model_ref}
                                    );
         if ($pfam_ref eq '0') {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Error while analysing $input_file for PFams.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error processing input files".
                           "\nError while analysing $input_file for PFams.\n\n";
            close ERRORLOG;
            #return (0); #return to main
            #keep analysis going - solve problem at re-run
         }
      }
   }

   #run TIGRFam analysis if selected
   if (${$args{auto_ini_ref}}{TIGRfam_selector} == 1) {
      #sort database order
      my @tigr_db = split/ ; /,${$args{auto_ini_ref}}{full_TIGRfam_db};
      @tigr_db = sort {$b cmp $a} @tigr_db;
      ${$args{auto_ini_ref}}{full_TIGRfam_db} = join (' ; ', @tigr_db);
      undef @tigr_db;

      #process each file
      foreach my $input_file (@{$args{input_list_ref}}) {
         #reset status frame
         ${$args{progress_bar}}->configure(-label=>" TIGRfam ");
         ${$args{main_window}}->update;

         my ($tigrfam_ref) = &TIGRfam_file(main_window      => $args{main_window},
                                           progress_bar     => $args{progress_bar},
                                           auto_ini_ref     => $args{auto_ini_ref},
                                           ini_ref          => $args{ini_ref},
                                           directory        => ${$args{ini_ref}}{input_files},
                                           input_file       => $input_file,
                                           databases        => ${$args{auto_ini_ref}}{full_TIGRfam_db},
                                           output_folder    => ${$args{ini_ref}}{TIGRfam_results},
                                           combined_orf_ref => $args{gene_model_ref}
                                          );
         if ($tigrfam_ref eq '0') {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Error while analysing $input_file for TIGRfams.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error processing input files".
                           "\nError while analysing $input_file for TIGRfams.\n\n";
            close ERRORLOG;
            #return (0); #return to main
            #keep analysis going - solve problem at re-run
         }
      }
   }

   #run tRNA analysis if selected
   if (${$args{auto_ini_ref}}{trna_selector} == 1) {
      #tRNAscan runs on whole sequence, hence multithreading is done on input files and not on gene models
      #process each input file is perfomerd in subroutine
      ${$args{progress_bar}}->configure(-label=>" tRNA ");
      ${$args{main_window}}->update;

      my ($trna_ref) = &trna(main_window      => $args{main_window},
                             progress_bar     => $args{progress_bar},
                             auto_ini_ref     => $args{auto_ini_ref},
                             ini_ref          => $args{ini_ref},
                             directory        => ${$args{ini_ref}}{input_files},
                             input_array_ref  => $args{input_list_ref},
                             input_type       => 'fasta'
                            );
      if ($trna_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while identifying tRNA structures.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nError while identifying tRNA structures.\n\n";
         close ERRORLOG;
         #return (0); #return to main
         #keep analysis going - solve problem at re-run
      }
   }

   #run SignalP analysis if selected
   if (${$args{auto_ini_ref}}{signalp_selector} == 1) {
      #process each file
      foreach my $input_file (@{$args{input_list_ref}}) {
         #reset status frame
         ${$args{progress_bar}}->configure(-label=>" SignalP ");
         ${$args{main_window}}->update;

         my ($signalp_ref) = &signalp(main_window      => $args{main_window},
                                      progress_bar     => $args{progress_bar},
                                      auto_ini_ref     => $args{auto_ini_ref},
                                      ini_ref          => $args{ini_ref},
                                      directory        => ${$args{ini_ref}}{input_files},
                                      input_file       => $input_file,
                                      combined_orf_ref => $args{gene_model_ref}
                                    );
         if ($signalp_ref eq '0') {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Error while analysing $input_file for Signal cleavage site.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error processing input files".
                           "\nError while analysing $input_file for Signal cleavage site.\n\n";
            close ERRORLOG;
            #return (0); #return to main
            #keep analysis going - solve problem at re-run
         }
      }
      #create summary file for ALL input files
      my ($signalp_ref) = &signalp_summary(main_window      => $args{main_window},
                                           progress_bar     => $args{progress_bar},
                                           auto_ini_ref     => $args{auto_ini_ref},
                                           ini_ref          => $args{ini_ref},

                                          );
      if ($signalp_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while creating SignalP summary file.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nError while creating SignalP summary file.\n\n";
         close ERRORLOG;
         #return (0); #return to main
         #keep analysis going - solve problem at re-run
      }
   }

   #run TMHMM analysis if selected
   if (${$args{auto_ini_ref}}{tmhmm_selector} == 1) {
      #process each file
      foreach my $input_file (@{$args{input_list_ref}}) {
         #reset status frame
         ${$args{progress_bar}}->configure(-label=>" TMHMM ");
         ${$args{main_window}}->update;

         my ($tmhmm_ref) = &tmhmm_file(main_window      => $args{main_window},
                                       progress_bar     => $args{progress_bar},
                                       auto_ini_ref     => $args{auto_ini_ref},
                                       ini_ref          => $args{ini_ref},
                                       directory        => ${$args{ini_ref}}{input_files},
                                       input_file       => $input_file,
                                       combined_orf_ref => $args{gene_model_ref}
                                      );
         if ($tmhmm_ref eq '0') {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Error while analysing $input_file for transmembrane helices.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error processing input files".
                           "\nError while analysing $input_file for transmembrane helices.\n\n";
            close ERRORLOG;
            #return (0); #return to main
            #keep analysis going - solve problem at re-run
         }
      }
   }

   #run Transterm analysis if selected
   if (${$args{auto_ini_ref}}{terminator_selector} == 1) {
      #Transterm runs on whole sequence, hence multithreading is done on input files and not on gene models
      #process each input file is perfomerd in subroutine
      ${$args{progress_bar}}->configure(-label=>" Transterm ");
      ${$args{main_window}}->update;

      my ($transterm_ref) = &transterm(main_window      => $args{main_window},
                                       progress_bar     => $args{progress_bar},
                                       auto_ini_ref     => $args{auto_ini_ref},
                                       ini_ref          => $args{ini_ref},
                                       directory        => ${$args{ini_ref}}{input_files},
                                       input_array_ref  => $args{input_list_ref},
                                       combined_orf_ref => $args{gene_model_ref},
                                       input_type       => 'fasta'
                                      );
      if ($transterm_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while identifying terminator structures.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nError while identifying terminator structures.\n\n";
         close ERRORLOG;
         #return (0); #return to main
         #keep analysis going - solve problem at re-run
      }
   }

   #run rRNA analysis if selected
   if (${$args{auto_ini_ref}}{rrna_selector} == 1) {
      #rRNA Blast runs on whole sequence, hence multithreading is done on input files and not on gene models
      #process each input file is perfomerd in subroutine
      ${$args{progress_bar}}->configure(-label=>" ribosomal RNA ");
      ${$args{main_window}}->update;

      #remove summary file (will be re-established during parsing)
      unlink "${$args{ini_ref}}{rrna_results}\/rRNA_summary.txt";

      my ($rrna_ref) = &rrna(main_window      => $args{main_window},
                             progress_bar     => $args{progress_bar},
                             auto_ini_ref     => $args{auto_ini_ref},
                             ini_ref          => $args{ini_ref},
                             directory        => ${$args{ini_ref}}{input_files},
                             input_array_ref  => $args{input_list_ref},
                             combined_orf_ref => $args{gene_model_ref},
                             input_type       => 'fasta'
                            );
      if ($rrna_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while identifying ribosomal RNA structures.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nError while identifying ribosomal RNA structures.\n\n";
         close ERRORLOG;
         #return (0); #return to main
         #keep analysis going - solve problem at re-run
      }
   }

   #run vector screen analysis if selected
   if (${$args{auto_ini_ref}}{vector_screen_selector} == 1) {
      #vector screen runs on whole sequence, hence multithreading is done on input files and not on gene models
      #process each input file is perfomerd in subroutine
      ${$args{progress_bar}}->configure(-label=>" Vector Screen ");
      ${$args{main_window}}->update;

      my ($vector_screen_ref) = &vector_screen(main_window      => $args{main_window},
                                               progress_bar     => $args{progress_bar},
                                               auto_ini_ref     => $args{auto_ini_ref},
                                               ini_ref          => $args{ini_ref},
                                               directory        => ${$args{ini_ref}}{input_files},
                                               input_list_ref  => $args{input_list_ref},
                                               combined_orf_ref => $args{gene_model_ref},
                                               input_type       => 'fasta'
                                              );
      if ($vector_screen_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while scrrening for vector sequences in input files.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nError while screening for vector sequences in input files.\n\n";
         close ERRORLOG;
         #return (0); #return to main
         #keep analysis going - solve problem at re-run
      }
   }

   #run CRISPR analysis if selected
   if (${$args{auto_ini_ref}}{CRISPR_selector} == 1) {
      foreach my $input_file (@{$args{input_list_ref}}) {
         #reset status frame
         ${$args{progress_bar}}->configure(-label=>" CRISPR ");
         ${$args{main_window}}->update;

         my ($CRISPR_status) = &CRISPR(main_window      => $args{main_window},
                                       progress_bar     => $args{progress_bar},
                                       auto_ini_ref     => $args{auto_ini_ref},
                                       ini_ref          => $args{ini_ref},
                                       directory        => ${$args{ini_ref}}{input_files},
                                       input_array_ref  => $args{input_list_ref},
                                      );
         if ($CRISPR_status eq '0') {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Error while finding CRISPRs for $input_file.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error processing input files".
                           "\nError while finding CRISPRs for $input_file.\n\n";
            close ERRORLOG;
            #return (0); #return to main
            #keep analysis going - solve problem at re-run
         }
      }
   }

   #run non-coding RNA analysis if selected
   if (${$args{auto_ini_ref}}{ncrna_selector} == 1) {
      #Infernal runs on whole sequence, hence multithreading is done on input files and not on gene models
      #process each input file is perfomerd in subroutine
      ${$args{progress_bar}}->configure(-label=>" non-coding RNA ");
      ${$args{main_window}}->update;

      my ($ncrna_ref) = &ncrna(main_window      => $args{main_window},
                               progress_bar     => $args{progress_bar},
                               auto_ini_ref     => $args{auto_ini_ref},
                               ini_ref          => $args{ini_ref},
                               directory        => ${$args{ini_ref}}{input_files},
                               input_array_ref  => $args{input_list_ref},
                               combined_orf_ref => $args{gene_model_ref},
                               input_type       => 'fasta'
                              );
      if ($ncrna_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while identifying non-coding RNA structures.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nError while identifying non-coding RNA structures.\n\n";
         close ERRORLOG;
         #return (0); #return to main
         #keep analysis going - solve problem at re-run
      }
   }

}

sub run_compiler {
   my %args = @_;
   my ($file_ref, @list, $label, $max_length, $orf_to_id, $id_to_orf, $max_count, $counter,
       @COGs, $COGcode_to_number, $COGcode_to_header, $COGcode_to_letter, $COGletter_to_family, $COGcode_to_phylo,
       $Pfamcode_to_descriptor, $Pfamcode_to_number, $Pfamcode_to_access, $Pfamcode_to_family,
       $IPaccess_to_code, $IPaccess_to_descriptor,
       $TIGRdomain, $TIGRrole_number, $TIGRrole_name, $TIGRrole_GO, $TIGRrole_GoClass,
       $rfam_model);

   #create status box
   &progress_bar_1(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'Compiling features',
                   label        => 'Compiling features'
                  );
   &show_pbar_1;

   #reading translation tables for original COG
   if (${$args{auto_ini_ref}}{COG_selector} == 1 && ${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
      &update_pbar_1(title        => 'Compiling features',
                     label        => 'Translating COG'
                    );
      &translate_COG2003(main_window         => $args{main_window},
                         progress_bar        => $args{progress_bar},
                         auto_ini_ref        => $args{auto_ini_ref},
                         ini_ref             => $args{ini_ref},
                         COGcode_to_number   => \$COGcode_to_number,
                         COGcode_to_header   => \$COGcode_to_header,
                         COGcode_to_letter   => \$COGcode_to_letter,
                         COGcode_to_phylo    => \$COGcode_to_phylo,
                         COGletter_to_family => \$COGletter_to_family
                        );

   }

   #reading translation tables for Pfam and Interpro
   if (${$args{auto_ini_ref}}{Pfam_selector} == 1) {
      &update_pbar_1(title        => 'Compiling features',
                     label        => 'Translating PFam'
                    );
      my $HMMERver = 3;
      if (${$args{auto_ini_ref}}{use_Pfam2} == 1) {
         $HMMERver = 2;
      } else {
         $HMMERver = 3;
      }
      &translate_pfam(main_window            => $args{main_window},
                      progress_bar           => $args{progress_bar},
                      auto_ini_ref           => $args{auto_ini_ref},
                      ini_ref                => $args{ini_ref},
                      Pfamcode_to_descriptor => \$Pfamcode_to_descriptor,
                      Pfamcode_to_number     => \$Pfamcode_to_number,
                      Pfamcode_to_access     => \$Pfamcode_to_access,
                      Pfamcode_to_family     => \$Pfamcode_to_family,
                      HMMERver               => $HMMERver
                     );

      &update_pbar_1(title        => 'Compiling features',
                     label        => 'Translating Interpro'
                    );
      &translate_interpro(main_window            => $args{main_window},
                          progress_bar           => $args{progress_bar},
                          auto_ini_ref           => $args{auto_ini_ref},
                          ini_ref                => $args{ini_ref},
                          IPaccess_to_code       => \$IPaccess_to_code,
                          IPaccess_to_descriptor => \$IPaccess_to_descriptor,
                         );
   }

   #reading TIGRfam translation if selected
   if (${$args{auto_ini_ref}}{TIGRfam_selector} == 1) {
      &update_pbar_1(title        => 'Compiling features',
                     label        => 'Translating TIGRfam'
                    );
      &translate_tigrfam(main_window            => $args{main_window},
                         progress_bar           => $args{progress_bar},
                         auto_ini_ref           => $args{auto_ini_ref},
                         ini_ref                => $args{ini_ref},
                         TIGRdomain             => \$TIGRdomain,
                         TIGRrole_number        => \$TIGRrole_number,
                         TIGRrole_name          => \$TIGRrole_name,
                         TIGRrole_GO            => \$TIGRrole_GO,
                         TIGRrole_GoClass       => \$TIGRrole_GoClass
                        );
   }

   #reading non-coding RNA descriptors if selected
   if (${$args{auto_ini_ref}}{ncrna_selector} == 1) {
      &update_pbar_1(title        => 'Compiling features',
                     label        => 'Translating Rfam models'
                    );
      &translate_rfam(main_window            => $args{main_window},
                      progress_bar           => $args{progress_bar},
                      auto_ini_ref           => $args{auto_ini_ref},
                      ini_ref                => $args{ini_ref},
                      rfam_model             => \$rfam_model
                     );
   }

   #test for Consolidated directory, preparation for sorting data if selected
   if (-d ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results') {
      &update_pbar_1(title        => 'Compiling features',
                     label        => 'Cleaning up Consolidated Results Folder'
                    );
      &cleanup(${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results');
   }

   #create initial directory structure if not already exists
   if (!-d ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results') {
      mkdir ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results' or die "Error creating sorted directory \n";
   }

   #setup Blast summary file if selected
   if (${$args{auto_ini_ref}}{blast_summary} == 1) {
      &update_pbar_1(title        => 'Compiling features',
                     label        => 'Creating Blast summary'
                    );
      &blast_summary(main_window            => $args{main_window},
                     progress_bar           => $args{progress_bar},
                     auto_ini_ref           => $args{auto_ini_ref},
                     ini_ref                => $args{ini_ref},
                     input_list_ref         => $args{input_list_ref}
                    );

   }

   #compile all results foreach input file
   $max_count = $#{$args{input_list_ref}} + 1;
   $counter = 0;
   foreach my $input_file (@{$args{input_list_ref}}) {
      my (%genbank, @feature_list); #local Genbank entries
                                    #these will be populated by individual parsing modules

      ${$args{progress_bar}}->configure(-label=>" Compiling results for $input_file");
      ${$args{main_window}}->update;
      $counter++;
      &update_pbar_1(title        => "Compiling features",
                     label        => "Compiling results for $input_file",
                     progress     => ($counter / $max_count) * 100,
                    );

      #map correct gene model entries for respective input file
      @list = grep /$input_file\___/, @{$args{gene_model_ref}};

      #create proper ORF numbers if necessary (sorted by left_bd) and associate with ID list
      $counter = 1;
      my @sorted =
                  map  $_->[0] =>
                  sort { $a->[1] <=> $b->[1] }
                  map  [ $_, m/^$input_file\___\d+___(\d+)___/ ]
                  => @list;
      foreach my $entry (@sorted) {
         my ($name, $id);
         'reset' =~ /reset/;
         $entry =~ m/^(.+?)___(\d+)___/;
         $name = $1;
         $id = $2;
         $orf_to_id->{$counter} = $entry;
         $id_to_orf->{$name.'_'.$id} = $counter;
         $counter++;
      }
      undef @sorted;
      undef @list;

      #run Blast parser for each file if selected
      if (${$args{auto_ini_ref}}{blast_selector} == 1) {
         &blast_parser(main_window      => $args{main_window},
                       progress_bar     => $args{progress_bar},
                       auto_ini_ref     => $args{auto_ini_ref},
                       ini_ref          => $args{ini_ref},
                       genbank_ref      => \%genbank,
                       feature_list_ref => \@feature_list,
                       orf_to_id        => $orf_to_id,
                       counter          => $counter
                      );
      } elsif (${$args{auto_ini_ref}}{blast_selector}  == 0 &&
               ${$args{auto_ini_ref}}{gene_model_only} == 1) {# &&
               #${$args{auto_ini_ref}}{create_genbank}  == 1) {
         &gene_model_parser(main_window      => $args{main_window},
                            progress_bar     => $args{progress_bar},
                            auto_ini_ref     => $args{auto_ini_ref},
                            ini_ref          => $args{ini_ref},
                            genbank_ref      => \%genbank,
                            feature_list_ref => \@feature_list,
                            orf_to_id        => $orf_to_id,
                            counter          => $counter,
                            input_file       => $input_file
                           );
      }

      #run COG parser for each file if selected
      if (${$args{auto_ini_ref}}{COG_selector} == 1) {
         &cog_parser(main_window          => $args{main_window},
                     progress_bar         => $args{progress_bar},
                     auto_ini_ref         => $args{auto_ini_ref},
                     ini_ref              => $args{ini_ref},
                     genbank_ref          => \%genbank,
                     feature_list_ref     => \@feature_list,
                     orf_to_id            => $orf_to_id,
                     counter              => $counter,
                     COGcode_to_number    => \$COGcode_to_number,
                     COGcode_to_header    => \$COGcode_to_header,
                     COGcode_to_letter    => \$COGcode_to_letter,
                     COGcode_to_phylo     => \$COGcode_to_phylo,
                     COGletter_to_family  => \$COGletter_to_family,
                     arCOGcode_to_COG     => \$arCOGcode_to_COG,
                     arCOGacc_to_org      => \$arCOGacc_to_org,
                     COG2008genomes       => \$COG2008genomes,
                     COG2008_to_def       => \$COG2008_to_def,
                     COG2008_to_acc       => \$COG2008_to_acc,
                     COG2014genomes       => \$COG2014genomes,
                     COG2014_to_def       => \$COG2014_to_def,
                     COG2014_to_acc       => \$COG2014_to_acc,
                     COG2014_to_refseq    => \$COG2014_to_refseq,
                     arCOGcode2014_to_COG => \$arCOGcode2014_to_COG,
                     arCOGacc2014_to_org  => \$arCOGacc2014_to_org,
                     POG2013_to_gi        => \$POG2013_to_gi,
                     COGletter_to_family  => \$COGletter_to_family,
                    );
      }

      #run PFam parser for each file if selected
      if (${$args{auto_ini_ref}}{Pfam_selector} == 1) {
         if (${$args{auto_ini_ref}}{use_Pfam2} == 1) {
            &pfam_parser(main_window      => $args{main_window},
                         progress_bar     => $args{progress_bar},
                         auto_ini_ref     => $args{auto_ini_ref},
                         ini_ref          => $args{ini_ref},
                         genbank_ref      => \%genbank,
                         feature_list_ref => \@feature_list,
                         orf_to_id        => $orf_to_id,
                         counter          => $counter,
                         Pfamcode_to_descriptor => \$Pfamcode_to_descriptor,
                         Pfamcode_to_number     => \$Pfamcode_to_number,
                         Pfamcode_to_access     => \$Pfamcode_to_access,
                         Pfamcode_to_family     => \$Pfamcode_to_family,
                         IPaccess_to_code       => \$IPaccess_to_code,
                         IPaccess_to_descriptor => \$IPaccess_to_descriptor
                        );
         } else {
            &pfam3_parser(main_window      => $args{main_window},
                          progress_bar     => $args{progress_bar},
                          auto_ini_ref     => $args{auto_ini_ref},
                          ini_ref          => $args{ini_ref},
                          genbank_ref      => \%genbank,
                          feature_list_ref => \@feature_list,
                          orf_to_id        => $orf_to_id,
                          counter          => $counter,
                          Pfamcode_to_descriptor => \$Pfamcode_to_descriptor,
                          Pfamcode_to_number     => \$Pfamcode_to_number,
                          Pfamcode_to_access     => \$Pfamcode_to_access,
                          Pfamcode_to_family     => \$Pfamcode_to_family,
                          IPaccess_to_code       => \$IPaccess_to_code,
                          IPaccess_to_descriptor => \$IPaccess_to_descriptor
                         );
         }
      }

      #run TIGRfam parser for each file if selected
      if (${$args{auto_ini_ref}}{TIGRfam_selector} == 1) {
         #Test version of TIGRfam databases. Version > 10 requires Hmmer3
         my $require_hmmer3 = 0;
         my @TIGRfam_db = split /\s+\;\s+/, ${$args{auto_ini_ref}}{TIGRfam_db};

         #test if hmmer3 is required (TIGRfam ver > 10)
         if ($TIGRfam_db[0] =~ m/TIGRFAMs_\d+/) {
            $TIGRfam_db[0] =~ m/TIGRFAMs_(\d+)/;
            if ($1 >=10) {
               $require_hmmer3 = 1;
            }
         }
         if ($require_hmmer3 == 0) {
            &tigrfam_parser(main_window      => $args{main_window},
                            progress_bar     => $args{progress_bar},
                            auto_ini_ref     => $args{auto_ini_ref},
                            ini_ref          => $args{ini_ref},
                            input_file       => $input_file,
                            genbank_ref      => \%genbank,
                            feature_list_ref => \@feature_list,
                            gene_model_ref   => $args{gene_model_ref},
                            orf_to_id        => $orf_to_id,
                            counter          => $counter,
                            TIGRdomain       => \$TIGRdomain,
                            TIGRrole_number  => \$TIGRrole_number,
                            TIGRrole_name    => \$TIGRrole_name,
                            TIGRrole_GO      => \$TIGRrole_GO,
                            TIGRrole_GoClass => \$TIGRrole_GoClass
                           );
         } else {
            &tigrfam3_parser(main_window      => $args{main_window},
                             progress_bar     => $args{progress_bar},
                             auto_ini_ref     => $args{auto_ini_ref},
                             ini_ref          => $args{ini_ref},
                             input_file       => $input_file,
                             genbank_ref      => \%genbank,
                             feature_list_ref => \@feature_list,
                             gene_model_ref   => $args{gene_model_ref},
                             orf_to_id        => $orf_to_id,
                             counter          => $counter,
                             TIGRdomain       => \$TIGRdomain,
                             TIGRrole_number  => \$TIGRrole_number,
                             TIGRrole_name    => \$TIGRrole_name,
                             TIGRrole_GO      => \$TIGRrole_GO,
                             TIGRrole_GoClass => \$TIGRrole_GoClass
                            );
         }
      }

      #run SignalP parser for each file if selected
      if (${$args{auto_ini_ref}}{signalp_selector} == 1) {
         &signalp_parser(main_window      => $args{main_window},
                         progress_bar     => $args{progress_bar},
                         auto_ini_ref     => $args{auto_ini_ref},
                         ini_ref          => $args{ini_ref},
                         genbank_ref      => \%genbank,
                         feature_list_ref => \@feature_list,
                         filename         => $input_file,
                         counter          => $counter
                         );
      }

      #run TMHMM parser for each file if selected
      if (${$args{auto_ini_ref}}{tmhmm_selector} == 1) {
         &tmhmm_parser(main_window      => $args{main_window},
                       progress_bar     => $args{progress_bar},
                       auto_ini_ref     => $args{auto_ini_ref},
                       ini_ref          => $args{ini_ref},
                       genbank_ref      => \%genbank,
                       feature_list_ref => \@feature_list,
                       filename         => $input_file,
                       counter          => $counter
                       );
      }

      #run CRISPR parser for each file if selected
      if (${$args{auto_ini_ref}}{CRISPR_selector} == 1) {
         &CRISPR_parser(main_window      => $args{main_window},
                        progress_bar     => $args{progress_bar},
                        auto_ini_ref     => $args{auto_ini_ref},
                        ini_ref          => $args{ini_ref},
                        genbank_ref      => \%genbank,
                        feature_list_ref => \@feature_list,
                        filename         => $input_file,
                        counter          => $counter
                        );
      }

      #run tRNAscan parser for each file if selected
      if (${$args{auto_ini_ref}}{trna_selector} == 1) {
         &trnascan_parser(main_window      => $args{main_window},
                          progress_bar     => $args{progress_bar},
                          auto_ini_ref     => $args{auto_ini_ref},
                          ini_ref          => $args{ini_ref},
                          genbank_ref      => \%genbank,
                          feature_list_ref => \@feature_list,
                          filename         => $input_file,
                          counter          => $counter
                          );
      }

      #run transterm parser for each file if selected
      if (${$args{auto_ini_ref}}{terminator_selector} == 1) {
         &transterm_parser(main_window      => $args{main_window},
                           progress_bar     => $args{progress_bar},
                           auto_ini_ref     => $args{auto_ini_ref},
                           ini_ref          => $args{ini_ref},
                           genbank_ref      => \%genbank,
                           feature_list_ref => \@feature_list,
                           filename         => $input_file,
                           counter          => $counter
                           );
      }

      #run ncRNA parser for each file if selected
      if (${$args{auto_ini_ref}}{ncrna_selector} == 1) {
         &ncrna_parser(main_window      => $args{main_window},
                       progress_bar     => $args{progress_bar},
                       auto_ini_ref     => $args{auto_ini_ref},
                       ini_ref          => $args{ini_ref},
                       genbank_ref      => \%genbank,
                       feature_list_ref => \@feature_list,
                       rfam_model       => \$rfam_model,
                       filename         => $input_file,
                       counter          => $counter
                      );
      }

      #run rRNA parser for each file if selected
      if (${$args{auto_ini_ref}}{rrna_selector} == 1) {
         &rrna_parser(main_window      => $args{main_window},
                      progress_bar     => $args{progress_bar},
                      auto_ini_ref     => $args{auto_ini_ref},
                      ini_ref          => $args{ini_ref},
                      genbank_ref      => \%genbank,
                      feature_list_ref => \@feature_list,
                      rfam_model       => \$rfam_model,
                      filename         => $input_file,
                      counter          => $counter
                      );
      }

      #run Vector screen parser for each file if selected
      if (${$args{auto_ini_ref}}{vector_screen_selector} == 1) {
         &vector_parser(main_window      => $args{main_window},
                        progress_bar     => $args{progress_bar},
                        auto_ini_ref     => $args{auto_ini_ref},
                        ini_ref          => $args{ini_ref},
                        genbank_ref      => \%genbank,
                        feature_list_ref => \@feature_list,
                        filename         => $input_file,
                        counter          => $counter
                        );
      }

      #run RBSfinder parser for each file if selected
      if (${$args{auto_ini_ref}}{runRBSfinder} == 1) {
         &rbs_parser(main_window      => $args{main_window},
                     progress_bar     => $args{progress_bar},
                     auto_ini_ref     => $args{auto_ini_ref},
                     ini_ref          => $args{ini_ref},
                     genbank_ref      => \%genbank,
                     feature_list_ref => \@feature_list,
                     filename         => $input_file,
                     counter          => $counter
                     );
      }

      #build Genbank files -> new module to create header, make standard header, and submit to gb_compiler
      my $local_file = $input_file;
      $local_file =~ s/_GAMOLAdna$//;
      &compile_gb(main_window      => $args{main_window},
                  progress_bar     => $args{progress_bar},
                  auto_ini_ref     => $args{auto_ini_ref},
                  ini_ref          => $args{ini_ref},
                  input_file       => $input_file,
                  genbank_ref      => \%genbank,
                  feature_list_ref => \@feature_list,
                  genbank_header   => $args{genbank_header}->{$local_file},
                  genbank_source   => $args{genbank_source}->{$local_file}
                 );

      #sort results into individual folders if selected
      if (${$args{auto_ini_ref}}{sort_results} == 1) {
         ${$args{progress_bar}}->configure(-label=>" Sorting results ");
         ${$args{main_window}}->update;
         &sort_results (main_window      => $args{main_window},
                        progress_bar     => $args{progress_bar},
                        auto_ini_ref     => $args{auto_ini_ref},
                        ini_ref          => $args{ini_ref},
                        input_file       => $input_file,
                        orf_to_id        => $orf_to_id,
                        id_to_orf        => $id_to_orf,
                        max_ORF_number   => $counter
                       );

         #copy Blast result file to Consolidated Results folder
         opendir READ, ${$args{ini_ref}}{blast_results};
         my @summary = grep /^Blast_Summary_/, readdir (READ);
         closedir READ;
         foreach my $file (@summary) {
            copy(${$args{ini_ref}}{blast_results}.'/'.$file, ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$file);
         }
         undef @summary;
      }

      #reset association hash
      undef $orf_to_id;
   }

   #free memory
   undef $COGcode_to_number;
   undef $COGcode_to_header;
   undef $COGcode_to_letter;
   undef $COGcode_to_phylo;
   undef $COGletter_to_family;
   undef $Pfamcode_to_descriptor;
   undef $Pfamcode_to_number;
   undef $Pfamcode_to_access;
   undef $Pfamcode_to_family;
   undef $IPaccess_to_code;
   undef $IPaccess_to_descriptor;
   undef $TIGRdomain;
   undef $TIGRrole_number;
   undef $TIGRrole_name;
   undef $TIGRrole_GO;
   undef $TIGRrole_GoClass;

   &hide_pbar_1;
   return($counter);
}

sub from_gb {
   my %args = @_;
   #processing each Genbank file individually
   my @corenames = ();
   my @combined_gm = ();

   #initialise progress bar
   &progress_bar_1(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'Processing Genbank files',
                   label        => 'GENBANK'
                  );
   &show_pbar_1;

   foreach my $gb_file (@{$args{input_list_ref}}) {
      #max number of Genbank entries
      my $max_count = $#{$args{input_list_ref}} + 1;
      my $label = "Processing Genbank file $gb_file";
      $progress = 0;

      #parse submitted Genbank file
      $progress++;
      &update_pbar_1(title     => 'Gene model',
                     label     => "Parsing Genbank file $gb_file",
                     progress  => ($progress / $max_count) * 100
                    );
      my ($header, $source, $core, $DNAseq, $start_ref,
          $feature_list_ref, $genbank_ref, $feature_counter) = &gb_parser (main_window   => $args{main_window},
                                                                           progress_bar  => $args{progress_bar},
                                                                           auto_ini_ref  => $args{auto_ini_ref},
                                                                           ini_ref       => $args{ini_ref},
                                                                           gb_file       => ${$args{ini_ref}}{input_files}.'/'.$gb_file
                                                                          );

      #save DNA sequence in separate file
      $progress++;
      &update_pbar_1(title     => 'Gene model',
                     label     => "Saving DNA sequence and moving genbank file to ${$args{ini_ref}}{move_msfasta}",
                     progress  => ($progress / $max_count) * 100,
                    );

      open WRITE, "+>${$args{ini_ref}}{input_files}".'/'.$gb_file.'_GAMOLAdna';
      print WRITE '>'.$gb_file."\n".$DNAseq;
      close WRITE;

      #building proper gene model
      $progress++;
      &update_pbar_1(title     => 'Gene model',
                     label     => "Building gene model for $gb_file",
                     progress  => ($progress / $max_count) * 100,
                    );
      my ($gene_model_ref) = &build_gene_model(main_window      => $args{main_window},
                                               progress_bar     => $args{progress_bar},
                                               auto_ini_ref     => $args{auto_ini_ref},
                                               ini_ref          => $args{ini_ref},
                                               gb_file          => $gb_file,
                                               genbank_ref      => $genbank_ref,
                                               feature_list_ref => $feature_list_ref
                                              );
      if ($gene_model_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error processing $gb_file.\nSkipping file for analysis",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nError processing $gb_file.\nSkipping file for analysis\n\n";
         close ERRORLOG;
         next;
      }

      #add genemodel to pooled gm array
      @combined_gm = (@combined_gm, @{$gene_model_ref});

      my %seen = ();
      @combined_gm = grep { ! $seen{$_} ++ } @combined_gm;

      #move Genbank file
      unless (-d ${$args{ini_ref}}{move_gb}) {
         mkdir(${$args{ini_ref}}{move_gb}, 0777) or do {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not create msfasta directory in  ${$args{ini_ref}}{input_files}",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error processing input files".
                           "\nCould not create msfasta directory in  ${$args{ini_ref}}{input_files}\n\n";
            close ERRORLOG;
            &hide_pbar_1;
            return (0); #return to main
         }
      }
      move (${$args{ini_ref}}{input_files}.'/'.$gb_file , ${$args{ini_ref}}{move_gb}.'/'.$gb_file);

      #add current corename to array for further processing
      push (@corenames, $gb_file.'_GAMOLAdna');
   }

   &hide_pbar_1;
   return (\@combined_gm, \@corenames);
}

sub build_gene_model {
   my %args = @_;
   my $preferred_key = 'CDS';
   my $max_feature = 0;

   #do we have CDS entries? if not choose gene entries - if not, error message
   #re-use existing gene model if updating Genbank files
   my @key = ();
   if (${$args{auto_ini_ref}}{update_genbank} == 1) {
      @key = grep{/\_CDS\_/i} @{${args{feature_list_ref}}};
      if ($#key < 0) {
         @key = grep{/\_gene\_/i} @{${args{feature_list_ref}}};
         if ($#key >= 0) {
            $preferred_key = 'gene';
         } else {
            ${$args{progress_bar}}->configure(-label=>"No CDS or gene key found in current Genbank file $args{gb_file}.\nGenerating new gene model.");
            ${$args{main_window}}->update;

            #generate new gene model from scratch for current genbank file
            my ($orf_list_ref) = &get_genemodel_fasta(main_window   => $args{main_window},
                                                      progress_bar  => $args{progress_bar},
                                                      auto_ini_ref  => $args{auto_ini_ref},
                                                      ini_ref       => $args{ini_ref},
                                                      input_file    => $args{gb_file}.'_GAMOLAdna',
                                                      genbank_new   => '1'
                                                      );
            if ($orf_list_ref eq '0') {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Error while generating and reading gene model for $args{gb_file}.",
                                                             -buttons => ['OK'],
                                                             -bitmap  => 'info');
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error processing input files".
                              "\nError while generating and reading gene model for $args{gb_file}.\n\n";
               close ERRORLOG;
               return (0); #return to main
            }
            return ($orf_list_ref);
         }
      }
   }
   #build new gene model for Genbank file is selected
   elsif (${$args{auto_ini_ref}}{create_genbank} == 1) {
      ${$args{progress_bar}}->configure(-label=>"Building new gene model for current Genbank file $args{gb_file}.\nGenerating new gene model.");
      ${$args{main_window}}->update;

      #generate new gene model from scratch for current genbank file
      my ($orf_list_ref) = &get_genemodel_fasta(main_window   => $args{main_window},
                                                progress_bar  => $args{progress_bar},
                                                auto_ini_ref  => $args{auto_ini_ref},
                                                ini_ref       => $args{ini_ref},
                                                input_file    => $args{gb_file}.'_GAMOLAdna',
                                                genbank_new   => '1'
                                                );
      if ($orf_list_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while generating and reading gene model for $args{gb_file}.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nError while generating and reading gene model for $args{gb_file}.\n\n";
         close ERRORLOG;
         return (0); #return to main
      }
      return ($orf_list_ref);
   }

   #sort start positions
   my @sorted =
      map  $_->[0] =>
      sort { $a->[1] <=> $b->[1] }
      map  [ $_, m/^(\d+)\_/ ]
      => @key;

   #defined highest features number
   foreach my $entry (@sorted) {
      'reset' =~ m/reset/;
      $entry =~ m/\d+_\d+_$preferred_key\_(\d+)/;
      my $ORF_number = $1;
      if ($ORF_number > $max_feature) {$max_feature = $ORF_number};
   }

   #build gene model and split joined features if necessary
   foreach my $entry (@sorted) {
      my ($ORF_number, $left_bd, $right_bd, $orientation, @boundaries);

      #get ORF number
      'reset' =~ m/reset/;
      $entry =~ m/\d+_\d+_$preferred_key\_(\d+)/i;
      $ORF_number = $1;

      unless (defined $ORF_number) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse feature boundaries from key -> $entry\.\nContent: ${$args{genbank_ref}}{$entry}",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Central Station: Error processing input files".
                        "\n1. Could not parse feature boundaries from key -> $entry\.\nContent: ${$args{genbank_ref}}{$entry}\n\n";
         close ERRORLOG;
         return (0);
      }

      #check orientation
      'reset' =~ m/reset/;
      if (${$args{genbank_ref}}{$entry} =~ m/^[^\n\r]+complement/) {
         $orientation = 'antisense';
      } else {
         $orientation = 'sense';
      }

      #catch boundary line from entry
      'reset' =~ m/reset/;
      ${$args{genbank_ref}}{$entry} =~ m/^([^\/]*?)\//s;
      my $boundaries = $1;
      unless (defined $boundaries && $boundaries =~ m/\S+/) {
         'reset' =~ m/reset/;
         ${$args{genbank_ref}}{$entry} =~ m/(.+)/s;
         $boundaries = $1;
      }
      @boundaries = ($boundaries =~ m/(\d+\D*?\d+)/g);

      if ($#boundaries < 0) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse feature boundaries from key -> $entry\.\nContent: ${$args{genbank_ref}}{$entry}",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Central Station: Error processing input files".
                        "\n2. Could not parse feature boundaries from key -> $entry\.\nContent: ${$args{genbank_ref}}{$entry}\nBoundaries content: $boundaries\.\n\n";
         close ERRORLOG;
         return (0);
      }
      #increase max_feature if joined features exist
      my $counter = 0;
      foreach my $boundary (@boundaries) {
         $counter++;
         my ($left_bd, $right_bd);
         'reset'=~m/reset/;
         $boundary =~ m/(\d+)\D*?(\d+)/;
         $left_bd  = $1;
         $right_bd = $2;
         unless (defined $right_bd && defined $left_bd) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse feature boundaries from key -> $entry\.\nContent: ${$args{genbank_ref}}{$entry}",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Central Station: Error processing input files".
                           "\n3. Could not parse feature boundaries from key -> $entry\.\nContent: ${$args{genbank_ref}}{$entry}\n\n";
            close ERRORLOG;
            return (undef,undef,0);
         }

         #parse and push into orf_list array
         if ($counter == 1) {
            if ($orientation eq 'sense') {
               push (@orf_list, $args{gb_file}.'_GAMOLAdna___'.$ORF_number.'___'.$left_bd.'___'.$right_bd.'___sense');
            } elsif ($orientation eq 'antisense') {
               push (@orf_list, $args{gb_file}.'_GAMOLAdna___'.$ORF_number.'___'.$left_bd.'___'.$right_bd.'___antisense');
            }
         } elsif ($counter > 1) {
            $max_feature++;
            if ($orientation eq 'sense') {
               push (@orf_list, $args{gb_file}.'_GAMOLAdna___'.$max_feature.'___'.$left_bd.'___'.$right_bd.'___sense');
            } elsif ($orientation eq 'antisense') {
               push (@orf_list, $args{gb_file}.'_GAMOLAdna___'.$max_feature.'___'.$left_bd.'___'.$right_bd.'___antisense');
            }
         }
      }
   }

   return (\@orf_list);
}

sub wanted {
   if ($_ =~ m/^$rename_filename.*?_$rename_id/) {
      push (@rename_files, $File::Find::name);
   }
}

sub sort_results {
   my %args = @_;
   my (@results, @files, $local_file, $suffix);

   #grab result files
   ${$args{progress_bar}}->configure(-label=>"Sorting result data");
   ${$args{main_window}}->update;

   #generate sorted hits
   $local_file = $args{input_file};
   $local_file =~ s/_GAMOLAdna$//;
   'reset' =~ m/reset/;
   $local_file =~ s/\.(gb|gbk)$//i;
   $suffix = $1;
   unless (defined $suffix) {$suffix = 'gb'};

   #clear directory first
   if (-d ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file) {
      ${$args{progress_bar}}->configure(-label=>"Removing old Consolidated data files");
      ${$args{main_window}}->update;
      &cleanup(${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file);
   }

   #local directory
   mkdir ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file;

   #Genbank file:
   copy(${$args{ini_ref}}{results}.'/'.$local_file.'.'.$suffix , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/'.$local_file.'.'.$suffix);

   #Blast
   ${$args{progress_bar}}->configure(-label=>"Sorting result data...Blast");
   ${$args{main_window}}->update;
   opendir RESULTS, ${$args{ini_ref}}{blast_results};
   @files = grep /^$args{input_file}\_\d+$/, readdir(RESULTS);
   closedir RESULTS;
   mkdir ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/Blast_database';
   foreach my $file (@files) {
      my ($db, $id);
      'reset' =~ m/reset/;
      $file =~ m/$args{input_file}\_(\d+)$/;
      $id = $1;
      copy(${$args{ini_ref}}{blast_results}.'/'.$file , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/Blast_database/'.$local_file.'.'.$suffix.'_'.$id);
   }

   #COG
   ${$args{progress_bar}}->configure(-label=>"Sorting result data...COG");
   ${$args{main_window}}->update;
   opendir RESULTS, ${$args{ini_ref}}{COG_results};
   @files = grep /^$args{input_file}\_COG_\d+$/, readdir(RESULTS);
   closedir RESULTS;
   mkdir ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/COG_database';
   foreach my $file (@files) {
       my ($db, $id);
      'reset' =~ m/reset/;
      $file =~ m/$args{input_file}(_.+)_(\d+)$/;
      $db = $1;
      $id = $2;
      copy(${$args{ini_ref}}{COG_results}.'/'.$file , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/COG_database/'.$local_file.'.'.$suffix.$db.'_'.$id);

   }

   #Pfam
   ${$args{progress_bar}}->configure(-label=>"Sorting result data...PFam");
   ${$args{main_window}}->update;
   opendir RESULTS, ${$args{ini_ref}}{pfam_results};
   @files = grep /^$args{input_file}\_/, readdir(RESULTS);
   closedir RESULTS;
   mkdir ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/PFam_database';
   foreach my $file (@files) {
      my ($db, $id);
      'reset' =~ m/reset/;
      $file =~ m/$args{input_file}(_.+)_(\d+)$/;
      $db = $1;
      $id = $2;
      copy(${$args{ini_ref}}{pfam_results}.'/'.$file , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/PFam_database/'.$local_file.'.'.$suffix.$db.'_'.$id);
   }

   #TIGRfam
   ${$args{progress_bar}}->configure(-label=>"Sorting result data...TIGRfam");
   ${$args{main_window}}->update;
   opendir RESULTS, ${$args{ini_ref}}{TIGRfam_results};
   @files = grep /^$args{input_file}\_/, readdir(RESULTS);
   closedir RESULTS;
   mkdir ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/TIGRfam_database';
   foreach my $file (@files) {
      my ($db, $id);
      'reset' =~ m/reset/;
      $file =~ m/$args{input_file}(_.+)_(\d+)$/;
      $db = $1;
      $id = $2;
      copy(${$args{ini_ref}}{TIGRfam_results}.'/'.$file , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/TIGRfam_database/'.$local_file.'.'.$suffix.$db.'_'.$id);
   }

   #SignalP
   ${$args{progress_bar}}->configure(-label=>"Sorting result data...SignalP directories");
   ${$args{main_window}}->update;
   opendir RESULTS, ${$args{ini_ref}}{signalp_results};
   @files = grep /^$args{input_file}\_SignalP_/, readdir(RESULTS);
   closedir RESULTS;
   mkdir ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/SignalP_database';
   #copy and rename directories
   foreach my $file (@files) {
      @rename_files = ();
      my ($db, $id);
      'reset' =~ m/reset/;
      $file =~ m/$args{input_file}(_.+)_(\d+)$/;
      $db = $1;
      $id = $2;
      #copy directories
      dircopy(${$args{ini_ref}}{signalp_results}.'/'.$file , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/SignalP_database/'.$local_file.'.'.$suffix.$db.'_'.$id);
      #find individual SignalP files
      ${$args{progress_bar}}->configure(-label=>"Sorting result data...SignalP individual files, id $id");
      ${$args{main_window}}->update;
      my @dir = (${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/SignalP_database/'.$local_file.'.'.$suffix.$db.'_'.$id);
      $rename_filename = $local_file;
      find (sub {push (@rename_files, $File::Find::name) if (-f && $_ =~ m/^$rename_filename/) }, ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/SignalP_database/'.$local_file.'.'.$suffix.$db.'_'.$id);
   }


   #Transterm
   ${$args{progress_bar}}->configure(-label=>"Sorting result data...Transterm");
   ${$args{main_window}}->update;
   opendir RESULTS, ${$args{ini_ref}}{transterm_results};
   @files = grep /^$args{input_file}\_transterm/, readdir(RESULTS);
   closedir RESULTS;
   mkdir ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/Transterm_database';
   foreach my $file (@files) {
      my ($extension);
      'reset' =~ m/reset/;
      $file =~ m/.+(_\S+)$/;
      $extension = $1;
      copy(${$args{ini_ref}}{transterm_results}.'/'.$file , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/Transterm_database/'.$local_file.'.'.$suffix.$extension);
   }

   #CRISPR
   ${$args{progress_bar}}->configure(-label=>"Sorting result data...CRISPR");
   ${$args{main_window}}->update;
   opendir RESULTS, ${$args{ini_ref}}{CRISPR_results};
   @files = grep /^$args{input_file}\.CRISPR/, readdir(RESULTS);
   closedir RESULTS;
   mkdir ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/CRISPR_database';
   foreach my $file (@files) {
      copy(${$args{ini_ref}}{CRISPR_results}.'/'.$file , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/CRISPR_database/'.$file);
   }

   #tRNAscan
   ${$args{progress_bar}}->configure(-label=>"Sorting result data...tRNAscan");
   ${$args{main_window}}->update;
   opendir RESULTS, ${$args{ini_ref}}{trnascan_results};
   @files = grep /$args{input_file}\_/, readdir(RESULTS);
   closedir RESULTS;
   mkdir ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/tRNAscan_database';
   foreach my $file (@files) {
      my ($extension);
      'reset' =~ m/reset/;
      $file =~ m/$args{input_file}(_.+)$/;
      $extension = $1;
      copy(${$args{ini_ref}}{trnascan_results}.'/'.$file , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/tRNAscan_database/'.$local_file.'.'.$suffix.$extension);
   }

   #tmhmm
   ${$args{progress_bar}}->configure(-label=>"Sorting result data...TMHMM");
   ${$args{main_window}}->update;
   mkdir ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/TMHMM_database';
   #copy directories
   dircopy(${$args{ini_ref}}{tmhmm_results}.'/TMHMM_'.$args{input_file} , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/TMHMM_database');

   #non-coding RNA
   ${$args{progress_bar}}->configure(-label=>"Sorting result data...non-coding RNAs");
   ${$args{main_window}}->update;
   opendir RESULTS, ${$args{ini_ref}}{rfam_results};
   @files = grep /^$args{input_file}(_chunk_\d+)?\.(blastn|RF\d+\.res|Rfam)$/, readdir(RESULTS);
   closedir RESULTS;
   mkdir ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/non_coding_RNA_database';
   foreach my $file (@files) {
      copy(${$args{ini_ref}}{rfam_results}.'/'.$file , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/non_coding_RNA_database/'.$file);
   }

   #ribosomal RNA
   ${$args{progress_bar}}->configure(-label=>"Sorting result data...ribosomal RNAs");
   ${$args{main_window}}->update;
   opendir RESULTS, ${$args{ini_ref}}{rrna_results};
   @files = grep /^$args{input_file}(\.|\-)/, readdir(RESULTS);
   closedir RESULTS;
   mkdir ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/rRNA_database';
   foreach my $file (@files) {
      copy(${$args{ini_ref}}{rrna_results}.'/'.$file , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/rRNA_database/'.$file);
   }

   #vector screen
   ${$args{progress_bar}}->configure(-label=>"Sorting result data...vector screen results");
   ${$args{main_window}}->update;
   opendir RESULTS, ${$args{ini_ref}}{vector_results};
   @files = grep /^$args{input_file}(\.|\-)/, readdir(RESULTS);
   closedir RESULTS;
   mkdir ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/vector_database';
   foreach my $file (@files) {
      copy(${$args{ini_ref}}{vector_results}.'/'.$file , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$local_file.'/vector_database/'.$file);
   }

   return (1);
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

sub msGenbank {
   my %args       = @_;
   my @temp_gb    = ();
   my @temp_fasta = ();
   #does msfasta folder exist? if not, create
   #use for msGenbank as well
   unless (-d ${$args{ini_ref}}{move_msfasta}) {
      mkdir(${$args{ini_ref}}{move_msfasta}, 0777) or do {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not create msfasta directory in  ${$args{ini_ref}}{input_files}",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info'
                                                       );
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nCould not create msfasta directory in  ${$args{ini_ref}}{input_files}\n\n";
         close ERRORLOG;
         &hide_pbar_1;
         return (0); #return to main
      }
   }

   #split msGenbank files into single Genbank entries and move resulting file to msfasta folder
   ${$args{progress_bar}}->configure(-label=>"Checking for msGenbank content");
   ${$args{main_window}}->update;

   foreach my $genbank_file (@{$args{input_list_ref}}) {
      #read file into memory
      my ($file_ref) = &slurp(main_window   => $args{main_window},
                              auto_ini_ref  => $args{auto_ini_ref},
                              ini_ref       => $args{ini_ref},
                              filename      => $genbank_file,
                              directory     => ${$args{ini_ref}}{input_files}
                             );
      if ($file_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not open input file $genbank_file",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info'
                                                       );
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nCould not open input file $genbank_file\n\n";
         close ERRORLOG;
         &hide_pbar_1;
         return (0); #return to main
      }
      #is it an msGenbank file?
      my $count = 0;
      #remove false msGBfiles with multpiple EOFs
      ${$file_ref}  =~ s/(\/\/\s*)+/\/\//gs;

      $count = (${$file_ref} =~ s/(\n\/\/)/$1/gi);
      if ($count > 1) {
         if (${$args{auto_ini_ref}}{concatenate_input_files} == 1) {
            #process msGenbank file
            #save as individal GB files and add to cluster file hash
            #convert into individual files instead
            ${$args{progress_bar}}->configure(-label=>"Processing msGenbank file $genbank_file");
            ${$args{main_window}}->update;
            #msGenbank file found -> split and move
            my $gb_cluster        = '';
            my @msgenbank         = ();
            my $new_cluster_entry = '';
            my $new_gb_name       = '';
            @msgenbank            = split(/\n\/\//,${$file_ref});
            foreach my $genbank (@msgenbank) {
               next unless ($genbank =~ m/\w+/);
               #catch locus
               $genbank      =~ s/^\s+//gs;
               $genbank      =~ m/^LOCUS\s+(\S+)\s/s;
               my $new_name  = $1;
               if ($new_name !~ /\w+/) {
                  my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                                -text    => "Could not grab one of the LOCUS tags from $genbank_file",
                                                                -buttons => ['OK'],
                                                                -bitmap  => 'info'
                                                                );
                  $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
                  $error_msg-> Show();
                  open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
                  print ERRORLOG "Error processing input files".
                                 "\nCould not grab one of the LOCUS tags from $genbank_file\n\n".
                                 "...$genbank ...";;
                  close ERRORLOG;
                  &hide_pbar_1;
                  return (0); #return to main
               }
               $new_name       =~ s/\W/_/g;
               $new_name       =~ s/_+$//;
               $new_gb_name = $genbank_file;
               $new_gb_name    =~ s/\.\w+$//; #remove extension
               $new_gb_name    =~ s/\W/_/g;
               $new_gb_name    =~ s/\_+$//;
               open WRITE, "+>${$args{ini_ref}}{input_files}".'/'.$new_gb_name.'_'.$new_name;
               print WRITE $genbank."\n\/\/";
               close WRITE;
               #update new input file array and cluster hash
               $gb_cluster .= $new_gb_name.'_'.$new_name."\n";
               $new_name    = "";
            }
            @msgenbank = ();

            #move msGenbank file to new directory
            move (${$args{ini_ref}}{input_files}.'/'.$genbank_file, ${$args{ini_ref}}{move_msfasta}.'/'.$genbank_file);

            #capture original msGenbank file to move back to input file folder
            push (@{$args{msgenbank_move_ref}}, $genbank_file);

            #create valid cluster entry for hash
            $gb_cluster =~ s/\n$//s;
            $new_cluster_entry = "File_name\:\t$new_gb_name\.conc\n".
                                 $gb_cluster;

            #add to concatenate hash
            ${$args{auto_ini_ref}}{concatenation_counter}++;
            ${$args{concatenate_clusters}}{${$args{auto_ini_ref}}{concatenation_counter}} = $new_cluster_entry;

            #concatenate now
            &concatenate(progress_bar         => $args{progress_bar},
                         main_window          => $args{main_window},
                         auto_ini_ref         => $args{auto_ini_ref},
                         ini_ref              => $args{ini_ref},
                         concatenate_clusters => $args{concatenate_clusters}
                        ); #concatenate_clusters hash contains concat_file name plus individual files

            #move source files to msfasta directory
            foreach my $entry (keys %{ $args{concatenate_clusters} }) {
               my @conc_files = split/\n/, $args{concatenate_clusters}->{$entry};
               undef $file_ref;
               foreach my $file (@conc_files) {
                  next if ($file =~ m/^File_name\:/);
                  if (-e ${$args{ini_ref}}{input_files}.'/'.$file) {
                     my ($status) = move (${$args{ini_ref}}{input_files}.'/'.$file , ${$args{ini_ref}}{move_msfasta}.'/'.$file);
                     if ($status eq '0') {
                        my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                                      -text    => "Error moving original file $file for concatenate cluster $entry.\nOriginal file will remain in file queue.",
                                                                      -buttons => ['OK'],
                                                                      -bitmap  => 'info'
                                                                     );
                        $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
                        $error_msg-> Show();
                        open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
                        print ERRORLOG "Error moving concatenation selections".
                                       "\nError moving original file $file for concatenate cluster $entry.\nOriginal file will remain in file queue.\n\n";
                        close ERRORLOG;
                     }
                  }
               }
            }

            #add to entry file list
            push (@temp_fasta, $new_gb_name.'.conc.cb');

            next;
         } elsif (${$args{auto_ini_ref}}{concatenate_input_files} == 0) {
            #convert into individual files instead
            my @msgenbank         = ();
            my $new_gb_name       = '';
            @msgenbank            = split(/\n\/\//,${$file_ref});
            foreach my $genbank (@msgenbank) {
               next unless ($genbank =~ m/\w+/);
               #catch locus
               $genbank      =~ m/^LOCUS\s+(\S+)\s/s;
               my $new_name  = $1;
               if ($new_name !~ /\w+/) {
                  my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                                -text    => "Could not grab one of the LOCUS tags from $genbank_file",
                                                                -buttons => ['OK'],
                                                                -bitmap  => 'info'
                                                                );
                  $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
                  $error_msg-> Show();
                  open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
                  print ERRORLOG "Error processing input files".
                                 "\nCould not grab one of the LOCUS tags from $genbank_file\n\n".
                                 "...$genbank ...";
                  close ERRORLOG;
                  &hide_pbar_1;
                  return (0); #return to main
               }
               $new_name     =~ s/\W/_/g;
               $new_name     =~ s/_+$//;
               $new_gb_name  = $genbank_file;
               $new_gb_name  =~ s/\.\w+$//; #remove extension
               $new_gb_name  =~ s/\W/_/g;
               $new_gb_name  =~ s/\_+$//;
               open WRITE, "+>${$args{ini_ref}}{input_files}".'/'.$new_gb_name.'_'.$new_name.'.gb';
               print WRITE $genbank."\n\/\/";
               close WRITE;
               #update new input file array
               push (@temp_gb, $new_gb_name.'_'.$new_name.'.gb');
               $new_name = '';
               #move msGenbank file to new directory
               move (${$args{ini_ref}}{input_files}.'/'.$genbank_file, ${$args{ini_ref}}{move_msfasta}.'/'.$genbank_file);
            }
            next;
         }
      }
      #update new input file array in case input is regular Genbank
      push (@temp_gb, $genbank_file);
   }

   #update fasta input array
   @{$args{input_list_ref}} = @temp_gb;
   push (@{$args{fasta_list_ref}}, @temp_fasta);
   @temp_gb    = ();
   @temp_fasta = ();

   return;
}

sub msfasta {
   my %args = @_;
   my @temp_fasta = ();

   #does msfasta folder exist? if not, create
   unless (-d ${$args{ini_ref}}{move_msfasta}) {
      mkdir(${$args{ini_ref}}{move_msfasta}, 0777) or do {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not create msfasta directory in  ${$args{ini_ref}}{input_files}",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info'
                                                       );
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nCould not create msfasta directory in  ${$args{ini_ref}}{input_files}\n\n";
         close ERRORLOG;
         &hide_pbar_1;
         return (0); #return to main
      }
   }

   #split msfasta files into single fasta entries and move resulting file to msfasta folder
   ${$args{progress_bar}}->configure(-label=>"Checking for msfasta content");
   ${$args{main_window}}->update;

   foreach my $fasta_file (@{$args{input_list_ref}}) {
      #read file into memory
      my ($file_ref) = &slurp(main_window   => $args{main_window},
                              auto_ini_ref  => $args{auto_ini_ref},
                              ini_ref       => $args{ini_ref},
                              filename      => $fasta_file,
                              directory     => ${$args{ini_ref}}{input_files}
                             );
      if ($file_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not open input file $fasta_file",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info'
                                                       );
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error processing input files".
                        "\nCould not open input file $fasta_file\n\n";
         close ERRORLOG;
         &hide_pbar_1;
         return (0); #return to main
      }
      #is it an msfasta file?
      my $count = 0;
      #my @mscount = split/\n\>/,${$file_ref};
      #$count = scalar @mscount;
      #undef @mscount;
      $count = (${$file_ref} =~ s/(\n\>)/$1/gi);
      if ($count >= 1) {
         #process msfasta file
         if (${$args{auto_ini_ref}}{concatenate_input_files} == 1) {
            &concatenate_individual_file (progress_bar         => $args{progress_bar},
                                          main_window          => $args{main_window},
                                          auto_ini_ref         => $args{auto_ini_ref},
                                          ini_ref              => $args{ini_ref},
                                          concatenate_file_ref => $file_ref,
                                          filename             => $fasta_file,
                                          concatenate_clusters => $args{concatenate_clusters}
                                         );
            #update input file list
            push (@temp_fasta, $fasta_file.'.cb');
            next;
         } elsif (${$args{auto_ini_ref}}{concatenate_input_files} == 0) {
            #convert into individual files instead
            ${$args{progress_bar}}->configure(-label=>"Processing msfasta file $fasta_file");
            ${$args{main_window}}->update;
            #msfasta file found -> split and move
            my @msfasta = ();
            @msfasta = split(/\n\>/,${$file_ref});
            foreach my $fasta (@msfasta) {
               #catch header
               $fasta =~ m/^\s*\>?\s*(.*?)\s+[^\n]*?\n(.+)/s;
               my $new_name = $1;
               my $seq = $2;
               $new_name =~ s/\W/_/g;
               $new_name =~ s/_+$//;
               #$fasta =~ s/^\s*\>//;
               $seq =~ s/\s//gs;
               open WRITE, "+>${$args{ini_ref}}{input_files}".'/'.$fasta_file.'_'.$new_name;
               print WRITE '>'.$new_name."\n".$seq;
               close WRITE;
               #update new input file array
               push (@temp_fasta, $fasta_file.'_'.$new_name);

               $new_name = "";
            }
            @msfasta = ();
            #move msfasta file to new directory
            move (${$args{ini_ref}}{input_files}.'/'.$fasta_file, ${$args{ini_ref}}{move_msfasta}.'/'.$fasta_file);
            next;
         }
      } else {
         #if not, check header anyway and modify if necessary
         unlink ${$args{ini_ref}}{input_files}.'/'.$fasta_file;
         open WRITE, "+>".${$args{ini_ref}}{input_files}.'/'.$fasta_file;
         ${$file_ref} =~ /^\>([^\n]+)\n(.+)/s;
         my $new_name = $1;
         my $seq = $2;
         $new_name =~ s/\W/_/g;
         $new_name =~ s/_+$//;
         $seq =~ s/\s//gs;
         print WRITE '>'.$new_name."\n".$seq;
         close WRITE;
      }

      #update new input file array in case input is regular fasta
      push (@temp_fasta, $fasta_file);
   }

   #update fasta input array
   @{$args{input_list_ref}} = @temp_fasta;
   @temp_fasta = ();

   return;
}

sub concatenate {
   my %args = @_;
   my ($total_cluster_number);

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'Concatenating selected file clusters',
                   label        => ''
                  );
   &show_pbar_2;

   #start iterating over all selected clusters
   #this creates an overall hash to store sequences and contig orders for individual files and clusters
   $total_cluster_number = keys %{$args{concatenate_clusters}};
   while (my($key, $file) = each %{$args{concatenate_clusters}} ) {
      my @conc_files        = ();
      my $cluster_file_name = '';
      my $skip_flag         = 0;

      #then grab individual files for this cluster
      @conc_files        = split/\n/, $file;

      #determine proper file name from hash
      'reset'            =~ m/reset/;
      $conc_files[0]     =~ m/^File_name\:\t(.+)/;
      $cluster_file_name = $1;
      if ($conc_files[0] =~ m/^File_name\:\t/) {
         shift @conc_files; #remove header line from list, got chosen file name
      }

      #does cluster file, concatenated file and contig order file exist? If yes, read files, put values to hash and skip
      if (-e ${$args{ini_ref}}{input_files}.'/'.$cluster_file_name.'.cb'                          &&
          -s ${$args{ini_ref}}{input_files}.'/'.$cluster_file_name.'.cb' > 0                      &&
          -e ${$args{ini_ref}}{input_files}.'/contig_order/'.$cluster_file_name.'.contig_order'     &&
          -s ${$args{ini_ref}}{input_files}.'/contig_order/'.$cluster_file_name.'.contig_order' > 0
         ) {
         #read ccontig_roder file and test if complete
         (my $file_ref) = slurp_cmd(auto_ini_ref  => $args{auto_ini_ref},
                                    ini_ref       => $args{ini_ref},
                                    filename      => $cluster_file_name.'.contig_order',
                                    directory     => ${$args{ini_ref}}{input_files}.'/contig_order'
                                   );
         if ($$file_ref =~ m/\/\/$/) {
            #set flag to ignore processing of conc_files_hash
            $skip_flag = 1;
         }
         undef $file_ref;
      }

      #test if each file is still present (in case of later re-runs), if not, ignore and delete from cluster
      my @temp = ();
      &update_pbar_2(label    => 'Testing file presence',
                     progress => ($key / $total_cluster_number) * 100,
                    );
      foreach my $entry (@conc_files) {
         unless (-e ${$args{ini_ref}}{input_files}.'/'.$entry  ||
                 -e ${$args{ini_ref}}{move_msfasta}.'/'.$entry ||
                 -e ${$args{ini_ref}}{move_gb}.'/'.$entry
                ) {
            open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
            print WRITE "Error in concatenating files.".
                        "\nError finding file $entry in directory ${$args{ini_ref}}{input_files}.\nFile will be ignored for concatenation\n\n";
            close WRITE;
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Error finding file $entry in directory ${$args{ini_ref}}{input_files}.\nFile will be ignored for concatenation",
                                                          -bitmap  => 'error',
                                                          -buttons => ['ok']);
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            next;
         }
         push (@temp, $entry) ;
      }
      @conc_files = @temp;
      undef @temp;

      #test file format for each file, put sequence and header into hash
      &update_pbar_2(label => 'Concatenating files');

      #files could be in 3 directories
      my @file_dir = (${$args{ini_ref}}{input_files}, ${$args{ini_ref}}{move_msfasta}, ${$args{ini_ref}}{move_gb});
      # $conc_files_hash this is where conc sequence and contig order are stored, key is cluster counter, seq+order as anon hash
      my $file_counter = 0;
      my $previous_pos = 1;
      my $scaf_counter = 1; #add artificial contig names for each scaffold
      foreach my $entry (@conc_files) {
         my $file_ref = 0;
         #are we having fasta or msfasta or genbank?
         my $count = 0; #check where the file resides if present
         while ($file_ref eq '0') {
            ($file_ref) = &slurp_cmd (main_window  => $args{main_window},
                                      progress_bar => $args{progress_bar},
                                      auto_ini_ref => $args{auto_ini_ref},
                                      ini_ref      => $args{ini_ref},
                                      directory    => $file_dir[$count],
                                      filename     => $entry,
                                      quiet        => 1
                                     );
            if ($file_ref eq '0') {
               $count++;
            }
            last if ($count > $#file_dir);
         }
         if ($file_ref eq '0' ) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Error reading sequence file $entry in cluster $key.\nSkipping",
                                                          -bitmap  => 'error',
                                                          -buttons => ['ok']
                                                          );
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error compiling Genbank files".
                           "\nError reading sequence file $entry in cluster $key.\nSkipping\n\n";
            close ERRORLOG;

            #abort if file has content, something must be wrong. Ignore and log if filesize is 0, no sequence in file
            if (-s "$file_dir[$count]\/$entry" > 0) {
               &hide_pbar_2;
               return (0);
            } else {
               next;
            }
         }

         #set the right key for the combined file cluster
         $file_counter++;

         #remove empty lines in the beginning
         $$file_ref =~ s/^\s+//s;

         #parse individual files and msfasta into common hash
         if ($$file_ref =~ /^\>/) { #msfasta format
            my ($new_fasta, $new_seq);
            my (@individual_seqs, $contig_order);
            $new_seq = "";
            my $c_length = 0;

            $$file_ref = "\n".$$file_ref;
            @individual_seqs = split /\n>/, $$file_ref;

            foreach my $sequence (@individual_seqs) {
               if ($sequence !~ /\w+/) {next};
               my $name = "";
               my $seq = "";
               $sequence = ">".$sequence;
               'reset' =~ m/reset/;
               $sequence =~ m/>(\S+)\s*[^\n]*?\n(.+)/s;
               $name = $1;
               $seq = $2;
               unless (defined $seq) {
                  my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                                -text    => "Cannot parse header and seq in file $entry.\nSkipping",
                                                                -bitmap  => 'error',
                                                                -buttons => ['ok']);
                  $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
                  $error_msg-> Show();
                  open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
                  print ERRORLOG "Error concatenating files".
                                 "\nCannot parse header and seq in file $entry.\nSkipping.\n$sequence\n\n";
                  close ERRORLOG;
                  next;
               }

               #clean-up sequence
               $seq =~ s/\s+//igs;
               $seq =~ s/[^a-z]+//igs;

               #substitute internal "N" boundaries with default stop sequence
               if (${$args{auto_ini_ref}}{replace_N} == 1) {
                  $seq =~ s/N{1,}/${$args{ini_ref}}{spacer}/ig;
               }

               #length of contig
               $c_length = length($seq);
               $name = $name."\t".$previous_pos."\t".($c_length + $previous_pos - 1);

               $previous_pos = $previous_pos + length(${$args{ini_ref}}{spacer}) + $c_length;

               $contig_order .= $name."\n";
               $new_seq .= $seq.${$args{ini_ref}}{spacer};
            }

            #replace multi spacer occurences
            $new_seq =~ s/${$args{ini_ref}}{spacer}${$args{ini_ref}}{spacer}+/${$args{ini_ref}}{spacer}/gs;
            #remove whitespace
            $new_seq =~ s/\s+//gs;
            $new_seq =~ s/[^a-z]+//igs;

            #this give contig order and concatenated sequence in 2 variables
            $contig_order .= "\n";

            #hash reference to build proper hash structure for all concatenated files
            $conc_files_hash-> { $key } -> { $file_counter } = { contig_order => $contig_order,
                                                                 sequence     => $new_seq,
                                                                 filename     => $entry,
                                                                 filedir      => $file_dir[$count],
                                                                 conc_file    => $cluster_file_name
                                                               };
            $new_seq      = '';
            $contig_order = '';
         } elsif ($$file_ref =~ /^locus/i) { #Genbank format
            my ($contig_order, $name);
            my ($header, $source, $core, $DNAseq);
            $header = 0;
            #send to gb
            my $count = 0; #check where the file resides if present
            while ($header eq '0') {
               ($header, $source, $core, $DNAseq) = &gb_parser (main_window   => $args{main_window},
                                                                progress_bar  => $args{progress_bar},
                                                                auto_ini_ref  => $args{auto_ini_ref},
                                                                ini_ref       => $args{ini_ref},
                                                                gb_file       => $file_dir[$count].'/'.$entry,
                                                                process_gb    => 'sequence',
                                                                quiet         => '1'
                                                               );
               $count++;
               last if ($count > $#file_dir);
            }
            if ($header eq '0') {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Error reading Genbank file $entry for concatenation. Skipping for analysis and moving file to msfasta folder",
                                                             -buttons => ['OK'],
                                                             -bitmap  => 'info'
                                                            );
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error reading Genbank file for conatenation".
                              "\nError reading Genbank file $entry for concatenation. Skipped from analysis and moved file to msfasta folder.\n$source\n\n";
               close ERRORLOG;
               undef $header;
               undef $source;
               undef $core;
               undef $DNAseq;

               #adjust file counter
               $file_counter--;

               #move to msfasta folder
               move (${$args{ini_ref}}{input_files}.'/'.$entry, ${$args{ini_ref}}{move_msfasta}.'/'.$entry);

               next;
            }
            if ($DNAseq !~ /\w+/) {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Error reading Genbank file $entry for concatenation, no DNA sequence found. Skipping for analysis and moving file to msfasta folder",
                                                             -buttons => ['OK'],
                                                             -bitmap  => 'info'
                                                            );
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error reading Genbank file for conatenation".
                              "\nError reading Genbank file $entry for concatenation, no DNA sequence found. Skipped from analysis and moved file to msfasta folder.\n$source\n\n";
               close ERRORLOG;
               undef $header;
               undef $source;
               undef $core;
               undef $DNAseq;

               #move to msfasta folder
               move (${$args{ini_ref}}{input_files}.'/'.$entry, ${$args{ini_ref}}{move_msfasta}.'/'.$entry);

               #adjust file counter
               $file_counter--;

               next;
            }

            #contig_order: retrieve Locus designator and definition from GB header
            'reset' =~ m/reset/;
            $header =~ m/LOCUS\s+(\S+?)\s.*?DEFINITION\s+(.+?)(ACCESSION|VERSION|KEYWORDS|SOURCE|FEATURES)/s;
            $name   = $1.'_'.$2;
            $name   =~ s/\n//gs;
            $name   =~ s/\s+/ /gs;
            $name   =~ s/\,\s*complete.+//s;
            $name   =~ s/\,\s*whole genome.+//s;

            #does GB already have contig boundaries?
            if ($core =~ m/fasta_record    \d+\.\.\d+/) {
               my @boundaries = ($core =~ m/fasta_record    (\d+\.\.\d+)/gs);
               foreach my $bd (@boundaries) {
                  $bd =~ m/(\d+)\.\.(\d+)/;
                  $contig_order .= $name."\t".$1."\t".$2."\n";
               }
               undef @boundaries;
            } elsif ($core =~ m/(gap             |assembly_gap    )\d+\.\.\d+/) {
               my $left_bd  = 1;
               my @gaps     = ();
               if ($core =~ m/gap             \d+\.\.\d+/) {
                  @gaps     = ($core =~ m/gap             (\d+\.\.\d+)/gs);
               } else {
                  @gaps     = ($core =~ m/assembly_gap    (\d+\.\.\d+)/gs);
               }
               my $temp_seq = '';
               foreach my $gap (@gaps) {
                  $gap =~ m/(\d+)\.\.(\d+)/;
                  my ($gap_left, $gap_right) = ($1, $2);

                  #define new name
                  my $c_name = $name;
                  $c_name    =~ s/\s+/_/gs;
                  $c_name    =~ s/\.$\_*//;

                  #grap subsequences
                  my $seq = substr ($DNAseq, $left_bd - 1, $gap_left - $left_bd);
                  #print "\n...$name ...$counter ...$left_bd ... $gap_left ...$gap_right ...\n...$seq ...\n";my$ttt=<STDIN>;
                  $temp_seq .= "\n\>$c_name\_$scaf_counter\n$seq";

                  #redefine new left boundary
                  $left_bd = $gap_right + 1;
                  $scaf_counter++;
               }

               #capture last entry
               {
                  #define new name
                  my $c_name = $name;
                  $c_name    =~ s/\s+/_/gs;
                  $c_name    =~ s/\.\_*$//;

                  #grap subsequences
                  my $seq = substr ($DNAseq, $left_bd - 1);
                  #print "\n...$name ...$counter ...$left_bd ... $gap_left ...$gap_right ...\n...$seq ...\n";my$ttt=<STDIN>;
                  $temp_seq .= "\n\>$c_name\_$scaf_counter\n$seq";

                  $scaf_counter++;
               }

               my @individual_seqs = split /\n>/, $temp_seq;
               my $new_seq         = '';
               #parse for contig boundaires and include non bleeding spacer
               foreach my $sequence (@individual_seqs) {
                  if ($sequence !~ /\w+/) {next};
                  my $name     = '';
                  my $seq      = '';
                  my $c_length = 0;
                  $sequence = '>'.$sequence;
                  'reset' =~ m/reset/;
                  $sequence =~ m/>(\S+)\s*[^\n]*?\n(.+)/s;
                  $name = $1;
                  $seq = $2;
                  unless (defined $seq) {
                     my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                                   -text    => "Cannot parse header and seq in file $entry.\nSkipping",
                                                                   -bitmap  => 'error',
                                                                   -buttons => ['ok']);
                     $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
                     $error_msg-> Show();
                     open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
                     print ERRORLOG "Error concatenating files".
                                    "\nCannot parse header and seq in file $entry.\nSkipping.\n$sequence\n\n";
                     close ERRORLOG;
                     next;
                  }

                  #clean-up sequence
                  $seq =~ s/\s+//igs;

                  #substitute internal "N" boundaries with default stop sequence
                  if (${$args{auto_ini_ref}}{replace_N} == 1) {
                     $seq =~ s/N{1,}/${$args{ini_ref}}{spacer}/ig;
                  }

                  #length of contig
                  $c_length = length($seq);
                  $name = $name."\t".$previous_pos."\t".($c_length + $previous_pos - 1);

                  $previous_pos = $previous_pos + length(${$args{ini_ref}}{spacer}) + $c_length;

                  $contig_order .= $name."\n";
                  $new_seq .= $seq.${$args{ini_ref}}{spacer};
               }

               #replace multi spacer occurances
               $new_seq =~ s/${$args{ini_ref}}{spacer}${$args{ini_ref}}{spacer}+/${$args{ini_ref}}{spacer}/gs;
               $new_seq =~ s/\s+//gs;
               $DNAseq = $new_seq;

               #this give contig order and concatenated seuqence in 2 variables
               $contig_order .= "\n";
               undef @gaps;
            } else {
               #my $end_pos   = $previous_pos + length(${$args{ini_ref}}{spacer}) + length($DNAseq);
               $name = $name."\t".$previous_pos."\t".(length($DNAseq) + $previous_pos - 1);
               $contig_order = $name."\n";
               $previous_pos = $previous_pos + length(${$args{ini_ref}}{spacer}) + length($DNAseq);
            }

            #hash reference to build proper hash structure for all concatenated files
            $conc_files_hash-> { $key } -> { $file_counter } = { contig_order => $contig_order,
                                                                 sequence     => $DNAseq,
                                                                 filename     => $entry,
                                                                 filedir      => $file_dir[$count],
                                                                 conc_file    => $cluster_file_name
                                                               };
            undef $header;
            undef $source;
            undef $core;
            undef $DNAseq;
         } else {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Format for input file $entry is not recognised. Skipping for concatenation.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info'
                                                         );
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error concatenating input files".
                           "\nFormat for input file $entry is not recognised. Skipping for concatenation.\n\n";
            close ERRORLOG;
         }
      }
   }

   &update_pbar_2(label    => 'Writing concatenated sequences');
   #make contig_order directory if it does not exists
   unless (-d ${$args{ini_ref}}{input_files}.'/contig_order') {
      mkdir (${$args{ini_ref}}{input_files}.'/contig_order', 0777);
   }

   #write concatenated sequences and contig order files
   foreach my $cluster (sort {$a<=>$b} keys %{ $conc_files_hash } ) {
      my ($conc_file_name, $conc_sequence, $conc_order);

      #does cluster file, concatenated file and contig order file exist? If yes, skip if re-using results; else delete files
      if (${$args{auto_ini_ref}}{reuse_results} == 1) {
         if (-e ${$args{ini_ref}}{input_files}.'/'.$conc_files_hash->{$cluster}->{'1'}->{conc_file}.'.cb'                          &&
             -s ${$args{ini_ref}}{input_files}.'/'.$conc_files_hash->{$cluster}->{'1'}->{conc_file}.'.cb' > 0                      &&
             -e ${$args{ini_ref}}{input_files}.'/contig_order/'.$conc_files_hash->{$cluster}->{'1'}->{conc_file}.'.contig_order'     &&
             -s ${$args{ini_ref}}{input_files}.'/contig_order/'.$conc_files_hash->{$cluster}->{'1'}->{conc_file}.'.contig_order' > 0 &&
             -e ${$args{auto_ini_ref}}{work_dir}.'/'.$cluster.'.cluster_files'                               &&
             -s ${$args{auto_ini_ref}}{work_dir}.'/'.$cluster.'.cluster_files' > 0
            ) {
            my ($file_ref) = &slurp(main_window   => $args{main_window},
                                    auto_ini_ref  => $args{auto_ini_ref},
                                    ini_ref       => $args{ini_ref},
                                    filename      => $conc_files_hash->{$cluster}->{'1'}->{conc_file}.'.contig_order',
                                    directory     => ${$args{ini_ref}}{input_files}.'/contig_order',
                                    quiet         => 1
                                   );
            #ignore and move to next cluster if contig_order file is complete (last file created)
            if ($$file_ref =~ m/\/\/$/) {
               undef $file_ref;
               next;
            }
            undef $file_ref;
         }
      } elsif (${$args{auto_ini_ref}}{reuse_results} == 0) {
         unlink ${$args{ini_ref}}{input_files}.'/'.$conc_files_hash->{$cluster}->{'1'}->{conc_file}.'.cb';
         unlink ${$args{ini_ref}}{input_files}.'/contig_order/'.$conc_files_hash->{$cluster}->{'1'}->{conc_file}.'.contig_order';
      }

      #create concatenated sequences for each cluster
      foreach my $counter (sort {$a<=>$b} keys %{ $conc_files_hash->{$cluster} } ) {
         my @temp = split /\n/, $conc_files_hash->{$cluster}->{$counter}->{contig_order};
         foreach my $entry (@temp) {
            next if ($entry !~ /\S+/);
            $conc_order .= $conc_files_hash->{$cluster}->{$counter}->{filename}."\t".$entry."\n";
         }
         undef @temp;

         $conc_sequence .= $conc_files_hash->{$cluster}->{$counter}->{sequence}.${$args{ini_ref}}{spacer};
      }
      #replace multi spacer occurances
      $conc_sequence =~ s/${$args{ini_ref}}{spacer}${$args{ini_ref}}{spacer}+/${$args{ini_ref}}{spacer}/gs;
      $conc_sequence =~ s/\s+//gs;
      #define filename for concatenated sequence
      'reset' =~ m/reset/;
      $conc_file_name = $conc_files_hash->{$cluster}->{'1'}->{conc_file};

      unless (defined $conc_file_name) {$conc_file_name = $cluster};

      #write sequence
      open  WRITE, "+>${$args{ini_ref}}{input_files}".'/'.$conc_file_name.'.cb';
      print WRITE '>'.$conc_file_name." Length=".length($conc_sequence)."\n".$conc_sequence."\n";
      close WRITE;

      #write contig order
      open  WRITE, "+>${$args{ini_ref}}{input_files}".'/contig_order/'.$conc_file_name.'.contig_order';
      print WRITE $conc_order;
      print WRITE "\n\/\/"; #make sure the file is complete
      close WRITE;

      #first, write files to concatenate into file
      open WRITE, "+>".${$args{auto_ini_ref}}{work_dir}.'/'.$cluster.'.cluster_files';
      print WRITE "File_name:\t".$conc_files_hash->{$cluster}->{'1'}->{conc_file}."\n";
      foreach my $counter (sort {$a<=>$b} keys %{ $conc_files_hash->{$cluster} } ) {
         print WRITE $conc_files_hash->{$cluster}->{$counter}->{filename}."\n";
      }
      close WRITE;

      undef $conc_sequence;
      undef $conc_order;
      undef $conc_file_name;
   }

   #input file list is read afterwards, so will automatically update
   &hide_pbar_2;
   return (1);
}

sub concatenate_individual_file {
   #this comes only from MSFASTA files.
   my %args = @_;
   my ($file);
   my (@entries);
   my $contig_order = '';
   my $local_count  = 0;
   my $chosen_name  = '';
   my $new_seq      = '';
   my $skip_flag    = 0;

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Concatenating individual file $args{filename}",
                   label        => ''
                  );
   &show_pbar_2;

   #predefine filenames already used;
   foreach my $cluster (sort {$a<=>$b} keys %{ $conc_files_hash } ) {
      unless (defined $local_name_hash{$conc_files_hash->{$cluster}->{'1'}->{conc_file}}) {
         $local_name_hash{$conc_files_hash->{$cluster}->{'1'}->{conc_file}} = 1;
      }
   }

   #if current filename is not yet used, initialise hash
   unless (defined $local_name_hash{$args{filename}}) {$local_name_hash{$args{filename}} = 0};

   #determine proper file name for cluster. By default use first entry, then iterate through rest of files, then use appendices
   my $total_count = keys %local_name_hash;

   #easy way, choose one of the selected files
   foreach my $entry (sort keys %local_name_hash) {
      $local_count++;
      $local_name_hash{$entry}++;
      if ($local_name_hash{$entry} == 1) {
         $chosen_name = $entry;
         last;
      }
   }
   #if already taken, increase local counter until unique name has been found.
   if ($local_count == $total_count) {
      $local_count = 0;
      while ($chosen_name !~ /\w+/) {
         foreach my $entry (sort keys %local_name_hash) {
            $local_count++;
            $local_name_hash{$args{filename}.'_'.$local_count}++;
            if ($local_name_hash{$args{filename}.'_'.$local_count} == 1) {
               $chosen_name = $args{filename}.'_'.$local_count;
               last;
            }
         }
      }
   }

   #does cluster file, concatenated file and contig order file exist? If yes, read files, put values to hash and skip
   my $conc_file_name = $chosen_name;
   if (-e ${$args{ini_ref}}{input_files}.'/'.$chosen_name.'.cb'                          &&
       -s ${$args{ini_ref}}{input_files}.'/'.$chosen_name.'.cb' > 0                      &&
       -e ${$args{ini_ref}}{input_files}.'/contig_order/'.$chosen_name.'.contig_order'     &&
       -s ${$args{ini_ref}}{input_files}.'/contig_order/'.$chosen_name.'.contig_order' > 0
      ) {
      #read ccontig_order file and test if complete
      (my $file_ref) = slurp_cmd(auto_ini_ref  => $args{auto_ini_ref},
                                 ini_ref       => $args{ini_ref},
                                 filename      => $chosen_name.'.contig_order',
                                 directory     => ${$args{ini_ref}}{input_files}.'/contig_order'
                                );
      if ($$file_ref =~ m/\/\/$/) {
         #set flag to ignore processing of conc_files_hash
         $skip_flag = 1;
      }
   }

   #increase cluster count
   $args{auto_ini_ref}{concatenation_counter}++;

   ${$args{concatenate_file_ref}} = "\n".${$args{concatenate_file_ref}};
   @entries = split /\n>/, ${$args{concatenate_file_ref}};
   my $c_length     = 0;
   my $previous_pos = 1;
   foreach my $entry (@entries) {
      if ($entry !~ /\w+/) {next};
      my $name = '';
      my $seq  = '';
      $entry   = '>'.$entry;
      'reset'  =~ m/reset/;
      $entry   =~ m/^>\s*(\S+)[^\n]*[\n](.+)/s;
      $name    = $1;
      $seq     = $2;

      unless (defined $seq) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Cannot parse header and seq in contig:\n$entry\n for file $args{filename}.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info'
                                                       );
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error concatenating individual input file ".
                        "\nCannot parse header and seq in contig:\n$entry\n for file $args{filename}.\n\n";
         close ERRORLOG;
         &hide_pbar_1;
         return (0); #return to main
      }

      #clean-up sequence
      $seq =~ s/\s+//igs;
      $seq =~ s/[^a-z]+//igs;

      #substitute internal "N" boundaries with default stop sequence
      if (${$args{auto_ini_ref}}{replace_N} == 1) {
         $seq =~ s/N{1,}/${$args{ini_ref}}{spacer}/ig;
      }

      #length of contig
      $c_length = length($seq);
      $name = $args{filename}."\t".$name."\t".$previous_pos."\t".($c_length + $previous_pos - 1);

      $previous_pos = $previous_pos + length(${$args{ini_ref}}{spacer}) + $c_length;

      $contig_order .= $name."\n";
      $new_seq .= $seq.${$args{ini_ref}}{spacer};
   }

   #replace multi spacer occurances
   $new_seq =~ s/${$args{ini_ref}}{spacer}${$args{ini_ref}}{spacer}+/${$args{ini_ref}}{spacer}/gs;
   $new_seq =~ s/\s+//gs;
   $new_seq =~ s/[^a-z]+//igs;

   #hash reference to build proper hash structure for all concatenated files
   $conc_files_hash-> { $args{auto_ini_ref}{concatenation_counter} } -> { '1' } = { contig_order => $contig_order,
                                                                                    sequence     => $new_seq,
                                                                                    filename     => $args{filename},
                                                                                    filedir      => ${$args{ini_ref}}{input_files},
                                                                                    conc_file    => $chosen_name
                                                                                  };

   #ignore if skip flag is set
   if ($skip_flag == 0) {
      &update_pbar_2(label    => "Writing concatenated sequences for file $args{filename}");

      #make contig_order directory if it does not exists
      unless (-d ${$args{ini_ref}}{input_files}.'/contig_order') {
         mkdir (${$args{ini_ref}}{input_files}.'/contig_order', 0777);
      }

      #write sequence
      open  WRITE, "+>${$args{ini_ref}}{input_files}".'/'.$chosen_name.'.cb';
      print WRITE '>'.$args{filename}." Length=".length($new_seq)."\n".$new_seq."\n";
      close WRITE;

      #write contig order
      &update_pbar_2(label    => "Writing contig order for file $chosen_name");
      open  WRITE, "+>${$args{ini_ref}}{input_files}".'/contig_order/'.$chosen_name.'.contig_order';
      print WRITE $contig_order;
      print WRITE "\n\/\/"; #make sure the file is complete
      close WRITE;

      #move original file
      &update_pbar_2(label    => "Moving original file $args{filename}");
      rename(${$args{ini_ref}}{input_files}.'/'.$args{filename}, ${$args{ini_ref}}{move_msfasta}.'/'.$args{filename});
   }

   &hide_pbar_2;
   return (1);
}

sub contig_boundaries {
   my %args = @_;
   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'Entering Contig Boundaries',
                   label        => ''
                  );
   &show_pbar_2;

   foreach my $cluster (sort {$a<=>$b} keys %{ $conc_files_hash} ) {
      my $gb_dir          = '';
      my $new_genbank     = '';
      my $colour          = 10;
      my $feature_counter = 0;
      my @contig_names    = ();
      my @order           = ();
      my $contig_counter  = 1;
      my $file_name       = $conc_files_hash->{$cluster}->{'1'}->{conc_file};
      &update_pbar_2(label    => "Entering contig boundaries for $file_name");

      #define GB directory
      if (${$args{auto_ini_ref}}{sort_results} == 1) {
         $gb_dir = ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$file_name.'.cb';
      } else {
         $gb_dir = ${$args{ini_ref}}{results};
      }

      #read and parse created Genbank file
      &update_pbar_2(label    => "Entering contig boundaries for $file_name\nParsing Genbank file");
      my ($header, $source, $core, $DNAseq, $start_ref,
          $feature_list_ref, $genbank_ref, undef) = &gb_parser (main_window   => $args{main_window},
                                                                progress_bar  => $args{progress_bar},
                                                                auto_ini_ref  => $args{auto_ini_ref},
                                                                ini_ref       => $args{ini_ref},
                                                                gb_file       => $gb_dir.'/'.$file_name.'.cb.gb',
                                                                DNA_cleanup   => 0
                                                               );
      if ($header eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not read Genbank file $file_name\.cb\.gb to enter contig boundaries. Skipping for concatenation.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info'
                                                      );
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error entering contig boundaries".
                        "\nCould not read Genbank file $file_name\.cb\.gb to enter contig boundaries.\n$source\nSkipping for concatenation.\n\n";
         close ERRORLOG;
         next;
      }

      #iterate
      &update_pbar_2(label    => "Entering contig boundaries for $file_name\nParsing boundaries");
      foreach my $counter (sort {$a<=>$b} keys %{ $conc_files_hash->{$cluster} } ) {
         #get current contig order
         my @boundaries = split /\n/, $conc_files_hash->{$cluster}->{$counter}->{'contig_order'};
         foreach my $entry (@boundaries) {
            my ($name, $left_bd, $right_bd);
            'reset' =~ m/reset/;
            $entry =~ m/([^\t]+)\t(\d+)\t(\d+)/;
            ($name, $left_bd, $right_bd) = ($1, $2, $3);
            unless (defined $right_bd) {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Cannot parse contig boundaries for $file_name.\nSkipping",
                                                             -bitmap  => 'error',
                                                             -buttons => ['ok']);
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error concatenating files".
                              "\nCannot parse contig boundaries for $file_name.\nSkipping\n$entry\n\n";
               close ERRORLOG;
               next;
            }
            push (@order, $left_bd."_".$right_bd."_".$name);
         }
      }

      &update_pbar_2(label    => "Entering contig boundaries for $file_name. Creating new features");
      #generate fasta_record keys
      foreach my $entry (@order) {
         my $value = "";
         $feature_counter++;
         'reset' =~ m/reset/;
         $entry =~ m/([^_]+)_([^_]+)_([^\n]+)/;
         my $start = $1;
         my $stop = $2;
         my $name = $3;
         #$name =~ m/([^\t]*?)\t/;
         #$name = $1;
         $name =~ s/\s//gs;

         $value = "fasta_record    $start\.\.$stop\n".
                  "                     \/product\=\"Contig: $name, Source: $file_name\"\n".
                  "                     \/label\=$name\n".
                  "                     /colour=$colour";
         $$genbank_ref{$start.'_'.$stop.'_fasta_record_'.$feature_counter} = $value;
         push (@$feature_list_ref, $start.'_'.$stop.'_fasta_record_'.$feature_counter);
         if ($colour == 10) {
             $colour = 11;
         } else {
            $colour = 10;
         }
      }

      #right format for genbank source?
      unless ($source =~ /^     /) {
         $source = '     '.$source;
      }

      #do we have feature/qualifier comment?
      unless ($header =~ /FEATURES             Location\/Qualifiers/ ||
              $source =~ /FEATURES             Location\/Qualifiers/) {
         $header .= "\n".'FEATURES             Location/Qualifiers';
      }

      #define key order
      my @gb_order = qw(RBS gene CDS COG_match PFAM_match TIGR_match ncRNA signalP TMM terminator CRISPR);

      #sort feature_list_ref via start position and genbank_key_order
      foreach my $key (@gb_order) {
         map { $_ =~ s/$key/$contig_counter\-$key/ } @$feature_list_ref;
         $contig_counter++;
      }

      #add generic number to those entries not listed in key order
      foreach my $entry (@$feature_list_ref) {
         $entry =~ s/^(\d+_\d+_)(\D+\_.+)/$1$contig_counter\-$2/;
      }

      my @sorted =
         map  $_->[0] =>
         sort { $a->[1] <=> $b->[1] ||
                $a->[2] <=> $b->[2] }
         map  [ $_, m/^(\d+)\_/, m/^\d+_\d+_(\d+)\-.*?_\d+/ ]
         => @$feature_list_ref;

      $contig_counter = 1;
      foreach my $entry (@sorted) {
         #remove genbank order key again
         $entry =~ s/_\d+\-(.+)$/\_$1/;

         #create access key for maintaining previous annotation
         my $short = $entry;
         $short =~ s/_\d+$//;

         #now add to GB file
         chomp($$genbank_ref{$entry});
         $$genbank_ref{$entry} =~ s/^\n//;

         #test if feature needs to be modified
         if ($entry =~ /~/) {
            $$genbank_ref{$entry} =~ s/^(\s*\w+?)~/$1_/;
         }

         if ($$genbank_ref{$entry} !~ /^     /) {
            $$genbank_ref{$entry} = "     ".$$genbank_ref{$entry};
         }
         $new_genbank .= "\n".$$genbank_ref{$entry};
      }
      $new_genbank =~ s/^\n//;

      #delete old Genbank file
      &update_pbar_2(label    => "Entering contig boundaries for $file_name. Replacing Genbank file");
      unlink $gb_dir.'/'.$file_name.'.cb.gb';

      #write new GB file
      open WRITE, '+>'.$gb_dir.'/'.$file_name.'.cb.gb';
      print WRITE $header.$source.$new_genbank."\n".$DNAseq;
      close WRITE;

      #release memmory
      undef $header;
      undef $source;
      undef $core;
      undef $DNAseq;
      undef $start_ref;
      undef $feature_list_ref;
      undef $genbank_ref;
      undef $feature_counter;

   }
   &hide_pbar_2;
   return (1);
}

sub move_cb_files {
   my %args = @_;
   foreach my $cluster (sort {$a<=>$b} keys %{ $conc_files_hash } ) {
      my $cluster_file_name = '';
      $cluster_file_name = $conc_files_hash->{$cluster}->{'1'}->{conc_file};

      #move conc.fasta files, contig order files and cluster files to results or consolidated folder
         #move original files back
      if (-e ${$args{ini_ref}}{move_msfasta}.'/'.$cluster_file_name) {
         move (${$args{ini_ref}}{move_msfasta}.'/'.$cluster_file_name , ${$args{ini_ref}}{input_files}.'/'.$cluster_file_name);
      }
      #move conc files
      if (${$args{auto_ini_ref}}{sort_results} == 1) {
         #move cluster file
         if (-e ${$args{auto_ini_ref}}{work_dir}.'/'.$cluster.'.cluster_files') {
            move (${$args{auto_ini_ref}}{work_dir}.'/'.$cluster.'.cluster_files' , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$cluster_file_name.'/'.$cluster.'.cluster_files');
         }
         #move concatenated fasta file
         if (-e ${$args{ini_ref}}{input_files}.'/'.$cluster_file_name.'.cb') {
            move (${$args{ini_ref}}{input_files}.'/'.$cluster_file_name.'.cb' , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$cluster_file_name.'.cb/'.$cluster_file_name.'.cb');
         }
         #move contig order file
         if (-e ${$args{ini_ref}}{input_files}.'/contig_order/'.$cluster_file_name.'.contig_order') {
            move (${$args{ini_ref}}{input_files}.'/contig_order/'.$cluster_file_name.'.contig_order' , ${$args{auto_ini_ref}}{work_dir}.'/Consolidated_results/'.$cluster_file_name.'.cb/'.$cluster_file_name.'.contig_order');
         }
      } else {
         #move cluster file
         if (-e ${$args{auto_ini_ref}}{work_dir}.'/'.$cluster.'.cluster_files') {
            move (${$args{auto_ini_ref}}{work_dir}.'/'.$cluster.'.cluster_files' , ${$args{ini_ref}}{results}.'/'.$cluster.'.cluster_files');
         }
         #move concatenated fasta file
         if (-e ${$args{ini_ref}}{input_files}.'/'.$cluster_file_name.'.cb') {
            move (${$args{ini_ref}}{input_files}.'/'.$cluster_file_name.'.cb' , ${$args{ini_ref}}{results}.'/'.$cluster_file_name.'.cb');
         }
         #move contig order file
         if (-e ${$args{ini_ref}}{input_files}.'/contig_order/'.$cluster_file_name.'.contig_order') {
            move (${$args{ini_ref}}{input_files}.'/contig_order/'.$cluster_file_name.'.contig_order' , ${$args{ini_ref}}{results}.'/'.$cluster_file_name.'.contig_order');
         }
      }
   }
   #move original msGenbank files back
   foreach my $entry (@{$args{msgenbank_move_ref}}) {
      move (${$args{ini_ref}}{move_msfasta}.'/'.$entry, ${$args{ini_ref}}{input_files}.'/'.$entry);
   }
   return;
}

1;
