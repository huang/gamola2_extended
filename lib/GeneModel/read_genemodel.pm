#!/opt/ActivePerl-5.8/bin/perl

#read external genemodel: read external gene model files, initial support for Glimmer and GFF
#input arguments: main_window, ini_ref, auto_ini_ref, input_array reference

#GFF reader:
#            orf_list reference (includes filename and ORF ID, left and right bd and orientation)
#internal gene model reader: same output, reads generated internal gene model into memory

#both models will test for and INCLUDE stop codons

package GeneModel::read_genemodel;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&read_external_gene_model &read_internal_gene_model &read_internal_gene_model_frame @orf_list);
use vars qw();

use Tk;
use initialise::read_me              qw(:DEFAULT);
use ProgrammeModules::sequence       qw(:DEFAULT);
use ProgrammeModules::genbank_parser qw(:DEFAULT);
use Cwd;

#local variables
my (%args, @seqfile, $jjj, @orf_list, $include_stop);

#initialise


sub read_external_gene_model {
   my %args = @_;
   my $gm_ref = 0;
   ############
   # for now, only GFF and internal gene model (with "combined" extension) formats are supported
   #########

   #check if gene model file with "GFF" extension exists
   if ($args{input_file} =~ /\.gff$/i || -e ${$args{auto_ini_ref}}{ext_gm_folder}.'/'.$args{input_file}.'.gff'
                                      || -e ${$args{auto_ini_ref}}{ext_gm_folder}.'/'.$args{input_file}.'.GFF') {
      ($gm_ref) = &read_gff(main_window   => $args{main_window},
                            progress_bar  => $args{progress_bar},
                            auto_ini_ref  => $args{auto_ini_ref},
                            ini_ref       => $args{ini_ref},
                            gm_dir        => ${$args{auto_ini_ref}}{ext_gm_folder},
                            input_file    => $args{input_file}
                           );
   } else {
      #else try internal gene model format
      ($gm_ref) = &read_internal_gene_model(main_window   => $args{main_window},
                                            progress_bar  => $args{progress_bar},
                                            auto_ini_ref  => $args{auto_ini_ref},
                                            ini_ref       => $args{ini_ref},
                                            gm_dir        => ${$args{auto_ini_ref}}{ext_gm_folder},
                                            input_file    => $args{input_file}
                                           );
   }

   if ($gm_ref eq 0) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Error while reading external gene model for $args{input_file}",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error reading gene models for file $args{input_file}".
                     "\nError while reading external gene model for $args{input_file}\n\n";
      close ERRORLOG;
      return (0); #return to main
   }

   return ($gm_ref);
}


sub read_internal_gene_model {
   my %args = (input_file => '',
               gm_dir     => ${$args{ini_ref}}{genemodel_output},
               @_);
   @orf_list = ();

   #reading internal gene model to memory
   ${$args{progress_bar}}->configure(-label=>"Importing internal gene model for $args{input_file}");
   ${$args{main_window}}->update;

   #make sure each input file has a gene model
   #chdir $args{gm_dir};
   if (!-e $args{gm_dir}                      .'/'.$args{input_file}.'.combined' &&
       !-e ${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.combined') {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Missing internal gene model file for $args{input_file}",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error reading gene models for file $args{input_file}".
                     "\nMissing internal gene model file for $args{input_file}\n\n";
      close ERRORLOG;
      ${$args{progress_bar}}->configure(-label=>"Missing internal gene model file for $args{input_file}");
      ${$args{main_window}}->update;
      return (0, undef);
   }
   #read gene model to memory for each input file
   ${$args{progress_bar}}->configure(-label=>"Analysing $args{input_file}");
   ${$args{main_window}}->update;

   undef $include_stop;
   #open gene model file
   if (-e $args{gm_dir}.'/'.$args{input_file}.'.combined') {
      #return with empty file notice
      if (-s $args{gm_dir}.'/'.$args{input_file}.'.combined' < 1) {
         my $dna_seq = '';
         my ($file_ref) = &slurp_cmd (directory   => ${$args{ini_ref}}{input_files},
                                      filename    => $args{input_file}
                                      );
         #parsing input file
         if (${$file_ref} =~ /^\>/) {
            ${$file_ref}  =~ m/^(\>[^\n]+)\n(.+)/s;
            $dna_seq      = $2;
         } else {
            $dna_seq = ${$file_ref};
         }
         push (@orf_list, $args{input_file}.'___1___1___3___sense');
         push (@orf_list, $args{input_file}.'___2___'. (length ($dna_seq) - 2) .'___'.length ($dna_seq).'___sense');
         push (@orf_list, $args{input_file}.'___3___1___3___antisense');
         push (@orf_list, $args{input_file}.'___4___'. (length ($dna_seq) - 2) .'___'.length ($dna_seq).'___antisense');
         #write dummy gene model for later updates
         #open WRITE, "+>".$args{gm_dir}.'/'.$args{input_file}.'.combined';
         #print WRITE ' 1  1  3  [+1    10.00'."\n".
         #            ' 2  '. (length ($dna_seq) - 2) .'  '. length ($dna_seq) .'  [+1    10.00'."\n";
         #close WRITE;

         undef $dna_seq;
         undef $file_ref;
         ${$args{progress_bar}}->configure(-label=>"Gene model import finished");
         ${$args{main_window}}->update;
         return (\@orf_list);
      }
      open READ, "<$args{gm_dir}\/$args{input_file}\.combined" ;
   }
   #open alternative location
   else {
      #return with empty file notice
      if (-s ${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.combined' < 1) {
         my $dna_seq = '';
         my ($file_ref) = &slurp_cmd (directory   => ${$args{ini_ref}}{input_files},
                                      filename    => $args{input_file}
                                      );
         #parsing input file
         if ($file_ref eq '0') {
            print "\nERROR reading input file $args{input_file} in directory ${$args{ini_ref}}{input_files}\n";<>;
         }

         if (${$file_ref} =~ /^\>/) {
            ${$file_ref}  =~ m/^(\>[^\n]+)\n(.+)/s;
            $dna_seq      = $2;
         } else {
            $dna_seq = ${$file_ref};
         }
         push (@orf_list, $args{input_file}.'___1___1___3___sense');
         push (@orf_list, $args{input_file}.'___2___'. (length ($dna_seq) - 2) .'___'. length ($dna_seq) .'___sense');
         push (@orf_list, $args{input_file}.'___3___1___3___antisense');
         push (@orf_list, $args{input_file}.'___4___'. (length ($dna_seq) - 2) .'___'.length ($dna_seq).'___antisense');
         #write dummy gene model for later updates
         #open WRITE, "+>".${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.combined';
         #print WRITE ' 1  1  3  [+1    10.00'."\n".
         #            ' 2  '. (length ($dna_seq) - 2) .'  '. length ($dna_seq) .'  [+1    10.00'."\n";
         #close WRITE;

         undef $dna_seq;
         undef $file_ref;
         ${$args{progress_bar}}->configure(-label=>"Gene model import finished");
         ${$args{main_window}}->update;
         return (\@orf_list);
      }
      open READ, "<${$args{ini_ref}}{genemodel_output}\/$args{input_file}\.combined" ;
   }
   while (<READ>) {
      if ($_ !~ /\w+/) {next}; #skip empty lines
      my ($ID, $left_bd, $right_bd, $orient, $length);
      'reset' =~ m/reset/;
      m/^\s*(\d+)\s+(\d+)\s+(\d+)\s+(.+)/;
      $ID = $1;
      $left_bd = $2;
      $right_bd = $3;
      $orient = $4;

      unless (defined $ID && defined $left_bd && defined $right_bd && defined $orient) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not read entry $_ in file $args{input_file}\.combined. Skipping entry",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error reading gene models for file $args{input_file}".
                        "\nCould not read entry $_ in file $args{input_file}\.combined. Skipping entry\n\n";
         close ERRORLOG;
         next;
      }

      $length = $right_bd - $left_bd;
      $orient =~ s/\s*\[(\-|\+).+/$1/;
      if ($orient eq '-') {
         ($left_bd, $right_bd) = ($right_bd, $left_bd);
      }

      #test for stop codon
      unless (defined $include_stop && $include_stop =~ /(yes|no)/) {
         $include_stop = &test_stop_codon(main_window   => $args{main_window},
                                          auto_ini_ref  => $args{auto_ini_ref},
                                          ini_ref       => $args{ini_ref},
                                          filename      => $args{input_file},
                                          left_bd       => $left_bd,
                                          right_bd      => $right_bd,
                                          orient        => $orient);
         if ($include_stop eq '0') {return (0)};
      }

      if ($orient eq '+') {
         if ($include_stop eq 'no') {
            $right_bd += 3;
         }
         #$gene_model{$input.'_'.$ID} = "$ID  $left_bd  $right_bd  \[+1  L\=$length  r\=1\]";
         push (@orf_list, $args{input_file}.'___'.$ID.'___'.$left_bd.'___'.$right_bd.'___sense');
      } elsif ($orient eq '-') {
         if ($include_stop eq 'no') {
            $left_bd -= 3;
         }
         #$gene_model{$input.'_'.$ID} = "$ID  $right_bd  $left_bd  \[-1  L\=$length  r\=1\]";
         push (@orf_list, $args{input_file}.'___'.$ID.'___'.$left_bd.'___'.$right_bd.'___antisense');
      } else {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Problems parsing internal gene model file $args{input_file} at line\n$_\nUnable to determine orientation of $orient",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error reading gene models for file $args{input_file}".
                        "\nProblems parsing internal gene model file $args{input_file} at line\n$_\nUnable to determine orientation of $orient\n\n";
         close ERRORLOG;
         return (0);
      }
   }
   close READ;
   ${$args{progress_bar}}->configure(-label=>"Gene model import finished");
   ${$args{main_window}}->update;
   return (\@orf_list);
}

sub read_internal_gene_model_frame {
   my %args = (input_file => '',
               @_);
   @orf_list = ();
   #reading internal gene model to memory
   ${$args{progress_bar}}->configure(-label=>"Importing internal gene model for $args{input_file}");
   ${$args{main_window}}->update;

   #make sure each input file has a gene model
   chdir ${$args{ini_ref}}{genemodel_output};
   unless (-e $args{input_file}.'.combined') {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Missing internal gene model file for $args{input_file}",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error reading gene models for file $args{input_file}".
                     "\nMissing internal gene model file for $args{input_file}\n\n";
      close ERRORLOG;
      ${$args{progress_bar}}->configure(-label=>"Missing internal gene model file for $args{input_file}");
      ${$args{main_window}}->update;
      return (0, undef);
   }

   #read gene model to memory for each input file
   ${$args{progress_bar}}->configure(-label=>"Analysing $args{input_file}");
   ${$args{main_window}}->update;
   undef $include_stop;
   open READ, "<$args{input_file}".'.combined' ;
   while (<READ>) {
      if ($_ !~ /\w+/) {next}; #skip empty lines
      my $ID = "";
      my $left_bd = "";
      my $right_bd = "";
      my $orient = "";
      my $frame = "";
      m/^\s*(\d+)\s+(\d+)\s+(\d+)\s+(.+)/;

      unless (defined $1 && defined $2 && defined $3 && defined $4) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not read entry $_ in file $args{input_file}\.combined. Skipping entry",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error reading gene models for file $args{input_file}".
                        "\nCould not read entry $_ in file $args{input_file}\.combined. Skipping entry\n\n";
         close ERRORLOG;
         next;
      }

      $ID = $1;
      $left_bd = $2;
      $right_bd = $3;
      $orient = $4;
      $orient =~ s/\s*\[([\-\+])(\d).+/$1/;
      $frame = $2;
      #test for stop codon
      unless (defined $include_stop && $include_stop =~ /(yes|no)/) {
         $include_stop = &test_stop_codon(main_window   => $args{main_window},
                                          auto_ini_ref  => $args{auto_ini_ref},
                                          ini_ref       => $args{ini_ref},
                                          filename      => $args{input_file},
                                          left_bd       => $left_bd,
                                          right_bd      => $right_bd,
                                          orient        => $orient);
         if ($include_stop eq '0') {
            return (0);
         }
      }

      if ($orient eq '+') {
         if ($include_stop eq 'no') {
            $right_bd += 3;
         }
         #$gene_model{$input.'_'.$ID} = "$ID  $left_bd  $right_bd  \[+1  L\=$length  r\=1\]";
         push (@orf_list, $args{input_file}.'___'.$ID.'___'.$left_bd.'___'.$right_bd.'___sense___'.$frame);
      } elsif ($orient eq '-') {
         if ($include_stop eq 'no') {
            $left_bd -= 3;
         }
         #$gene_model{$input.'_'.$ID} = "$ID  $right_bd  $left_bd  \[-1  L\=$length  r\=1\]";
         push (@orf_list, $args{input_file}.'___'.$ID.'___'.$left_bd.'___'.$right_bd.'___antisense___'.$frame);
      } else {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Problems parsing internal gene model file $args{input_file} at line\n$_\nUnable to determine orientation of $orient",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error reading gene models for file $args{input_file}".
                        "\nProblems parsing internal gene model file $args{input_file} at line\n$_\nUnable to determine orientation of $orient\n\n";
         close ERRORLOG;
         return (0);
      }
   }
   close READ;
   ${$args{progress_bar}}->configure(-label=>"Gene model import finished");
   ${$args{main_window}}->update;

   return (\@orf_list);
}

sub read_gff {
   my %args = (input_file => '',
               @_);
   my (%seen, $joined, $counter);

   #convert gff gene model files into internal Glimmer compatible format
   ${$args{progress_bar}}->configure(-label=>"Importing GFF gene models");
   ${$args{main_window}}->update;

   #make sure each input file has a GFF gene model
   chdir ${$args{auto_ini_ref}}{ext_gm_folder};
   $args{input_file} =~ s/\.gff$//;

   unless (-e $args{input_file}.'.gff') {
      ${$args{progress_bar}}->configure(-label=>"Missing GFF gene model file for $args{input_file}");
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Missing GFF gene model file for $args{input_file} in folder ${$args{auto_ini_ref}}{ext_gm_folder}\.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error reading gene models for file $args{input_file}".
                     "\nMissing GFF gene model file for $args{input_file} in folder ${$args{auto_ini_ref}}{ext_gm_folder}\.\n\n";
      close ERRORLOG;
      ${$args{main_window}}->update;
      return (0);
   }

   #convert each GFF file to Glimmer compatible format
   ${$args{progress_bar}}->configure(-label=>"Analysing $args{input_file}");
   ${$args{main_window}}->update;
   $counter = 1;
   undef $include_stop;
   open READ, "<$args{input_file}".'.gff';
   while (<READ>) {
      if ($_ !~ /\w+/) {next}; #skip empty lines
      my ($feature, $left_bd, $right_bd, $orient, $length);
      'reset' =~ m/reset/;
      $_ =~ m/[^\t]+\t[^\t]+\t([^\t]+)\t(\d+)\t(\d+)\t[^\t]+\t([^\t])+\t/;
      $feature = $1;
      $left_bd = $2;
      $right_bd = $3;
      $orient = $4;

      #more columns?
      if ($left_bd =~ /\D+/) {
         'reset' =~ m/reset/;
         $_ =~ m/[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t(\d+)\t(\d+)\t[^\t]+\t([^\t])+\t/;
         $feature = $1;
         $left_bd = $2;
         $right_bd = $3;
         $orient = $4;
      }

      unless (defined $orient) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse provided GFF gene model fileno $args{input_file}\.gff.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error reading gene models for file $args{input_file}".
                        "\nCould not parse provided GFF gene model fileno $args{input_file}\.gff.\n\n";
         close ERRORLOG;
         return (0);
      }

      $length = $right_bd - $left_bd;
      $orient =~ s/\s*\[([\-\+])(\d).+/$1/;

      #test for stop codon
      unless (defined $include_stop && $include_stop =~ /(yes|no)/) {
         $include_stop = &test_stop_codon(main_window   => $args{main_window},
                                          auto_ini_ref  => $args{auto_ini_ref},
                                          ini_ref       => $args{ini_ref},
                                          filename      => $args{input_file},
                                          left_bd       => $left_bd,
                                          right_bd      => $right_bd,
                                          orient        => $orient);
         if ($include_stop eq '0') {return (0)};
      }

      #are there exons?
      if ($feature eq 'mRNA') {
         $joined =~ s/\,$//;
         push (@{$args{auto_ini_ref}{exons}}, $joined);
         $joined = '';
         $include_stop = 'yes'; #assume gene model is correct as is
      } else {
         $joined .= $left_bd.'..'.$right_bd.',';
      }

      #set gene boundary
      if ($feature ne 'mRNA') {
         $seen{$left_bd.'_'.$right_bd.'_'.$orient}++;
      }
      next if ($seen{$left_bd.'_'.$right_bd.'_'.$orient} > 1);

      #gene model has to have CDS
      unless ($feature =~ /CDS/i) {next};
      unless ($feature =~ /\w+/ && $left_bd =~ /\w+/ && $right_bd =~ /\w+/ && $orient =~ /.+/) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Problems parsing GFF file $args{input_file} at line\n$_",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error reading gene models for file $args{input_file}".
                        "\nProblems parsing GFF file $args{input_file} at line\n$_\n\n";
         close ERRORLOG;
         return (0);
      }
      if ($orient eq '+') {
         if ($include_stop eq 'no') {
            $right_bd += 3;
         }
         #$gene_model{$input.'_'.$ID} = "$counter  $left_bd  $right_bd  \[+1  L\=$length  r\=1\]";
         push (@orf_list, $args{input_file}.'___'.$counter.'___'.$left_bd.'___'.$right_bd.'___sense');
      } elsif ($orient eq '-') {
         if ($include_stop eq 'no') {
            $left_bd -= 3;
         }
         #$gene_model{$input.'_'.$ID} = "$counter  $right_bd  $left_bd  \[-1  L\=$length  r\=1\]";
         push (@orf_list, $args{input_file}.'___'.$counter.'___'.$left_bd.'___'.$right_bd.'___antisense');
      } else {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Problems parsing GFF file $args{input_file} at line\n$_\nUnable to determine orientation of $orient",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error reading gene models for file $args{input_file}".
                        "\nProblems parsing GFF file $args{input_file} at line\n$_\nUnable to determine orientation of $orient\n\n";
         close ERRORLOG;
         return (0);
      }
      $counter++;
   }

   #get last exons
   $joined =~ s/\,$//;
   push (@{$args{auto_ini_ref}{exons}}, $joined);
   undef $joined;
   undef %seen;

   close READ;
   ${$args{progress_bar}}->configure(-label=>"Gene model import finished");
   ${$args{main_window}}->update;
   return (\@orf_list);
}

sub test_stop_codon {
   my %args = @_;
   my $format = '';
   my $seq_ref;
   my $test_seq;

   #test if input file is GB or fasta
   open READFILE, '<'.${$args{ini_ref}}{input_files}.'/'.$args{filename} or do {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Could not open file $args{filename} in directory ${$args{ini_ref}}{input_files}",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error reading file $args{filename}".
                     "\nCould not open file $args{filename} in directory ${$args{ini_ref}}{input_files}\n\n";
      close ERRORLOG;
      return (0);
   };
   while (<READFILE>) {
      if ($_ =~ /^\>/) {
         $format = 'fasta'; #fasta format
      } elsif ($_ =~ /^locus/i) {
         $format = 'gb'; #Genbank format
      } else {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Format for input file $args{filename} is not recognised.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error reading file $args{filename}".
                        "\nFormat for input file $args{filename} is not recognised.\n\n";
         close ERRORLOG;
         return (0);
      }
      last;
   }
   close READFILE;

   #get nt sequence
   if ($format eq 'fasta') {
      ($seq_ref) = &slurp (main_window => $args{main_window},
                           directory   => ${args{ini_ref}}{input_files},
                           filename    => $args{filename}
                          );
      if ($seq_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not read input file $args{filename}.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error reading file $args{filename}".
                        "\nCould not read input file $args{filename}.\n\n";
         close ERRORLOG;
         return (0);
      }
      $$seq_ref =~ s/^\>[^\n]+?\n//; #eliminate header
      $$seq_ref =~ s/\s//gs; #eliminate whitespaces
      $test_seq = $$seq_ref; #transfer to local var
   } elsif ($format eq 'gb') {
      (undef, undef, undef, $test_seq, undef, undef, undef, undef) = gb_parser (main_window => $args{main_window},
                                                                                gb_file   => ${args{ini_ref}}{input_files}.'/'.$args{filename},
                                                                                process_gb => 'sequence'
                                                                               );
   }

   #test for presence of stop position
   if ($args{orient} =~ /\+/) {
      my $test = substr($test_seq, ($args{right_bd} - 4), 3);
      if ($test eq /(TGA|TAG|TAA)/) {
            return ('yes');
      }
   } elsif ($args{orient} =~ /\-/) {
      my $test = substr($test_seq, ($args{left_bd} - 1), 3);
      if ($test eq /(TCA|CTA|TTA)/) {
            return ('yes');
      }
   }
   undef $test_seq;
   return ('no'); #no stop codon included
}

1;