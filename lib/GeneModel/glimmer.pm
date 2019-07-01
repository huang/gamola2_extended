#!/opt/ActivePerl-5.8/bin/perl

#Glimmer 2 and 3 program calls
#input arguments: progress_bar, directory, filename, ini_ref, auto_ini_ref


package GeneModel::glimmer;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&glimmer);
use vars qw();

use Basics::MesgBox;
use initialise::read_me qw(:DEFAULT);

#local variables
my (%args, $tl);




sub glimmer {
   my %args = @_;
   my $dna_seq;
   my $glimmer_coord="";
   my $str = "";
   my $header = "";
   my $res = "";

format FASTA =
^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<~~
$dna_seq
.

   #generating status message
   ${$args{progress_bar}}->configure(-label=>"Generating Glimmer gene model for $args{input_file}.");
   ${$args{main_window}}->update;


   #test for training file and/or create a new one if selected
   if (!-e ${$args{auto_ini_ref}}{selected_gl_model} || ${$args{auto_ini_ref}}{make_training_file} == 1) {
      &make_training_glimmer(main_window   => $args{main_window},
                             progress_bar  => $args{progress_bar},
                             auto_ini_ref  => $args{auto_ini_ref},
                             ini_ref       => $args{ini_ref},
                             input_file    => $args{input_file},
                             directory     => ${$args{ini_ref}}{input_files}
                            );

   }

   #reading input file
   my ($file_ref) = &slurp(main_window => \$args{progress_bar},
                           directory   => $args{directory},
                           filename    => $args{input_file}
                          );
   #parsing input file
   if (${$file_ref} =~ /^\>/) {
      ${$file_ref} =~ m/^(\>[^\n]+)\n(.+)/s;
      $header = $1;
      $dna_seq = $2;
   } else {
      $header = '>'.$args{input_file};
      $dna_seq = ${$file_ref};
   }

   #cleanup and format DNA sequence
   $dna_seq =~ s/\s//gs;
   $dna_seq =~ s/[^acgtmkrybvdhwsn]/n/igs;
   $dna_seq = lc($dna_seq);

   open  FASTA, '>', \$str or die;
   write FASTA;
   close FASTA;

   #MODIFIED
   #open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
   #print ERRORLOG "Error tempfile: ${$args{ini_ref}}{tempfile}\n\n";
   #close ERRORLOG;

   #write tempfile for glimmer
   open  WRITE, "+>".${$args{ini_ref}}{tempfile} or return (0);
   print WRITE $header."\n".$str."\n";
   close WRITE;

   #prepare Glimmer run
   $glimmer_coord = $args{input_file}.'.glcoord';

   #adapt to Glimmer version
   if (${$args{auto_ini_ref}}{use_glimmer2} == 1) {
      ##stringent GLIMMER setting
      $res = `${$args{ini_ref}}{glprog2}/glimmer2 ${$args{ini_ref}}{tempfile} ${$args{auto_ini_ref}}{selected_gl_model} -l -f`;

      if ($res =~ /\w+/) {
         open WRITE, '+>'.${$args{ini_ref}}{genemodel_output}.'/'.$glimmer_coord or do {
            my $error_msg = ${$args{main_window}}->Dialog(-title  => 'Error',
                                                          -text   => "Error1: ${$args{ini_ref}}{glprog2}/glimmer2 ${$args{ini_ref}}{tempfile} ${$args{auto_ini_ref}}{selected_gl_model} -l -f: $res! while creating Glimmer output file $glimmer_coord".
                                                                     "\nin directory ${$args{ini_ref}}{genemodel_output}",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg->Show();
            open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error while creating Glimmer output file $glimmer_coord".
                           "\nin directory ${$args{ini_ref}}{genemodel_output}\n\n";
            close ERRORLOG;
         };
         print WRITE $res;
         close WRITE;
      } else {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error2 creating Glimmer gene model for $args{input_file}.\nCheck Error.log for detailed information.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error creating Glimmer gene model for file $args{input_file}".
                        "\nUsing gene model ${$args{auto_ini_ref}}{selected_gl_model} on Glimmer2\n\n";
         close ERRORLOG;

         #create empty dummy file
         open  WRITE, "+> ${$args{ini_ref}}{genemodel_output}/$glimmer_coord\.mod";
         print WRITE ' ';
         close WRITE;

         #remove temp file
         unlink ${$args{ini_ref}}{tempfile};

         return (1);
      }

      #reformat Glimmer 2 output
      &glimmer2_format(main_window   => $args{main_window},
                       auto_ini_ref  => $args{auto_ini_ref},
                       ini_ref       => $args{ini_ref},
                       directory     => ${$args{ini_ref}}{genemodel_output},
                       filename      => $glimmer_coord
                      );
   } elsif (${$args{auto_ini_ref}}{use_glimmer3} == 1) {
      #define translation table
      'reset' =~ m/reset/;
      ${$args{auto_ini_ref}}{translation_table} =~ m/\s*(\d+)\:/;
      my $translation_code = $1;
      unless(defined $translation_code) {$translation_code = 1};
      `${$args{ini_ref}}{glprog3}/glimmer3 -z $translation_code -g 100 -l ${$args{ini_ref}}{tempfile} ${$args{auto_ini_ref}}{selected_gl_model} ${$args{ini_ref}}{genemodel_output}/$glimmer_coord 2>&1`;

      unless (-e ${$args{ini_ref}}{genemodel_output}.'/'.$glimmer_coord.'.predict') {
        my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                      -text    => "Error: creating Glimmer gene model for $args{input_file}.\nProbably wrong model file ${$args{auto_ini_ref}}{selected_gl_model}.\nCheck Error.log for detailed information.",
                                                      -buttons => ['OK'],
                                                      -bitmap  => 'info');
        $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
        $error_msg-> Show();
        open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
        print ERRORLOG "Error3 creating Glimmer gene model for file $args{input_file}".
                       "\nUsing gene model ${$args{auto_ini_ref}}{selected_gl_model} on Glimmer3: ${$args{ini_ref}}{glprog3}/glimmer3 -z $translation_code -g 100 -l ${$args{ini_ref}}{tempfile} ${$args{auto_ini_ref}}{selected_gl_model} ${$args{ini_ref}}{genemodel_output}/$glimmer_coord 2>&1\n\n";
        close ERRORLOG;

        #create empty dummy file
        open WRITE, "+> ${$args{ini_ref}}{genemodel_output}/$glimmer_coord\.mod";
        print WRITE ' ';
        close WRITE;

        return (1);
      }

      #reformat Glimmer 3 output
      &glimmer3_format(main_window   => $args{main_window},
                       auto_ini_ref  => $args{auto_ini_ref},
                       ini_ref       => $args{ini_ref},
                       directory     => ${$args{ini_ref}}{genemodel_output},
                       filename      => $glimmer_coord
                      );
   }
   unlink ${$args{ini_ref}}{tempfile};
   return (1);
}

sub glimmer2_format {
   my %args = @_;
   my $glimmer_mod = "";

   $glimmer_mod = $args{directory}.'/'.$args{filename}.'.mod';

   my ($gene_model_ref) = &slurp(main_window   => $args{main_window},
                                 directory     => $args{directory},
                                 filename      => $args{filename}
                                );


   #test if putative genes were predicted, if not, create empty file and return
   if (${$gene_model_ref} =~/Putative Genes:/) {
      ${$gene_model_ref} =~ s/^.*?Putative Genes://s;
   } else {
      open WRITE, "+>$glimmer_mod";
      close WRITE;
      return;
   }

   #make the ORF numbers sequentially
   my $increase_orf_number = 1;
   my @sq_array = ();
   @sq_array = split/\n/,${$gene_model_ref};
   open WRITE, "+>$glimmer_mod";
   foreach (@sq_array) {
      if ($_ !~ /\w+/) {next};
      s/(\s+)\d+(\s+)/$1$increase_orf_number$2/;
      $increase_orf_number++;
      print WRITE $_."\n";
   }
   close WRITE;
}

sub glimmer3_format {
   my %args = @_;
   my $glimmer_mod = "";
   my $temporder = "";
   $glimmer_mod = $args{directory}.'/'.$args{filename}.'.mod';

   open GLOUT, $args{directory}.'/'.$args{filename}.'.predict';
   #convert to Glimmer 2 output
   while (<GLOUT>) {
      #skip header
      if ($_ !~ /^orf/) {next};
      #remove new ORF labeling
      $_ =~ s/orf[0]*/    /;
      #add syntax
      $_ =~ s/\+(\d)/\[\+$1/;
      $_ =~ s/\-(\d)/\[\-$1/;
      #reduce stop position by 3 to match Glimmer2 output
      #for forward or reverse strand
      if ($_ =~ /\[\+/) {
         m/(\s+\d+\s+\d+\s+)(\d+)(.+)/;
         my $new_stp = "";
         $new_stp = $2-3;
         $_ = $1.$new_stp.$3."\n";
      } elsif ($_ =~ /\[\-/) {
         m/(\s+\d+\s+)(\d+)(\s+)(\d+)(.+)/;
         my $new_stp = "";
         my $new_stt = "";
         $new_stp = $2;
         $new_stt = $4+3;
         $_ = $1.$new_stp.$3.$new_stt.$5."\n";
      }

      $temporder.=$_;
   }
   close GLOUT;

   open WRITE, "+>$glimmer_mod";
   print WRITE $temporder."\n";
   close WRITE;
}

sub make_training_glimmer {
   my %args = @_;

   #find long orfs
   ${$args{progress_bar}}->configure(-label=>"Finding long ORFs for $args{input_file}.");
   ${$args{main_window}}->update;
   if (${$args{auto_ini_ref}}{use_glimmer2} == 1) {
      `${$args{ini_ref}}{glprog2}/long-orfs $args{directory}/$args{input_file} > ${$args{ini_ref}}{glprog2}/long_orf.txt`;
   } elsif (${$args{auto_ini_ref}}{use_glimmer3} == 1) {
      `${$args{ini_ref}}{glprog3}/long-orfs $args{directory}/$args{input_file}   ${$args{ini_ref}}{glprog3}/long_orf.txt 2>&1`;
   }

   #extract ORF sequences
   ${$args{progress_bar}}->configure(-label=>"Extracing long ORF sequences for $args{input_file}.");
   ${$args{main_window}}->update;
   if (${$args{auto_ini_ref}}{use_glimmer2} == 1) {
      `${$args{ini_ref}}{glprog2}/extract $args{directory}/$args{input_file} ${$args{ini_ref}}{glprog2}/long_orf.txt > ${$args{ini_ref}}{glprog2}/long_orf.seq`;
      unlink ${$args{ini_ref}}{glprog2}.'/long_orf.txt';
   } elsif (${$args{auto_ini_ref}}{use_glimmer3} == 1) {
      `${$args{ini_ref}}{glprog3}/extract $args{directory}/$args{input_file} ${$args{ini_ref}}{glprog3}/long_orf.txt > ${$args{ini_ref}}{glprog3}/long_orf.seq 2>&1`;
      unlink ${$args{ini_ref}}{glprog3}.'/long_orf.txt';
   }

   #build ICM
   ${$args{progress_bar}}->configure(-label=>"Building Glimmer ICM for $args{input_file}.");
   ${$args{main_window}}->update;

   'reset' =~ m/reset/;
   ${$args{auto_ini_ref}}{translation_table} =~ m/\s*(\d+)\:/;
   my $translation_code = $1;
   unless(defined $translation_code) {$translation_code = 1};

   if (${$args{auto_ini_ref}}{use_glimmer2} == 1) {
      chdir ${$args{ini_ref}}{glprog2};
      `./build-icm < ${$args{ini_ref}}{glprog2}/long_orf.seq > ${$args{ini_ref}}{glimmer_model}/$args{input_file}\.model2`;
      chdir ${$args{auto_ini_ref}}{work_dir};
      unlink ${$args{ini_ref}}{glprog2}.'/long_orf.seq';
   } elsif (${$args{auto_ini_ref}}{use_glimmer3} == 1) {
      `${$args{ini_ref}}{glprog3}/build-icm -r -z $translation_code ${$args{ini_ref}}{glimmer_model}/$args{input_file}\.model3 < ${$args{ini_ref}}{glprog3}/long_orf.seq 2>&1`;
      unlink ${$args{ini_ref}}{glprog3}.'/long_orf.seq';
   }

   #set new gene model file:
   if (${$args{auto_ini_ref}}{use_glimmer2} == 1) {
      ${$args{auto_ini_ref}}{selected_gl_model} = ${$args{ini_ref}}{glimmer_model}.'/'.$args{input_file}.'.model2';
      ${$args{auto_ini_ref}}{gl_short_file}     = $args{input_file}.'.model2';
   } elsif (${$args{auto_ini_ref}}{use_glimmer3} == 1) {
      ${$args{auto_ini_ref}}{selected_gl_model} = ${$args{ini_ref}}{glimmer_model}.'/'.$args{input_file}.'.model3';
      ${$args{auto_ini_ref}}{gl_short_file}     = $args{input_file}.'.model3';
   }

   return;

}


1;
