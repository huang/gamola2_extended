#!/opt/ActivePerl-5.8/bin/perl5.8.8

#Critica program calls
#input arguments: progress_bar, directory, filename, ini_ref, auto_ini_ref
#output: one concatenated blastn result file in input_directory "xxx.blast"

package GeneModel::critica_call;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&run_critica);
use vars qw();

#use Basics::MesgBox;
use initialise::read_me qw(:DEFAULT);
use GeneModel::blast_contigs;
use GeneModel::make_blastpairs;
use GeneModel::setup;
use GeneModel::iterate_critica;

#local variables
my (%args, $tl);


sub run_critica {
   my %args = @_;

   #add scripts and binary location to environmental path in case they are not already in there
   unless ($ENV{'PATH'} =~ /${$args{auto_ini_ref}}{critica_bin}/ && $ENV{'PATH'} =~ /${$args{auto_ini_ref}}{critica_scripts}/) {
      my ($status) = &check_critica_binaries(main_window  => $args{main_window},
                                             ini_ref      => $args{ini_ref},
                                             auto_ini_ref => $args{auto_ini_ref}
                                            );
      if ($status == 0) {
         return (0);
         #create empty paths
         #${$args{auto_ini_ref}}{critica_bin} = "";
         #${$args{auto_ini_ref}}{critica_scripts} = "";
      }
      $ENV{'PATH'} = ${$args{auto_ini_ref}}{critica_bin}.':'.${$args{auto_ini_ref}}{critica_scripts}.':'.$ENV{'PATH'};
   }

   #start Critica with blast contigs
   ${$args{progress_bar}}->configure(-label=>"Starting Critica for $args{filename}");
   ${$args{main_window}}->update;
   &blastcontigs(main_window   => $args{main_window},
                 progress_bar  => $args{progress_bar},
                 auto_ini_ref  => $args{auto_ini_ref},
                 ini_ref       => $args{ini_ref},
                 filename      => $args{filename},
                 directory     => ${$args{ini_ref}}{input_files},
                 blast         => "${$args{ini_ref}}{blast_executables}/blastall -p blastn -g F -e 1e-4 -a ${$args{auto_ini_ref}}{CPU}",
                );
    #open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
    #print ERRORLOG "Error blastall: ${$args{ini_ref}}{blast_executables}/blastall -p blastn -g F -e 1e-4 -a ${$args{auto_ini_ref}}{CPU}.".
    #                 "\n\n";  
    #close ERRORLOG;

   #make blast pairs
   &blastpair_critica(main_window   => $args{main_window},
                      progress_bar  => $args{progress_bar},
                      auto_ini_ref  => $args{auto_ini_ref},
                      ini_ref       => $args{ini_ref},
                      filename      => $args{filename},
                      directory     => ${$args{ini_ref}}{input_files}
                     );

   #test if critica created blast pairs - if not, abort with error message
   unless (-s ${$args{ini_ref}}{input_files}.'/'.$args{filename}.'.blast.pairs' > 0) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Error creating Blast pairs in Critica for $args{filename}",
                                                    -buttons => ['OK'],
                                                    -bitmap    => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error creating Blast pairs in Critica for $args{filename}".
                     "\n\n";
      close ERRORLOG;

      &delete_temp_files(main_window   => $args{main_window},
                         progress_bar  => $args{progress_bar},
                         auto_ini_ref  => $args{auto_ini_ref},
                         ini_ref       => $args{ini_ref},
                         filename      => $args{filename},
                         directory     => ${$args{ini_ref}}{input_files}
                        );
      return (0);
   }

   #run "scanblastpairs" on results
   ${$args{progress_bar}}->configure(-label=>"Scanning Critica Blast pairs for $args{filename}");
   ${$args{main_window}}->update;
   `scanblastpairs ${$args{ini_ref}}{input_files}/$args{filename} ${$args{ini_ref}}{input_files}/$args{filename}.blast.pairs ${$args{ini_ref}}{input_files}/$args{filename}.triplets`;

   #iterate over critica
   ${$args{progress_bar}}->configure(-label=>"Iterating Critica for $args{filename}");
   ${$args{main_window}}->update;
   &iterate_critica(main_window   => $args{main_window},
                    progress_bar  => $args{progress_bar},
                    auto_ini_ref  => $args{auto_ini_ref},
                    ini_ref       => $args{ini_ref},
                    filename      => $args{filename},
                    directory     => ${$args{ini_ref}}{input_files},
                    iterations     => 3,
                    fraction       => 0.8,
                    iterations     => 3,
                    addlong        => 0,
                    nosd           => 0,
                    noprom         => 1,
                    startiteration => 0,
                    dicodonscores  => "",
                    initscores     => "",
                    options        => "",
                    fixedoptions   => "",
                    gencode        => "",
                   );
   #delete critica temp files; only keep final CDS model
   &delete_temp_files(main_window   => $args{main_window},
                      progress_bar  => $args{progress_bar},
                      auto_ini_ref  => $args{auto_ini_ref},
                      ini_ref       => $args{ini_ref},
                      filename      => $args{filename},
                      directory     => ${$args{ini_ref}}{input_files}
                     );

   #modify final CDS model
   &modify_final_model(main_window   => $args{main_window},
                       progress_bar  => $args{progress_bar},
                       auto_ini_ref  => $args{auto_ini_ref},
                       ini_ref       => $args{ini_ref},
                       filename      => $args{filename},
                       directory     => ${$args{ini_ref}}{genemodel_output}
                      );

   return(1);

}

sub delete_temp_files {
   my %args = @_;

   #delete tempfiles in input directory
   unlink ${$args{ini_ref}}{input_files}.'/'.$args{filename}.'.blast';
   unlink ${$args{ini_ref}}{input_files}.'/'.$args{filename}.'.blast.pairs';
   unlink ${$args{ini_ref}}{input_files}.'/'.$args{filename}.'.triplets';

   #delete tempfiles in gene model directory
   foreach (my $i=0;$i<${$args{auto_ini_ref}}{iterations};$i++) {
      unlink ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.$i.'.cds';
      unlink ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.$i.'.crit';
      unlink ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.$i.'.dicod';
      unlink ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.$i.'.init';
      unlink ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.$i.'.reject';
      unlink ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.$i.'.sdscores';
      unlink ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.$i.'.mesg';
   }
   #unlink temp files for final interation
   unlink ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.${$args{auto_ini_ref}}{iterations}.'.crit';
   unlink ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.${$args{auto_ini_ref}}{iterations}.'.dicod';
   unlink ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.${$args{auto_ini_ref}}{iterations}.'.init';
   unlink ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.${$args{auto_ini_ref}}{iterations}.'.reject';
   unlink ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.${$args{auto_ini_ref}}{iterations}.'.sdscores';
   unlink ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.${$args{auto_ini_ref}}{iterations}.'.mesg';
}

sub modify_final_model {
   my %args = @_;
   my @critica_array;
   my $orientation;

   #Glimmer like Critica gene model
   open WRITE, "+>".$args{directory}.'/'.$args{filename}.'.critica.mod';

   #read model into memory
   my ($gene_model_ref) = &slurp(main_window   => $args{main_window},
                                 auto_ini_ref  => $args{auto_ini_ref},
                                 ini_ref       => $args{ini_ref},
                                 filename      => $args{filename}.${$args{auto_ini_ref}}{iterations}.'.cds',
                                 directory     => $args{directory}
                                );
   if (${$gene_model_ref} eq '0' || ${$gene_model_ref} !~ /\w+/) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Empty Critica gene model for $args{filename}.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error generating Critica Model for file $args{filename}".
                     "\nEmpty Critica gene model for $args{filename}.\n\n";
      close ERRORLOG;
      close WRITE;
      return (0); #return to main
   }

   @critica_array = split/\n/,${$gene_model_ref};
   my $counter = 1;
   foreach my $critica (@critica_array) {
      my ($start, $stop, $pvalue, $length, $rbs_distance, $rbs_pattern);
      $critica =~ m/\s*[^\s]+\s+(\d+)\s+(\d+)\s+([^\s]+)\s+\d+\s+\d+\s+\d+\s+[^\s]+\s+\w+\s+[^\s]+\s+(\d+)\s+([^\s]+)/;
      $start = $1;
      $stop = $2;
      $pvalue = $3;
      $rbs_distance = $4;
      $rbs_pattern = $5;

      $length = abs($start - $stop);
      if ($start < $stop) {
         $stop -= 3;
         #RBS present?
         if ($rbs_pattern ne '-') {
            $orientation = '[+1 L='.$length.' r='.$pvalue.']'."\tsense__".($start-$rbs_distance).'__'.($start-$rbs_distance+length($rbs_pattern)).'__'.$rbs_pattern;
         } else {
            $orientation = '[+1 L='.$length.' r='.$pvalue.']';
         }
      }
      #antisense direction?
      elsif ($start > $stop) {
         $stop += 3;
         #RBS present?
         if ($rbs_pattern ne '-') {
            $orientation = '[-1 L='.$length.' r='.$pvalue.']'."\tantisense__".($start+$rbs_distance).'__'.($start+$rbs_distance+length($rbs_pattern)).'__'.$rbs_pattern;
         } else {
            $orientation = '[-1 L='.$length.' r='.$pvalue.']';
         }
      }
      print WRITE ' '.$counter.'    '.$start.'    '.$stop.'  '.$orientation."\n";
      $counter++;
      $rbs_pattern = "";
   }
   close WRITE;
   @critica_array = ();
   return(1);

}

1;