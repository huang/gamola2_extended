#!/opt/ActivePerl-5.8/bin/perl

# Script to automate iterative runs of CRITICA
#input arguments: progress_bar, directory, filename, ini_ref, auto_ini_ref
#output: one gene model result file in input_directory "xxx.critica_genemodel"

package GeneModel::iterate_critica;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&iterate_critica);
use vars qw();

use Tk;
use Basics::MesgBox;
use Proc::Background;
use initialise::read_me qw(:DEFAULT);
use GeneModel::blast_contigs;
use GeneModel::setup;
use Critica;

#local variables
my (%args, $scriptname, $tripletsfile, $seqfile, $outputprefix,
    $totalinits, $i, $options, $scratch, $criticaoutput, $crit,
    $mesg, $cdsfile, $rejectfile, $sdscores, $promscores, $dicodonscores,
    $initscores, $err, $fixedoptions, $gencode, $addlong, $iterations,
    $nosd, $noprom, $diflag, $fraction, $startiteration);

sub iterate_critica {
   my %args = @_;
   my $timescale = '80';
   $fixedoptions = "";
   $gencode = "";
   $options = "";
   $err = 0;
   $startiteration = 0;
   $fraction=0.8;
   $addlong=0;
   $nosd=0;
   $noprom=1;

   $scriptname=getScriptName();

   &process_options(auto_ini_ref  => $args{auto_ini_ref});

   $tripletsfile=$args{directory}.'/'.$args{filename}.'.triplets';
   $seqfile=$args{directory}.'/'.$args{filename};
   $outputprefix=${$args{ini_ref}}{genemodel_output}.'/'.$args{filename};

   if(!-e($seqfile)) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "$scriptname: No such file as $seqfile!",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error iterating Critica for file $args{input_file}".
                     "\n$scriptname: No such file as $seqfile!\n\n";
      close ERRORLOG;
      return (0, "Error in Critica iteration.");
   }
   if(!-e($tripletsfile)) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "$scriptname: No such file as $tripletsfile!",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error iterating Critica for file $args{input_file}".
                     "\n$scriptname: No such file as $tripletsfile!\n\n";
      close ERRORLOG;
      return (0, "Error in Critica iteration.");
   }

   $totalinits=$outputprefix.".init";


   for($i=$startiteration;$i<=${args{auto_ini_ref}}{iterations};$i++) {
     $options.=$fixedoptions;
     $scratch=$outputprefix.".tmp";
     $criticaoutput=$outputprefix.$i;
     $crit=$criticaoutput.".crit";
     $mesg=$criticaoutput.".mesg";
     $cdsfile=$criticaoutput.".cds";
     $rejectfile=$criticaoutput.".reject";
     $sdscores=$criticaoutput.".sdscores";
     $promscores=$criticaoutput.".promscores";
     $dicodonscores=$criticaoutput.".dicod";
     $initscores=$criticaoutput.".init";

     ${$args{progress_bar}}->configure(-label=>"Executing Critica for $args{filename} in iteration $i");
     ${$args{main_window}}->update;
     my $local_start_time = time();
     my $status = Proc::Background->new('perl', ${$args{auto_ini_ref}}{work_dir}."\/lib\/Basics\/progress_bar.pl", "Critica", "Executing Critica for $args{filename} in iteration $i", $timescale);
     $err=`critica $gencode $options $seqfile $tripletsfile > $crit`;

     #$err=system("critica $seqfile $tripletsfile $gencode $options > $crit 2> $mesg");
     $timescale = (time() - $local_start_time) * 6;
     if ($timescale == 0) {$timescale = 2};
     $status->die;
     if ($err) {
        my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                      -text    => "Critica died with error message ". $err/256 ." on file $args{filename}",
                                                      -buttons => ['OK'],
                                                      -bitmap  => 'info');
        $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
        $error_msg-> Show();
        open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
        print ERRORLOG "Error iterating Critica for file $args{input_file}".
                       "\nCritica died with error message ". $err/256 ." on file $args{filename}\n\n";
        close ERRORLOG;
        return (0, "Error in Critica iteration.");
     }

     bestInits($crit,$scratch);
     if (($i==0) && ($addlong>0)) {
       ${$args{progress_bar}}->configure(-label=>"Executing Addlongorfs for $args{filename}");
       ${$args{main_window}}->update;
       $status = Proc::Background->new('perl', ${$args{auto_ini_ref}}{work_dir}."\/lib\/Basics\/progress_bar.pl", "Addlongorfs", "Executing Addlongorfs for $args{filename}");
       $err=system("addlongorfs -orf-aa-length=$addlong $gencode $scratch $seqfile >> $scratch");
       $status->die;
       if ($err) {
          my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                        -text    => "Addlongorfs died with error message ". $err/256 ." on file $args{filename}",
                                                        -buttons => ['OK'],
                                                        -bitmap  => 'info');
          $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
          $error_msg-> Show();
          open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
          print ERRORLOG "Error iterating Critica for file $args{input_file}".
                         "Addlongorfs died with error message ". $err/256 ." on file $args{filename}\n\n";
          close ERRORLOG;
          return (0, "Error in Critica iteration.");
       }
     }
     ${$args{progress_bar}}->configure(-label=>"Executing Removeoverlaps for $args{filename}");
     ${$args{main_window}}->update;
     $status = Proc::Background->new('perl', ${$args{auto_ini_ref}}{work_dir}."\/lib\/Basics\/progress_bar.pl", "Removeoverlaps", "Executing Removeoverlaps for $args{filename}");
     $err=system("removeoverlaps $seqfile $scratch 20 > $cdsfile 2>$rejectfile");
     $status->die;
     if ($err) {
        my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                      -text    => "Removeoverlaps died with error message ". $err/256 ." on file $args{filename}",
                                                      -buttons => ['OK'],
                                                      -bitmap  => 'info');
        $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
        $error_msg-> Show();
        open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
        print ERRORLOG "Error iterating Critica for file $args{input_file}".
                       "Removeoverlaps died with error message ". $err/256 ." on file $args{filename}\n\n";
        close ERRORLOG;
        return (0, "Error in Critica iteration.");
     }
     cdsSort($cdsfile);
     if ($i<${args{auto_ini_ref}}{iterations}) {
       if (!$nosd) {
         computeSDScores($crit,$sdscores);
       }
       if (!$noprom) {
         computePromScores($crit,$promscores);
       }
       computeInitScores($cdsfile,$seqfile,$initscores);
       if ($i==0) {
         $diflag="-fraction-coding=$fraction";
       }
       else {
         $diflag="";
       }
       ${$args{progress_bar}}->configure(-label=>"Executing Dicodontable for $args{filename}");
       ${$args{main_window}}->update;

       $status = Proc::Background->new('perl', ${$args{auto_ini_ref}}{work_dir}."\/lib\/Basics\/progress_bar.pl", "Dicodontable", "Executing Dicodontable for $args{filename}");
       $err=system("dicodontable $diflag $cdsfile $seqfile > $dicodonscores 2>> $mesg");
       $status->die;
       if ($err) {
          my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                        -text    => "Dicodontable died with error message ". $err/256 ." on file $args{filename}",
                                                        -buttons => ['OK'],
                                                        -bitmap  => 'info');
          $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
          $error_msg-> Show();
          open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
          print ERRORLOG "Error iterating Critica for file $args{input_file}".
                         "Dicodontable died with error message ". $err/256 ." on file $args{filename}\n\n";
          close ERRORLOG;
          return (0, "Error in Critica iteration.");
       }
       $options="-dicodon-scores=$dicodonscores -init-scores=$initscores";
       if (!$nosd) {
         $options.=" -sd-scores=$sdscores";
       }
       if (!$noprom) {
         $options.=" -prom-scores=$promscores";
       }

     }
   }
   $criticaoutput=$outputprefix.$i;
   $cdsfile=$criticaoutput.".cds";
   unlink $scratch;
}

sub usage {
  printf("Reference: Badger, Jonathan H. and Gary J. Olsen. CRITICA:\n");
  printf("Coding Region Identification Tool Invoking Comparative Analysis.\n");
  printf("Molecular Biology and Evolution 16: 512-524 (1999)\n\n");

  printf("$scriptname: [options] output-name seq-file triplets-scores\n");
  printf("Valid $scriptname options:\n");
  printf("\t-iterations=number (default: 3)\n");
  printf("\t-scoring-matrix=file\n");
  printf("\t-initial-dicod=file\n");
  printf("\t-initial-init=file\n");
  printf("\t-initial-sd=file\n");
  printf("\t-initial-prom=file\n");
  printf("\t-no-sdscores\n");
  printf("\t-prom-find\n");
  printf("\t-threshold=value\n");
  printf("\t-alpha=value\n");
  printf("\t-fraction-coding=value\n");
  printf("\t-add-longorfs=length\n");
  printf("\t-genetic-code=number\n");
  printf("\t-strict-threshold\n");
  printf("\t-frameshift-threshold=value\n");
  die("\t-quick-stats\n");
}

sub process_options {
   my %args = @_;
   if (${args{auto_ini_ref}}{scoring_matrix_state} == 1) {
      $fixedoptions.=' -scoring-matrix='.$args{scoring_matrix};
   }
   if (${args{auto_ini_ref}}{dicodon_scores_state} == 1) {
      $startiteration=1;
      $options.=' -dicodon-scores='.${args{auto_ini_ref}}{dicodon_scores_file};
   }
   if (${args{auto_ini_ref}}{init_scores_state} == 1) {
      $startiteration=1;
      $options.=' -init-scores='.${args{auto_ini_ref}}{init_scores_file};
   }
   if (${args{auto_ini_ref}}{sd_scores_state} == 1) {
      $startiteration=1;
      $options.=' -sd-scores='.${args{auto_ini_ref}}{sd_scores_file};
   }
   if (${args{auto_ini_ref}}{prom_scores_state} == 1) {
      $startiteration=1;
      $options.=' -prom-scores='.${args{auto_ini_ref}}{prom_scores_file};
   }
   if (${args{auto_ini_ref}}{no_sdscores_state} == 1) {
      $nosd=1;
      $fixedoptions.=' -no-sdscores';
   }
   if (${args{auto_ini_ref}}{prom_find_state} == 1) {
      $noprom=0;
      $fixedoptions.=' -prom-find';
   }
   if (${args{auto_ini_ref}}{threshold_state} == 1) {
      $fixedoptions.=' -threshold='.${args{auto_ini_ref}}{threshold};
   }
   if (${args{auto_ini_ref}}{alpha_state} == 1) {
      $fixedoptions.=' -alpha='.${args{auto_ini_ref}}{alpha};
   }
   if (${args{auto_ini_ref}}{quick_stats_state} == 1) {
      $fixedoptions.=" -quick-stats";
   }
   if (${args{auto_ini_ref}}{frameshift_threshold_state} == 1) {
      $fixedoptions.=' -frameshift-threshold='.${args{auto_ini_ref}}{frameshift_threshold};
   }
   #if (${args{auto_ini_ref}}{genetic_code_state} == 1) {
      $gencode='-genetic-code='.${args{auto_ini_ref}}{genetic_code};
   #}
   if (${args{auto_ini_ref}}{strict_threshold_state} == 1) {
      $fixedoptions.=" -strict-threshold";
   }
}

1;