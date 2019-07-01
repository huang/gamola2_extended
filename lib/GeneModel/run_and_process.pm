#!/opt/ActivePerl-5.8/bin/perl

#run and process internal gene models:
#input arguments: main_window, ini_ref, auto_ini_ref,


package GeneModel::run_and_process;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&make_genemodel);
use vars qw();

use initialise::read_me              qw(:DEFAULT);
use ProgrammeModules::sequence       qw(:DEFAULT);
use ProgrammeModules::genbank_parser qw(:DEFAULT);
use GeneModel::glimmer               qw(:DEFAULT);
use GeneModel::prodigal              qw(:DEFAULT);
use GeneModel::critica_call          qw(:DEFAULT);
use GeneModel::igblast               qw(:DEFAULT);
use GeneModel::rbsfinder             qw(:DEFAULT);
use Cwd;
use File::Copy;

#local variables
my (%args, $status, $input_list_ref, @fasta_input, @genbank_input);

sub make_genemodel {
   my %args = @_;

   ${$args{progress_bar}}->configure(-label=>"Determining gene models for $args{input_file}");
   ${$args{main_window}}->update;
   #run Glimmer if selected
   if (${$args{auto_ini_ref}}{use_glimmer2} == 1 || ${$args{auto_ini_ref}}{use_glimmer3} == 1) {
      ${$args{progress_bar}}->configure(-label=>"Creating Glimmer gene model for $args{input_file}");
      ${$args{main_window}}->update;
      my ($status) = &glimmer(main_window   => $args{main_window},
                              progress_bar  => $args{progress_bar},
                              auto_ini_ref  => $args{auto_ini_ref},
                              ini_ref       => $args{ini_ref},
                              input_file    => $args{input_file},
                              directory     => ${$args{ini_ref}}{input_files}
                             );
      #open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      #print ERRORLOG "Error_status: ($status)";
      #close ERRORLOG;

      if ($status == 0) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while generating Glimmer gene model for ${$args{ini_ref}}{input_files},
						        $args{auto_ini_ref}, $args{ini_ref}, $args{input_file}. Aborting analysis",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error generating Glimmer Model for file $args{input_file}".
                        "\nError while generating Glimmer gene model for $args{input_file}. Aborting analysis\n\n";
         close ERRORLOG;
         return (0); #return to main
      }
   }
   
   #run Prodigal if selected
   if (${$args{auto_ini_ref}}{use_prodigal} == 1) {
      ${$args{progress_bar}}->configure(-label=>"Creating Prodigal gene model for $args{input_file}");
      ${$args{main_window}}->update;
      my ($status) = &prodigal(main_window   => $args{main_window},
                               progress_bar  => $args{progress_bar},
                               auto_ini_ref  => $args{auto_ini_ref},
                               ini_ref       => $args{ini_ref},
                               input_file    => $args{input_file},
                               directory     => ${$args{ini_ref}}{input_files}
                              );
      if ($status == 0) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while generating Prodigal gene model for $args{input_file}. Aborting analysis",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error generating Prodigal Model for file $args{input_file}".
                        "\nError while generating Prodigal gene model for $args{input_file}. Aborting analysis\n\n";
         close ERRORLOG;
         return (0); #return to main
      }
   }
   
   #run Critica if selected
   if (${$args{auto_ini_ref}}{use_critica} == 1) {
      ${$args{progress_bar}}->configure(-label=>"Creating Critica gene model for $args{input_file}");
      ${$args{main_window}}->update;

      my ($status) = &run_critica(main_window   => $args{main_window},
                                  progress_bar  => $args{progress_bar},
                                  auto_ini_ref  => $args{auto_ini_ref},
                                  ini_ref       => $args{ini_ref},
                                  filename      => $args{input_file},
                                  directory     => ${$args{ini_ref}}{input_files}
                                 );
      if ($status == 0) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while generating Critica gene model for $args{input_file}. Aborting analysis",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error generating Critica Model for file $args{input_file}".
                        "\nError while generating Critica gene model for $args{input_file}. Aborting analysis\n\n";
         close ERRORLOG;
         return (0); #return to main
      }
   }

   #merge gene models if .combined file doesn't exists
   #unless (-e ${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.combined') {
   #   if (-e ${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.glcoord.mod' && -e ${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.critica.mod') {
         ${$args{progress_bar}}->configure(-label=>"Merging gene models for $args{input_file}");
         ${$args{main_window}}->update;
         &combine_gene_models(main_window   => $args{main_window},
                              progress_bar  => $args{progress_bar},
                              auto_ini_ref  => $args{auto_ini_ref},
                              ini_ref       => $args{ini_ref},
                              filename      => $args{input_file},
                              directory     => ${$args{ini_ref}}{input_files}
                             );
   #   }
   #}

   #create a "combined" file in case it doesn't exist
   #if ((${$args{auto_ini_ref}}{use_critica} == 0 && ${$args{auto_ini_ref}}{use_glimmer2} == 0 && ${$args{auto_ini_ref}}{use_glimmer3} == 1) ||
   #    (${$args{auto_ini_ref}}{use_critica} == 0 && ${$args{auto_ini_ref}}{use_glimmer2} == 1 && ${$args{auto_ini_ref}}{use_glimmer3} == 0) ||
   #    (${$args{auto_ini_ref}}{use_critica} == 1 && ${$args{auto_ini_ref}}{use_glimmer2} == 0 && ${$args{auto_ini_ref}}{use_glimmer3} == 0) ) {
   #   if (-e (${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.critica.mod')) {
   #      copy(${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.critica.mod', ${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.combined');
   #   }
   #   if (-e (${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.glcoord.mod')) {
   #      copy(${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.glcoord.mod', ${$args{ini_ref}}{genemodel_output}.'/'.$args{input_file}.'.combined');
   #   }
   #

   #run intergenic Blast if selected
   #if option is selected to remove small ORFs from IG determination, remove from .combined gene model as well!
   if (${$args{auto_ini_ref}}{runintergenicblast} == 1) {
      ${$args{progress_bar}}->configure(-label=>"Checking intergenic regions for $args{input_file}");
      ${$args{main_window}}->update;
      my ($status) = &run_igblast(main_window   => $args{main_window},
                                   progress_bar  => $args{progress_bar},
                                   auto_ini_ref  => $args{auto_ini_ref},
                                   ini_ref       => $args{ini_ref},
                                   input_file    => $args{input_file},
                                   directory     => ${$args{ini_ref}}{input_files}
                                  );
      if ($status == 0) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while checking Intergenic Regions for $args{input_file}. Aborting analysis",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error running Intergenic Regions for file $args{input_file}".
                        "\nError while running RBSfinder for $args{input_file}. Aborting analysis\n\n";
         close ERRORLOG;
         return (0); #return to main
      }
   }

   #run RBS finder if selected
   if (${$args{auto_ini_ref}}{runRBSfinder} == 1) {
      ${$args{progress_bar}}->configure(-label=>"Running RBSfinder for $args{input_file}");
      ${$args{main_window}}->update;
      #run 2 iterations of RBSfinder, RBSoutput is reference to RBSfinder genemodel
      my ($RBSoutput_ref) = &runRBSfinder(main_window   => $args{main_window},
                                          progress_bar  => $args{progress_bar},
                                          auto_ini_ref  => $args{auto_ini_ref},
                                          ini_ref       => $args{ini_ref},
                                          filename      => $args{input_file},
                                          directory     => ${$args{ini_ref}}{input_files}
                                         );
      if (${$RBSoutput_ref} eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while running RBSfinder for $args{input_file}. Aborting analysis",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error running RBSfinder for file $args{input_file}".
                        "\nError while running RBSfinder for $args{input_file}. Aborting analysis\n\n";
         close ERRORLOG;
         return (0); #return to main
      }
      #modify combined genemodel
      ${$args{progress_bar}}->configure(-label=>"Merging RBSfinder into gene model for $args{input_file}");
      ${$args{main_window}}->update;
      &RBS2genemodel(main_window   => $args{main_window},
                     progress_bar  => $args{progress_bar},
                     auto_ini_ref  => $args{auto_ini_ref},
                     ini_ref       => $args{ini_ref},
                     filename      => $args{input_file},
                     directory     => ${$args{ini_ref}}{input_files},
                     RBSoutput     => $RBSoutput_ref
                   );
   }
   ${$args{progress_bar}}->configure(-label=>"Successfully created gene model for $args{input_file}");
   ${$args{main_window}}->update;

   return(1);
}

sub combine_gene_models {
   my %args = @_;
   #my ($glimmer_file_ref, $critica_file_ref, $prodigal_file_ref, @glimmer_gm, @critica_gm, @prodigal_gm, %glimmer_left_bd, %glimmer_right_bd, %glimmer_orientation,
   #    %critica_left_bd, %critica_right_bd, %critica_orientation, %critica_rbs_distance, %critica_rbs_pattern, $max_gl_orf, $cr_orf,
   #    %prodigal_left_bd, %prodigal_right_bd, %prodigal_orientation, %prodigal_rbs_distance, %prodigal_rbs_pattern, $max_gl_orf, $prodigal_orf);
   #my ($glimmer, $critica, $prodigal, $delete_critica, $delete_glimmer, $delete_prodigal);
   my ($combined_gm);
   ${$args{progress_bar}}->configure(-label=>"Merging gene models for $args{filename}");
   ${$args{main_window}}->update;

   #read each gene model individually, then merge into the combined model
   #read Glimmer
   if ((${$args{auto_ini_ref}}{use_glimmer2} == 1 || ${$args{auto_ini_ref}}{use_glimmer3} == 1) &&
       -e ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.'.glcoord.mod') {
      #read gene model
      my ($gene_model_ref) = &read_gene_model(main_window   => $args{main_window},
                                              progress_bar  => $args{progress_bar},
                                              auto_ini_ref  => $args{auto_ini_ref},
                                              ini_ref       => $args{ini_ref},
                                              model_name    => 'Glimmer',
                                              filename      => $args{filename}.'.glcoord.mod',
                                              directory     => ${$args{ini_ref}}{genemodel_output}
                                             );
      #check status, if OK, continue with merging into combined gene model
      if ($gene_model_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing Glimmer gene model for $args{filename}",
                                                       -buttons => ['OK'],
                                                       -bitmap    => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing Glimmer gene model for $args{filename}".
                        "\n\n";
         close ERRORLOG;
      } else {
         #merge with combined GM
         &merge_gene_model(main_window   => $args{main_window},
                           progress_bar  => $args{progress_bar},
                           auto_ini_ref  => $args{auto_ini_ref},
                           ini_ref       => $args{ini_ref},
                           filename      => $args{filename},
                           model_name    => 'Glimmer',
                           gm_ref        => $gene_model_ref,
                           combined_gm   => \$combined_gm
                          );
      }
   }
   #read Critica
   if (${$args{auto_ini_ref}}{use_critica} == 1 &&
       -e ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.'.critica.mod') {
      #read gene model
      my ($gene_model_ref) = &read_gene_model(main_window   => $args{main_window},
                                              progress_bar  => $args{progress_bar},
                                              auto_ini_ref  => $args{auto_ini_ref},
                                              ini_ref       => $args{ini_ref},
                                              model_name    => 'Critica',
                                              filename      => $args{filename}.'.critica.mod',
                                              directory     => ${$args{ini_ref}}{genemodel_output}
                                             );
      #check status, if OK, continue with merging into combined gene model
      if ($gene_model_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing Critica gene model for $args{input_file}",
                                                       -buttons => ['OK'],
                                                       -bitmap    => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing Critica gene model for $args{input_file}".
                        "\n\n";
         close ERRORLOG;
      } else {
         #merge with combined GM
         &merge_gene_model(main_window   => $args{main_window},
                           progress_bar  => $args{progress_bar},
                           auto_ini_ref  => $args{auto_ini_ref},
                           ini_ref       => $args{ini_ref},
                           filename      => $args{filename},
                           model_name    => 'Critica',
                           gm_ref        => $gene_model_ref,
                           combined_gm   => \$combined_gm
                          );
      }
   }

   #read Prodigal
   if (${$args{auto_ini_ref}}{use_prodigal} == 1 &&
       -e ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.'.prodigal.coord.mod') {
      #read gene model
      my ($gene_model_ref) = &read_gene_model(main_window   => $args{main_window},
                                              progress_bar  => $args{progress_bar},
                                              auto_ini_ref  => $args{auto_ini_ref},
                                              ini_ref       => $args{ini_ref},
                                              model_name    => 'Prodigal',
                                              filename      => $args{filename}.'.prodigal.coord.mod',
                                              directory     => ${$args{ini_ref}}{genemodel_output}
                                             );
      #check status, if OK, continue with merging into combined gene model
      if ($gene_model_ref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error parsing Prodigal gene model for $args{input_file}",
                                                       -buttons => ['OK'],
                                                       -bitmap    => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing Prodigal gene model for $args{input_file}".
                        "\n\n";
         close ERRORLOG;
      } else {
         #merge with combined GM
         &merge_gene_model(main_window   => $args{main_window},
                           progress_bar  => $args{progress_bar},
                           auto_ini_ref  => $args{auto_ini_ref},
                           ini_ref       => $args{ini_ref},
                           filename      => $args{filename},
                           model_name    => 'Prodigal',
                           gm_ref        => $gene_model_ref,
                           combined_gm   => \$combined_gm
                          );
      }
   }

   #write combined gene model
   open WRITE, "+>".${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.'.combined';
   my $new_count = 1;
   foreach my $left_bd (sort {$a <=> $b} keys %{$combined_gm}) {
      foreach my $right_bd (keys %{$combined_gm->{$left_bd}}) {
         #sense orientation
         if ($combined_gm->{$left_bd}->{$right_bd} =~ /^\[\+\d/) {
            print WRITE ' '.$new_count.'  '.$left_bd.'  '.$right_bd.'  '.$combined_gm->{$left_bd}->{$right_bd}."\n";
         }
         #antisense orientation
         elsif ($combined_gm->{$left_bd}->{$right_bd} =~ /^\[\-\d/) {
            print WRITE ' '.$new_count.'  '.$right_bd.'  '.$left_bd.'  '.$combined_gm->{$left_bd}->{$right_bd}."\n";
         }
         $new_count++;
      }
   }
   close WRITE;
   undef $combined_gm;
}

sub merge_gene_model {
   my %args = @_;
   my ($left_bd, $right_bd, $delete, $delete_combined);

   #map boundaries from one Gene caller onto the combined model
   ${$args{progress_bar}}->configure(-label=>"Simplifying gene model for $args{filename}");
   ${$args{main_window}}->update;

   #test if more than one gene caller was selected. If not, skip merging
   my $number_of_genemodels = ${$args{auto_ini_ref}}{use_glimmer2} +
                              ${$args{auto_ini_ref}}{use_glimmer3} +
                              ${$args{auto_ini_ref}}{use_prodigal} +
                              ${$args{auto_ini_ref}}{use_critica}  +
                              ${$args{auto_ini_ref}}{runintergenicblast};

   if ($number_of_genemodels > 1) {
      #get boundaries from current gene model ORF
      foreach my $left_bd (keys %{${$args{gm_ref}}}) {
         ${$args{progress_bar}}->configure(-label=>"Simplifying gene model for $args{filename} at position $left_bd");
         ${$args{main_window}}->update;
         #could have more than one shared left_bd!
         foreach my $right_bd (keys %{${$args{gm_ref}}->{$left_bd}}) {
            #get boundaries from combined gene model ORF
            foreach my $combined_left_bd (keys %{${$args{combined_gm}}}) {
               foreach my $combined_right_bd (keys %{${$args{combined_gm}}->{$combined_left_bd}}) {
                  #skip if ORFs are not intersecting
                  next if ($combined_right_bd <= $left_bd || $combined_left_bd >= $right_bd);

                  #delete if identical
                  if ($left_bd == $combined_left_bd && $right_bd == $combined_right_bd) {
                     $delete->{$left_bd} = $right_bd;
                     next;
                  }

                  #delete if same right bd and different left bd.
                  if (defined ${$args{combined_gm}}->{$combined_left_bd}->{$right_bd}) {
                     if    ($left_bd >  $combined_left_bd) {$delete->{$left_bd} = $right_bd}
                     elsif ($left_bd <= $combined_left_bd) {$delete_combined->{$combined_left_bd} = $combined_right_bd};
                     next;
                  }

                  #delete if same left bd
                  if (defined ${$args{combined_gm}}->{$left_bd}) {
                     if    ($right_bd <  $combined_right_bd) {$delete->{$left_bd} = $right_bd}
                     elsif ($right_bd >= $combined_right_bd) {$delete_combined->{$combined_left_bd} = $combined_right_bd};
                     next;
                  }
               }
            }
         }
      }

      #modify hashes
      foreach my $combined_left_bd (keys %{$delete_combined}) {
         delete ${$args{combined_gm}}->{$combined_left_bd}->{ $delete_combined->{$combined_left_bd} };
         my $number = keys %{${$args{combined_gm}}->{$combined_left_bd}};
         if ($number == 0) { delete ${$args{combined_gm}}->{$combined_left_bd} };
      }
      foreach my $left_bd (keys %{$delete}) {
         delete  ${$args{gm_ref}}->{$left_bd}->{ $delete->{$left_bd} };
         my $number =  keys %{${$args{gm_ref}}->{$left_bd}};
         if ($number == 0) { delete ${$args{gm_ref}}->{$left_bd} };
      }
      undef $delete;
      undef $delete_combined;
   }

   #merge remaining orfs between Gene model and combined gene model
   ${$args{progress_bar}}->configure(-label=>"Merging gene models for $args{filename}");
   ${$args{main_window}}->update;
   foreach my $left_bd (keys %{${$args{gm_ref}}}) {
      foreach my $right_bd (keys %{${$args{gm_ref}}->{$left_bd}}) {
         ${$args{combined_gm}}->{$left_bd}->{$right_bd} = ${$args{gm_ref}}->{$left_bd}->{$right_bd};
      }
   }
}

sub read_gene_model {
   my %args = @_;
   my ($file_ref, @gene_model, $gene_model);

   ${$args{progress_bar}}->configure(-label=>"Reading $args{model_name} gene models for $args{filename}");
   ${$args{main_window}}->update;
   ($file_ref) = &slurp(main_window   => $args{main_window},
                        auto_ini_ref  => $args{auto_ini_ref},
                        ini_ref       => $args{ini_ref},
                        filename      => $args{filename},
                        directory     => $args{directory}
                       );
   if ($file_ref eq '0') {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Could not open $args{model_name} gene model file $args{filename}",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error parsing $args{model_name} results for file $args{filename}".
                     "\nCould not open $args{model_name} gene model file $args{filename}\n\n";
      close ERRORLOG;
      return (0); #return to main
   }

   @gene_model = split/\n/,${$file_ref};
   #generate boundary hashes
   foreach my $orf (@gene_model) {
      my ($id, $start, $stop, $rest);
      'reset' =~ m/reset/;
      $orf =~ m/^\s*(\d+)\s+(\d+)\s+(\d+)\s+(.+)/;
      $id    = $1;
      $start = $2;
      $stop  = $3;
      $rest  = $4;

      unless (defined $id && defined $start && defined $stop && defined $rest) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not read $args{model_name} entry $orf. Skipping entry",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing $args{model_name} results for file $args{filename}".
                        "\nCould not read $args{model_name} entry $orf. Skipping entry\n\n";
         close ERRORLOG;
         next;
      }

      #add to left and right boundary hash
      if ($rest =~ /\[\-\d/) {
         $gene_model->{$stop}->{$start} = $rest;
      } else {
         $gene_model->{$start}->{$stop} = $rest;
      }
   }
   return (\$gene_model);
}

sub runRBSfinder {
   my %args = @_;
   my $RBSoutput_ref = 0;
   #my $timescale = '80';

   ${$args{progress_bar}}->configure(-label=>"Running RBSfinder in first iteration for $args{filename}");
   ${$args{main_window}}->update;

   #read gene model file
   my $gene_coord_file_ref =slurp(main_window   => $args{main_window},
                                  auto_ini_ref  => $args{auto_ini_ref},
                                  ini_ref       => $args{ini_ref},
                                  filename      => $args{filename}.'.combined',
                                  directory     => ${$args{ini_ref}}{genemodel_output}
                                 );

   #run 1st iteration of RBSfinder
   $RBSoutput_ref = RBSfinder(main_window    => $args{main_window},
                              progress_bar   => $args{progress_bar},
                              auto_ini_ref   => $args{auto_ini_ref},
                              ini_ref        => $args{ini_ref},
                              filename       => $args{filename},
                              directory      => ${$args{ini_ref}}{input_files},
                              gene_coord_ref => $gene_coord_file_ref
                             );
   #$timescale = (time() - $local_start_time) * 6;
   #if ($timescale == 0) {$timescale = 2};
   #$status->die;
   if (${$RBSoutput_ref} eq '0') {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Error while running 1st RBSfinder for $args{filename}. Aborting analysis",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error running RBSfinder for file $args{filename}".
                     "\nError in 1st iteration\n\n";
      close ERRORLOG;
      return (0); #return to main
   }

   #reformat for 2nd iteration
   my @temp = split /\n/,${$RBSoutput_ref};
   unlink ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.'.rbsfinder';
   my $RBS_to_glimmer = "";
   foreach my $pos (@temp) {
      if ($pos !~ /\w+/) {next};
      my ($orf, $left_bd, $right_bd, $rbs);
      $pos =~ m/\s*(\d+)\s+(\d+)\s+(\d+)\s+([^\s]+)\s+/;
      $orf = $1;
      $left_bd = $2;
      $right_bd = $3;
      $rbs = $4;

      if ($left_bd < $right_bd) { #sense orientation
         my $distance = $right_bd - $left_bd;
         $RBS_to_glimmer .= "    ".$orf."    ".$left_bd."    ".$right_bd.'    [+1 L='.$distance.' r=-1]'."\n";
      } elsif ($left_bd > $right_bd) { #antisense orientation
         my $distance = $right_bd - $left_bd;
         $RBS_to_glimmer .= "    ".$orf."    ".$left_bd."    ".$right_bd.'    [-1 L='.$distance.' r=-1]'."\n";
      }
   }

   ${$args{progress_bar}}->configure(-label=>"Running RBSfinder in second iteration for $args{filename}");
   ${$args{main_window}}->update;

   $RBSoutput_ref = RBSfinder(main_window    => $args{main_window},
                              progress_bar   => $args{progress_bar},
                              auto_ini_ref   => $args{auto_ini_ref},
                              ini_ref        => $args{ini_ref},
                              filename       => $args{filename},
                              directory      => ${$args{ini_ref}}{input_files},
                              gene_coord_ref => \$RBS_to_glimmer
                             );
   #$status->die;

   if (${$RBSoutput_ref} eq '0') {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Error while running 2nd RBSfinder for $args{filename}. Aborting analysis",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error running RBSfinder for file $args{filename}".
                     "\nError in 2nd iteration\n\n";
      close ERRORLOG;
      return (0); #return to main
   }
   open WRITE, "+>".${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.'.rbsfinder';
   print WRITE ${$RBSoutput_ref};
   close WRITE;

   return ($RBSoutput_ref);
}

sub RBS2genemodel {
   my %args = @_;
   #my %glimmer_prob = ();
   my %glimmer_rbs = ();

   #reading combined gene model coord file
   my $gene_coord_file_ref =slurp(main_window   => $args{main_window},
                                  auto_ini_ref  => $args{auto_ini_ref},
                                  ini_ref       => $args{ini_ref},
                                  filename      => $args{filename}.'.combined',
                                  directory     => ${$args{ini_ref}}{genemodel_output}
                                 );
   my @Glimmerentries = ();
   @Glimmerentries = split /\n/, ${$gene_coord_file_ref};
   #catch probabilities and RBS pattern and enter into hash
   foreach (@Glimmerentries) {
      m/^\s*(\d+).*?(\[.+)/;
      #$glimmer_prob{$1} = $2;
      $glimmer_rbs{$1}  = $2;
      unless (defined $glimmer_rbs{$1} && $glimmer_rbs{$1} =~ /\w+/) {$glimmer_rbs{$1} = ""};
   }

   #get all RBS entries from RBSfinder
   my @RBSentries = ();
   @RBSentries = split /\n/, ${$args{RBSoutput}};

   #push into appropriate arrays and hashes
   my @new_gene_model = ();
   foreach (@RBSentries) {
      my $ORF_name    = "";
      my $new_start   = "";
      my $stop        = "";
      my $RBS_pattern = "";
      my $RBS_pos     = "";
      my $old_start   = "";
      my $move        = "false";
      my $complement  = "false";
      if ($_ !~ /^\s*\d+/) {next}; #skip non-entry lines
      #catch relevant start and stop positions
      m/\s*(\d+)\s+(\d+)\s+(\d+)\s+([ACGTXN\-]+)\s+([0-9ACGTN\-]+)\D+\d+\D+(\d+)/;
      $ORF_name    = $1;
      $new_start   = $2;
      $stop        = $3;
      $RBS_pattern = $4;
      $RBS_pos     = $5;
      $old_start   = $6;

      #skip ORFs without RBS
      unless (defined ($RBS_pos)) {
         my $tl = ${$args{main_window}}->Toplevel(-title => 'Error');
         $tl->Label(-text => "$_ is not a valid entry. Skipping") -> pack();
         $tl->raise;
         $tl->update;
         $tl->after(3000);
         $tl->destroy;
         open WRITE, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print WRITE "\nRBSfinder: $_ is not a valid entry. Skipping in file $args{filename}\n";
         close WRITE;
         next;
      }

      #unchanged start position
      if ($RBS_pos =~ /\-\-\-/) {
         #add new gene model to array
         my $temp = "";
         if ($new_start > $stop) { #complement ORF
            my $distance = $new_start - $stop;
            $glimmer_rbs{$ORF_name} =~ s/L\=\d+/L\=$distance/;
            $temp = "    $ORF_name  $new_start  $stop  $glimmer_rbs{$ORF_name}\n";
         } else {
            my $distance = $stop - $new_start;
            $glimmer_rbs{$ORF_name} =~ s/L\=\d+/L\=$distance/;
            $temp = "    $ORF_name  $new_start  $stop  $glimmer_rbs{$ORF_name}\n";
         }
         push (@new_gene_model, $temp);
      } elsif ($RBS_pos =~ /\d+/) { #modified start position
         if ($RBS_pos == 0) {next};

         my $RBS_stop = $RBS_pos+length($RBS_pattern);
         my $RBS_start = $RBS_pos;

         #check if start position has been moved
         if ($new_start != $old_start) {
            $move = "true";
         }

         #check if start < stop and modify accordingly
         if ($new_start > $stop) {
            $RBS_start = $RBS_pos-length($RBS_pattern);
            $RBS_stop = $RBS_pos;
            $complement = "true";
         }

         #add new gene model to array
         my $temp = "";
         if ($new_start > $stop) { #complement ORF
            my $distance = $new_start - $stop;
            $glimmer_rbs{$ORF_name} =~ s/L\=\d+/L\=$distance/;
            $RBS_start = $RBS_start+1;
            $temp = "    $ORF_name  $new_start  $stop  $glimmer_rbs{$ORF_name}\tantisense__$RBS_start\_\_$RBS_stop\_\_$RBS_pattern\n";
         } else {
            my $distance = $stop - $new_start;
            $glimmer_rbs{$ORF_name} =~ s/L\=\d+/L\=$distance/;
            $RBS_stop = $RBS_stop-1;
            $temp = "    $ORF_name  $new_start  $stop  $glimmer_rbs{$ORF_name}\tsense__$RBS_start\_\_$RBS_stop\_\_$RBS_pattern\n";
         }
         push (@new_gene_model, $temp);
      } else { #error
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Not a valid option for RBSfinder in line\n$_.\nSkipping.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "\nRBSfinder: Not a valid option for RBSfinder in line\n$_.\n Skipping in file $args{filename}\n\n";
         close ERRORLOG;
         next;
      }
   }

   #modify Glimmer coord file
   unlink ${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.'.combined';
   open WRITE, "+>".${$args{ini_ref}}{genemodel_output}.'/'.$args{filename}.'.combined';
   foreach (@new_gene_model) {
      print WRITE $_;
   }
   close WRITE;
   return (1);

}

1;
