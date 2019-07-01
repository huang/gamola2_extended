#!/opt/ActivePerl-5.8/bin/perl5.8.8

#initialise and modify variables

package initialise::init;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;

$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&read_ini &save_ini &mod_vars %ini %ini_comment %auto_ini $display_parameters $auto_display_parameters $read_file %args) ;
use vars qw(%ini %ini_comment %auto_ini $display_parameters $auto_display_parameters $read_file);

use initialise::read_me;
use initialise::recompile qw(:DEFAULT);
use GeneModel::setup;
use File::Find;
use File::Copy;
use File::Path;

#local variables
my (%args, $save_file, $nw, $function, $key, $value, $file, $comment, @types,
    %ini_gadget, %ini_gadget_1, %ini_gadget_browse, $setup_tl, $setup_tw, $setup_frame, $setup_browse,
    $user, $auto, $local_databases);

@types = (["Config Files", '.default', 'TEXT'], ["All Files", "*"] );


sub read_ini {
   my %args = @_;

   if (defined ${$args{progress_bar}}) {
      ${$args{progress_bar}}->configure(-label=>'Reading parameters');
   }
   ${$args{main_window}}->update;

   #reset display paramters
   @{$args{display_parameters}} = ();
   @{$args{auto_display_parameters}} = ();

   #read ini file
   unless (-e $args{file}) {
      $args{file} = ${$args{main_window}}->getOpenFile(-filetypes   => \@types,
                                                       -initialfile => $args{file});
      unless (defined $args{file} && -e $args{file}) {return(0)};
   };

   {
      local( $/, *ENTRY ) ;
      open( ENTRY, $args{file} ) or die "\nCant read file $args{file}\n";
      $file = <ENTRY>;
      close ENTRY;
   }

   #processing
   $file =~ m/(.*?)####################(.+)/gs;
   $user = $1;
   $auto = $2;
   unless (defined $user && defined $auto) {
      die "\ncannot untangle user and auto setting in $args{file}\n";
   }
   $user = "\n".$user;

   #generate user hash
   'reset' =~ m/reset/;
   while ($user =~ m/\n\s*?\#([^\n\r]*?)[\n\r]\s*(\S*?)\s+([^\n\r]+)/gs) {
      $comment = $1;
      $function = $2;
      $value = $3;
      unless (defined ($comment) && defined ($function) && defined ($value)) {
         print "\nCan't process entry $_ \n...$comment ...$function ...$value ...\n";
         die;
      }

      #clear pathway
      $value =~ s/\s*$//;
      $value =~ s/[\/\\]$//;

      ${$args{ini_ref}}{$function} = $value;
      ${$args{ini_comment_ref}}{$function} = $comment;
      #generate ordered hash for display
      push (@{$args{display_parameters}}, $function);

      undef $comment;
      undef $function;
      undef $value;
   }

   #generate auto_ini
   $auto =~ s/# Don't modify below this point\s*[\n\r]+//s;
   $auto = "\n".$auto;
   my @order;
   while ($auto =~ m/\n+\s*(\S*?)\s+([^\n\r]+)/gs) {
      my ($feature, $value);
      $feature = $1;
      $value = $2;
      if ($feature =~ /genbank_header/) {
         my @header = split/ ; /, $value;
         foreach my $entry (@header) {
            my ($key, $value);
            'reset' =~ /reset/;
            $entry =~ m/(\w+)\s+(.+)/;
            $key = $1;
            $value = $2;
            unless (defined $value) {
               ${$args{main_window}}->messageBox(-title   => 'Error parsing initial setup',
                                                 -message => "Error parsing entry $entry",
                                                 -type    => 'Error'
                                                );
            }
            ${$args{auto_ini_ref}}{$key} = $value;
            #resubstitute new_lines
            if (${$args{auto_ini_ref}}{$key} =~ /__newline__/) {
               ${$args{auto_ini_ref}}{$key} =~ s/__newline__/\n/gs;
            }

            #generate order array
            push (@order, $key);

         }
         ${$args{auto_ini_ref}}{'gb_order'} = \@order;
         next;
      }
      ${$args{auto_ini_ref}}{$feature} = $value;
      push (@{$args{auto_display_parameters}}, $feature);
   }

   #add current working directory to auto_ini
   if (defined $args{home_dir}) {
      ${$args{auto_ini_ref}}{work_dir} = $args{home_dir};
   }

   #configure blast filter for empty string
   if (${$args{auto_ini_ref}}{filter_blast} =~ /^---$/) {
      ${$args{auto_ini_ref}}{filter_blast} = "";
   }

   #configure SOAP::lite for proxy
   if (${$args{ini_ref}}{proxy} =~ /\w+/ ) {
      $ENV{HTTP_proxy} = ${$args{ini_ref}}{proxy};
   }

   #configure critica environment variables
   if (${$args{auto_ini_ref}}{scoring_matrix_file} =~ /^---$/) {
      ${$args{auto_ini_ref}}{scoring_matrix_file} = "";
   }
   if (${$args{auto_ini_ref}}{dicodon_scores_file} =~ /^---$/) {
      ${$args{auto_ini_ref}}{dicodon_scores_file} = "";
   }
   if (${$args{auto_ini_ref}}{init_scores_file} =~ /^---$/) {
      ${$args{auto_ini_ref}}{init_scores_file} = "";
   }
   if (${$args{auto_ini_ref}}{sd_scores_file} =~ /^---$/) {
      ${$args{auto_ini_ref}}{sd_scores_file} = "";
   }
   if (${$args{auto_ini_ref}}{prom_scores_file} =~ /^---$/) {
      ${$args{auto_ini_ref}}{prom_scores_file} = "";
   }

   #add critica paths to environment
   if ($ENV{'PATH'} !~ /${$args{auto_ini_ref}}{critica_bin}/) {
      $ENV{'PATH'} = ${$args{auto_ini_ref}}{critica_bin}.':'.${$args{auto_ini_ref}}{critica_scripts}.':'.$ENV{'PATH'};
   }

   #add tRNAscan paths to environment
   if ($ENV{'PATH'} !~ /${$args{ini_ref}}{trnascan_dir}/) {
      $ENV{'PATH'} = ${$args{ini_ref}}{trnascan_dir}.':'.$ENV{'PATH'};
   }

   #update internal and external gene model choices
   if (defined $args{path_widget}) {
      ${$args{path_widget}}->configure(-text => $auto_ini{ext_gm_folder},
                                      );
   }

   ${$args{main_window}}->update;

   return (1);
}

sub save_ini {
   my %args = @_;
   chdir $auto_ini{work_dir};
   $args{file} = ${$args{main_window}}->getSaveFile(-filetypes        => \@types,
                                                    -initialfile      => $args{file},
                                                    -defaultextension => '.default',
                                                    -title            => 'Saving GAMOLA default configuration file'
                                                   );
   unless (defined $args{file}) {return};

   open WRITE, "+>$args{file}" or die "\nCan't access $args{file} for writing";

   #write general parameters
   foreach $function (@{$args{display_parameters}}) {
      print WRITE "\#".${$args{ini_comment_ref}}{$function}."\n";
      print WRITE $function."\t".${$args{ini_ref}}{$function}."\n\n";
   }

   #configure blast filter for empty string
   if (${$args{auto_ini_ref}}{filter_blast} =~ /^\s*$/) {
      ${$args{auto_ini_ref}}{filter_blast} = '---';
   }
   #configure critica paramters for empty string
   if (${$args{auto_ini_ref}}{scoring_matrix_file} =~ /^\s*$/) {
      ${$args{auto_ini_ref}}{scoring_matrix_file} = '---';
   }
   if (${$args{auto_ini_ref}}{dicodon_scores_file} =~ /^\s*$/) {
      ${$args{auto_ini_ref}}{dicodon_scores_file} = '---';
   }
   if (${$args{auto_ini_ref}}{init_scores_file} =~ /^\s*$/) {
      ${$args{auto_ini_ref}}{init_scores_file} = '---';
   }
   if (${$args{auto_ini_ref}}{sd_scores_file} =~ /^\s*$/) {
      ${$args{auto_ini_ref}}{sd_scores_file} = '---';
   }
   if (${$args{auto_ini_ref}}{prom_scores_file} =~ /^\s*$/) {
      ${$args{auto_ini_ref}}{prom_scores_file} = '---';
   }

   #write auto parameters
   print WRITE "\n\n####################\n# Don't modify below this point\n\n";
   foreach $function (@{$args{auto_display_parameters}}) {
      print WRITE $function."\t".${$args{auto_ini_ref}}{$function}."\n";
   }

   #consolidate default Genbank header
   my $header = '';
   foreach my $value (@{${$args{auto_ini_ref}}{'gb_order'}}) {
      #subtitute new_lines
      if (${$args{auto_ini_ref}}{$value} =~ /\n\r?/) {
         ${$args{auto_ini_ref}}{$value} =~ s/\n\r?/__newline__/gs;
      }
      ${$args{auto_ini_ref}}{$value} =~ s/__newline__$//;
      #remove leading spaces
      ${$args{auto_ini_ref}}{$value} =~ s/^\s+//;

      $header .= $value.'  '.${$args{auto_ini_ref}}{$value}.' ; ';
   }
   $header =~ s/ ; $//;
   print WRITE 'genbank_header'."\t".$header;
   close WRITE;
   undef $header;
   return;
}

sub mod_vars {
   my %args = @_;
   if ($args{file} =~ /\w+/) {
      &setup (main_window             => $args{main_window},
              file                    => $args{file},
              ini_ref                 => $args{ini_ref},
              auto_ini_ref            => $args{auto_ini_ref},
              ini_comment_ref         => $args{ini_comment_ref},
              display_parameters      => $args{display_parameters},
              auto_display_parameters => $args{auto_display_parameters}
             );
   }
}

sub setup {
   my %args      = @_;
   my $max_width = 20;
   if (!Tk::Exists ($setup_tl)) {
      &setup_tl (main_window             => $args{main_window},
                 file                    => $args{file},
                 ini_ref                 => $args{ini_ref},
                 auto_ini_ref            => $args{auto_ini_ref},
                 ini_comment_ref         => $args{ini_comment_ref},
                 display_parameters      => $args{display_parameters},
                 auto_display_parameters => $args{auto_display_parameters}
                );
   } else {
      $setup_tl->state('normal');
   }

   #push values of ini hash into Text Widget
   foreach $function (@{$args{display_parameters}}) {
      #define required max width
      if ($max_width < length(${$args{ini_ref}}{$function})) { $max_width = length(${$args{ini_ref}}{$function})};

      if (! Tk::Exists($ini_gadget{$function})) {
         #add label
         $ini_gadget{$function} = $setup_tw->Label(-text    => "${$args{ini_comment_ref}}{$function}",
                                                   -justify => "right",
                                                   -relief  => 'groove',
                                                   -width   => 60);
         if ($function =~ /(cutoff|standard|internet|internal|structural)/) {
            $ini_gadget{$function}->configure(-background => 'grey');
         }
         if (${$args{ini_ref}}{$function} =~ /[\/]+/) {
            $ini_gadget_browse{$function} = $setup_tw->Button(-bitmap  => '@'.${$args{auto_ini_ref}}{work_dir}.'/lib/initialise/tree.xbm',
                                                              -command => sub {
                                                                                my $dir = &dir_select(main_window => $args{main_window},
                                                                                                      label       => 'Select location of Blast databases'
                                                                                                     );
                                                                                if (defined $dir && $dir =~ /\w+/) {
                                                                                   ${$args{ini_ref}}{$function} = $dir;
                                                                                   $setup_tw->update;
                                                                                }
                                                                               }
                                                                 );
         } else {
            $ini_gadget_browse{$function} = $setup_tw->Label(-bitmap => '@'.${$args{auto_ini_ref}}{work_dir}.'/lib/initialise/empty.xbm'); #dummy
         }
         $setup_tw     -> windowCreate('end',-window => $ini_gadget{$function});
         $setup_tw     -> windowCreate('end',-window => $ini_gadget_browse{$function});
      } else {
         $setup_tw ->update;
      }
      #add values
      if (! Tk::Exists($ini_gadget_1{$function})) {
         $ini_gadget_1{$function} = $setup_tw->Entry(-textvariable => \${$args{ini_ref}}{$function});
         if ($function =~ /(cutoff|standard|internet|internal|structural)/) {
            $ini_gadget_1{$function}->configure(-background => 'grey');
         }
         $setup_tw->windowCreate('end', -window => $ini_gadget_1{$function});
      } else {
         $ini_gadget_1{$function} -> update;
         $setup_tw -> update;
         next;
      }

      $setup_tw->insert('end', "\n");
   }

   #set appropriate width
   foreach $function (@{$args{display_parameters}}) {
      $ini_gadget_1{$function}->configure(-width => $max_width);
   }

}

sub setup_tl {
   my %args = @_;

   $setup_tl= ${$args{main_window}}->Toplevel(-title => 'Variable Setup');
   $setup_tl->minsize(300,300);
   $setup_tl->geometry('-200+500');
   $setup_tl->optionAdd('*BorderWidth' => 1);
   $setup_tl->focus;
   $setup_frame = $setup_tl->Frame->pack (-side => 'bottom');
   $setup_frame -> Button (-text => "Exit",
                           -command => sub { $setup_tl -> state('withdraw'); return; }) -> pack (-side => 'left');
   $setup_frame -> Button (-text => "Save",
                           -command => sub { &save_ini (main_window             => $args{main_window},
                                                        file                    => $args{file},
                                                        ini_ref                 => $args{ini_ref},
                                                        auto_ini_ref            => $args{auto_ini_ref},
                                                        ini_comment_ref         => $args{ini_comment_ref},
                                                        display_parameters      => $args{display_parameters},
                                                        auto_display_parameters => $args{auto_display_parameters}
                                                       );
                                             $setup_tl ->state('withdraw');
                                             return;
                                           }) -> pack (-side => 'left');
   $setup_frame -> Button (-text => "Generate Standard Setup",
                           -command => sub { &standard_setup (main_window             => $args{main_window},
                                                              file                    => $args{file},
                                                              ini_ref                 => $args{ini_ref},
                                                              auto_ini_ref            => $args{auto_ini_ref},
                                                              ini_comment_ref         => $args{ini_comment_ref},
                                                              display_parameters      => $args{display_parameters},
                                                              auto_display_parameters => $args{auto_display_parameters});
                                           }) -> pack (-side => 'left');

   $setup_tw = $setup_tl ->Scrolled("Text", -width => 120,
                                            -wrap => 'none',
                                            -scrollbars => 'sw'
                                            ) -> pack(-expand => 1, -fill => 'both');
}

sub standard_setup {
   my %args = @_;
   my (%location, $analyses, %root_dirs, %secondary_dirs);

   $root_dirs{'Results'}             = '/Results';
   $root_dirs{'Input sequences'}     = '/Input_sequences';
   $root_dirs{'Programmes'}          = '/Programmes/';

   #create two default root directories
   while (my ($key, $value) = each(%root_dirs)) {
      unless (-d ${$args{auto_ini_ref}}{work_dir}.$value) {
         mkdir (${$args{auto_ini_ref}}{work_dir}.$value, 0777) or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot create $key directory",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return;
         };
      }
   }
   unless (-d ${$args{auto_ini_ref}}{work_dir}.'/Programmes/GeneModelPrediction') {
      mkdir (${$args{auto_ini_ref}}{work_dir}.'/Programmes/GeneModelPrediction', 0777) or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create GeneModelPrediction directory",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return;
      };
   }
   ${$args{ini_ref}}{'input_files'} = ${$args{auto_ini_ref}}{work_dir}.'/Input_sequences';

   #default hash for directories
   $location{'signalp_dir'}             = '/Programmes/signalp';
   $location{'trnascan_dir'}            = '/Programmes/tRNA_Scan_SE';
   $location{'tmhmm_dir'}               = '/Programmes/tmhmm';
   $location{'transterm_dir'}           = '/Programmes/TransTerm';
   $location{'CRISPR_dir'}              = '/Programmes/CRISPR';
   $location{'CRISPR_results'}          = '/Results/CRISPR_results';
   $location{'infernal_dir'}            = '/Programmes/noncodingRNA';
   $location{'blast_executables'}       = '/Programmes/Blast';
   $location{'blast_plus_executables'}  = '/Programmes/Blast_plus';
   $location{'pfam_executable'}         = '/Programmes/hmmer';
   $location{'pfam3_executable'}        = '/Programmes/hmmer3';
   $location{'glimmer_model'}           = '/Programmes/GeneModelPrediction/Glimmer_models';
   $location{'glprog2'}                 = '/Programmes/GeneModelPrediction/Glimmer2';
   $location{'glprog3'}                 = '/Programmes/GeneModelPrediction/Glimmer3';
   $location{'critica'}                 = '/Programmes/GeneModelPrediction/critica';
   $location{'prodigal'}                = '/Programmes/GeneModelPrediction/Prodigal';
   $location{'ig_results'}              = '/Results/IG_Blast_results';
   $location{'blast_results'}           = '/Results/Blast_results';
   $location{'COG_results'}             = '/Results/COG_results';
   $location{'TIGRfam_results'}         = '/Results/TIGRfam_results';
   $location{'session_comparison'}      = '/Results/Session_comparison';
   $location{'trnascan_results'}        = '/Results/tRNAscan_results';
   $location{'rfam_results'}            = '/Results/non_coding_RNA_results';
   $location{'signalp_results'}         = '/Results/SignalP_results';
   $location{'transterm_results'}       = '/Results/TransTerm_results';
   $location{'results'}                 = '/Results/Genbank_annotation';
   $location{'genemodel_output'}        = '/Results/gene_models';
   $location{'pfam_results'}            = '/Results/PFam_results';
   $location{'annotation_transfer_dir'} = '/Results/Annotation_transfer';
   $location{'tmhmm_results'}           = '/Results/tmhmm_results';
   $location{'move_msfasta'}            = '/Input_sequences/ms_fasta_files';
   $location{'move_gb'}                 = '/Input_sequences/ms_fasta_files';
   $location{'rrna_results'}            = '/Results/rRNA_results';
   $location{'vector_results'}          = '/Results/vector_results';
   $location{'Sequin_submission'}       = '/Results/Sequin_submission';
   #$location{'tracker_dir'}             = '/Results/Tracker';

   while (($key, $value) = each(%location)) {
      unless (-d ${$args{auto_ini_ref}}{work_dir}.$value) {
         mkdir (${$args{auto_ini_ref}}{work_dir}.$value, 0777) or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot create $key directory",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return;
         };
      }
      ${$args{ini_ref}}{$key} = ${$args{auto_ini_ref}}{work_dir}.$value;
   }


   #tempfile
   ${$args{ini_ref}}{tempfile}     = ${$args{auto_ini_ref}}{work_dir}.'/temp';

   #define new locations for databases
   $analyses->{'Blast databases'}                 = {(path => 'blast_db_path',
                                                      name => 'Blast_db')};
   $analyses->{'COG databases'}                   = {(path => 'COG_db_path',
                                                      name => 'COG_db')};
   $analyses->{'PFam databases'}                  = {(path => 'pfam_db_path',
                                                      name => 'Pfam_db')};
   $analyses->{'Pfam descriptor'}                 = {(path => 'Pfam_descriptor',
                                                      name => 'Pfam_db')};
   $analyses->{'InterPro descriptors'}            = {(path => 'Interpro_descriptor',
                                                      name => 'Interpro_db')};
   $analyses->{'TIGRfam databases'}               = {(path => 'TIGRfam_db_path',
                                                      name => 'TIGRfam_db')};
   $analyses->{'Gene Ontology descriptor'}        = {(path => 'GO_db_path',
                                                      name => 'TIGRfam_db')};
   $analyses->{'Non-coding RNA (Rfam) databases'} = {(path => 'rfam_db_path',
                                                      name => 'Rfam_db')};
   $analyses->{'Ribosomal RNA database'}          = {(path => 'rrna_db_path',
                                                      name => 'rRNA_db')};

   $setup_tl->state('withdrawn');

   #create a set directories in default location or let user select
   my $dbs = ${$args{main_window}}->MesgBox(-title     => 'Database directories',
                                            -text      => "Choose your own directories [Custom] or use default locations [Default]",
                                            -icon      => 'QUESTION',
                                            -buttons   => ['Custom', 'Default'],
                                            -defbutton => 'Default',
                                            );
   my $db_choice = $dbs->Show;

   if ($db_choice eq 'Custom') {
      foreach my $key (keys (%{$analyses})) {
         ${$args{ini_ref}}{$analyses->{$key}->{'path'}} = &dir_select(main_window => $args{main_window},
                                                                      label       => "Select location of $key"
                                                                     );
      }
   } elsif ($db_choice eq 'Default') {
      foreach my $key (keys (%{$analyses})) {
         unless (-d ${$args{ini_ref}}{$analyses->{$key}->{'path'}}) {
            ${$args{ini_ref}}{$analyses->{$key}->{'path'}} = ${$args{auto_ini_ref}}{work_dir}.'/'.$analyses->{$key}->{'name'};
            ${$args{ini_ref}}{$analyses->{$key}->{'path'}} =~ s/(${$args{auto_ini_ref}}{work_dir})+/${$args{auto_ini_ref}}{work_dir}/; #remove extra instances that may creep up
            unless (-d ${$args{ini_ref}}{$analyses->{$key}->{'path'}}) { #second check in case there is a duplicate directory (i.e. Pfam_db)
               mkdir (${$args{ini_ref}}{$analyses->{$key}->{'path'}}, 0777) or do {
                  ${$args{main_window}}->messageBox(-title   => 'Error',
                                                    -message => "Cannot create standard $key directory: ...${$args{ini_ref}}{$analyses->{$key}->{'path'}}...",
                                                    -icon    => 'error',
                                                    -type    => 'ok');
                  return;
               }
            }
         }
      }
   }

   $setup_tl->state('normal');
   #set default vector_db path to Blast db path
   ${$args{ini_ref}}{'vector_db_path'}      = ${$args{ini_ref}}{'blast_db_path'};
   ${$args{auto_ini_ref}}{'full_vector_db'} = ${$args{ini_ref}}{'blast_db_path'}.'/UniVec';

   #setup subdirectories for COG databases
   my @COG_dbs = qw(COG2003 COG2008 COG2014 arCOG arCOG2014 POG2013);
   foreach my $COG_db (@COG_dbs) {
      unless (-d ${$args{ini_ref}}{'COG_db_path'}.'/'.$COG_db) {
         mkdir (${$args{ini_ref}}{'COG_db_path'}.'/'.$COG_db, 0777) or do {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Cannot create standard COG directory: ${$args{ini_ref}}{'COG_db_path'}\/$COG_db.",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return;
         };
      }
   }

   #check for TIGR info separately
   unless (-d ${$args{ini_ref}}{TIGRfam_db_path}.'/TIGRinfo') {
      mkdir (${$args{ini_ref}}{TIGRfam_db_path}.'/TIGRinfo', 0777)  or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create TIGRfam info directory in ${$args{ini_ref}}{TIGRfam_db_path}",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return;
      };
   }
   ${$args{ini_ref}}{TIGRfam_info} = ${$args{ini_ref}}{TIGRfam_db_path}.'/TIGRinfo';

   #check for non-coding RNA CM dataset separately
   unless (-d ${$args{ini_ref}}{rfam_db_path}.'/cm') {
      mkdir (${$args{ini_ref}}{rfam_db_path}.'/cm', 0777)  or do {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Cannot create RFam cm dataset directory in ${$args{ini_ref}}{rfam_db_path}",
                                           -icon    => 'error',
                                           -type    => 'ok');
         return;
      };
   }
   ${$args{ini_ref}}{rfam_cm_path} = ${$args{ini_ref}}{rfam_db_path}.'/cm';

   #define auto_ini_ref directories
   ${$args{auto_ini_ref}}{'selected_gl_model'} = ${$args{auto_ini_ref}}{work_dir}.'/Programmes/GeneModelPrediction/Glimmer_models/E.coli.model2';
   ${$args{auto_ini_ref}}{'ext_gm_folder'}     = ${$args{auto_ini_ref}}{work_dir}.'/Results/gene_models';
   ${$args{auto_ini_ref}}{'full_critica_db'}   = ${$args{ini_ref}}{'blast_db_path'}.'/nt';
   ${$args{auto_ini_ref}}{'critica_bin'}       = ${$args{ini_ref}}{'critica'}.'/bin';
   ${$args{auto_ini_ref}}{'intergenicblastdb'} = ${$args{ini_ref}}{'blast_db_path'}.'/nr';
   ${$args{auto_ini_ref}}{'full_blast_db'}     = ${$args{ini_ref}}{'blast_db_path'}.'/nr';
   ${$args{auto_ini_ref}}{'full_COG_db'}       = ${$args{ini_ref}}{'COG_db_path'}.'/COG2003/COG';
   ${$args{auto_ini_ref}}{'full_Pfam_db'}      = ${$args{ini_ref}}{'pfam_db_path'}.'/Pfam-A.hmm';
   ${$args{auto_ini_ref}}{'full_TIGRfam_db'}   = ${$args{ini_ref}}{'TIGRfam_db_path'}.'/TIGRFAMs_15.0_HMM.LIB';

   ${$args{auto_ini_ref}}{'gl_short_file'}       = 'E.coli.model2';
   ${$args{auto_ini_ref}}{'selected_critica_db'} = 'nt';
   ${$args{auto_ini_ref}}{'blast_db'}            = 'nr';
   ${$args{auto_ini_ref}}{'COG_db'}              = 'COG';
   ${$args{auto_ini_ref}}{'Pfam_db'}             = 'Pfam-A.hmm';
   ${$args{auto_ini_ref}}{'TIGRfam_db'}          = 'TIGRFAMs_15.0_HMM.LIB';

   #set correct Blast executable directory
   find (sub {${$args{ini_ref}}{'blast_executables'} = $File::Find::dir if (/^blastall$/i)}, ${$args{ini_ref}}{'blast_executables'} );

   #set correct Blast executable directory
   find (sub {${$args{ini_ref}}{'blast_plus_executables'} = $File::Find::dir if (/^blastp$/i)}, ${$args{ini_ref}}{'blast_plus_executables'} );

   #move database files from the separate 'databases.rar' file if present.
   find (sub {$local_databases = $File::Find::dir if (/^databases.rar$/i)}, ${$args{auto_ini_ref}}{work_dir} );
   #make a temp directory
   if (defined $local_databases && -e $local_databases.'/databases.rar') {
      #test if all databases are present
      unless (-e ${$args{ini_ref}}{'pfam_db_path'}.'/relnotes.txt'        &&
              -e ${$args{ini_ref}}{'Interpro_descriptor'}.'/interpro.txt' &&
              -e ${$args{ini_ref}}{'rfam_db_path'}.'/rfam.txt'            &&
              -e ${$args{ini_ref}}{'TIGRfam_db_path'}.'/TIGR_ROLE_NAMES') {

         ${$args{main_window}}->messageBox(-title   => 'Extracting local databases',
                                           -message => "Found the local databases file. Extraction will take a few minutes and will begin after clicking on 'OK'\nBe patient :-)",
                                           -icon    => 'info',
                                           -type    => 'ok');

         mkdir ${$args{auto_ini_ref}}{work_dir}.'/tmp_databases';
         chdir ${$args{auto_ini_ref}}{work_dir}.'/tmp_databases';
         `unrar x -o+ $local_databases/databases.rar`;

         #grab all directory names
         opendir INPUT, ${$args{auto_ini_ref}}{work_dir}.'/tmp_databases';
         my @tmp_db_dirs = grep !/^\./, readdir(INPUT);
         closedir INPUT;

         #distribute directory contents to respective databse folders
         foreach my $db_type (@tmp_db_dirs) {
            chdir ${$args{auto_ini_ref}}{work_dir}.'/tmp_databases/'.$db_type;
            opendir INPUT, ${$args{auto_ini_ref}}{work_dir}.'/tmp_databases/'.$db_type;
            my @tmp_db_files = grep !/^\./, readdir(INPUT);
            closedir INPUT;
            #blast database files
            if ($db_type =~ m/^Swissprot/) {
               foreach my $file (@tmp_db_files) {
                  unlink ${$args{ini_ref}}{'blast_db_path'}.'/'.$file;
                  move($file, ${$args{ini_ref}}{'blast_db_path'}.'/'.$file);
               }
            }
            #PFam database files
            if ($db_type =~ m/^PFam/) {
               foreach my $file (@tmp_db_files) {
                  unlink ${$args{ini_ref}}{'pfam_db_path'}.'/'.$file;
                  move($file, ${$args{ini_ref}}{'pfam_db_path'}.'/'.$file);
               }
            }
            #TIGRfam database files
            if ($db_type =~ m/^TIGRfam/) {
               foreach my $file (@tmp_db_files) {
                  next if (-d $file);
                  unlink ${$args{ini_ref}}{'TIGRfam_db_path'}.'/'.$file;
                  move($file, ${$args{ini_ref}}{'TIGRfam_db_path'}.'/'.$file);
               }
               #move TIGRinfo files
               opendir INPUT, ${$args{auto_ini_ref}}{work_dir}.'/tmp_databases/'.$db_type.'/TIGRinfo';
               my @tmp_tigr_files = grep !/^\./, readdir(INPUT);
               closedir INPUT;
               unless (-d ${$args{ini_ref}}{'TIGRfam_db_path'}.'/TIGRinfo') {
                  mkdir ${$args{ini_ref}}{'TIGRfam_db_path'}.'/TIGRinfo';
               }
               foreach my $file (@tmp_tigr_files) {
                  unlink ${$args{ini_ref}}{'TIGRfam_db_path'}.'/TIGRinfo/'.$file;
                  move(${$args{auto_ini_ref}}{work_dir}.'/tmp_databases/'.$db_type.'/TIGRinfo/'.$file, ${$args{ini_ref}}{'TIGRfam_db_path'}.'/TIGRinfo/'.$file);
               }
            }
            #RFam database files
            if ($db_type =~ m/^RFam/) {
               foreach my $file (@tmp_db_files) {
                  unlink ${$args{ini_ref}}{'rfam_db_path'}.'/'.$file;
                  move($file, ${$args{ini_ref}}{'rfam_db_path'}.'/'.$file);
               }
            }
            #Interpro database files
            if ($db_type =~ m/^Interpro/) {
               foreach my $file (@tmp_db_files) {
                  unlink ${$args{ini_ref}}{'Interpro_descriptor'}.'/'.$file;
                  move($file, ${$args{ini_ref}}{'Interpro_descriptor'}.'/'.$file);
               }
            }
         }
         chdir ${$args{auto_ini_ref}}{work_dir};
         rmtree(${$args{auto_ini_ref}}{work_dir}.'/tmp_databases');
      }
   }


   #compile programmes
   &recompile(progress_bar  => $args{progress_bar},
              main_window   => $args{main_window},
              auto_ini_ref  => $args{auto_ini_ref},
              ini_ref       => $args{ini_ref},
              force_compile => 1
             );

   #save new default ini file
   chdir ${$args{auto_ini_ref}}{work_dir};
   &save_ini (main_window             => $args{main_window},
              file                    => ${$args{auto_ini_ref}}{work_dir}.'/Gamola.default',
              ini_ref                 => $args{ini_ref},
              auto_ini_ref            => $args{auto_ini_ref},
              ini_comment_ref         => $args{ini_comment_ref},
              display_parameters      => $args{display_parameters},
              auto_display_parameters => $args{auto_display_parameters});

   $setup_tl->state('withdrawn');

   &setup    (main_window             => $args{main_window},
              file                    => ${$args{auto_ini_ref}}{work_dir}.'/Gamola.default',
              ini_ref                 => $args{ini_ref},
              auto_ini_ref            => $args{auto_ini_ref},
              ini_comment_ref         => $args{ini_comment_ref},
              display_parameters      => $args{display_parameters},
              auto_display_parameters => $args{auto_display_parameters});

   #remove last saved to establish default values
   unlink ${$args{auto_ini_ref}}{work_dir}.'/last.saved';

   ${$args{main_window}}->messageBox(-title   => 'Generated default setup',
                                     -message => "Most directories have been setup now. Modify remaining values to your specifications",
                                     -icon    => 'info',
                                     -type    => 'ok');

   return;
}


1;