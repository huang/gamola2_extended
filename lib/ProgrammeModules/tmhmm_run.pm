#!/opt/ActivePerl-5.8/bin/perl

#tmhmm program calls and reformatting of results
#input arguments: main_window, progress_bar, auto_ini_ref, ini_ref,
#                 directory, input_array_ref, output_folder, input_type

package ProgrammeModules::tmhmm_run;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&tmhmm_file);
use vars qw();

use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use ProgrammeModules::tmhmm    qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);
use Cwd;

#local variables
my (%args, @tmhmm_tasks);

sub tmhmm_file {
   my %args = @_;
   my @list = ();
   my @tmhmm_results = ();
   my %exist_tmhmm_file = ();
   my $progress = 0;
   @tmhmm_tasks = ();

   #setup parameters for tmhmm
   my $opt_basedir = ${$args{ini_ref}}{tmhmm_dir};	#MC
   my $opt_d = 0;          # DEBUGGING
   my $opt_workdir = ".";  # Working dir.
   my $opt_wwwdir = ".";   # The place where the www server looks for files
                           # (The www name for the working dir)
   my $opt_serverhome = ".";
   my $opt_html = 0;       # Produce HTML output
   my $opt_short = 0;      # Short output format
   my $opt_plot = 1;       # Produce graphics
   my $opt_v1 = 0;         # Use old model (version 1)

   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "TMHMM analysis",
                   label        => 'TMHMM analysis'
                  );
   &show_pbar_2;

   #slurp up sequence
   my ($seq_ref) = &slurp(main_window  => $args{main_window},
                          progress_bar => $args{progress_bar},
                          auto_ini_ref => $args{auto_ini_ref},
                          ini_ref      => $args{ini_ref},
                          directory    => ${$args{ini_ref}}{input_files},
                          filename     => $args{input_file}
                         );
   if ($seq_ref eq '0') {
      &hide_pbar_2;
      return (0);
   }

   #clean up sequence
   ${$seq_ref} =~ s/^\>[^\n]*?\n//;
   ${$seq_ref} =~ s/\s//gs;

   #map correct entries for respective input file
   @list = grep /$args{input_file}\__/, @{$args{combined_orf_ref}};

   #define working dir for tmhmm
   my $work_dir = ${$args{ini_ref}}{tmhmm_results}.'/TMHMM_'.$args{input_file};

   #remove existing entries if re-use results
   if (${$args{auto_ini_ref}}{reuse_results} == 1) {
      my (@jobs, @TMHMM_results, %exist_tmhmm_file);
      &update_pbar_2(title        => 'TMHMM analysis',
                     label        => "Testing for existing results",
                     progress     => 1,
                    );

      #read existing results
      opendir SEQINPUT, $work_dir;
      @TMHMM_results = grep /_tmhmm$/, readdir(SEQINPUT);
      closedir SEQINPUT;
      #create hash
      foreach (@TMHMM_results) {
         $exist_tmhmm_file{$_} = '1';
      }
      #iterate over gene model
      foreach my $entry (@list) {
         my ($ID, $left_bd, $right_bd, $orientation);
         'reset' =~ m/reset/;
         $entry =~ m/$args{input_file}___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
         $ID = $1;
         $left_bd = $2;
         $right_bd = $3;
         $orientation = $4;
         #skip Blast if result file exists and re-use is activated
         next if (exists $exist_tmhmm_file{$ID.'_'.$left_bd.'_'.$right_bd.'_'.$orientation.'_tmhmm'} && -s $work_dir.'/'.$ID.'_'.$left_bd.'_'.$right_bd.'_'.$orientation.'_tmhmm' > 0);
         #add if required
         push (@jobs, $entry);
      }
      @list = @jobs;
      undef @jobs;
      undef @TMHMM_results;
      undef %exist_tmhmm_file;
   }

   #max number of fasta entries
   my $max_count = $#list + 1;
   if ($max_count < 1) {$max_count = 1};

   #iterate over array
   foreach my $entry (@list) {
      #update progress bar
      $progress++;
      if (($progress % 100) == 0) {
         &update_pbar_2(title        => "TMHMM analysis",
                        label        => "Reading entry for $args{input_file}",
                        progress     => ($progress / $max_count) * 100,
                       );
      }

      #get boundary parameters
      'reset' =~ m/reset/;
      my ($ID, $left_bd, $right_bd, $orientation);
      $entry =~ m/$args{input_file}\___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
      $ID = $1;
      $left_bd = $2;
      $right_bd = $3;
      $orientation = $4;
      my $query_seq_nt = "";

      unless (defined $orientation && $orientation =~ /(sense|antisense)/) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Error parsing gene model for $args{input_file}",
                                           -type    => 'OK',
                                           -icon    => 'info');
         &hide_pbar_2;
         return (0); #return to main
      }

      #extract nt sequence from ORF
      $query_seq_nt = substr(${$seq_ref}, ($left_bd - 1), ($right_bd - $left_bd + 1));

      #translate nt sequence to aa if selected
      my $query_seq_aa = '';
      $query_seq_aa = nt2aa(main_window  => $args{main_window},
                            progress_bar => $args{progress_bar},
                            auto_ini_ref => $args{auto_ini_ref},
                            ini_ref      => $args{ini_ref},
                            module       => 'TMHMM',
                            orientation  => $orientation,
                            sequence     => $query_seq_nt,
                            filename     => $args{input_file},
                            left_bd      => $left_bd,
                            right_bd     => $right_bd
                           );

      #eliminate stop codon for tmhmm
      $query_seq_aa =~ s/\*/X/gs;
      $query_seq_aa =~ s/X$//;
      #create array for new tmhmm tasks
      push (@tmhmm_tasks, $args{input_file}.'___'.$ID.'___'.$left_bd.'___'.$right_bd.'___'.$orientation.'___'.$query_seq_aa);
   }

   #skip exisiting results? if yes, read results from appropriate results directory
   ${$args{progress_bar}}->configure(-label=>" TMHMM ");
   ${$args{main_window}}->update;
   if (${$args{auto_ini_ref}}{reuse_results} == 1) {
      if (-d $work_dir) {
         opendir SEQINPUT, $work_dir;
         @tmhmm_results = grep /_tmhmm$/, readdir(SEQINPUT);

         #create hash
         %exist_tmhmm_file = ();
         foreach (@tmhmm_results) {
            $exist_tmhmm_file{$_} = '1';
         }
         closedir SEQINPUT;
      } else {
         mkdir "$work_dir",0777;
      }
   } else {
      if (-d $work_dir) {
         system "rm -R $work_dir";
         mkdir "$work_dir",0777;
      } else {
         mkdir "$work_dir",0777;
      }
   }

   #start multithreading PFAM hmmer
   for (my $i=1; $i<=$#tmhmm_tasks+2; $i = $i+${$args{auto_ini_ref}}{CPU}) {
      my $count = $i + ${$args{auto_ini_ref}}{CPU};
      if ($count > $#tmhmm_tasks+2) {$count = $#tmhmm_tasks+2};

      &update_pbar_2(title        => "TMHMM analysis",
                     label        => "TMHMM analysis for $args{input_file}",
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
            #grab entries
            my $sequence = "";
            my $left_bd = "";
            my $right_bd = "";
            my $orientation = "";
            my $ORF_number = "";
            $tmhmm_tasks[$j-1] =~ m/$args{input_file}___(\d+)___(\d+)___(\d+)___([a-z]+)___(.+)/;
            $ORF_number = $1;
            $left_bd = $2;
            $right_bd = $3;
            $orientation = $4;
            $sequence = $5;

            #run tmhmm for selected feature
            #generate fasta header
            my $header = "";
            $header = $ORF_number.'_'.$left_bd.'_'.$right_bd.'_'.$orientation;

            #write input file
            open WRITE, "+>$work_dir/$header";
            print WRITE '>'.$header."\n".$sequence;
            close WRITE;
            #call tmhmm
            my ($run_tmhmm) = &new($opt_basedir, $opt_d, $opt_workdir, $opt_wwwdir, $opt_serverhome, $opt_html, $opt_short, $opt_plot, $opt_v1, $header, $work_dir);

            #unlink seq file
            unlink "$work_dir/$header";

            #write results
            open WRITE, "+>$work_dir/$header\_tmhmm";
            print WRITE $header."\n------\n".$run_tmhmm;
            close WRITE;

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
   undef %exist_tmhmm_file;

   &hide_pbar_2;
   return (1);
}


1;