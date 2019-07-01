#!/opt/ActivePerl-5.8/bin/perl

#Transterm program calls and reformatting of results
#input arguments: main_window, progress_bar, auto_ini_ref, ini_ref,
#                 directory, input_array_ref, output_folder, input_type

package ProgrammeModules::transterm;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&transterm);
use vars qw();

use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);
use Cwd;
use File::Copy;

#local variables
my (%args, $newdir);

sub transterm {
   my %args = @_;
   my (@list);

   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Transterm analysis",
                   label        => 'Transterm analysis'
                  );
   &show_pbar_2;

   #remove existing entries if re-use results
   if (${$args{auto_ini_ref}}{reuse_results} == 1) {
      my (@transterm_results, %exist_transterm_file);
      &update_pbar_2(title        => 'Transterm analysis',
                     label        => "Testing for existing results",
                     progress     => 1,
                    );

      #read existing results
      opendir SEQINPUT, ${$args{ini_ref}}{transterm_results};
      @transterm_results = grep /^.*_transterm$/, readdir(SEQINPUT);
      closedir SEQINPUT;
      #create hash
      foreach (@transterm_results) {
         $exist_transterm_file{$_} = '1';
      }
      #iterate over gene model
      foreach my $entry (@{ $args{input_array_ref} }) {
         #skip if result file exists and re-use is activated
         next if (exists $exist_transterm_file{$entry.'_transterm'} && -s ${$args{ini_ref}}{transterm_results}.'/'.$entry.'_transterm' > 0);
         #add if required
         push (@list, $entry);
      }
      undef @transterm_results;
      undef %exist_transterm_file;
   } else {
      @list = @{ $args{input_array_ref} };
   }

   #max number of fasta entries
   my $max_count = $#list + 1;
   if ($max_count < 1) {$max_count = 1};

   #start multithreading transterm
   for (my $i=0; $i<=$#list + 1; $i = $i+${$args{auto_ini_ref}}{CPU}) {
      my $count = $i + ${$args{auto_ini_ref}}{CPU};
      if ($count > $#list + 1) {$count = $#list + 1};

      &update_pbar_2(title        => "Transterm analysis",
                     label        => "Transterm analysis $count of $max_count",
                     progress     => ($count / $max_count) * 100,
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
            my (@job_list, @terminators, $gene_model, $transterm);

            #map correct entries for respective input file
            @job_list = grep /${$args{input_array_ref}}[$j]\__/, @{$args{combined_orf_ref}};

            #change header to match
            my ($seq_ref) = &slurp_cmd(auto_ini_ref => $args{auto_ini_ref},
                                       ini_ref      => $args{ini_ref},
                                       directory    => ${$args{ini_ref}}{input_files},
                                       filename     => $list[$j]
                                      );
            if ($seq_ref eq '0') {
              #exit safely
               CORE::exit();
            }

            ${$seq_ref} =~ s/^>.*?\n/\>${$args{input_array_ref}}[$j] \n/;
            open WRITE, "+>${$args{ini_ref}}{transterm_dir}\/${$args{input_array_ref}}[$j]";
            print WRITE ${$seq_ref};
            close WRITE;
            undef $seq_ref;

            #generate Genemodel file for Transterm
            foreach my $entry (@job_list) {
               #get boundary parameters
               my ($ID, $left_bd, $right_bd, $orientation);
               'reset' =~ m/reset/;
               $entry =~ m/$args{input_file}\___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
               $ID = $1;
               $left_bd = $2;
               $right_bd = $3;
               $orientation = $4;

               unless (defined $orientation && $orientation =~ /(sense|antisense)/) {
                  print "Error in Transterm analysis with file $list[$j]";
                  open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
                  print WRITE "\nError in Transterm analysis with file $list[$j]\n\n";
                  close WRITE;
                  #exit safely
                  CORE::exit();
               }

               if ($orientation eq "antisense") {
                  $gene_model .= $ID."\t".$right_bd."\t".$left_bd."\t".$list[$j]."\n";
               } else {
                  $gene_model .= $ID."\t".$left_bd."\t".$right_bd."\t".$list[$j]."\n";
               }
            }
            open WRITE, "+>${$args{ini_ref}}{transterm_dir}\/$list[$j]\.coords" or do {
               print "Error in Transterm analysis: cannot write dna coordinate for file $list[$j]";
               open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
               print WRITE "\nError in Transterm analysis: cannot write dna coordinate for file $list[$j]\n\n";
               close WRITE;
               #exit safely
               CORE::exit();
            };
            print WRITE $gene_model;
            close WRITE;

            #call transterm2
            my $coord = ${$args{ini_ref}}{transterm_dir}.'/'.$list[$j].'.coords';
            my $input = ${$args{ini_ref}}{transterm_dir}.'/'.$list[$j];
            my $data  = ${$args{ini_ref}}{transterm_dir}.'/expterm.dat';
            my $terminator = ${$args{ini_ref}}{transterm_dir}.'/transterm';
            $transterm = `$terminator -p $data $input $coord`;

            #write results file if there are any generated
            if ($transterm =~ /\w+/) {
               open WRITE, "+>${$args{ini_ref}}{transterm_results}/$list[$j]\_transterm" or do {
                  print "Error writing Transterm results for file $list[$j]";
                  open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
                  print WRITE "\nError writing Transterm results for file $list[$j]\n\n";
                  close WRITE;
                  #exit safely
                  CORE::exit();
               };
               print WRITE $transterm;
               close WRITE;
            }

            #delete temp files
            unlink ${$args{ini_ref}}{transterm_dir}.'/'.$list[$j].'.coords';
            unlink ${$args{ini_ref}}{transterm_dir}.'/'.$list[$j];

            #exit safely
            CORE::exit();
         } else {
            die "Couldn\’t fork\: $!\n";
         }

      }
      #wait
      foreach (@childs) {
         waitpid($_, 0);
      }
   }
   &hide_pbar_2;
   return (1);
}


1;