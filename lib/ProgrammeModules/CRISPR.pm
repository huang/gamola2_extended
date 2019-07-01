#!/opt/ActivePerl-5.8/bin/perl

#CRT program calls and reformatting of results
#input arguments: main_window, progress_bar, auto_ini_ref, ini_ref,
#                 directory, input_array_ref, output_folder, input_type

package ProgrammeModules::CRISPR;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&CRISPR);
use vars qw();

use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);
use Cwd;
use File::Copy;

#local variables
my (%args, $newdir);

sub CRISPR {
   my %args = @_;
   my (@list);

   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "CRISPR analysis",
                   label        => 'CRISPR analysis'
                  );
   &show_pbar_2;

   #reuse existing entries if re-use results
   if (${$args{auto_ini_ref}}{reuse_results} == 1) {
      my (@CRISPR_results, %exist_CRISPR_file);
      &update_pbar_2(title        => 'CRISPR analysis',
                     label        => "Testing for existing results",
                     progress     => 1,
                    );

      #read existing results
      opendir SEQINPUT, ${$args{ini_ref}}{CRISPR_results};
      @CRISPR_results = grep /\.CRISPR$/, readdir(SEQINPUT);
      closedir SEQINPUT;
      #create hash
      foreach (@CRISPR_results) {
         $exist_CRISPR_file{$_} = '1';
      }
      #iterate over input files
      foreach my $entry (@{ $args{input_array_ref} }) {
         #skip if result file exists and re-use is activated
         next if (exists $exist_CRISPR_file{$entry.'.CRISPR'} && -s ${$args{ini_ref}}{CRISPR_results}.'/'.$entry.'.CRISPR' > 0);
         #add if required
         push (@list, $entry);
      }
      undef @CRISPR_results;
      undef %exist_CRISPR_file;
   } else {
      @list = @{ $args{input_array_ref} };
   }

   #max number of fasta entries
   my $max_count = $#list + 1;
   if ($max_count < 1) {$max_count = 1};

   #start multithreading CRISPR
   for (my $i=0; $i<=$#list + 1; $i = $i+${$args{auto_ini_ref}}{CPU}) {
      my $count = $i + ${$args{auto_ini_ref}}{CPU};
      if ($count > $#list + 1) {$count = $#list + 1};

      &update_pbar_2(title        => "CRISPR analysis",
                     label        => "CRISPR analysis $count of $max_count",
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
            my (@CRISPRs, $CRISPR);

            #define input and putput files
            my $exec   = ${$args{ini_ref}}{CRISPR_dir}.'/'.'CRT_java.jar';
            my $input  = ${$args{ini_ref}}{input_files}.'/'.$list[$j];
            my $output = ${$args{ini_ref}}{CRISPR_results}.'/'.$list[$j].'.CRISPR';

            #call CRT
            `java -cp $exec crt -minNR ${$args{auto_ini_ref}}{CRISPR_minNR} -minRL ${$args{auto_ini_ref}}{CRISPR_minRL} -maxRL ${$args{auto_ini_ref}}{CRISPR_maxRL} -minSL ${$args{auto_ini_ref}}{CRISPR_minSL} -maxSL ${$args{auto_ini_ref}}{CRISPR_maxSL} $input $output` ;

            #test if output file was written:
            unless (-e $output && -s $output > 0) {
               print "Error in CRISPR analysis: cannot write output for file $list[$j]";
               open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
               print WRITE "\nError in CRISPR analysis: cannot write output for file $list[$j]\n\n";
               close WRITE;
               #exit safely
               CORE::exit();
            };

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
   &hide_pbar_2;
   return (1);
}


1;