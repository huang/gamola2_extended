#!/opt/ActivePerl-5.8/bin/perl

#vector_screen to identify vector sequences in input, uses Blast program calls
#input arguments: main_window, progress_bar, auto_ini_ref, ini_ref,
#                 directory, input_array_ref, output_folder, input_type
#feature key: Vector_mtach

package ProgrammeModules::vector_screen;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&vector_screen);
use vars qw();

use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);
use Cwd;
use File::Copy;
use File::Find;

#local variables
my (%args, $newdir);

sub vector_screen {
   my %args = @_;
   my @templates = ();
   my @joblist   = ();
   my $max_count = @{$args{input_list_ref}};
   my $count     = 1;

   #spacer sequence
   my $spacer     = ${$args{ini_ref}}{spacer};
   my $inv_spacer = orient_sequence(seq         => ${$args{ini_ref}}{spacer},
                                    orientation => 'antisense');

   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'Screening for vector sequences',
                   label        => 'vector Sequences'
                  );
   &show_pbar_2;

   #break up GAMOLA linked contigs, this will increase hit rate for long concatenated sequences
   foreach my $entry (@{ $args{input_list_ref} }) {
      my $file_ref;
      my $local_position = 1;
      my @contig_order   = ();

      &update_pbar_2(title        => "Screening for vector sequences",
                     label        => "Breaking up concatenated sequences",
                     progress     => ($count / $max_count) * 100,
                    );
      $count++;

      $file_ref = slurp(main_window  => $args{main_window},
                        progress_bar => $args{progress_bar},
                        auto_ini_ref => $args{auto_ini_ref},
                        ini_ref      => $args{ini_ref},
                        directory    => $args{directory},
                        filename     => $entry
                       );
      #remove fasta header
      $$file_ref =~ s/^>[^\n]+?\n//;
      $$file_ref =~ s/\s//gs;

      #capture individual contigs is sequence
      @contig_order = split/$spacer/i, $$file_ref;
      #only one entry? try inverse primer
      if ($#contig_order < 1) {
         @contig_order = split/$inv_spacer/i, $$file_ref;
      }
      #create ordered array with all relevant data
      foreach my $contig (@contig_order) {
         my ($left_bd, $right_bd);
         $left_bd        = $local_position;
         $right_bd       = $left_bd + length($contig) - 1;
         #newlocal_position
         $local_position = $right_bd + length($spacer) + 1;
         push (@joblist, $entry.'___'.$left_bd.'___'.$right_bd.'___'.$contig);
      }
   }

   #generate list of jobs to do
   foreach my $entry (@joblist) {
      #test if re-use
      if (${$args{auto_ini_ref}}{reuse_results} == 1) {
         $entry =~ m/(.+?)___(\d+)___/;
         my ($name, $left_bd) = ($1, $2);
         next if (-e ${$args{ini_ref}}{vector_results}.'/'.$name.'___'.$left_bd.'.blastn' && -s ${$args{ini_ref}}{vector_results}.'/'.$name.'___'.$left_bd.'.blastn' > 0);
         #add if required
         push (@templates, $entry);
      }
      #add always if no re-use
      elsif (${$args{auto_ini_ref}}{reuse_results} == 0) {
         push (@templates, $entry);
      }
   }

   #max number of fasta entries
   $max_count = @templates;
   if ($max_count < 1) {$max_count = 1};

   #start multithreading BlastN round for vector screen
   for (my $i=0; $i<=$#templates+1; $i = $i+${$args{auto_ini_ref}}{CPU}) {
      my $count = $i + ${$args{auto_ini_ref}}{CPU};
      if ($count > $#templates+1) {$count = $#templates+1};

      &update_pbar_2(title        => "Screening for vector sequences",
                     label        => "Vector screen:  $count of $max_count",
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
            my (@list, $blastresult, $blast_hits_ref);

            #parse entry
            $templates[$j] =~ m/(.+?)___(\d+)___(\d+)___(.+)/;
            my ($name, $left_bd, $right_bd, $seq) = ($1, $2, $3, $4);

            #create nice fasta format
            $seq = '>'.$name.' Left_bd: '.$left_bd.' Right_bd: '.$right_bd."\n".$seq;

            #blastN against nt sequence
            if (open RESULT, "-|") { # original process
               local $/;
               $blastresult = <RESULT>;
            } else { # child
               if (open STDIN, "-|") { # child
                  #exec "./blastall -p blastn -d ${$args{ini_ref}}{rrna_db_path}/$db -n T -I T -g T -m 9 -W 56 -F F -b 1000000 -v 1000000";

                  if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
                     my ($dust);
                     chdir ${$args{ini_ref}}{blast_executables};
                     if (${$args{auto_ini_ref}}{vector_dust} == 1) {
                        $dust = 'T';
                     } else {
                        $dust = 'F';
                     }
                     exec "./blastall -p blastn -q ${$args{auto_ini_ref}}{vector_penalty} -r ${$args{auto_ini_ref}}{vector_reward} -G ${$args{auto_ini_ref}}{vector_gapopen} -E ${$args{auto_ini_ref}}{vector_gapextend} -F $dust -e ${$args{auto_ini_ref}}{vector_evalue} -Y ${$args{auto_ini_ref}}{vector_xdrop} -d ${$args{auto_ini_ref}}{full_vector_db}";
                  } else {
                     chdir ${$args{ini_ref}}{blast_plus_executables};
                     my ($dust, $soft_mask);
                     if (${$args{auto_ini_ref}}{vector_dust} == 1) {
                        $dust = 'yes';
                     } else {
                        $dust = 'no';
                     }
                     if (${$args{auto_ini_ref}}{vector_soft_masking} == 1) {
                        $soft_mask = 'true';
                     } else {
                        $soft_mask = 'false';
                     }
                     exec "./blastn -db ${$args{auto_ini_ref}}{full_vector_db} -penalty ${$args{auto_ini_ref}}{vector_penalty} -gapopen ${$args{auto_ini_ref}}{vector_gapopen} -gapextend ${$args{auto_ini_ref}}{vector_gapextend} -dust $dust -soft_masking $soft_mask -evalue ${$args{auto_ini_ref}}{vector_evalue} -xdrop_ungap ${$args{auto_ini_ref}}{vector_xdrop}";
                  }
                  die "Cannot exec: $!";
               } else { # grandchild
                  print $seq;
                  CORE::exit();
               }
            }
            chdir ${$args{auto_ini_ref}}{work_dir};

            # write Blast output file
            open WRITE, "+>${$args{ini_ref}}{vector_results}/$name\___$left_bd\.blastn" or do {
               open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
               print WRITE "Error in ribosomal RNA analysis".
                           "\nError writing BlastN result file $name\___$left_bd\.blastn to directory ${$args{ini_ref}}{vector_results}\n\n";
               close WRITE;
               print "\nVector screening analysis: Error writing Blast result file $name\___$left_bd\.blastn to directory ${$args{ini_ref}}{vector_results}\n";
               CORE::exit();
            };
            print WRITE $blastresult;
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

   #parse blast results
   &update_pbar_2(title        => "Vector screen analysis",
                  label        => "Parsing all results",
                  progress     => 1,
                 );
   &hide_pbar_2;
   return (1);
}


1;