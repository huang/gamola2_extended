#!/opt/ActivePerl-5.8/bin/perl -w
# Script to run BLASTN on a FASTA files containing split contigs,
# creating the sequence fragments to be BLASTed if needed

#input: main window, filename, blast_options

package GeneModel::blast_contigs;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&blastcontigs);
use vars qw();

#use Basics::MesgBox;
use Basics::progress_bar qw(:DEFAULT);
use initialise::read_me  qw(:DEFAULT);
use Critica;

my ($blastn, $x, $parameters, $db, @flags, $dirname, $f, @output, @files);

sub blastcontigs {
   my %args = @_;
   ${$args{progress_bar}}->Label(-text=>"Generating Critica Blast fragments") ->grid (-row => 0, -column => 0, -sticky => 'ew');
   ${$args{main_window}}->update;
   $dirname=getPrefixName($args{filename})."dir";
   if (!-e("$dirname")) {
      &createBlastFragments($args{directory}.'/'.$dirname, $args{directory}.'/'.$args{filename}, 3000, 100);
   }
   opendir(DIR,$args{directory}.'/'.$dirname) || die("Can't find $args{directory}".'/'."$dirname! Aborting\n");
   @files=readdir(DIR);
   closedir(DIR);

   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Critica analysis",
                   label        => 'Critica Blast analysis'
                  );
   &show_pbar_2;

   #start multithreading Blast
   for (my $i=1; $i<=$#files+2; $i += ${$args{auto_ini_ref}}{CPU}) {
      my $count = $i + ${$args{auto_ini_ref}}{CPU};
      if ($count > $#files+2) {$count = $#files+2};

      ${$args{progress_bar}}->Label(-text=>"Generating Critica Blast fragments for files $i to $count out of ".($#files +1)) ->grid (-row => 0, -column => 0, -sticky => 'ew');
      ${$args{main_window}}->update;
      &update_pbar_2(title    => 'Critica analysis',
                     label    => "Analysing Critica fragment files $i to $count out of ".($#files +1),
                     progress => (($count / $#files) * 100)
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
            my $output = '';
            if ($files[$j-1] =~ /^[\.]+$/) {CORE::exit()};
            $f=$args{directory}.'/'.$dirname.'/'.$files[$j-1];
            if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
               chdir ${$args{ini_ref}}{blast_executables};
               $output = `${$args{ini_ref}}{blast_executables}/blastall -p blastn -g F -e 1e-4 -d ${$args{auto_ini_ref}}{full_critica_db} -i $f`; #child
               chdir ${$args{auto_ini_ref}}{work_dir};
            } elsif (${$args{auto_ini_ref}}{blast_plus} == 1) {
               chdir ${$args{ini_ref}}{blast_plus_executables};
               $output = `${$args{ini_ref}}{blast_plus_executables}/blastn -ungapped -evalue 1e-4 -db ${$args{auto_ini_ref}}{full_critica_db} -query $f`;
               chdir ${$args{auto_ini_ref}}{work_dir};
            }
            #open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            #print ERRORLOG "Error generating Critica Model for file $args{filename}".
            #         "\nError while opening directory '.$args{directory}.'/'.$dirname: ${$args{ini_ref}}{blast_plus_executables}/blastn -ungapped -evalue 1e-4 -db ${$args{auto_ini_ref}}{full_critica_db} -query $f!\n\n";
            #close ERRORLOG;
            
            open RESULT, "+>$f".'_blastn' or die "\nCouldn't create file $f\n";
            print RESULT $output;
            close RESULT;

            #exit safely
            CORE::exit();
         } else {
            die "Couldn\'t fork\: $!\n";
            Tk::exit;
         }
      }
      #wait
      foreach (@childs) {
         waitpid($_, 0);
      }
   }
   undef @files;

   #read all results and combine to original output
   opendir BLAST, $args{directory}.'/'.$dirname or do {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Error while opening directory '.$args{directory}.'/'.$dirname",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error generating Critica Model for file $args{filename}".
                     "\nError while opening directory '.$args{directory}.'/'.$dirname\n\n";
      close ERRORLOG;
      return (0); #return to main
   };
   my @results = grep /_blastn$/, readdir(BLAST);
   my $file;

   #open tempfile for blast-results
   open WRITE, '+>'.$args{directory}.'/'.$args{filename}.'.blast';
   #open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
   #print ERRORLOG "Error open file\n\n";
   #close ERRORLOG;
   
   foreach my $entry (@results) {
      {
         local ( $/, *ENTRY ) ;
         open ENTRY, "<".$args{directory}.'/'.$dirname.'/'.$entry;
         $file = <ENTRY>;
         close ENTRY;
      }

      print WRITE $file;

      #open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      #print ERRORLOG "Error before s/_blastn $entry\n\n";
      #close ERRORLOG;
      unlink $args{directory}.'/'.$dirname.'/'.$entry;
      #$entry =~ s/_blastn$//;
      #open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      #print ERRORLOG "Error after s/_blastn $entry\n\n";
      #close ERRORLOG;
      
      #unlink $args{directory}.'/'.$dirname.'/'.$entry;
   }

   #rmdir($args{directory}.'/'.$dirname);
   close WRITE;
   undef @results;
   &hide_pbar_2;
   return (1);
}
1;