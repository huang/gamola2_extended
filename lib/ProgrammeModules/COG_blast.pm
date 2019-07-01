#!/opt/ActivePerl-5.8/bin/perl

#COG program calls
#input arguments: progress_bar, nt_seq, left_bd, right_bd, orientation, ini_ref, auto_ini_ref


package ProgrammeModules::COG_blast;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&COG_file);
use vars qw();

use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);
use Cwd;

#local variables
my (%args, $ret_seq, $psi_aa_db, $psi_nt_db, $blastresult, $COG_db_name);

sub COG_blast_seq {
   my %args = (aa_input_file => 'temp.aa', @_);
   my $blastresult = "";

   #reformat sequence into fasta format
   my ($seq_ref_nt) = &seq2fasta (main_window  => \$args{main_window},
                                  progress_bar => \$args{progress_bar},
                                  auto_ini_ref => $args{auto_ini_ref},
                                  ini_ref      => $args{ini_ref},
                                  sequence     => $args{sequence_nt},
                                  orientation  => $args{orientation}
                                 );
   my ($seq_ref_aa) = &seq2fasta (main_window  => \$args{main_window},
                                  progress_bar => \$args{progress_bar},
                                  auto_ini_ref => $args{auto_ini_ref},
                                  ini_ref      => $args{ini_ref},
                                  sequence     => $args{sequence_aa},
                                  orientation  => $args{orientation}
                                 );

   #get current working directory
   my $curdir = "";
   $curdir = getcwd();

   #change to Blast directory
   if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
      chdir ${$args{ini_ref}}{blast_executables};
   } else {
      chdir ${$args{ini_ref}}{blast_plus_executables};
   }

   if (open RESULT, "-|") { # original process
      local $/;
      $blastresult = <RESULT>;
   } else { # child
      if (open STDIN, "-|") { # child
         if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
            exec "./blastall -a ${$args{auto_ini_ref}}{CPU} -p blastp -d ${$args{auto_ini_ref}}{full_COG_db}"; #child
         } else {
            exec "./blastp -num_threads ${$args{auto_ini_ref}}{CPU} -db ${$args{auto_ini_ref}}{full_COG_db}"; #child
         }
         die "Cannot exec: $!";
      } else { # grandchild
         print ${$seq_ref_aa};
         CORE::exit;
      }
   }

   #change back to previous working dir
   chdir $curdir;

   #define COG database
   if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
      $COG_db_name = 'Original 2003 COG database';
   } elsif (${$args{auto_ini_ref}}{COG_type} eq 'arCOG') {
      $COG_db_name = 'Archaeal COG database';
   } elsif (${$args{auto_ini_ref}}{COG_type} eq 'arCOG2014') {
      $COG_db_name = 'Updated archaeal 2014 COG database';
   } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2008') {
      $COG_db_name = 'Updated 2008 COG database';
   } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2014') {
      $COG_db_name = 'Updated 2014 COG database';
   } elsif (${$args{auto_ini_ref}}{COG_type} eq 'POG2013') {
      $COG_db_name = 'Original POG (Phage COG) database';
   } else {
      ${$args{progress_bar}}->configure(-label=>"Error defining database label for database ${$args{auto_ini_ref}}{full_COG_db}");
      ${$args{main_window}}->update;
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in COG Blast analysis".
                  "\nError defining database label for database ${$args{auto_ini_ref}}{full_COG_db}\n\n";
      close WRITE;
   }

   #build COG header
   my $COG_header = "COG analysis: $COG_db_name\n".
                    "Name of the input sequence:\t$args{input_file}\n\n".
                    "DNA sequence:\t$args{sequence_nt}\n".
                    "Deduced aminoacid sequence:\t$args{sequence_aa}\n\n".
                    "Length in aminoacids:\t".length($args{sequence_aa})."\n".
                    "Gene model summary:\tORF-designation $args{orf_id}\tLeft boundary $args{left_bd}\tRight boundary $args{right_bd}\tOrientation $args{orientation}\n";

   #valid result?
   unless ($blastresult =~ /\*\*\s+No\s+hits\s+found\s+\*\*/igs ||
           $blastresult =~ /Identities/                         ||
           $blastresult =~ /\n\n\s*No COG Hits found/igs           ||
           $blastresult =~ /\n\s*\*\*\*\*\* No hits found/igs)             {
      ${$args{progress_bar}}->configure(-label=>"Error blasting for COG results for $args{input_file} using database ${$args{auto_ini_ref}}{full_COG_db}");
      ${$args{main_window}}->update;
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in COG Blast analysis".
                  "\nError blasting for COG results  for $args{input_file} using database ${$args{auto_ini_ref}}{full_COG_db}\n$COG_header\n\n";
      close WRITE;
      #print "No COG results found for $args{input_file} using database ${$args{auto_ini_ref}}{full_COG_db} and Sequence: ${$seq_ref_nt} \n\n..$blastresult ...";<STDIN>;
      return (0);
   };

   #reformat blast results
   ($blastresult) = &reformat_COG_blast_res(main_window          => $args{main_window},
                                            progress_bar         => $args{progress_bar},
                                            auto_ini_ref         => $args{auto_ini_ref},
                                            ini_ref              => $args{ini_ref},
                                            blast_result         => $blastresult,
                                            input_file           => $args{input_file},
                                            orf_id               => $args{orf_id},
                                            COGcode_to_number    => $args{COGcode_to_number},
                                            COGcode_to_header    => $args{COGcode_to_header},
                                            COGcode_to_letter    => $args{COGcode_to_letter},
                                            COGcode_to_phylo     => $args{COGcode_to_phylo},
                                            COGletter_to_family  => $args{COGletter_to_family},
                                            COG2008genomes       => $args{COG2008genomes},
                                            COG2008_to_def       => $args{COG2008_to_def},
                                            COG2008_to_acc       => $args{COG2008_to_acc},
                                            arCOGcode_to_COG     => $args{arCOGcode_to_COG},
                                            arCOGacc_to_org      => $args{arCOGacc_to_org},
                                            COG2014genomes       => $args{COG2014genomes},
                                            COG2014_to_def       => $args{COG2014_to_def},
                                            COG2014_to_acc       => $args{COG2014_to_acc},
                                            COG2014_to_refseq    => $args{COG2014_to_refseq},
                                            arCOGcode2014_to_COG => $args{arCOGcode2014_to_COG},
                                            arCOGacc2014_to_org  => $args{arCOGacc2014_to_org},
                                            POG2013_to_gi        => $args{POG2013_to_gi},
                                            COGletter_to_family  => $args{COGletter_to_family},
                                           );

   return ($COG_header.$blastresult);
}

sub COG_blast_threaded {
   my %args = (aa_input_file => 'temp.aa', @_);
   my $blastresult = "";

   #reformat sequence into fasta format
   my ($seq_ref_nt) = &seq2fasta (sequence     => $args{sequence_nt},
                                  orientation  => $args{orientation}
                                 );
   my ($seq_ref_aa) = &seq2fasta (sequence     => $args{sequence_aa},
                                  orientation  => $args{orientation}
                                 );

   #get current working directory
   my $curdir = "";
   $curdir = getcwd();

   #change to Blast directory
   if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
      chdir ${$args{ini_ref}}{blast_executables};
   } else {
      chdir ${$args{ini_ref}}{blast_plus_executables};
   }

   if (open RESULT, "-|") { # original process
      local $/;
      $blastresult = <RESULT>;
   } else { # child
      if (open STDIN, "-|") { # child
         if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
            exec "./blastall -p blastp -d ${$args{auto_ini_ref}}{full_COG_db}"; #child
         } else {
            exec "./blastp -db ${$args{auto_ini_ref}}{full_COG_db}"; #child
         }
         die "Cannot exec: $!";
      } else { # grandchild
         print ${$seq_ref_aa};
         CORE::exit;
      }
   }

   #change back to previous working dir
   chdir $curdir;

   #define COG database
   if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
      $COG_db_name = 'Original 2003 COG database';
   } elsif (${$args{auto_ini_ref}}{COG_type} eq 'arCOG') {
      $COG_db_name = 'Archaeal COG database';
   } elsif (${$args{auto_ini_ref}}{COG_type} eq 'arCOG2014') {
      $COG_db_name = 'Updated archaeal 2014 COG database';
   } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2008') {
      $COG_db_name = 'Updated 2008 COG database';
   } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2014') {
      $COG_db_name = 'Updated 2014 COG database';
   } elsif (${$args{auto_ini_ref}}{COG_type} eq 'POG2013') {
      $COG_db_name = 'Original POG (Phage COG) database';
   } else {
      ${$args{progress_bar}}->configure(-label=>"Error defining database label for database ${$args{auto_ini_ref}}{full_COG_db}");
      ${$args{main_window}}->update;
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in COG Blast analysis".
                  "\nError defining database label for database ${$args{auto_ini_ref}}{full_COG_db}\n\n";
      close WRITE;
   }

   #build COG header
   my $COG_header = "COG analysis: $COG_db_name\n".
                    "Name of the input sequence:\t$args{input_file}\n\n".
                    "DNA sequence:\t$args{sequence_nt}\n".
                    "Deduced aminoacid sequence:\t$args{sequence_aa}\n\n".
                    "Length in aminoacids:\t".length($args{sequence_aa})."\n".
                    "Gene model summary:\tORF-designation $args{orf_id}\tLeft boundary $args{left_bd}\tRight boundary $args{right_bd}\tOrientation $args{orientation}\n";

   #valid result?
   unless ($blastresult =~ /\n\s*\*\*\*\*\* No hits found \*\*\*\*\*/igs ||
           $blastresult =~ /Identities/igs                              ||
           $blastresult =~ /\n\nNo COG Hits found/igs) {
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in COG Blast analysis".
                  "\nError in COG result for $args{input_file} using database ${$args{auto_ini_ref}}{full_COG_db}\n$COG_header\n$blastresult\n\n";
      close WRITE;
      #print "No COG results found for $args{input_file} using database ${$args{auto_ini_ref}}{full_COG_db} and Sequence: ${$seq_ref_nt}";
      return (0);
   };

   #reformat blast results
   ($blastresult) = &reformat_COG_threaded(main_window          => $args{main_window},
                                           progress_bar         => $args{progress_bar},
                                           auto_ini_ref         => $args{auto_ini_ref},
                                           ini_ref              => $args{ini_ref},
                                           blast_result         => $blastresult,
                                           input_file           => $args{input_file},
                                           orf_id               => $args{orf_id},
                                           COGcode_to_number    => $args{COGcode_to_number},
                                           COGcode_to_header    => $args{COGcode_to_header},
                                           COGcode_to_letter    => $args{COGcode_to_letter},
                                           COGcode_to_phylo     => $args{COGcode_to_phylo},
                                           COGletter_to_family  => $args{COGletter_to_family},
                                           COG2008genomes       => $args{COG2008genomes},
                                           COG2008_to_def       => $args{COG2008_to_def},
                                           COG2008_to_acc       => $args{COG2008_to_acc},
                                           arCOGcode_to_COG     => $args{arCOGcode_to_COG},
                                           arCOGacc_to_org      => $args{arCOGacc_to_org},
                                           COG2014genomes       => $args{COG2014genomes},
                                           COG2014_to_def       => $args{COG2014_to_def},
                                           COG2014_to_acc       => $args{COG2014_to_acc},
                                           COG2014_to_refseq    => $args{COG2014_to_refseq},
                                           arCOGcode2014_to_COG => $args{arCOGcode2014_to_COG},
                                           arCOGacc2014_to_org  => $args{arCOGacc2014_to_org},
                                           POG2013_to_gi        => $args{POG2013_to_gi},
                                           COGletter_to_family  => $args{COGletter_to_family},
                                         );

   return ($COG_header.$blastresult);
}

sub reformat_COG_blast_res {
   my %args                = @_;
   my $line1               = '';
   my $orig_blast_overview = '';
   my $alignments          = '';
   my @acc                 = ();
   my $line_length         = 0;
   my $acc_length          = 0;
   my $genome_length       = 0;
   my $annotation_length   = 0;
   my $class_length        = 0;
   my $subj_length         = 0;
   my $new_overview        = '';

   #catch original blast overview
   $args{blast_result} =~ m/Sequences producing significant alignments:                      \(bits\) Value\s*\n(.*?)\n\n\>/s;
   $orig_blast_overview = $1;
   unless (defined $orig_blast_overview) {
      $args{blast_result} =~ m/Sequences producing significant alignments:\s+\(Bits\)\s+Value\s*\n(.*?)\n\n\>/is;
      $orig_blast_overview = $1;
   }

   unless (defined ($orig_blast_overview)) {
      ${$args{progress_bar}}->configure(-label=>"Could not grab Blast overview from current aa sequence".
                                                "\nNo significant Blast results were generated. Skipping");
      ${$args{main_window}}->update;
      return ("\n\nNo COG Hits found\n");
   }

   if ($orig_blast_overview !~ /\w+/) {
      ${$args{progress_bar}}->configure(-label=>"No significant Blast results were generated. Skipping");
      ${$args{main_window}}->update;
      return ("\n\nNo COG Hits found\n");
   }

   $orig_blast_overview = "\n".$orig_blast_overview."\n";

   #catch aligments
   $args{blast_result} =~ m/\n\n\>(.*?)Lambda\s+K\s+H/s;
   $alignments = '>'.$1;
   unless ($alignments =~ /\w{2,}/) {
      ${$args{progress_bar}}->configure(-label=>"Could not grab Blast alignments from current aa sequence");
      ${$args{main_window}}->update;
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in COG Blast analysis".
                  "\nCould not grab Blast alignments from current aa sequence for $args{input_file}\n\n";
      close WRITE;
   }

   #build array for accession numbers from blast overview
   if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
      @acc = ($orig_blast_overview =~ m/\n(\S*?) +\d+ +[^\n]+/sg);
      unless (scalar @acc >=0) {
         @acc = ($orig_blast_overview =~ m/lcl\|(\S*?)\s+/sg);
      }
      unless ((scalar @acc) >= 1) {
         @acc = ($orig_blast_overview =~ m/\n(\S+)\s+/sg);
      }
   } elsif (${$args{auto_ini_ref}}{COG_type} =~ m/^(arCOG|arCOG2014|COG2008|COG2014)$/) {
      @acc = ($orig_blast_overview =~ m/\n(\w+\|\S+)\|? +/sg);
      unless ((scalar @acc) >= 1) {
         @acc = ($orig_blast_overview =~ m/\n(\S+)\s+/sg);
      }
   } elsif (${$args{auto_ini_ref}}{COG_type} =~ m/^POG2013$/) {
      #@acc = ($orig_blast_overview =~ m/\n.*?ref\|([^\|]+)\|/sg); #needs 'ref' acc in all cases
      #unless ((scalar @acc) >= 1) {
         @acc = ($orig_blast_overview =~ m/\n(\S+)\s+/sg);
      #}
   }
   @acc = splice(@acc,0,${$args{auto_ini_ref}}{blast_entry_number}); #restrict number of entries to set value in $key
   #built new Blast overview
   #find out maximum text length

   #first, max acc name length
   foreach (@acc) {
      my $acc_temp    = '';
      my $temp        = '';
      $acc_temp       = "\>$_";

      my $nometa  = quotemeta ($_);
      'reset'     =~ m/reset/;
      $alignments =~ m/\>($nometa.*?)\n\s+Score \= /s;
      $temp = $1;
      unless (defined $temp) {
         $alignments =~ m/\>\w+\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }
      #POG?
      unless (defined $temp) {
         $alignments =~ m/\>lcl.+? ref\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }

      #determine maximum length of entry line
      $temp =~ s/\n?\s+Length\s*=\s*.+$//;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $temp =~ s/\s+$//g;

      #determine max length of acc
      if (length($acc_temp) > $acc_length) {$acc_length  = length($acc_temp)};
   }

   #next, max genome name length
   foreach (@acc) {
      my $acc_temp    = '';
      my $temp        = '';
      $acc_temp       = "\>$_";

      my $nometa  = quotemeta ($_);
      'reset'     =~ m/reset/;
      $alignments =~ m/\>($nometa.*?)\n\s+Score \= /s;
      $temp = $1;
      unless (defined $temp) {
         $alignments =~ m/\>\w+\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }
      #POG?
      unless (defined $temp) {
         $alignments =~ m/\>lcl.+? ref\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }

      #determine maximum length of entry line
      $temp =~ s/\n?\s+Length\s*=\s*.+$//;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $temp =~ s/\s+$//g;

      #determine genome & annotation name length for COG 2003
      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
         if (length(${$args{COGcode_to_phylo}}->{$_}) > $genome_length) {
            $genome_length = length(${$args{COGcode_to_phylo}}->{$_});
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2008') {
         $_ =~ s/^\w+\|//;
         if (length(${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$_}->{'genome_abbr'}}) > $genome_length) {
            $genome_length = length(${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$_}->{'genome_abbr'}});
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2014') {
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         if (length(${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'genome_abbr'}) > $genome_length) {
            $genome_length = length(${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'genome_abbr'});
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'arCOG2014') { #arCOG doesn't have genome name info
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         if (length(${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'genome'}) > $genome_length) {
            $genome_length = length(${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'genome'});
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'POG2013') {
         if (length(${$args{POG2013_to_gi}}->{$_}->{'phage'}) > $genome_length) {
            $genome_length = length(${$args{POG2013_to_gi}}->{$_}->{'phage'});
         }
      }
   }

   #next, max class name length
   foreach (@acc) {
      my $acc_temp    = '';
      my $temp        = '';
      $acc_temp       = "\>$_";

      my $nometa  = quotemeta ($_);
      'reset'     =~ m/reset/;
      $alignments =~ m/\>($nometa.*?)\n\s+Score \= /s;
      $temp = $1;
      unless (defined $temp) {
         $alignments =~ m/\>\w+\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }
      #POG?
      unless (defined $temp) {
         $alignments =~ m/\>lcl.+? ref\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }

      #determine maximum length of entry line
      $temp =~ s/\n?\s+Length\s*=\s*.+$//;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $temp =~ s/\s+$//g;

      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2008') {
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]//g;

         if (length(${$args{COG2008_to_acc}}->{$COG_protID}->{'COG2008class'}) > $class_length) {
            $class_length = length(${$args{COG2008_to_acc}}->{$COG_protID}->{'COG2008class'});
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2014') {
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         if (length(${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}.', '.${$args{COGletter_to_family}}->{${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}}) > $class_length) {
            $class_length = length(${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}.', '.${$args{COGletter_to_family}}->{${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}});
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'arCOG') {
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         if (length(${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'COGclass'}) > $class_length) {
            $class_length = length(${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'COGclass'});
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'arCOG2014') {
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         if (length(${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGclass'}.', '.${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'}) > $class_length) {
            $class_length = length(${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGclass'}.', '.${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'});
         }
      }
   }

   #next, max annotation length
   foreach (@acc) {
      my $acc_temp    = '';
      my $temp        = '';
      $acc_temp       = "\>$_";

      my $nometa  = quotemeta ($_);
      'reset'     =~ m/reset/;
      $alignments =~ m/\>($nometa.*?)\n\s+Score \= /s;
      $temp = $1;
      unless (defined $temp) {
         $alignments =~ m/\>\w+\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }
      #POG?
      unless (defined $temp) {
         $alignments =~ m/\>lcl.+? ref\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }

      #determine maximum length of entry line
      $temp =~ s/\n?\s+Length\s*=\s*.+$//;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $temp =~ s/\s+$//g;

      #determine genome & annotation name length for COG 2003
      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
         my $additional_COG_info = 'COG2003: '.      ${$args{COGcode_to_number}}->{$_}.
                                  '; Letter: '.      ${$args{COGcode_to_letter}}->{$_}.
                                  '; Genome: '.      ${$args{COGcode_to_phylo}}->{$_}.
                                  '; Annotation: '.  ${$args{COGcode_to_header}}->{$_}.
                                  ' '
                                 ;
         if (length($additional_COG_info) > $annotation_length) {
            $annotation_length = length($additional_COG_info);
         }
      }

      #determine genome & annotation name length for COG 2008
      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2008') {
         $_ =~ s/^\w+\|//;
         my $additional_COG_info = '';
         #COG2008 defined?
         if (${$args{COG2008_to_acc}}->{$_}->{'COG2008code'} =~ m/\w+/) {
            $additional_COG_info = ' COG2008: '.   ${$args{COG2008_to_acc}}->{$_}->{'COG2008code'}.
                                   '; Class: '.     ${$args{COG2008_to_acc}}->{$_}->{'COG2008class'}.
                                   '; Genome: '.    ${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$_}->{'genome_abbr'}}." " x ($genome_length - length(${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$_}->{'genome_abbr'}})).
                                   '; Annotation: '.${$args{COG2008_to_def}}->{${$args{COG2008_to_acc}}->{$_}->{'COG2008code'}}->{'annotation'}.
                                   ' '
                                   ;
            if (${$args{COG2008_to_acc}}->{$_}->{'COG2008code'} eq 'NULL') {
               $additional_COG_info = ' COG2008: not def'.
                                      '; Class: not def'.
                                      '; Genome: '.    ${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$_}->{'genome_abbr'}}." " x ($genome_length - length(${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$_}->{'genome_abbr'}})).
                                      '; Annotation: '.${$args{COG2008_to_def}}->{${$args{COG2008_to_acc}}->{$_}->{'COG2008code'}}->{'annotation'}.
                                      ' '
                                      ;
            }
         } else {
            $additional_COG_info = ' COG2008: Not def; Class: X; Genome: Not defined'." " x ($genome_length - 11).'; Annotation: Not available ';
         }
         if (length($additional_COG_info) > $annotation_length) {
            $annotation_length = length($additional_COG_info);
         }
      }

      #determine annotation length for COG2014
      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2014') {
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;

         my $additional_COG_info = ' COG2014: '.    ${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}.
                                   '; Class: '.     ${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}.', '.${$args{COGletter_to_family}}->{${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}}. ' ' x ($class_length - length(${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}) - 2 - length(${$args{COGletter_to_family}}->{${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}})).
                                   '; Genome: '.    ${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'genome_abbr'}.' ' x ($genome_length - length(${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'genome_abbr'})).
                                   '; Annotation: '.${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'annotation'}.
                                   ' '
                                   ;
         if (length($additional_COG_info) > $annotation_length) {
            $annotation_length = length($additional_COG_info);
         }
      }

      #determine annotation length for arCOG
      if (${$args{auto_ini_ref}}{COG_type} eq 'arCOG') {
         my $additional_COG_info = '';
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//;

         #arCOG defined?
         if (${$args{arCOGacc_to_org}} ->{$COG_protID}->{'arCOG'} =~ m/\w+/) {
            $additional_COG_info = ' arCOG: '  .    ${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}.
                                   '; Class: '  .   ${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'COGclass'}.
                                   '; Cluster: '.   ${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'COGcluster'}.
                                   '; Annotation: '.${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'arCOGannotation'}.
                                   ' '
                                   ;
         } else {
            $additional_COG_info = ' arCOG: Not def   ; Class: X; Cluster: Not def'.' Annotation: Not available ';
         }
         if (length($additional_COG_info) > $annotation_length) {
            $annotation_length = length($additional_COG_info);
         }
      }

      #determine info length for arCOG2014
      if (${$args{auto_ini_ref}}{COG_type} eq 'arCOG2014') {
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         my $additional_COG_info = ' arCOG2014: '.  ${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}.
                                   '; Class: '.     ${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGclass'}.', '.
                                                    ${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'}.
                                                    ' ' x ($class_length - length(${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGclass'}) - 2 - length(${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'})).
                                   '; Genome: '.    ${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'genome'}.' ' x ($genome_length - length(${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'genome'})).
                                   '; Annotation: '.${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'}.
                                   ' '
                                   ;
         if (length($additional_COG_info) > $annotation_length) {
            $annotation_length = length($additional_COG_info);
         }
      }

      #determine info length for POG2013
      if (${$args{auto_ini_ref}}{COG_type} eq 'POG2013') {
         my $additional_COG_info = '';
         #POG2013 defined?
         if (${$args{POG2013_to_gi}}->{$_}->{'class'} =~ m/\w+/) {
            $additional_COG_info = ' POG2013: '.    ${$args{POG2013_to_gi}}->{$_}->{'class'}.
                                   '; Phage: '.     ${$args{POG2013_to_gi}}->{$_}->{'phage'}.' ' x ($genome_length - length(${$args{POG2013_to_gi}}->{$_}->{'phage'})).
                                   '; Annotation: '.${$args{POG2013_to_gi}}->{$_}->{'annotation'}.
                                   ' '
                                   ;
         } else {
            $additional_COG_info = ' POG2013: Not def; Phage: Not defined'." " x ($genome_length - 11).'; Annotation: Not available ';
         }

         if (length($additional_COG_info) > $annotation_length) {
            $annotation_length = length($additional_COG_info);
         }
      }
   }

   #generate new overview
   $orig_blast_overview .= "\n";
   foreach (@acc) {
      my $acc_temp            = '';
      my $temp                = '';
      my $subj_length         = '';
      my $COG_acc             = '';
      my $additional_COG_info = '';

      $COG_acc = $_;
      $COG_acc =~ s/^\S+?\|//;
      $COG_acc =~ s/[\|\s]//g;

      #define proper acc code for COG dbs
      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
         $additional_COG_info = 'COG2003: '.      ${$args{COGcode_to_number}}->{$COG_acc}.
                               '; Letter: '.      ${$args{COGcode_to_letter}}->{$COG_acc}.
                               '; Genome: '.      ${$args{COGcode_to_phylo}}->{$COG_acc}.
                               '; Annotation: '.  ${$args{COGcode_to_header}}->{$COG_acc};
      }

      if (${$args{auto_ini_ref}}{COG_type} eq 'arCOG') {
         #arCOG defined?
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         if (${$args{arCOGacc_to_org}} ->{$COG_protID}->{'arCOG'} =~ m/\w+/) {
            $additional_COG_info = ' arCOG: '  .${$args{arCOGacc_to_org}} ->{$COG_protID}->{'arCOG'}.
                                   '; Class: '  .${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'COGclass'}.
                                   '; Cluster: '.${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'COGcluster'}.
                                   '; Annotation: '.${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'arCOGannotation'}.
                                   ' '
                                  ;
         } else {
            $additional_COG_info = ' arCOG: Not def   ; Class: X; Cluster: Not def ; Annotation: Not available';
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'arCOG2014') {
         #arCOG2014 defined?
         my $COG_protID = $_; #ACC is different here
         $COG_protID    =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;

         if (${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'} =~ m/\w+/) {
            $additional_COG_info = ' arCOG2014: '.  ${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}.
                                   '; Class: '.     ${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGclass'}.', '.
                                                    ${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'}.
                                                    ' ' x ($class_length - length(${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGclass'}) - 2 - length(${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'})).
                                   '; Genome: '.    ${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'genome'}.' ' x ($genome_length - length(${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'genome'})).
                                   '; Annotation: '.${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'}.
                                   ' '
                                   ;
         } else {
            $additional_COG_info = ' arCOG2014: Not def   ; Class: X; Genome: Not defined'." " x ($genome_length - 11).'; Annotation: Not available';
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2008') {
         #COG2008 defined?
         if (${$args{COG2008_to_acc}}->{$COG_acc}->{'COG2008code'} =~ m/\w+/) {
            $additional_COG_info = ' COG2008: '.   ${$args{COG2008_to_acc}}->{$COG_acc}->{'COG2008code'}.
                                   '; Class: '.     ${$args{COG2008_to_acc}}->{$COG_acc}->{'COG2008class'}.
                                   '; Genome: '.    ${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$COG_acc}->{'genome_abbr'}}." " x ($genome_length - length(${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$COG_acc}->{'genome_abbr'}})).
                                   '; Annotation: '.${$args{COG2008_to_def}}->{${$args{COG2008_to_acc}}->{$COG_acc}->{'COG2008code'}}->{'annotation'}.
                                   ' '
                                   ;
            if (${$args{COG2008_to_acc}}->{$COG_acc}->{'COG2008code'} eq 'NULL') {
               $additional_COG_info = ' COG2008: not def'.
                                      '; Class: not def'.
                                      '; Genome: '.    ${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$COG_acc}->{'genome_abbr'}}." " x ($genome_length - length(${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$COG_acc}->{'genome_abbr'}})).
                                      '; Annotation: '.${$args{COG2008_to_def}}->{${$args{COG2008_to_acc}}->{$COG_acc}->{'COG2008code'}}->{'annotation'}.
                                      ' '
                                      ;
            }
         } else {
            $additional_COG_info = ' COG2008: Not def; Class: X; Genome: Not defined'." " x ($genome_length - 11).'; Annotation: Not available';
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2014') {
         #COG2014 defined?
         my $COG_protID = $_; #ACC is different here
         $COG_protID    =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         if (${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'} =~ m/\w+/) {
            $additional_COG_info = ' COG2014: '.    ${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}.
                                   '; Class: '.     ${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}.', '.${$args{COGletter_to_family}}->{${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}}. ' ' x ($class_length - length(${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}) - 2 - length(${$args{COGletter_to_family}}->{${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}})).
                                   '; Genome: '.    ${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'genome_abbr'}.' ' x ($genome_length - length(${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'genome_abbr'})).
                                   '; Annotation: '.${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'annotation'}.
                                   ' '
                                   ;
         } else {
            $additional_COG_info = ' COG2014: Not def; Class: X; Genome: Not defined'." " x ($genome_length - 11).'; Annotation: Not available';
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'POG2013') {
         #POG2013 defined?
         if (${$args{POG2013_to_gi}}->{$_}->{'class'} =~ m/\w+/) {
            $additional_COG_info = ' POG2013: '.    ${$args{POG2013_to_gi}}->{$_}->{'class'}.
                                   '; Phage: '.     ${$args{POG2013_to_gi}}->{$_}->{'phage'}.' ' x ($genome_length - length(${$args{POG2013_to_gi}}->{$_}->{'phage'})).
                                   '; Annotation: '.${$args{POG2013_to_gi}}->{$_}->{'annotation'}.
                                   ' '
                                   ;
         } else {
            $additional_COG_info = ' POG2013: Not def; Phage: Not defined'." " x ($genome_length - 11).'; Annotation: Not available ';
         }
      }

      $acc_temp = "\>$_";
      my $nometa = quotemeta ($_);
      'reset' =~ m/reset/;
      $alignments =~ m/\>$nometa(.*?)\n\s+Score = /s; #get new overview line
      $temp = $1;
      unless (defined $temp) {
         $alignments =~ m/\>lcl\|$nometa(.*?)\n\s+Score = /s; #get new overview line
         $temp = $1;
      }
      #POG?
      unless (defined $temp) {
         $alignments =~ m/\>lcl.+? ref\|$nometa\|\s+(.*?)\n\s+Score \= /s;
         $temp = $1;
      }
      next unless (defined $temp);

      $temp =~ s/\s*(Length\s*=\s*.+)$//;
      $subj_length = $1;
      $subj_length =~ s/\s+//g;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $orig_blast_overview =~ m/$nometa[^\n]*?(\d+)\s+([^ ]+)\s*\n/;
      my $score = $1;
      my $evalue = $2;
      if ($evalue =~ /^e/) {$evalue = "1".$evalue};

      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
         $_ =~ s/^\w+\|//;
         $new_overview .= $acc_temp." " x ($acc_length - length($acc_temp)).$additional_COG_info."~" x ($annotation_length - length($additional_COG_info) + 2)." ".$subj_length."  Score=".$score."  Expect=".$evalue."\n";
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2008') {
         $_ =~ s/^\w+\|//;
         $new_overview .= $acc_temp." " x ($acc_length - length($acc_temp)).$additional_COG_info."~" x ($annotation_length - length($additional_COG_info) + 2)." ".$subj_length."  Score=".$score."  Expect=".$evalue."\n";
      } elsif (${$args{auto_ini_ref}}{COG_type} =~ m/(COG2014|arCOG2014|arCOG)/) {
         my $COG_protID = $_; #ACC is different here and is protein ID
         $COG_protID    =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         $new_overview .= '>'.$COG_protID." " x ($acc_length - length($COG_protID)).$additional_COG_info."~" x ($annotation_length - length($additional_COG_info) + 2)." ".$subj_length."  Score=".$score."  Expect=".$evalue."\n";
      } elsif (${$args{auto_ini_ref}}{COG_type} =~ m/POG2013/) {
         $new_overview .= '>'.$_." " x ($acc_length - length($_)).$additional_COG_info."~" x ($annotation_length - length($additional_COG_info) + 2)." ".$subj_length."  Score=".$score."  Expect=".$evalue."\n";
      } else {
         $_ =~ s/^\w+\|//;
         $new_overview .= $acc_temp." " x ($acc_length - length($acc_temp)).$additional_COG_info.$temp."~" x ($line_length - length($temp)+5)." ".$subj_length."  Score=".$score."  Expect=".$evalue."\n";
      }
   }

   #attach the alignment to new overview
   $line1 = "\nBlast overview\n\n".$new_overview."\n\n\nComplete list of Blast results\n\n\n".$alignments;
   return ($line1);

}

sub reformat_COG_threaded {
   my %args                  = @_;
   my $line1                 = '';
   my $orig_blast_overview   = '';
   my $alignments            = '';
   my @acc                   = ();
   my $line_length           = 0;
   my $acc_length            = 0;
   my $genome_length         = 0;
   my $annotation_length     = 0;
   my $class_length          = 0;
   my $subj_length           = 0;
   my $new_overview          = '';

   #catch original blast overview
   $args{blast_result} =~ m/Sequences producing significant alignments:                      \(bits\) Value\s*\n(.*?)\n\n\>/is;
   $orig_blast_overview = $1;
   unless (defined $orig_blast_overview) {
      $args{blast_result} =~ m/Sequences producing significant alignments:\s+\(Bits\)\s+Value\s*\n(.*?)\n\n\>/is;
      $orig_blast_overview = $1;
   }

   unless (defined ($orig_blast_overview)) {
      return ("\n\nNo COG Hits found\n");
   }

   if ($orig_blast_overview !~ /\w+/) {
      return ("\n\nNo COG Hits found\n");
   }

   $orig_blast_overview = "\n".$orig_blast_overview."\n";

   #catch aligments
   $args{blast_result} =~ m/\n\n\>(.+)Lambda\s+K\s+H/s;
   $alignments = '>'.$1;
   unless ($alignments =~ /\S{2,}/) {
      print "Could not grab Blast alignments from current aa sequence";
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in COG Blast analysis".
                  "\nCould not grab Blast alignments from current aa sequence for $args{input_file}\n\n";
      close WRITE;
   }

   #build array for accession numbers from blast overview
   if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
      @acc = ($orig_blast_overview =~ m/\n(\S*?) +\d+ +[^\n]+/sg);
      unless ((scalar @acc) >= 1) {
         @acc = ($orig_blast_overview =~ m/lcl\|(\S*?)\s+/sg);
      }
      unless ((scalar @acc) >= 1) {
         @acc = ($orig_blast_overview =~ m/\n(\S+)\s+/sg);
      }
   } elsif (${$args{auto_ini_ref}}{COG_type} =~ m/^(arCOG|arCOG2014|COG2008|COG2014)$/) {
      @acc = ($orig_blast_overview =~ m/\n(\w+\|\S+)\|? +/sg);
      unless ((scalar @acc) >= 1) {
         @acc = ($orig_blast_overview =~ m/\n(\S+)\s+/sg);
      }
   } elsif (${$args{auto_ini_ref}}{COG_type} =~ m/^POG2013$/) {
      @acc = ($orig_blast_overview =~ m/\n.*?ref\|([^\|]+)\|/sg); #needs 'ref' acc in all cases
      unless ((scalar @acc) >= 1) {
         @acc = ($orig_blast_overview =~ m/\n(\S+)\s+/sg);
      }
   }
   @acc = splice(@acc,0,${$args{auto_ini_ref}}{blast_entry_number}); #restrict number of entries to set value in $key

   #built new Blast overview
   #find out maximum text length

   #first, max acc name length
   foreach (@acc) {
      my $acc_temp    = '';
      my $temp        = '';
      $acc_temp       = "\>$_";

      my $nometa  = quotemeta ($_);
      'reset'     =~ m/reset/;
      $alignments =~ m/\>($nometa.*?)\n\s+Score \= /s;
      $temp = $1;
      unless (defined $temp) {
         $alignments =~ m/\>\w+\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }
      #POG?
      unless (defined $temp) {
         $alignments =~ m/\>lcl.+? ref\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }

      #error message for standard COG only, all other COGs may have too many hits to be listed in alignments
      unless ($temp =~ /\w+/) {
         if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
            print "Could not grab overview line from current entry $_\n in file $args{input_file} for ORF id $args{orf_id}\nSearch term: $nometa...\n";
            open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
            print WRITE "Error in COG Blast analysis".
                        "\nCould not grab overview line from current entry $_\n$args{blast_result}\n\n";
            close WRITE;
         } else {
            next;
         }
      }

      #determine maximum length of entry line
      $temp =~ s/\n?\s+Length\s*=\s*.+$//;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $temp =~ s/\s+$//g;

      #determine max length of acc
      if (length($acc_temp) > $acc_length) {$acc_length  = length($acc_temp)};
   }

   #next, max genome name length
   foreach (@acc) {
      my $acc_temp    = '';
      my $temp        = '';
      $acc_temp       = "\>$_";

      my $nometa  = quotemeta ($_);
      'reset'     =~ m/reset/;
      $alignments =~ m/\>($nometa.*?)\n\s+Score \= /s;
      $temp = $1;
      unless (defined $temp) {
         $alignments =~ m/\>\w+\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }
      #POG?
      unless (defined $temp) {
         $alignments =~ m/\>lcl.+? ref\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }

      #error message for standard COG only, all other COGs may have too many hits to be listed in alignments
      unless ($temp =~ /\w+/) {
         if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
            print "Could not grab overview line from current entry $_\n in file $args{input_file} for ORF id $args{orf_id}\nSearch term: $nometa...\n";
            open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
            print WRITE "Error in COG Blast analysis".
                        "\nCould not grab overview line from current entry $_\n$args{blast_result}\n\n";
            close WRITE;
         } else {
            next;
         }
      }

      #determine maximum length of entry line
      $temp =~ s/\n?\s+Length\s*=\s*.+$//;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $temp =~ s/\s+$//g;

      #determine genome & annotation name length for COG 2003
      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
         if (length(${$args{COGcode_to_phylo}}->{$_}) > $genome_length) {
            $genome_length = length(${$args{COGcode_to_phylo}}->{$_});
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2008') {
         $_ =~ s/^\w+\|//;
         if (length(${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$_}->{'genome_abbr'}}) > $genome_length) {
            $genome_length = length(${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$_}->{'genome_abbr'}});
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2014') {
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         if (length(${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'genome_abbr'}) > $genome_length) {
            $genome_length = length(${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'genome_abbr'});
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'arCOG2014') { #arCOG doesn't have genome name info
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         if (length(${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'genome'}) > $genome_length) {
            $genome_length = length(${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'genome'});
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'POG2013') {
         if (length(${$args{POG2013_to_gi}}->{$_}->{'phage'}) > $genome_length) {
            $genome_length = length(${$args{POG2013_to_gi}}->{$_}->{'phage'});
         }
      }
   }

   #next, max class name length
   foreach (@acc) {
      my $acc_temp    = '';
      my $temp        = '';
      $acc_temp       = "\>$_";

      my $nometa  = quotemeta ($_);
      'reset'     =~ m/reset/;
      $alignments =~ m/\>($nometa.*?)\n\s+Score \= /s;
      $temp = $1;
      unless (defined $temp) {
         $alignments =~ m/\>\w+\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }
      #POG?
      unless (defined $temp) {
         $alignments =~ m/\>lcl.+? ref\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }

      #error message for standard COG only, all other COGs may have too many hits to be listed in alignments
      unless ($temp =~ /\w+/) {
         if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
            print "Could not grab overview line from current entry $_\n in file $args{input_file} for ORF id $args{orf_id}\nSearch term: $nometa...\n";
            open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
            print WRITE "Error in COG Blast analysis".
                        "\nCould not grab overview line from current entry $_\n$args{blast_result}\n\n";
            close WRITE;
         } else {
            next;
         }
      }

      #determine maximum length of entry line
      $temp =~ s/\n?\s+Length\s*=\s*.+$//;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $temp =~ s/\s+$//g;

      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2008') {
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         if (length(${$args{COG2008_to_acc}}->{$COG_protID}->{'COG2008class'}) > $class_length) {
            $class_length = length(${$args{COG2008_to_acc}}->{$COG_protID}->{'COG2008class'});
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2014') {
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         if (length(${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}.', '.${$args{COGletter_to_family}}->{${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}}) > $class_length) {
            $class_length = length(${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}.', '.${$args{COGletter_to_family}}->{${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}});
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'arCOG') {
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         if (length(${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'COGclass'}) > $class_length) {
            $class_length = length(${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'COGclass'});
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'arCOG2014') {
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         if (length(${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGclass'}.', '.${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'}) > $class_length) {
            $class_length = length(${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGclass'}.', '.${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'});
         }
      }
   }

   #next, max annotation length
   foreach (@acc) {
      my $acc_temp    = '';
      my $temp        = '';
      $acc_temp       = "\>$_";

      my $nometa  = quotemeta ($_);
      'reset'     =~ m/reset/;
      $alignments =~ m/\>($nometa.*?)\n\s+Score \= /s;
      $temp = $1;
      unless (defined $temp) {
         $alignments =~ m/\>\w+\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }
      #POG?
      unless (defined $temp) {
         $alignments =~ m/\>lcl.+? ref\|($nometa.*?)\n\s+Score \= /s;
         $temp = $1;
      }

      #error message for standard COG only, all other COGs may have too many hits to be listed in alignments
      unless ($temp =~ /\w+/) {
         if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
            print "Could not grab overview line from current entry $_\n in file $args{input_file} for ORF id $args{orf_id}\nSearch term: $nometa...\n";
            open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
            print WRITE "Error in COG Blast analysis".
                        "\nCould not grab overview line from current entry $_\n$args{blast_result}\n\n";
            close WRITE;
         } else {
            next;
         }
      }

      #determine maximum length of entry line
      $temp =~ s/\n?\s+Length\s*=\s*.+$//;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $temp =~ s/\s+$//g;

      #determine genome & annotation name length for COG 2003
      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
         my $additional_COG_info = 'COG2003: '.      ${$args{COGcode_to_number}}->{$_}.
                                  '; Letter: '.      ${$args{COGcode_to_letter}}->{$_}.
                                  '; Genome: '.      ${$args{COGcode_to_phylo}}->{$_}.
                                  '; Annotation: '.  ${$args{COGcode_to_header}}->{$_}.
                                  ' '
                                 ;
         if (length($additional_COG_info) > $annotation_length) {
            $annotation_length = length($additional_COG_info);
         }
      }

      #determine genome & annotation name length for COG 2008
      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2008') {
         $_ =~ s/^\w+\|//;
         my $additional_COG_info = '';
         #COG2008 defined?
         if (${$args{COG2008_to_acc}}->{$_}->{'COG2008code'} =~ m/\w+/) {
            $additional_COG_info = ' COG2008: '.   ${$args{COG2008_to_acc}}->{$_}->{'COG2008code'}.
                                   '; Class: '.     ${$args{COG2008_to_acc}}->{$_}->{'COG2008class'}.
                                   '; Genome: '.    ${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$_}->{'genome_abbr'}}." " x ($genome_length - length(${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$_}->{'genome_abbr'}})).
                                   '; Annotation: '.${$args{COG2008_to_def}}->{${$args{COG2008_to_acc}}->{$_}->{'COG2008code'}}->{'annotation'}.
                                   ' '
                                   ;
            if (${$args{COG2008_to_acc}}->{$_}->{'COG2008code'} eq 'NULL') {
               $additional_COG_info = ' COG2008: not def'.
                                      '; Class: not def'.
                                      '; Genome: '.    ${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$_}->{'genome_abbr'}}." " x ($genome_length - length(${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$_}->{'genome_abbr'}})).
                                      '; Annotation: '.${$args{COG2008_to_def}}->{${$args{COG2008_to_acc}}->{$_}->{'COG2008code'}}->{'annotation'}.
                                      ' '
                                      ;
            }
         } else {
            $additional_COG_info = ' COG2008: Not defined; Class: X; Genome: Not defined'." " x ($genome_length - 11).'; Annotation: Not available ';
         }
         if (length($additional_COG_info) > $annotation_length) {
            $annotation_length = length($additional_COG_info);
         }
      }

      #determine info length for COG 2014
      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2014') {
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         my $additional_COG_info = ' COG2014: '.    ${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}.
                                   '; Class: '.     ${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}.', '.${$args{COGletter_to_family}}->{${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}}. ' ' x ($class_length - length(${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}) - 2 - length(${$args{COGletter_to_family}}->{${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}})).
                                   '; Genome: '.    ${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'genome_abbr'}.' ' x ($genome_length - length(${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'genome_abbr'})).
                                   '; Annotation: '.${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'annotation'}.
                                   ' '
                                   ;
         if (length($additional_COG_info) > $annotation_length) {
            $annotation_length = length($additional_COG_info);
         }
      }

      #determine info length for arCOG
      if (${$args{auto_ini_ref}}{COG_type} eq 'arCOG') {
         my $additional_COG_info = '';
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         #arCOG defined?
         if (${$args{arCOGacc_to_org}} ->{$COG_protID}->{'arCOG'} =~ m/\w+/) {
            $additional_COG_info = ' arCOG: '  .    ${$args{arCOGacc_to_org}} ->{$COG_protID}->{'arCOG'}.
                                   '; Class: '  .   ${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'COGclass'}.
                                   '; Cluster: '.   ${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'COGcluster'}.
                                   '; Annotation: '.${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'arCOGannotation'}.
                                   ' '
                                   ;
         } else {
            $additional_COG_info = ' arCOG: Not def   ; Class: X; Cluster: Not def'.' Annotation: Not available ';
         }
         if (length($additional_COG_info) > $annotation_length) {
            $annotation_length = length($additional_COG_info);
         }
      }

      #determine info length for arCOG2014
      if (${$args{auto_ini_ref}}{COG_type} eq 'arCOG2014') {
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         my $additional_COG_info = ' arCOG2014: '.  ${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}.
                                   '; Class: '.     ${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGclass'}.', '.
                                                    ${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'}.
                                                    ' ' x ($class_length - length(${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGclass'}) - 2 - length(${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'})).
                                   '; Genome: '.    ${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'genome'}.' ' x ($genome_length - length(${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'genome'})).
                                   '; Annotation: '.${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'}.
                                   ' '
                                   ;
         if (length($additional_COG_info) > $annotation_length) {
            $annotation_length = length($additional_COG_info);
         }
      }

      #determine info length for POG2013
      if (${$args{auto_ini_ref}}{COG_type} eq 'POG2013') {
         my $additional_COG_info = '';
         #POG2013 defined?
         if (${$args{POG2013_to_gi}}->{$_}->{'class'} =~ m/\w+/) {
            $additional_COG_info = ' POG2013: '.    ${$args{POG2013_to_gi}}->{$_}->{'class'}.
                                   '; Phage: '.     ${$args{POG2013_to_gi}}->{$_}->{'phage'}.' ' x ($genome_length - length(${$args{POG2013_to_gi}}->{$_}->{'phage'})).
                                   '; Annotation: '.${$args{POG2013_to_gi}}->{$_}->{'annotation'}.
                                   ' '
                                   ;
         } else {
            $additional_COG_info = ' POG2013: Not def; Phage: Not defined'." " x ($genome_length - 11).'; Annotation: Not available ';
         }

         if (length($additional_COG_info) > $annotation_length) {
            $annotation_length = length($additional_COG_info);
         }
      }
   }

   #generate new overview
   $orig_blast_overview .= "\n";
   foreach (@acc) {
      my $acc_temp            = '';
      my $temp                = '';
      my $subj_length         = '';
      my $COG_acc             = '';
      my $additional_COG_info = '';

      $COG_acc = $_;
      $COG_acc =~ s/^\S+?\|//;
      $COG_acc =~ s/[\|\s]//g;

      #define proper acc code for COG dbs
      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
         $additional_COG_info = 'COG2003: '.      ${$args{COGcode_to_number}}->{$COG_acc}.
                               '; Letter: '.      ${$args{COGcode_to_letter}}->{$COG_acc}.
                               '; Genome: '.      ${$args{COGcode_to_phylo}}->{$COG_acc}.
                               '; Annotation: '.  ${$args{COGcode_to_header}}->{$COG_acc}.
                               ' '
                               ;
      }

      if (${$args{auto_ini_ref}}{COG_type} eq 'arCOG') {
         #arCOG defined?
         my $COG_protID = $_;
         $COG_protID =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         if (${$args{arCOGacc_to_org}} ->{$COG_protID}->{'arCOG'} =~ m/\w+/) {
            $additional_COG_info = ' arCOG: '  .     ${$args{arCOGacc_to_org}} ->{$COG_protID}->{'arCOG'}.
                                   '; Class: '  .    ${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'COGclass'}.
                                   '; Cluster: '.    ${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'COGcluster'}.
                                   '; Annotation: '. ${$args{arCOGcode_to_COG}}->{${$args{arCOGacc_to_org}}->{$COG_protID}->{'arCOG'}}->{'arCOGannotation'}.
                                   ' '
                                   ;
         } else {
            $additional_COG_info = ' arCOG: Not def   ; Class: X; Cluster: Not def'.' Annotation: Not available ';
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'arCOG2014') {
         #arCOG2014 defined?
         my $COG_protID = $_; #ACC is different here
         $COG_protID    =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;

         if (${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'} =~ m/\w+/) {
            $additional_COG_info = ' arCOG2014: '.  ${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}.
                                   '; Class: '.     ${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGclass'}.', '.
                                                    ${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'}.
                                                    ' ' x ($class_length - length(${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGclass'}) - 2 - length(${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'})).
                                   '; Genome: '.    ${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'genome'}.' ' x ($genome_length - length(${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'genome'})).
                                   '; Annotation: '.${$args{arCOGcode2014_to_COG}}->{${$args{arCOGacc2014_to_org}}->{$COG_protID}->{'arCOG2014'}}->{'arCOGannotation'}.
                                   ' '
                                   ;
         } else {
            $additional_COG_info = ' arCOG2014: Not def   ; Class: X; Genome: Not defined'." " x ($genome_length - 11).'; Annotation: Not available ';
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2008') {
         #COG2008 defined?
         if (${$args{COG2008_to_acc}}->{$COG_acc}->{'COG2008code'} =~ m/\w+/) {
            $additional_COG_info = ' COG2008: '.   ${$args{COG2008_to_acc}}->{$COG_acc}->{'COG2008code'}.
                                   '; Class: '.     ${$args{COG2008_to_acc}}->{$COG_acc}->{'COG2008class'}.
                                   '; Genome: '.    ${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$COG_acc}->{'genome_abbr'}}." " x ($genome_length - length(${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$COG_acc}->{'genome_abbr'}})).
                                   '; Annotation: '.${$args{COG2008_to_def}}->{${$args{COG2008_to_acc}}->{$COG_acc}->{'COG2008code'}}->{'annotation'}.
                                   ' '
                                   ;
            if (${$args{COG2008_to_acc}}->{$COG_acc}->{'COG2008code'} eq 'NULL') {
               $additional_COG_info = ' COG2008: not def'.
                                      '; Class: not def'.
                                      '; Genome: '.    ${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$COG_acc}->{'genome_abbr'}}." " x ($genome_length - length(${$args{COG2008genomes}}->{${$args{COG2008_to_acc}}->{$COG_acc}->{'genome_abbr'}})).
                                      '; Annotation: '.${$args{COG2008_to_def}}->{${$args{COG2008_to_acc}}->{$COG_acc}->{'COG2008code'}}->{'annotation'}.
                                      ' '
                                      ;
            }
         } else {
            $additional_COG_info = ' COG2008: Not def; Class: X; Genome: Not defined'." " x ($genome_length - 11).'; Annotation: Not available ';
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2014') {
         #COG2014 defined?
         my $COG_protID = $_; #ACC is different here
         $COG_protID    =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         my $COG_code_annotation = ${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}; #make sure this is one letter only
         $COG_code_annotation    =~ s/^(.).+/$1/;
         if (${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'} =~ m/\w+/) {
            $additional_COG_info = ' COG2014: '.    ${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}.
                                   '; Class: '.     ${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}.', '.
                                                    ${$args{COGletter_to_family}}->{$COG_code_annotation }.
                                                    ' ' x ($class_length - length(${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'class'}) - 2 - length(${$args{COGletter_to_family}}->{$COG_code_annotation})).
                                   '; Genome: '.    ${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'genome_abbr'}.' ' x ($genome_length - length(${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'genome_abbr'})).
                                   '; Annotation: '.${$args{COG2014_to_def}}->{${$args{COG2014_to_acc}}->{${$args{COG2014_to_refseq}}->{$COG_protID}}->{'COG2014code'}}->{'annotation'}.
                                   ' '
                                   ;
         } else {
            $additional_COG_info = ' COG2014: Not def; Class: X; Genome: Not defined'." " x ($genome_length - 11).'; Annotation: Not available ';
         }
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'POG2013') {
         #POG2013 defined?
         if (${$args{POG2013_to_gi}}->{$_}->{'class'} =~ m/\w+/) {
            $additional_COG_info = ' POG2013: '.    ${$args{POG2013_to_gi}}->{$_}->{'class'}.
                                   '; Phage: '.     ${$args{POG2013_to_gi}}->{$_}->{'phage'}.' ' x ($genome_length - length(${$args{POG2013_to_gi}}->{$_}->{'phage'})).
                                   '; Annotation: '.${$args{POG2013_to_gi}}->{$_}->{'annotation'}.
                                   ' '
                                   ;
         } else {
            $additional_COG_info = ' POG2013: Not def; Phage: Not defined'." " x ($genome_length - 11).'; Annotation: Not available ';
         }
      }

      $acc_temp   = "\>$_";
      my $nometa  = quotemeta ($_);
      'reset'     =~ m/reset/;
      $alignments =~ m/\>$nometa(.*?)\n\s+Score = /s; #get new overview line
      $temp       = $1;
      unless (defined $temp) {
         $alignments =~ m/\>\w+\|$nometa(.*?)\n\s+Score = /s; #get new overview line
         $temp       = $1;
      }
      #POG?
      unless (defined $temp) {
         $alignments =~ m/\>lcl.+? ref\|$nometa\|\s+(.*?)\n\s+Score \= /s;
         $temp = $1;
      }
      next unless (defined $temp);

      $temp =~ s/\s*(Length\s*=\s*.+)$//;
      $subj_length = $1;
      $subj_length =~ s/\s+//g;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $orig_blast_overview =~ m/$nometa[^\n]*?(\d+)\s+([^ ]+)\s*\n/;
      my $score = $1;
      my $evalue = $2;
      if ($evalue =~ /^e/) {$evalue = "1".$evalue};

      if (${$args{auto_ini_ref}}{COG_type} eq 'COG2003') {
         $_ =~ s/^\w+\|//;
         $new_overview .= $acc_temp." " x ($acc_length - length($acc_temp)).$additional_COG_info."~" x ($annotation_length - length($additional_COG_info) + 2)." ".$subj_length."  Score=".$score."  Expect=".$evalue."\n";
      } elsif (${$args{auto_ini_ref}}{COG_type} eq 'COG2008') {
         $_ =~ s/^\w+\|//;
         $new_overview .= $acc_temp." " x ($acc_length - length($acc_temp)).$additional_COG_info."~" x ($annotation_length - length($additional_COG_info) + 2)." ".$subj_length."  Score=".$score."  Expect=".$evalue."\n";
      } elsif (${$args{auto_ini_ref}}{COG_type} =~ m/(COG2014|arCOG2014|arCOG)/) {
         my $COG_protID = $_; #ACC is different here and is protein ID
         $COG_protID    =~ s/^\S+?\|//;
         $COG_protID =~ s/\.\d+[\|\s]*$//g;
         $new_overview .= '>'.$COG_protID." " x ($acc_length - length($COG_protID)).$additional_COG_info."~" x ($annotation_length - length($additional_COG_info) + 2)." ".$subj_length."  Score=".$score."  Expect=".$evalue."\n";
      } elsif (${$args{auto_ini_ref}}{COG_type} =~ m/POG2013/) {
         $new_overview .= '>'.$_." " x ($acc_length - length($_)).$additional_COG_info."~" x ($annotation_length - length($additional_COG_info) + 2)." ".$subj_length."  Score=".$score."  Expect=".$evalue."\n";
      } else {
         $_ =~ s/^\w+\|//;
         $new_overview .= $acc_temp." " x ($acc_length - length($acc_temp)).$additional_COG_info.$temp."~" x ($line_length - length($temp)+5)." ".$subj_length."  Score=".$score."  Expect=".$evalue."\n";
      }
   }

   #attach the alignment to new overview
   $line1 = "\nBlast overview\n\n".$new_overview."\n\n\nComplete list of Blast results\n\n\n".$alignments;
   return ($line1);

}

sub COG_file {
   my %args = (retry_counter  => 0,
               @_
              ); #limit the numbers of Blast retries in case of mising results
   my (@list, @COG_tasks);
   my $progress = 0;

   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'COG analysis',
                   label        => 'COG analysis'
                  );
   &show_pbar_2;

   #skip exisiting results?
   ${$args{progress_bar}}->configure(-label=>" COG ");
   ${$args{main_window}}->update;

   #map correct entries for respective input file
   @list = grep /$args{input_file}\___/, @{$args{combined_orf_ref}};

   #remove existing entries if re-use results
   if (${$args{auto_ini_ref}}{reuse_results} == 1) {
      my (@jobs, @COG_results, %exist_cog_file);
      &update_pbar_2(title        => 'COG analysis',
                     label        => "Testing for existing results",
                     progress     => 1,
                    );
      #read existing COG results
      opendir SEQINPUT, ${$args{ini_ref}}{COG_results};
      @COG_results = grep /$args{input_file}_COG_\d+$/, readdir(SEQINPUT);
      #create hash
      %exist_cog_file = ();
      foreach (@COG_results) {
         $exist_cog_file{$_} = '1';
      }

      #iterate through gene model
      foreach my $entry (@list) {
         my ($ID);
         $entry =~ m/$args{input_file}\___(\d+)___/;
         $ID    = $1;

         #skip COG if result file exists and re-use is activated
         next if (exists $exist_cog_file{$args{input_file}.'_COG_'.$ID} && -s ${$args{ini_ref}}{COG_results}.'/'.$args{input_file}.'_COG_'.$ID > 0);
         #add if required
         push (@jobs, $entry);
      }
      @list = @jobs;
      undef @jobs;
      undef @COG_results;
      undef %exist_cog_file;
   }

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

   #get boundary parameters
   my $ID = '';
   my $left_bd = '';
   my $right_bd = '';
   my $orientation = '';

   #max number of fasta entries
   my $max_count = $#list + 1;

   #run threaded or CPU clustered COG
   if (${$args{auto_ini_ref}}{COG_cluster} == 1) {
      #iterate over array
      foreach my $entry (@list) {
         #update progress bar
         $progress++;
         &update_pbar_2(title        => 'COG analysis',
                        label        => "Blasting $args{input_file}, $progress of $max_count",
                        progress     => ($progress / $max_count) * 100,
                       );
         ${$args{main_window}}->update;

         ($ID, $left_bd, $right_bd, $orientation) = "";
         $entry =~ m/$args{input_file}\___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
         $ID = $1;
         $left_bd = $2;
         $right_bd = $3;
         $orientation = $4;
         my $query_seq_nt = "";

         unless (defined $orientation && $orientation =~ /(sense|antisense)/) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Error parsing gene model for $args{input_file}",
                                                          -bitmap  => 'error',
                                                          -buttons => ['ok']);
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
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
                               module       => 'COG Blast',
                               orientation  => $orientation,
                               sequence     => $query_seq_nt,
                               filename     => $args{input_file},
                               left_bd      => $left_bd,
                               right_bd     => $right_bd
                              );

         my ($blastresult) = &COG_blast_seq(main_window          => $args{main_window},
                                            progress_bar         => $args{progress_bar},
                                            auto_ini_ref         => $args{auto_ini_ref},
                                            ini_ref              => $args{ini_ref},
                                            sequence_nt          => $query_seq_nt,
                                            sequence_aa          => $query_seq_aa,
                                            left_bd              => $left_bd,
                                            right_bd             => $right_bd,
                                            orientation          => $orientation,
                                            input_file           => $args{input_file},
                                            orf_id               => $ID,
                                            COGcode_to_number    => $args{COGcode_to_number},
                                            COGcode_to_header    => $args{COGcode_to_header},
                                            COGcode_to_letter    => $args{COGcode_to_letter},
                                            COGcode_to_phylo     => $args{COGcode_to_phylo},
                                            COGletter_to_family  => $args{COGletter_to_family},
                                            COG2008genomes       => $args{COG2008genomes},
                                            COG2008_to_def       => $args{COG2008_to_def},
                                            COG2008_to_acc       => $args{COG2008_to_acc},
                                            arCOGcode_to_COG     => $args{arCOGcode_to_COG},
                                            arCOGacc_to_org      => $args{arCOGacc_to_org},
                                            COG2014genomes       => $args{COG2014genomes},
                                            COG2014_to_def       => $args{COG2014_to_def},
                                            COG2014_to_acc       => $args{COG2014_to_acc},
                                            COG2014_to_refseq    => $args{COG2014_to_refseq},
                                            arCOGcode2014_to_COG => $args{arCOGcode2014_to_COG},
                                            arCOGacc2014_to_org  => $args{arCOGacc2014_to_org},
                                            POG2013_to_gi        => $args{POG2013_to_gi},
                                            COGletter_to_family  => $args{COGletter_to_family},
                                           );
         if ($blastresult eq '0') {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Error blasting ORF $ID in sequence $args{input_file} with COG database",
                                                          -bitmap  => 'error',
                                                          -buttons => ['ok']);
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            &hide_pbar_2;
            return (0); #return to main
         }
         #write Blast results
         open WRITE, "+>${$args{ini_ref}}{COG_results}\/$args{input_file}"."_COG_$ID";
         print WRITE $blastresult;
         close WRITE;
      }
   } elsif (${$args{auto_ini_ref}}{COG_threaded} == 1) {
      #iterate over array
      foreach my $entry (@list) {
         #update progress bar
         $progress++;
         if (($progress % 100) == 0) {
            &update_pbar_2(title        => 'COG analysis',
                           label        => "Reading gene model for $args{input_file}",
                           progress     => ($progress / $max_count) * 100,
                          );
            ${$args{main_window}}->update;
         }

         ($ID, $left_bd, $right_bd, $orientation) = "";
         $entry =~ m/$args{input_file}\___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
         $ID = $1;
         $left_bd = $2;
         $right_bd = $3;
         $orientation = $4;
         my $query_seq_nt = "";

         unless (defined $orientation && $orientation =~ /(sense|antisense)/) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Error parsing gene model for $args{input_file}",
                                                          -bitmap  => 'error',
                                                          -buttons => ['ok']);
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
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
                               module       => 'COG Blast',
                               orientation  => $orientation,
                               sequence     => $query_seq_nt,
                               filename     => $args{input_file},
                               left_bd      => $left_bd,
                               right_bd     => $right_bd
                              );
         #create array for new PFam tasks
         push (@COG_tasks, $args{input_file}.'___'.$ID.'___'.$left_bd.'___'.$right_bd.'___'.$orientation.'___'.$query_seq_aa.'___'.$query_seq_nt);
      }

      #start multithreading COG
      for (my $i=1; $i<=$#COG_tasks+2; $i = $i+${$args{auto_ini_ref}}{CPU}) {
         my $count = $i + ${$args{auto_ini_ref}}{CPU};
         if ($count > $#COG_tasks+2) {$count = $#COG_tasks+2};

         &update_pbar_2(title    => 'COG analysis',
                        label    => "COG analysis for $args{input_file}, $count of $max_count",
                        progress => ($count / ($#COG_tasks + 2)) * 100,
                       );
         ${$args{main_window}}->update;

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
               my $aa_sequence = "";
               my $nt_sequence = "";
               my $left_bd     = "";
               my $right_bd    = "";
               my $orientation = "";
               my $ORF_number  = "";
               $COG_tasks[$j-1] =~ m/$args{input_file}___(\d+)___(\d+)___(\d+)___(.+)___(.+)___(.+)/;
               $ORF_number  = $1;
               $left_bd     = $2;
               $right_bd    = $3;
               $orientation = $4;
               $aa_sequence = $5;
               $nt_sequence = $6;

               my ($blastresult) = &COG_blast_threaded(auto_ini_ref         => $args{auto_ini_ref},
                                                       ini_ref              => $args{ini_ref},
                                                       sequence_nt          => $nt_sequence,
                                                       sequence_aa          => $aa_sequence,
                                                       left_bd              => $left_bd,
                                                       right_bd             => $right_bd,
                                                       orientation          => $orientation,
                                                       input_file           => $args{input_file},
                                                       orf_id               => $ORF_number,
                                                       COGcode_to_number    => $args{COGcode_to_number},
                                                       COGcode_to_header    => $args{COGcode_to_header},
                                                       COGcode_to_letter    => $args{COGcode_to_letter},
                                                       COGcode_to_phylo     => $args{COGcode_to_phylo},
                                                       COGletter_to_family  => $args{COGletter_to_family},
                                                       COG2008genomes       => $args{COG2008genomes},
                                                       COG2008_to_def       => $args{COG2008_to_def},
                                                       COG2008_to_acc       => $args{COG2008_to_acc},
                                                       arCOGcode_to_COG     => $args{arCOGcode_to_COG},
                                                       arCOGacc_to_org      => $args{arCOGacc_to_org},
                                                       COG2014genomes       => $args{COG2014genomes},
                                                       COG2014_to_def       => $args{COG2014_to_def},
                                                       COG2014_to_acc       => $args{COG2014_to_acc},
                                                       COG2014_to_refseq    => $args{COG2014_to_refseq},
                                                       arCOGcode2014_to_COG => $args{arCOGcode2014_to_COG},
                                                       arCOGacc2014_to_org  => $args{arCOGacc2014_to_org},
                                                       POG2013_to_gi        => $args{POG2013_to_gi},
                                                       COGletter_to_family  => $args{COGletter_to_family},
                                                      );
               if ($blastresult eq '0') {
                  #print "Error in COG analysis. Error blasting ORF $ORF_number in sequence $args{input_file} with COG database";
                  open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
                  print WRITE "\nError in COG analysis. Error blasting ORF $ORF_number in sequence $args{input_file} with COG database\n\n";
                  close WRITE;
                  #exit safely
                  CORE::exit();
               }
               #write Blast results
               open WRITE, "+>${$args{ini_ref}}{COG_results}\/$args{input_file}"."_COG_$ORF_number";
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
   }
   undef @COG_tasks;

   #are really all COG results done? better check
   {
      my @Blast_results;
      &update_pbar_2(title        => 'COG Blasting',
                     label        => "Testing for complete results",
                     progress     => 1,
                    );

      #read existing results
      opendir SEQINPUT, ${$args{ini_ref}}{COG_results};
      @Blast_results = grep /^$args{input_file}_\d+$/, readdir(SEQINPUT);
      closedir SEQINPUT;
      #map correct entries for respective input file
      @list = grep /$args{input_file}\___/, @{$args{combined_orf_ref}};

      if ($#Blast_results ne $#list && $args{retry_counter} <= 3) {
         $args{retry_counter}++;
         &hide_pbar_2;
         undef @list;
         undef @Blast_results;
         &COG_file(main_window          => $args{main_window},
                   progress_bar         => $args{progress_bar},
                   auto_ini_ref         => $args{auto_ini_ref},
                   ini_ref              => $args{ini_ref},
                   directory            => ${$args{ini_ref}}{input_files},
                   input_file           => $args{input_file},
                   combined_orf_ref     => $args{gene_model_ref},
                   retry_counter        => $args{retry_counter},
                   COGcode_to_number    => $args{COGcode_to_number},
                   COGcode_to_header    => $args{COGcode_to_header},
                   COGcode_to_letter    => $args{COGcode_to_letter},
                   COGcode_to_phylo     => $args{COGcode_to_phylo},
                   COGletter_to_family  => $args{COGletter_to_family},
                   COG2008genomes       => $args{COG2008genomes},
                   COG2008_to_def       => $args{COG2008_to_def},
                   COG2008_to_acc       => $args{COG2008_to_acc},
                   arCOGcode_to_COG     => $args{arCOGcode_to_COG},
                   arCOGacc_to_org      => $args{arCOGacc_to_org},
                   COG2014genomes       => $args{COG2014genomes},
                   COG2014_to_def       => $args{COG2014_to_def},
                   COG2014_to_acc       => $args{COG2014_to_acc},
                   arCOGcode2014_to_COG => $args{arCOGcode2014_to_COG},
                   arCOGacc2014_to_org  => $args{arCOGacc2014_to_org},
                   POG2013_to_gi        => $args{POG2013_to_gi},
                   COGletter_to_family  => $args{COGletter_to_family},
                  );
      }
      undef @Blast_results;
   }
   undef @list;

   ${$args{progress_bar}}->configure(-label=>"Finished COG analysis");
   ${$args{main_window}}->update;
   &hide_pbar_2;
   ${$args{main_window}}->update;
}


1;
