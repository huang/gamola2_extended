#!/opt/ActivePerl-5.8/bin/perl
#Blast program calls and reformatting of results
#input arguments: progress_bar, nt_seq, left_bd, right_bd, orientation, ini_ref, auto_ini_ref


package ProgrammeModules::blast;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&blast_ig &blast_ig_clustered &blast_seq &blast_file &blast_hash);
use vars qw();

use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);
use db_selection::databases    qw(:DEFAULT);
use Cwd;

#local variables
my (%args, $ret_seq, $psi_aa_db, $psi_nt_db, $blastresult);

sub blast_seq {
   my %args = (aa_input_file => 'temp.aa', @_);
   my $blastresult = "";

   #reformat sequence into fasta format
   my ($seq_ref_nt) = &seq2fasta (main_window  => \$args{main_window},
                                  progress_bar => \$args{progress_bar},
                                  auto_ini_ref => $args{auto_ini_ref},
                                  ini_ref      => $args{ini_ref},
                                  sequence     => $args{sequence_nt},
                                  orientation  => $args{orientation},
                                  header       => $args{input_file}.'_'.$args{orf_id}
                                 );
   my ($seq_ref_aa) = &seq2fasta (main_window  => \$args{main_window},
                                  progress_bar => \$args{progress_bar},
                                  auto_ini_ref => $args{auto_ini_ref},
                                  ini_ref      => $args{ini_ref},
                                  sequence     => $args{sequence_aa},
                                  orientation  => $args{orientation},
                                  header       => $args{input_file}.'_'.$args{orf_id}
                                 );
   #get current working directory
   my $curdir = '';
   $curdir = getcwd();

   #convert to correct Blast type 'BlastN', 'BlastX', 'tBlastX', 'BlastP', 'gappedBlastP', 'PSI-Blast'
   my $blast_type;
   if (${$args{auto_ini_ref}}{blast_type}      =~ /^blastn$/i) {
      $blast_type = 'blastn';
   } elsif (${$args{auto_ini_ref}}{blast_type} =~ /^blastx$/i) {
      $blast_type = 'blastx';
   } elsif (${$args{auto_ini_ref}}{blast_type} =~ /^tblastx$/i) {
      $blast_type = 'tblastx';
   } elsif (${$args{auto_ini_ref}}{blast_type} =~ /^blastp$/i) {
      $blast_type = 'blastp';
   } elsif (${$args{auto_ini_ref}}{blast_type} =~ /^gappedblastp$/i) {
      $blast_type = 'blastp';
   } elsif (${$args{auto_ini_ref}}{blast_type} =~ /^PSI-Blast$/i) {
      $blast_type = 'psiblast';
   }

   if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
      #change to Blast directory
      chdir ${$args{ini_ref}}{blast_executables};
      if (${$args{auto_ini_ref}}{blast_type} =~ /(blastn|tblastx|blastx)/i) {
         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./blastall -a ${$args{auto_ini_ref}}{CPU} -p $blast_type -d ${$args{auto_ini_ref}}{full_blast_db}"; #child
               die "Cannot exec: $!";
            } else { # grandchild
               print ${$seq_ref_nt};
               CORE::exit;
            }
         }
      } elsif (${$args{auto_ini_ref}}{blast_type} =~ /(blastp)/i) {
         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./blastall -a ${$args{auto_ini_ref}}{CPU} -p $blast_type -d ${$args{auto_ini_ref}}{full_blast_db}"; #child
               die "Cannot exec: $!";
            } else { # grandchild
               print ${$seq_ref_aa};
               CORE::exit;
            }
         }
      } elsif (${$args{auto_ini_ref}}{blast_type} =~ /psi-blast/i) {
         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./blastpgp -d ${$args{auto_ini_ref}}{full_blast_db} -j 6";
               die "Cannot exec: $!";
            } else { # grandchild
               print ${$seq_ref_aa};
               CORE::exit;
            }
         }
      } elsif (${$args{auto_ini_ref}}{blast_type} =~ /gappedblastp/i) {
         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./blastpgp -a ${$args{auto_ini_ref}}{CPU} -d ${$args{auto_ini_ref}}{full_blast_db} -j 1"; # child
               die "Cannot exec: $!";
            } else { # grandchild
               print ${$seq_ref_aa};
               CORE::exit;
            }
         }
      }
   }
   elsif (${$args{auto_ini_ref}}{blast_plus} == 1) {
      #change to Blast_plus directory
      chdir ${$args{ini_ref}}{blast_plus_executables};

      if (${$args{auto_ini_ref}}{blast_type} =~ /(blastn|tblastx|blastx)/i) {
         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./$blast_type -num_threads ${$args{auto_ini_ref}}{CPU} -db ${$args{auto_ini_ref}}{full_blast_db}"; #child
               die "Cannot exec: $!";
            } else { # grandchild
               print ${$seq_ref_nt};
               CORE::exit;
            }
         }
      } elsif (${$args{auto_ini_ref}}{blast_type} =~ /(blastp|gappedblastp)/i) {

         #build options
         my $more_options = '';
         #set ungapped blastp
         if (${$args{auto_ini_ref}}{blast_type} =~ /^blastp$/i) {
            $more_options .= '-ungapped -comp_based_stats F';
         }

         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./$blast_type $more_options -num_threads ${$args{auto_ini_ref}}{CPU} -db ${$args{auto_ini_ref}}{full_blast_db}"; #child
               die "Cannot exec: $!";
            } else { # grandchild
               print ${$seq_ref_aa};
               CORE::exit;
            }
         }
      } elsif (${$args{auto_ini_ref}}{blast_type} =~ /psi-blast/i) {
         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./psiblast -num_threads ${$args{auto_ini_ref}}{CPU} -db ${$args{auto_ini_ref}}{full_blast_db} -num_iterations 6 -comp_based_stats 2"; #child
               die "Cannot exec: $!";
            } else { # grandchild
               print ${$seq_ref_aa};
               CORE::exit;
            }
         }
      } else {
         print "\nBlast type ${$args{auto_ini_ref}}{blast_type} for Blast plus not yet supported. Coming soon\n";
         CORE::exit;
      }
   }

   #change back to previous working dir
   chdir $curdir;

   #valid result?
   unless ($blastresult =~ /\*\*\s+No\s+hits\s+found\s+\*\*/is || $blastresult =~ /Identities/) {
      ${$args{progress_bar}}->configure(-label =>"Error in Blast analysis using database ${$args{auto_ini_ref}}{blast_db}");
      ${$args{main_window}}->update;
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "\nError in Blast analysis using database ${$args{auto_ini_ref}}{full_blast_db} and Sequence:\n ${$seq_ref_nt}\nBlastresult:\n$blastresult\n\n";
      close WRITE;
      return (0);
   };

   #build Blast header
   my $blast_header = "Performed Blast flavour:\t${$args{auto_ini_ref}}{blast_type}\n".
                      "Name of the input DNA-sequence:\t$args{input_file}\n\n".
                      "DNA sequence:\t$args{sequence_nt}\n".
                      "Deduced aminoacid sequence:\t$args{sequence_aa}\n\n".
                      "Length in aminoacids:\t".length($args{sequence_aa})."\n".
                      "Gene model summary:\tORF-designation $args{orf_id}\tLeft boundary $args{left_bd}\tRight boundary $args{right_bd}\tOrientation $args{orientation}\n";

   #reformat blast results
   ($blastresult) = &reformat_blast_res(main_window  => $args{main_window},
                                        progress_bar => $args{progress_bar},
                                        auto_ini_ref => $args{auto_ini_ref},
                                        ini_ref      => $args{ini_ref},
                                        blast_result => $blastresult,
                                       );

   return ($blast_header.$blastresult);
}

sub blast_threaded {
   my %args = (aa_input_file => 'temp.aa', @_);
   my $blastresult = "";

   #reformat sequence into fasta format
   my ($seq_ref_nt) = &seq2fasta (sequence     => $args{sequence_nt},
                                  orientation  => $args{orientation},
                                  header       => $args{input_file}.'_'.$args{orf_id}
                                 );
   my ($seq_ref_aa) = &seq2fasta (sequence     => $args{sequence_aa},
                                  orientation  => $args{orientation},
                                  header       => $args{input_file}.'_'.$args{orf_id}
                                 );

   #get current working directory
   my $curdir = '';
   $curdir    = getcwd();

   #convert to correct Blast type 'BlastN', 'BlastX', 'tBlastX', 'BlastP', 'gappedBlastP', 'PSI-Blast'
   my $blast_type;
   if (${$args{auto_ini_ref}}{blast_type} =~ /^blastn/i) {
      $blast_type = 'blastn';
   } elsif (${$args{auto_ini_ref}}{blast_type} =~ /^blastx/i) {
      $blast_type = 'blastx';
   } elsif (${$args{auto_ini_ref}}{blast_type} =~ /^tblastx/i) {
      $blast_type = 'tblastx';
   } elsif (${$args{auto_ini_ref}}{blast_type} =~ /^blastp/i) {
      $blast_type = 'blastp';
   } elsif (${$args{auto_ini_ref}}{blast_type} =~ /^gappedblastp$/i) {
      $blast_type = 'blastp';
   } elsif (${$args{auto_ini_ref}}{blast_type} =~ /^PSI-Blast$/i) {
      $blast_type = 'psiblast';
   }

   if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
      #change to Blast directory
      chdir ${$args{ini_ref}}{blast_executables};

      if (${$args{auto_ini_ref}}{blast_type} =~ /(blastn|tblastx)/i) {
         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./blastall -p ${$args{auto_ini_ref}}{blast_type} -d ${$args{auto_ini_ref}}{full_blast_db}"; #child
               die "Cannot exec: $!";
            } else { # grandchild
               print ${$seq_ref_nt};
               CORE::exit;
            }
         }
      } elsif (${$args{auto_ini_ref}}{blast_type} =~ /blastp/i) {
         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./blastall -p $blast_type -d ${$args{auto_ini_ref}}{full_blast_db}"; #child
               die "Cannot exec: $!";
            } else { # grandchild
               print ${$seq_ref_aa};
               CORE::exit;
            }
         }
      } elsif (${$args{auto_ini_ref}}{blast_type} =~ /blastx/i) {
         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./blastall -p ${$args{auto_ini_ref}}{blast_type} -d ${$args{auto_ini_ref}}{full_blast_db}"; #child
               die "Cannot exec: $!";
            } else { # grandchild
               print ${$seq_ref_nt};
               CORE::exit;
            }
         }
      } elsif (${$args{auto_ini_ref}}{blast_type} =~ /PSI-Blast/i) {
         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./blastpgp -d ${$args{auto_ini_ref}}{full_blast_db} -j 5";
               die "Cannot exec: $!";
            } else { # grandchild
               print ${$seq_ref_aa};
               CORE::exit;
            }
         }
      } elsif (${$args{auto_ini_ref}}{blast_type} =~ /gappedblastp/i) {
         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./blastpgp -d ${$args{auto_ini_ref}}{full_blast_db} -j 1"; # child
               die "Cannot exec: $!";
            } else { # grandchild
               print ${$seq_ref_aa};
               CORE::exit;
            }
         }
      }
   } elsif (${$args{auto_ini_ref}}{blast_plus} == 1) {
      #change to Blast_plus directory
      chdir ${$args{ini_ref}}{blast_plus_executables};

      if (${$args{auto_ini_ref}}{blast_type} =~ /(blastn|tblastx|blastx)/i) {
         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./$blast_type -db ${$args{auto_ini_ref}}{full_blast_db}"; #child
               die "Cannot exec: $!";
            } else { # grandchild
               print ${$seq_ref_nt};
               CORE::exit;
            }
         }
      } elsif (${$args{auto_ini_ref}}{blast_type} =~ /(blastp|gappedblastp)/i) {
         #build options
         my $more_options = '';
         #set ungapped blastp
         if (${$args{auto_ini_ref}}{blast_type} =~ /^blastp$/i) {
            $more_options = '-ungapped -comp_based_stats F';
         }

         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./$blast_type $more_options -db ${$args{auto_ini_ref}}{full_blast_db}"; #child
               die "Cannot exec: $!";
            } else { # grandchild
               print ${$seq_ref_aa};
               CORE::exit;
            }
         }
      } elsif (${$args{auto_ini_ref}}{blast_type} =~ /PSI-Blast/i) {
         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./$blast_type -db ${$args{auto_ini_ref}}{full_blast_db} -num_iterations 6 -comp_based_stats 2"; #child
               die "Cannot exec: $!";
            } else { # grandchild
               print ${$seq_ref_aa};
               CORE::exit;
            }
         }
      } else {
         print "\nBlast type ${$args{auto_ini_ref}}{blast_type} for Blast plus not yet supported. Coming soon\n";
         CORE::exit;
      }
   }

   #change back to previous working dir
   chdir $curdir;

   #valid result?
   unless ($blastresult =~ /\*\s*No\s*hits\s*found\s*\*/i || $blastresult =~ /Identities/) {
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in Blast analysis".
                  "\nNo Blast results found for $args{input_file} using database ${$args{auto_ini_ref}}{full_blast_db}\n\n";
      close WRITE;
      print "\nNo Blast results found for $args{input_file} using database ${$args{auto_ini_ref}}{full_blast_db}\n";
      return (0);
   };

   #build Blast header
   my $blast_header = "Performed Blast flavour:\t${$args{auto_ini_ref}}{blast_type}\n".
                      "Name of the input DNA-sequence:\t$args{input_file}\n\n".
                      "DNA sequence:\t$args{sequence_nt}\n".
                      "Deduced aminoacid sequence:\t$args{sequence_aa}\n\n".
                      "Length in aminoacids:\t".length($args{sequence_aa})."\n".
                      "Gene model summary:\tORF-designation $args{orf_id}\tLeft boundary $args{left_bd}\tRight boundary $args{right_bd}\tOrientation $args{orientation}\n";

   #reformat blast results
   ($blastresult) = &reformat_blast_threaded(auto_ini_ref => $args{auto_ini_ref},
                                             ini_ref      => $args{ini_ref},
                                             blast_result => $blastresult,
                                             header       => $args{input_file}.'_'.$args{orf_id}
                                            );

   return ($blast_header.$blastresult);
}

sub blast_hash_seq {
   my %args = (threaded => 0,
               @_);
   my $blastresult     = '';
   my $reformatted_res = '';
   my $CPU             = ${$args{auto_ini_ref}}{CPU}; #set no of CPUs for clustered Blast

   #test if threaded Blast
   if ($args{threaded} == 1) {
      $CPU = 1;
   }

   #get current working directory
   my $curdir = '';
   $curdir    = getcwd();

   if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
      #change to Blast directory
      chdir ${$args{ini_ref}}{blast_executables};

      if (open RESULT, "-|") { # original process
         local $/;
         $blastresult = <RESULT>;
      } else { # child
         if (open STDIN, "-|") { # child
            exec "./blastall -a $CPU -p $args{blast_flavour} -d $args{blast_db} -B $args{entry_count}"; #child
            die "Cannot exec: $!";
         } else { # grandchild
            print ${$args{entry_ref}};
            CORE::exit;
         }
      }
   } elsif (${$args{auto_ini_ref}}{blast_plus} == 1) {
      #change to Blast_plus directory
      chdir ${$args{ini_ref}}{blast_plus_executables};
      my $blast_type;
      if ($args{blast_flavour} =~ /^blastn/i) {
         $blast_type = 'blastn';
      } elsif ($args{blast_flavour} =~ /^blastx/i) {
         $blast_type = 'blastx';
      } elsif ($args{blast_flavour} =~ /^tblastx/i) {
         $blast_type = 'tblastx';
      }

      if (open RESULT, "-|") { # original process
         local $/;
         $blastresult = <RESULT>;
      } else { # child
         if (open STDIN, "-|") { # child
            exec "./$blast_type -db $args{blast_db} ";#-num_threads $CPU"; #child
            die "Cannot exec: $!";
         } else { # grandchild
            print ${$args{entry_ref}};
            CORE::exit;
         }
      }
   }

   #change back to previous working dir
   chdir $curdir;

   #valid result?
   unless ($blastresult =~ /\*\*\*\*\* No hits found \*\*\*\*\*/ || $blastresult =~ /Identities/) {
      #${$args{progress_bar}}->configure(-label =>"Error in Blast analysis using database ${$args{auto_ini_ref}}{blast_db}");
      #${$args{main_window}}->update;
      print "\nError in Blast analysis using database ${$args{auto_ini_ref}}{blast_db}\n";
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "\nError in Meta Blast analysis using database $args{blast_db} and Sequence:\n ${$args{entry_ref}}\nBlastresult:\n$blastresult...\n\n";
      close WRITE;
      return (0);
   };

   #reformat blast results
   my @Blast_res = split/\n\nReference:/,$blastresult;
   foreach my $entry (@Blast_res) {
      next unless ($entry =~ /\nQuery/s);
      my @individual_hits = split/\nQuery\=/, $entry;
      foreach my $hit (@individual_hits) {
         next if ($hit =~ m/Stephen F\. Altschul/); #ignore header section
         $reformatted_res .= &reformat_blast_hash(main_window         => $args{main_window},
                                                  progress_bar        => $args{progress_bar},
                                                  auto_ini_ref        => $args{auto_ini_ref},
                                                  ini_ref             => $args{ini_ref},
                                                  blast_result        => "\nQuery\=".$hit,
                                                  blast_flavour       => $args{blast_flavour},
                                                  blast_db            => $args{blast_db},
                                                  blast_display_limit => $args{blast_display_limit},
                                                 );
      }
   }
   undef $blastresult;

   return ($reformatted_res);
}

sub reformat_blast_res {
   my %args = @_;
   my $line1="";
   my $orig_blast_overview = "";
   my $alignments = "";
   my @acc = ();
   my $line_length = 0;
   my $acc_length = 0;
   my $subj_length = 0;
   my $new_overview = "";
   my ($score, $evalue);

   #catch original blast overview
   $args{blast_result} =~ m/Sequences producing significant alignments:\s+\(Bits\)\s+Value\s*N?\s*\n(.*?)\n\n\>/is;
   $orig_blast_overview = $1;
   unless (defined ($orig_blast_overview)) {
      ${$args{progress_bar}}->configure(-label => "Could not grab Blast overview from current aa sequence".
                                                  "\nNo significant Blast results were generated. Skipping");
      ${$args{main_window}}->update;sleep(1);
      return ("\n\nNo Blast Hits found\n");
   }
   if ($orig_blast_overview !~ /\w+/) {
      ${$args{progress_bar}}->configure(-label =>"No significant Blast results were generated. Skipping");
      ${$args{main_window}}->update;sleep(1);
      return ("\n\nNo Blast Hits found\n");
   }
   $orig_blast_overview = "\n".$orig_blast_overview;

   #catch aligments
   $args{blast_result} =~ m/\n\n(\>.*?)Lambda\s+K\s+H/s;
   $alignments = "\n".$1;
   unless ($alignments =~ /\w+/) {
      ${$args{progress_bar}}->configure(-label =>"Could not grab Blast alignments from current aa sequence");
      ${$args{main_window}}->update;sleep(1);
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in Blast analysis".
                  "\nCould not grab Blast alignments from results\n$args{blast_result}\n\n";
      close WRITE;
   }

   #test for correct format. There is a bug at NCBI nr that has the odd old '|' format appearing
   #to avoid a false recognition, the number of '|'s must be at least twice the number of entries.
   #using the alignment instead of overview to count the number of '>'
   my @alignment_entries  = ($alignments =~ m/\n\>/gs);
   my @vertical_seperator = ($alignments =~ m/\n\>\S+\|/gs);
   my $old_format         = 0; #assume new format
   if ((scalar @vertical_seperator) >= (scalar @alignment_entries)) {
       $old_format = 1;
   }

   #build array for accession numbers from blast overview
   if ($orig_blast_overview =~ m/\n[^\|]+?\|[^\|]*?\|/ && $old_format == 1) {
      @acc = ($orig_blast_overview =~ m/\n[^\|]+?\|([^\|]*?)\|/mg);
   } else {
      #check for new Blast format
      @acc = ($orig_blast_overview =~ m/\n(\S+)\s/mg); #this will take the new acc
   }
   @acc = splice(@acc,0,${$args{auto_ini_ref}}{blast_entry_number}); #restrict number of entries to set value in $key

   #built new Blast overview
   #find out maximum text length
   foreach (@acc) {
      next unless ($_ =~ m/\w+/); #removes entries without a proper acc number. Temporary fix until parser is re-written somehow
      my ($acc_temp, $temp);

      if ($alignments =~ m/^\s*\>\w+?\|.+\|/ && $old_format == 1) {
         $alignments =~ m/\n(\>\w+?\|$_\|)(.*?)\n\s+Score = /s;
         ($acc_temp, $temp) = ($1, $2);
      } else {
         $alignments =~ m/\n(\>$_)\s+(.*?)\n\s+Score = /s;
         ($acc_temp, $temp) = ($1, $2);
      }

      #restrict annotation to first set of entries
      if ($temp =~ m/\n\s+\w+\|/ && $old_format == 1) {
         $temp =~ s/\n\s+\w+\|.+Length\s*\=/   Length \=/s;
      } else {
         $temp =~ s/\n\s.+Length\s*\=/   Length \=/s;
      }

      unless ($temp =~ /\w+/) {
         ${$args{progress_bar}}->configure(-label =>"Could not grab overview line from current entry $_");
         ${$args{main_window}}->update;sleep(1);
         open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
         print WRITE "Error in Blast analysis".
                     "\nCould not grab overview line from current entry $_\nOriginal Bast results:\n$args{blast_result}\n\n";
         close WRITE;
      }
      #determine maximum length of entry line
      $temp =~ s/\n?\s+Length\s*=\s*.+$//;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $temp =~ s/RecName\://gs;
      $temp =~ s/Full\=//gs;
      if (length($temp) > $line_length) {$line_length = length($temp)};
      #determine max length of acc
      if (length($acc_temp) > $acc_length) {$acc_length = length($acc_temp)};
   }
   #generate new overview
   $orig_blast_overview .= "\n";
   foreach (@acc) {
      next unless ($_ =~ m/\w+/); #removes entries without a proper acc number. Temporary fix until parser is re-written somehow
      my ($acc_temp, $temp);
      #get new overview line
      if ($alignments =~ m/\>\w+?\|.+\|/ && $old_format == 1) {
         $alignments =~ m/\n(\>\w+?\|$_\|)(.*?)\n\s+Score = /s;
         ($acc_temp, $temp) = ($1, $2);
      } else {
         $alignments =~ m/\n(\>$_)\s+(.*?)\n\s+Score = /s;
         ($acc_temp, $temp) = ($1, $2);
      }

      #restrict annotation to first set of entries
      if ($temp =~ m/\n\s+\w+\|/ && $old_format == 1) {
         $temp =~ s/\n\s+\w+\|.+Length\s*\=/   Length \=/s;
      } else {
         $temp =~ s/\n\s.+Length\s*\=/   Length \=/s;
      }

      'reset' =~ m/reset/;
      $temp =~ s/\n?\s+(Length\s*=\s*.+)$//;
      $subj_length = $1;
      $subj_length =~ s/\s+//g;
      $temp =~ s/^\s*//;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $temp =~ s/RecName\://gs;
      $temp =~ s/Full\=//gs;
      'reset' =~ m/reset/;
      my $entry = '';
      if ($orig_blast_overview =~ m/(\w+\|$_\|/ && $old_format == 1) {
         $orig_blast_overview =~ m/(\w+\|$_\|[^\n]+?\n)/;
         $entry = $1;
         'reset' =~ m/reset/;
         $entry =~ m/.{68}\s+(\d[\d\.]+)\s+([\de\.\-]+).*\n/;
         ($score, $evalue) = ($1, $2);
      } else {
         $orig_blast_overview =~ m/($_.+?\n)/;
         $entry = $1;
         'reset' =~ m/reset/;
         $entry =~ m/([\d\.])+\s+([\d\.\-e]+)\s*$/;
         ($score, $evalue) = ($1, $2);
      }

      unless (defined $evalue) {
         print "\nBlast res: Could not grab Blast scores from overview for acc $_ in $args{header} \n";
         open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
         print WRITE "Error in Blast analysis".
                     "\nCould not grab Blast scores from overview for acc $_ in $args{header}\nOriginal Bast results:\n$args{blast_result}\n\n";
         close WRITE;
      }
      if ($evalue =~ /^e/) {$evalue = "1".$evalue};
      $new_overview .= $acc_temp." " x ($acc_length - length($acc_temp) + 1).$temp."~" x ($line_length - length($temp))." ".$subj_length." Score=".$score." Expect=".$evalue."\n";
   }

   #attach the alignment to new overview
   $line1 = "\nBlast overview\n\n".$new_overview."\n\n\nComplete list of Blast results\n\n\n".$alignments;
   return ($line1);

}

sub reformat_blast_threaded {
   my %args                = @_;
   my $line1               = '';
   my $orig_blast_overview = '';
   my $alignments          = '';
   my @acc                 = ();
   my $line_length         = 0;
   my $acc_length          = 0;
   my $subj_length         = 0;
   my $new_overview        = '';

   #catch original blast overview
   $args{blast_result} =~ m/Sequences producing significant alignments:\s+\(Bits\)\s+Value\s*N?\s*?\n(.*?)\n\n\>/is;
   $orig_blast_overview = $1;
   unless (defined ($orig_blast_overview)) {
      return ("\n\nNo Blast Hits found\n");
   }
   if ($orig_blast_overview !~ /\w+/) {
      return ("\n\nNo Blast Hits found\n");
   }
   $orig_blast_overview = "\n".$orig_blast_overview;

   #catch aligments
   $args{blast_result} =~ m/\n\n(\>.*?)Lambda\s+K\s+H/s;
   $alignments = "\n".$1;
   unless ($alignments =~ /\w+/) {
      print "\nCould not grab Blast alignments from current aa sequence\n";
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in Blast analysis".
                  "\nCould not grab Blast alignments from results\n$args{blast_result}\n\n";
      close WRITE;
   }

   #test for correct format. There is a bug at NCBI nr that has the odd old '|' format appearing
   #to avoid a false recognition, the number of '|'s must be at least twice the number of entries.
   #using the alignment instead of overview to count the number of '>'
   my @alignment_entries  = ($alignments =~ m/\n\>/gs);
   my @vertical_seperator = ($alignments =~ m/\n\>\S+\|/gs);
   my $old_format         = 0; #assume new format
   if ((scalar @vertical_seperator) >= (scalar @alignment_entries)) {
       $old_format = 1;
   }

   #build array for accession numbers from blast overview
   if ($orig_blast_overview =~ m/\n[^\|]+?\|[^\|]*?\|/ && $old_format == 1) {
      @acc = ($orig_blast_overview =~ m/\n[^\|]+?\|([^\|]*?)\|/mg);
   } else {
      #check for new Blast format
      @acc = ($orig_blast_overview =~ m/\n(\S+)\s/mg); #this will take the new acc
   }
   @acc = splice(@acc,0,${$args{auto_ini_ref}}{blast_entry_number}); #restrict number of entries to set value in $key

   #built new Blast overview
   #find out maximum text length
   #my $temp = 0;
   foreach (@acc) {
      #$temp++;
      #print "\n...$temp...".scalar @acc;
      next unless ($_ =~ m/\w+/); #removes entries without a proper acc number. Temporary fix until parser is re-written somehow
      my ($acc_temp, $temp);
      'reset' =~ m/reset/;
      if ($alignments =~ m/\n\>\w+?\|.+\|/ && $old_format == 1) {
         $alignments =~ m/\n(\>\w+?\|$_\|)(.*?)\n\s+Score = /s;
         ($acc_temp, $temp) = ($1, $2);
      } else {
         $alignments =~ m/\n(\>$_)\s+(.*?)\n\s+Score = /s;
         ($acc_temp, $temp) = ($1, $2);
      }

      #restrict annotation to first set of entries
      if ($temp =~ m/\n\s+\w+\|/ && $old_format == 1) {
         $temp =~ s/\n\s+\w+\|.+Length\s*\=/   Length \=/s;
      } else {
         $temp =~ s/\n\s.+Length\s*\=/   Length \=/s;
      }

      unless ($temp =~ /\w+/) {
         print "\nCould not grab overview line from current entry $_\n";
         open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
         print WRITE "Error in Blast analysis".
                     "\nCould not grab overview line from current entry $_\nOriginal Bast results:\n$args{blast_result}\n\n";
         close WRITE;
      }
      #determine maximum length of entry line
      $temp =~ s/\n?\s+Length\s*=\s*.+$//;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $temp =~ s/RecName\://gs;
      $temp =~ s/Full\=//gs;
      if (length($temp) > $line_length) {$line_length = length($temp)};
      #determine max length of acc
      if (length($acc_temp) > $acc_length) {$acc_length = length($acc_temp)};
   }

   #generate new overview
   $orig_blast_overview .= "\n";
   foreach (@acc) {
      next unless ($_ =~ m/\w+/); #removes entries without a proper acc number. Temporary fix until parser is re-written somehow
      my ($acc_temp, $temp, $score, $evalue);
      #get new overview line
      if ($alignments =~ m/\n\>\w+?\|.+\|/ && $old_format == 1) {
         $alignments =~ m/\n(\>\w+?\|$_\|)(.*?)\n\s+Score = /s;
         ($acc_temp, $temp) = ($1, $2);
      } else {
         $alignments =~ m/\n(\>$_)\s+(.*?)\n\s+Score = /s;
         ($acc_temp, $temp) = ($1, $2);
      }

      #restrict annotation to first set of entries
      if ($temp =~ m/\n\s+\w+\|/ && $old_format == 1) {
         $temp =~ s/\n\s+\w+\|.+Length\s*\=/   Length \=/s;
      } else {
         $temp =~ s/\n\s.+Length\s*\=/   Length \=/s;
      }

      'reset' =~ m/reset/;
      $temp =~ s/\n?\s+(Length\s*=\s*.+)$//;
      $subj_length = $1;
      $subj_length =~ s/\s+//g;
      $temp =~ s/^\s*//;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $temp =~ s/RecName\://gs;
      $temp =~ s/Full\=//gs;
      'reset' =~ m/reset/;
      #$orig_blast_overview =~ m/\w+\|$_\|[^\n]+?\s+(\d[\d\.]+)\s+([\de\.\-]+).*\n/;
      #($score, $evalue) = ($1, $2);
      my $entry = '';
      if ($orig_blast_overview =~ m/\w+\|.+?\|/ && $old_format == 1) {
         $orig_blast_overview =~ m/(\w+\|$_\|[^\n]+?\n)/;
         $entry = $1;
         'reset' =~ m/reset/;
         $entry =~ m/.{68}\s+(\d[\d\.]+)\s+([\de\.\-]+).*\n/;
         ($score, $evalue) = ($1, $2);
      } else {
         $orig_blast_overview =~ m/($_.+?\n)/;
         $entry = $1;
         'reset' =~ m/reset/;
         $entry =~ m/([\d\.])+\s+([\d\.\-e]+)\s*$/;
         ($score, $evalue) = ($1, $2);
      }

      unless (defined $evalue) {
         print "\nBlast threaded: Could not grab Blast scores from overview for acc $_ in $args{header} \n";
         open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
         print WRITE "Error in Blast analysis".
                     "\nCould not grab Blast scores from overview for acc $_ in $args{header}\nOriginal Bast results:\n$args{blast_result}\n\n";
         close WRITE;
      }
      if ($evalue =~ /^e/) {$evalue = "1".$evalue};
      $new_overview .= $acc_temp." " x ($acc_length - length($acc_temp) + 1).$temp."~" x ($line_length - length($temp))." ".$subj_length." Score=".$score." Expect=".$evalue."\n";
   }

   #attach the alignment to new overview
   $line1 = "\nBlast overview\n\n".$new_overview."\n\n\nComplete list of Blast results\n\n\n".$alignments;
   return ($line1);

}

sub reformat_blast_hash {
   my %args = @_;
   my ($line1, $orig_blast_overview, $alignments, $new_overview, $query, $score, $evalue);
   my @acc = ();
   my $line_length = 0;
   my $acc_length = 0;
   my $subj_length = 0;

   #get query name
   'reset' =~ m/reset/;
   $args{blast_result} =~ m/\nQuery\=\s(.*?)\nDatabase/s;
   $query = $1;
   unless (defined $query) {
      'reset' =~ m/reset/;
      $args{blast_result} =~ m/\nQuery\=\s(.*?)\n+Length/s;
      $query = $1;
   }
   unless (defined $query) {
      print "\nCould not grab Query name from current entry\n";
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in Meta Blast analysis".
         "\nCould not grab Query name from current entry $args{blast_result}\n\n";
      close WRITE;
      return;
   }
   $query =~ s/\(\d+\s+letters\)//s;
   $query =~ s/\n//gs;

   #return if not hits found
   if ($args{blast_result} =~ /\*\*\*\*\* No hits found \*\*\*\*\*/) {
      return ("\nBlast Flavour:\t$args{blast_flavour}\nBlast db:\t$args{blast_db}\nQuery name:\t$query\nBlast overview\n>No Blast Hits found\n\n\n");
   }

   #catch original blast overview
   'reset' =~ m/reset/;
   $args{blast_result} =~ m/Sequences producing significant alignments\:\s+\(bits\)\s+Value\s*N?\s*\n(.*?)\n\n\>/is;
   $orig_blast_overview = $1;

   unless (defined ($orig_blast_overview)) {
      return ("\nBlast Flavour:\t$args{blast_flavour}\nBlast db:\t$args{blast_db}\nQuery name:\t$query\nBlast overview\n>No Blast Hits found\n\n\n");
   }
   if ($orig_blast_overview !~ /\w+/) {
      return ("\nBlast Flavour:\t$args{blast_flavour}\nBlast db:\t$args{blast_db}\nQuery name:\t$query\nBlast overview\n>No Blast Hits found\n\n\n");
   }
   $orig_blast_overview = "\n".$orig_blast_overview;

   #catch aligments
   $args{blast_result} =~ m/\n\n\>(.+)/s;
   $alignments = '>'.$1;
   unless ($alignments =~ /\w+/) {
      print "\nCould not grab Blast alignments from current entry $args{blast_result}\n";
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in Blast analysis".
                  "\nCould not grab Blast alignments from results\n$args{blast_result}\n\n";
      close WRITE;

   }
   #remove trailing summary
   $alignments =~ s/\n\s+Database\:.*?\nLambda.+$//s;
   $alignments =~ s/\nLambda.+$//s;
   $alignments = "\n".$alignments;

   #test for correct format. There is a bug at NCBI nr that has the odd old '|' format appearing
   #to avoid a false recognition, the number of '|'s must be at least twice the number of entries.
   #using the alignment instead of overview to count the number of '>'
   my @alignment_entries  = ($alignments =~ m/\n\>/gs);
   my @vertical_seperator = ($alignments =~ m/\n\>\S+\|/gs);
   my $old_format         = 0; #assume new format
   if ((scalar @vertical_seperator) >= (scalar @alignment_entries)) {
       $old_format = 1;
   }

   #build array for accession numbers from blast overview
   if ($orig_blast_overview =~ m/\n[^\|]+?\|[^\|]*?\|/ && $old_format == 1) {
      @acc = ($orig_blast_overview =~ m/\n[^\|]+?\|\s*([^\|]+?)\|/mg);
   } else {
      #check for new Blast format
      @acc = ($orig_blast_overview =~ m/\n(\S+)\s/mg); #this will take the new acc
   }
   @acc = splice(@acc,0, $args{blast_display_limit},); #restrict number of entries to set value in the 'limit Blast display' option in the metagenome analysis

   #built new Blast overview
   #find out maximum text length
   foreach (@acc) {
      next unless ($_ =~ m/\w+/); #removes entries without a proper acc number. Temporary fix until parser is re-written somehow
      my ($acc_temp, $temp);
      if ($alignments =~ m/\>\w+?\|.+\|/ && $old_format == 1) {
         $alignments =~ m/(\>\w+?\|\s*$_\|)(.*?)\n\s+Score = /s;
         ($acc_temp, $temp) = ($1, $2);
      } else {
         $alignments =~ m/(\>$_)\s+(.*?)\n\s+Score = /s;
         ($acc_temp, $temp) = ($1, $2);
      }

      #restrict annotation to first set of entries
      if ($temp =~ m/\n\s+\w+\|/ && $old_format == 1) {
         $temp =~ s/\n\s+\w+\|.+Length\s*\=/   Length \=/s;
      } else {
         $temp =~ s/\n.+Length\s*\=/   Length \=/s;
      }

      unless ($temp =~ /\w+/) {
         print "\nCould not grab overview line from current entry $_\n";
         open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
         print WRITE "Error in Blast analysis".
                     "\nCould not grab overview line from current entry $_\n\n";
         close WRITE;
      }
      #determine maximum length of entry line
      $temp =~ s/\n?\s+Length\s*=\s*.+$//;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $temp =~ s/RecName\://gs;
      $temp =~ s/Full\=//gs;
      if (length($temp) > $line_length) {$line_length = length($temp)};
      #determine max length of acc
      if (length($acc_temp) > $acc_length) {$acc_length = length($acc_temp)};
      undef $temp;
      undef $acc_temp;
   }
   #generate new overview
   $orig_blast_overview .= "\n";
   foreach (@acc) {
      next unless ($_ =~ m/\w+/); #removes entries without a proper acc number. Temporary fix until parser is re-written somehow
      my ($acc_temp, $temp, $entry);
      #get new overview line
      if ($alignments =~ m/\>\w+?\|.+\|/ && $old_format == 1) {
         $alignments =~ m/(\>\w+?\|\s*$_\|)(.*?)\n\s+Score = /s;
         ($acc_temp, $temp) = ($1, $2);
      } else {
         $alignments =~ m/(\>$_)\s+(.*?)\n\s+Score = /s;
         ($acc_temp, $temp) = ($1, $2);
      }

      #restrict annotation to first set of entries
      if ($temp =~ m/\n\s+\w+\|/ && $old_format == 1) {
         $temp =~ s/\n\s+\w+\|.+Length\s*\=/   Length \=/s;
      } else {
         $temp =~ s/\n.+Length\s*\=/   Length \=/s;
      }

      'reset' =~ m/reset/;
      $temp =~ s/\n?\s+(Length\s*=\s*.+)$//;
      $subj_length = $1;
      $subj_length =~ s/\s+//g;
      $temp =~ s/^\s*//;
      $temp =~ s/\n//g;
      $temp =~ s/\s+/ /g;
      $temp =~ s/RecName\://gs;
      $temp =~ s/Full\=//gs;
      'reset' =~ m/reset/;
      $entry = '';
      if ($orig_blast_overview =~ m/\w+\|\s*$_\|/ && $old_format == 1) {
         $orig_blast_overview =~ m/(\w+\|\s*$_\|[^\n]+?)\n/;
         $entry = $1;
         'reset' =~ m/reset/;
         $entry =~ m/\s+(\d[\d\.]+)\s+([\de\.\-]+)\s*$/;
         ($score, $evalue) = ($1, $2);
      } else {
         $orig_blast_overview =~ m/($_.+?\n)/;
         $entry = $1;
         'reset' =~ m/reset/;
         $entry =~ m/([\d\.])+\s+([\d\.\-e]+)\s*$/;
         ($score, $evalue) = ($1, $2);
      }

      #$orig_blast_overview =~ m/\w+\|$_\|[^\n]+?\s+([\d\.]+)\s+([\de\.\-]+).*\n/;
      #$score, $evalue) = ($1, $2);
      unless (defined $evalue) {
         print "\nBlast hash: Could not grab Blast scores from overview for acc $_ in $entry.\n";
         open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
         print WRITE "Error in Blast analysis".
                     "\nCould not grab Blast scores from overview for acc $_ in $entry\n...$orig_blast_overview...\n\n";
         close WRITE;
      }
      if ($evalue =~ /^e/) {$evalue = "1".$evalue};
      $new_overview .= $acc_temp." " x ($acc_length - length($acc_temp) + 1).$temp."~" x ($line_length - length($temp))." ".$subj_length." Score=".$score." Expect=".$evalue."\n";
      undef $temp;
      undef $acc_temp;
      undef $entry;
      undef $subj_length;
   }

   #attach the alignment to new overview
   $line1 = "\nBlast Flavour:\t$args{blast_flavour}\nBlast db:\t$args{blast_db}\nQuery name:\t$query\nBlast overview\n\n".$new_overview."\n\n\n";#Complete list of Blast results\n\n\n".$alignments;
   return ($line1);

}

sub blast_file {
   my %args = (critica_select => 0,
               retry_counter  => 0,
               @_
              ); #limit the numbers of Blast retries in case of mising results
   my (@list, @Blast_tasks);
   my $progress = 0;

   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'Blasting',
                   label        => 'Blasting'
                  );
   &show_pbar_2;

   #before Blast analysis is started, better check if selected DB file is still accessible
   unless (-e ${$args{auto_ini_ref}}{full_blast_db} ||
           -e ${$args{auto_ini_ref}}{full_blast_db}.'.pal' ||
           -e ${$args{auto_ini_ref}}{full_blast_db}.'.nal' ||
           -e ${$args{auto_ini_ref}}{full_blast_db}.'.nhr' ||
           -e ${$args{auto_ini_ref}}{full_blast_db}.'.phr') {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Error acessing selected Blast database ${$args{auto_ini_ref}}{full_blast_db}.\nPlease select again",
                                                    -bitmap  => 'error',
                                                    -buttons => ['ok']);
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();

      &select_blast_db(main_window    => $args{main_window},
                       ini_ref        => $args{ini_ref},
                       auto_ini_ref   => $args{auto_ini_ref},
                       critica_select => $args{critica_select}
                      );
   }

   #map correct gene model entries for respective input file
   @list = grep /$args{input_file}\___/, @{$args{combined_orf_ref}};

   #remove existing entries if re-use results
   if (${$args{auto_ini_ref}}{reuse_results} == 1) {
      my (@jobs, @Blast_results, %exist_blast_file);
      &update_pbar_2(title        => 'Blasting',
                     label        => "Testing for existing results",
                     progress     => 1,
                    );

      #read existing results
      opendir SEQINPUT, ${$args{ini_ref}}{blast_results};
      @Blast_results = grep /^$args{input_file}\_\d+$/, readdir(SEQINPUT);
      closedir SEQINPUT;
      #create hash
      foreach (@Blast_results) {
         $exist_blast_file{$_} = '1';
      }
      #iterate over gene model
      foreach my $entry (@list) {
         my $ID  = '';
         'reset' =~ m/reset/;
         $entry  =~ m/$args{input_file}\___(\d+)___/;
         $ID     = $1;

         #skip Blast if result file exists and re-use is activated
         next if (exists $exist_blast_file{$args{input_file}.'_'.$ID} && -s ${$args{ini_ref}}{blast_results}.'/'.$args{input_file}.'_'.$ID > 0);
         #add if required
         push (@jobs, $entry);
      }
      @list = @jobs;
      undef @jobs;
      undef @Blast_results;
      undef %exist_blast_file;
   }

   #max number of fasta entries
   my $max_count = $#list + 1;

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

   #run threaded or CPU clustered Blast
   if (${$args{auto_ini_ref}}{blast_cluster} == 1) {
      #iterate over array
      foreach my $entry (@list) {
         my ($ID, $left_bd, $right_bd, $orientation);
         ${$args{progress_bar}}->configure(-label=>" ");
         ${$args{main_window}}->update;
         #update progress bar
         &update_pbar_2(title        => 'Blasting with clustered Blast',
                        label        => "Blasting $args{input_file}, $progress of $max_count",
                        progress     => ($progress / $max_count) * 100,
                       );
         $progress++;

         $ID = $left_bd = $right_bd = $orientation = '';
         $entry =~ m/$args{input_file}\___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
         $ID              = $1;
         $left_bd         = $2;
         $right_bd        = $3;
         $orientation     = $4;
         my $query_seq_nt = '';

         unless (defined $orientation && $orientation =~ /(sense|antisense)/) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Error parsing gene model $entry for $args{input_file}",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0); #return to main
         }

         #extract nt sequence from ORF
         $query_seq_nt = substr(${$seq_ref}, ($left_bd - 1), ($right_bd - $left_bd + 1));

         #translate nt sequence to aa if selected
         my $query_seq_aa = '';
         if (${$args{auto_ini_ref}}{blast_type} =~ /(blastp|gappedblastp|psi-blast)/i) {
            $query_seq_aa = nt2aa(main_window  => $args{main_window},
                                  progress_bar => $args{progress_bar},
                                  auto_ini_ref => $args{auto_ini_ref},
                                  ini_ref      => $args{ini_ref},
                                  module       => 'Blast',
                                  orientation  => $orientation,
                                  sequence     => $query_seq_nt,
                                  filename     => $args{input_file},
                                  left_bd      => $left_bd,
                                  right_bd     => $right_bd
                                 );
            return (0) if ($query_seq_aa eq '0');

            #save temp aa file for PSI blast
            if (${$args{auto_ini_ref}}{blast_type} =~ /psi-blast/i) {
               open WRITE, "+>${$args{ini_ref}}{blast_executables}\/temp.aa";
               print WRITE "\>$entry\n$query_seq_aa";
               close WRITE;
            }
         }

         my ($blastresult) = &blast_seq(main_window  => $args{main_window},
                                        progress_bar => $args{progress_bar},
                                        auto_ini_ref => $args{auto_ini_ref},
                                        ini_ref      => $args{ini_ref},
                                        sequence_nt  => $query_seq_nt,
                                        sequence_aa  => $query_seq_aa,
                                        left_bd      => $left_bd,
                                        right_bd     => $right_bd,
                                        orientation  => $orientation,
                                        input_file   => $args{input_file},
                                        orf_id       => $ID
                                       );
         if ($blastresult eq '0') {
            ${$args{main_window}}->messageBox(-title => 'Error',
                                              -message => "Error blasting ORF $ID in sequence $args{input_file}",
                                              -type => 'OK',
                                              -icon => 'info');
            &hide_pbar_2;
            return (0); #return to main
         }
         #write Blast results
         open WRITE, "+>${$args{ini_ref}}{blast_results}\/$args{input_file}\_$ID";
         print WRITE $blastresult;
         close WRITE;
      }
   } elsif (${$args{auto_ini_ref}}{blast_threaded} == 1) {
      #iterate over array
      foreach my $entry (@list) {
         my ($ID, $left_bd, $right_bd, $orientation);
         #update progress bar
         $progress++;
         if (($progress % 100) == 0) {
            &update_pbar_2(title        => 'Blast analysis with multi-threaded Blast',
                           label        => "Reading gene model for $args{input_file}",
                           progress     => ($progress / $max_count) * 100,
                          );
            ${$args{main_window}}->update;
         }

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
                               module       => 'Blast',
                               orientation  => $orientation,
                               sequence     => $query_seq_nt,
                               filename     => $args{input_file},
                               left_bd      => $left_bd,
                               right_bd     => $right_bd
                              );
         #create array for new Blast tasks
         push (@Blast_tasks, $args{input_file}.'___'.$ID.'___'.$left_bd.'___'.$right_bd.'___'.$orientation.'___'.$query_seq_aa.'___'.$query_seq_nt);
      }

      #start multithreading Blast
      for (my $i=1; $i<=$#Blast_tasks+2; $i = $i+${$args{auto_ini_ref}}{CPU}) {
         my $count = $i + ${$args{auto_ini_ref}}{CPU};
         if ($count > $#Blast_tasks+2) {$count = $#Blast_tasks+2};

         &update_pbar_2(title    => 'Blast analysis',
                        label    => "Blast analysis for $args{input_file}, ". ($count - 1) ." of $max_count",
                        progress => ($count / ($#Blast_tasks + 2)) * 100,
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
               $Blast_tasks[$j-1] =~ m/$args{input_file}___(\d+)___(\d+)___(\d+)___(.+)___(.+)___(.+)/;
               $ORF_number  = $1;
               $left_bd     = $2;
               $right_bd    = $3;
               $orientation = $4;
               $aa_sequence = $5;
               $nt_sequence = $6;

               my ($blastresult) = &blast_threaded(auto_ini_ref => $args{auto_ini_ref},
                                                   ini_ref      => $args{ini_ref},
                                                   sequence_nt  => $nt_sequence,
                                                   sequence_aa  => $aa_sequence,
                                                   left_bd      => $left_bd,
                                                   right_bd     => $right_bd,
                                                   orientation  => $orientation,
                                                   input_file   => $args{input_file},
                                                   orf_id       => $ORF_number
                                                  );
               if ($blastresult eq '0') {
                  print "Error in Blast analysis. Error blasting ORF $ORF_number in sequence $args{input_file} with Blast database";
                  open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
                  print WRITE "\nError in Blast analysis. Error blasting ORF $ORF_number in sequence $args{input_file} with Blast database\n\n";
                  close WRITE;
                  #exit safely
                  CORE::exit();
               }
               #write Blast results
               open WRITE, "+>${$args{ini_ref}}{blast_results}\/$args{input_file}"."_$ORF_number";
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

   #are really all Blast results done? better check
   {
      #my (@jobs, @Blast_results, %exist_blast_file);
      my @Blast_results = ();
      &update_pbar_2(title        => 'Blasting',
                     label        => "Testing for complete results",
                     progress     => 1,
                    );

      #read existing results
      opendir SEQINPUT, ${$args{ini_ref}}{blast_results};
      @Blast_results = grep /^$args{input_file}\_\d+$/, readdir(SEQINPUT);
      closedir SEQINPUT;
      #map correct gene model entries for respective input file
      @list = grep /$args{input_file}\___/, @{$args{combined_orf_ref}};

      if ($#Blast_results ne $#list && $args{retry_counter} <= 3) {
         $args{retry_counter}++;
         &hide_pbar_2;
         #print "\nIn $args{input_file}\: number of gene model ORFs $#list\, number of found Blast hits: $#Blast_results\n";
         undef @list;
         undef @Blast_results;
         #undef %exist_blast_file;
         &blast_file(main_window      => $args{main_window},
                     progress_bar     => $args{progress_bar},
                     auto_ini_ref     => $args{auto_ini_ref},
                     ini_ref          => $args{ini_ref},
                     directory        => ${$args{ini_ref}}{input_files},
                     input_file       => $args{input_file},
                     combined_orf_ref => $args{combined_orf_ref},
                     retry_counter    => $args{retry_counter}
                    );
      }
      undef @list;
      undef @Blast_results;
      #undef %exist_blast_file;
   }

   ${$args{progress_bar}}->configure(-label=>"Finished Blast analysis");
   ${$args{main_window}}->update;
   &hide_pbar_2;
   ${$args{main_window}}->update;
   return (1);
}

sub blast_hash {
   #used for metagenome analysis, requires header-seq hash as input
   my %args = @_;
   my (@list, @Blast_tasks, $counter, %seen, %catch, $max_count, @seen);
   my $progress = 0;
   $counter = 1;

   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'Blasting hash table',
                   label        => 'Blasting'
                  );
   &show_pbar_2;

   #remove existing entries if re-use results
   if (${$args{auto_ini_ref}}{reuse_results} == 1) {
      my (@temp_results); # @entries);
      &update_pbar_2(title        => 'Blasting',
                     label        => "Testing for existing results...Reading",
                     progress     => 1,
                    );

      #read existing results
      opendir SEQINPUT, ${$args{ini_ref}}{blast_results};
      @temp_results = grep /^metagenome\.$args{blast_flavour}(\_\d+)?$/, readdir(SEQINPUT);
      my $counter = 1;
      foreach my $entry (@temp_results) {
         #unlink if filesize 0
         if (-s ${$args{ini_ref}}{blast_results}.'/'.$entry == 0) {
            unlink ${$args{ini_ref}}{blast_results}.'/'.$entry;
            next;
         }

         #read results file line by line insted
         open READ, '<'.${$args{ini_ref}}{blast_results}.'/'.$entry;
         my $count = 1;
         while (<READ>) {
            if ($count % 100000 == 0) {
               &update_pbar_2(title        => "Blasting hash table",
                              label        => "Reading exisiting results entry $counter",
                              progress     => ($count / 10000000) * 100,
                             );
               if ($count >= 10000000 ) {$count = 1};
            }
            $count++;

            if ($_ =~ m/Query name\:/) {
               'reset' =~ m/reset/;
               m/Query name\:\s+(.+)\n/s;
               my $header = $1;
               $header =~ s/\s+$//;
               $header =~ s/^\s+//;
               $seen{$counter} = $header;
               $catch{$header} = 1;
               push (@seen, $header);
               $counter++;
            }
         }
         close READ;
      }
      undef @temp_results;
   }

   &update_pbar_2(title        => "Blasting hash table",
                  label        => "Finished exisiting result entries",
                  progress     => 100,
                 );

   #cluster reads into new array
   my $combined;
   $counter = 1;
   my $job_count = 1;
   $max_count = (scalar (keys %{ $args{seq_hash} })) + 2;
   foreach my $ID (keys %{ $args{seq_hash} }) {
      if ($job_count % 100 == 0) {
         &update_pbar_2(title        => 'Blasting',
                        label        => "Building Blast job list",
                        progress     => ($job_count / $max_count) * 100,
                       );
      }
      $job_count++;

      #skip if result already exists
      next if (${$args{auto_ini_ref}}{reuse_results} == 1 && exists $catch{$args{seq_hash}{$ID}->{header}}); #{

      #merge individual reads into cluster
      if ($counter < $args{cluster_seq}) {
         $combined .= '>'.$args{seq_hash}{$ID}->{header}."\n".$args{seq_hash}{$ID}->{seq}."\n\n";
         $counter++;
         next;
      }
      #reset counter for next cluster
      $counter = 1;

      #push new cluster into array
      if ($combined =~ /\w+/) {
         push (@Blast_tasks, $combined);
      }

      #setup combined for next run
      $combined = '>'.$args{seq_hash}{$ID}->{header}."\n".$args{seq_hash}{$ID}->{seq}."\n\n";
   }
   undef $job_count;
   undef %catch;
   undef @seen;
   undef %seen;

   #get last entry
   if ($combined =~ /\S+/) {
      push (@Blast_tasks, $combined);
      undef $combined;
   }

   &update_pbar_2(title        => "Blasting hash table",
                  label        => "Finished preparing Blast job list",
                  progress     => 100,
                 );

   #max number of fasta entries
   $max_count = $#Blast_tasks + 2;

   #run threaded or CPU clustered BLAST
   if ($args{cluster_blast} == 1) {
      #iterate over hash
      foreach my $entry (@Blast_tasks) {
         my ($ID, $left_bd, $right_bd, $orientation);
         my ($short_header, $entry_count);

         #define short header for progress bar
         $entry =~ m/^>(\S+?)\s/;
         $short_header = $1;
         unless (defined $short_header) {$short_header = ''};
         $short_header = substr($short_header, 0, 25);

         #determine how many entries are in current entry
         $entry_count = ($entry =~ tr/\>//);

         ${$args{progress_bar}}->configure(-label=>" ");
         ${$args{main_window}}->update;
         #update progress bar
         $progress++;
         &update_pbar_2(title        => 'Blasting',
                        label        => "Blasting $short_header, $progress of $max_count",
                        progress     => ($progress / $max_count) * 100,
                       );

         my ($blastresult) = &blast_hash_seq(auto_ini_ref        => $args{auto_ini_ref},
                                             ini_ref             => $args{ini_ref},
                                             entry_ref           => \$entry,
                                             entry_count         => $entry_count,
                                             blast_flavour       => $args{blast_flavour},
                                             blast_db            => $args{blast_db},
                                             blast_display_limit => $args{blast_display_limit},
                                            );
         if ($blastresult eq '0') {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Error blasting $entry in sequence $args{input_file}",
                                              -type    => 'OK',
                                              -icon    => 'info');
            &hide_pbar_2;
            return (0); #return to main
         }
         #write Blast results
         open WRITE, ">>".${$args{ini_ref}}{blast_results}.'/metagenome.'.$args{blast_flavour}.'_1';
         print WRITE $blastresult."\n\n";
         close WRITE;
      }
   } elsif (${$args{auto_ini_ref}}{blast_threaded} == 1) {
      my $max_count = $#Blast_tasks + 2;
      #start multithreading Blast
      for (my $i=1; $i<=$#Blast_tasks+2; $i = $i+${$args{auto_ini_ref}}{CPU}) {
         my ($short_header);
         my $count = $i + ${$args{auto_ini_ref}}{CPU};
         if ($count > $#Blast_tasks+2) {$count = $#Blast_tasks+2};

         #define short header for progress bar
         $Blast_tasks[$i - 1] =~ m/^>(\S+?)\s/;
         $short_header = $1;
         unless (defined $short_header) {$short_header = ''};
         $short_header = substr($short_header, 0, 25);

         &update_pbar_2(title    => 'Meta Blast analysis',
                        label    => "Blast analysis for $short_header, $count of $max_count",
                        progress => ($count / $max_count) * 100,
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
               #determine how many entries are in current entry
               my $entry_count = ($Blast_tasks[$j - 1] =~ tr/\>//);
               # child
               my ($blastresult) = &blast_hash_seq(auto_ini_ref        => $args{auto_ini_ref},
                                                   ini_ref             => $args{ini_ref},
                                                   entry_ref           => \$Blast_tasks[$j - 1],
                                                   entry_count         => $entry_count,
                                                   blast_flavour       => $args{blast_flavour},
                                                   blast_db            => $args{blast_db},
                                                   blast_display_limit => $args{blast_display_limit},
                                                   threaded            => '1'
                                                  );
               if ($blastresult eq '0') {
                  #define short header for progress bar
                  $Blast_tasks[$i - 1] =~ m/^>(\S+?)\s/;
                  my $indicator = $1;
                  unless (defined $indicator) {$indicator = ''};
                  print "Error in Blast analysis. Error blasting entry $indicator with Blast database $args{blast_db}";
                  open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
                  print WRITE "\nError in Meta Blast analysis. \nError blasting entry $indicator with Blast database $args{blast_db}\n\n";
                  close WRITE;
                  #exit safely
                  CORE::exit();
               }
               #write Blast results
               if (${$args{auto_ini_ref}}{reuse_results} == 1) {
                  open WRITE, ">>".${$args{ini_ref}}{blast_results}.'/metagenome.'.$args{blast_flavour}.'_'.$j;
               } else {
                  open WRITE, "+>".${$args{ini_ref}}{blast_results}.'/metagenome.'.$args{blast_flavour}.'_'.$j;
               }
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

   #merge all temp results
   opendir SEQINPUT, ${$args{ini_ref}}{blast_results};
   my $temp;
   my @temp_results = grep /^metagenome\.$args{blast_flavour}(\_\d+)?$/, readdir(SEQINPUT);
   $counter         = 1;
   $max_count       = $#temp_results + 2;
   closedir SEQINPUT;

   if (${$args{auto_ini_ref}}{reuse_results} == 1) {
      open WRITE, ">>".${$args{ini_ref}}{blast_results}.'/metagenome.'.$args{blast_flavour};
   } else {
      unlink ${$args{ini_ref}}{blast_results}.'/metagenome.'.$args{blast_flavour};
      open WRITE, ">>".${$args{ini_ref}}{blast_results}.'/metagenome.'.$args{blast_flavour};
   }

   foreach my $entry (@temp_results) {
      next if ($entry =~ m/^metagenome\.$args{blast_flavour}$/); #skip original metagenome summary file
      &update_pbar_2(title    => 'Meta Blast analysis',
                     label    => "Merging individual result file $entry, $counter of $max_count",
                     progress => ($counter / $max_count) * 100,
                    );
      $counter++;

      open READ, '<'.${$args{ini_ref}}{blast_results}.'/'.$entry;
      while (<READ>) {
         print WRITE $_;
      }
      close READ;

      if ($entry =~ /_\d+$/) {
         unlink ${$args{ini_ref}}{blast_results}.'/'.$entry;
      }
   }
   close WRITE;
   undef @temp_results;

   &hide_pbar_2;
   ${$args{main_window}}->update;
   return;
}

sub blast_ig {
   my %args = @_;
   my ($blastresult, $aa_seq, $aa_ref, $header, $nt_seq);

   #parse header and sequence, EXPECTED TO BE IN FASTA
   $args{sequence_nt} =~ m/>([^\n]*?)\n(.+)/s;
   $header = $1;
   $nt_seq = $2;
   $nt_seq =~ s/[^actgnACTGN]//gs;

   #translate to aa sequence
   ($aa_seq) = &nt2aa_unchecked(
                                sequence     => $nt_seq
                               );

   #format into fasta
   ($aa_ref) = &seq2fasta(header   => $header,
                          sequence => $aa_seq
                         );

   #get current working directory
   my $curdir = "";
   $curdir = getcwd();

   #change to Blast directory
   if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
      chdir ${$args{ini_ref}}{blast_executables};
      $blastresult = "";
      if (open RESULT, "-|") { # original process
         local $/;
         $blastresult = <RESULT>;
      } else { # child
         if (open STDIN, "-|") { # child
            exec "./blastpgp -d ${$args{auto_ini_ref}}{intergenicblastdb} -b 5"; #child
            die "Cannot exec: $!";
         } else { # grandchild
            print $$aa_ref;
            CORE::exit;
         }
      }
   } elsif (${$args{auto_ini_ref}}{blast_plus} == 1) {
      chdir ${$args{ini_ref}}{blast_plus_executables};
      $blastresult = "";
      if (open RESULT, "-|") { # original process
         local $/;
         $blastresult = <RESULT>;
      } else { # child
         if (open STDIN, "-|") { # child
            exec "./blastp -db ${$args{auto_ini_ref}}{intergenicblastdb} -num_descriptions 5 -num_alignments 5"; #child
            die "Cannot exec: $!";
         } else { # grandchild
            print $$aa_ref;
            CORE::exit;
         }
      }
   }
   #change back to previous working dir
   chdir $curdir;

   #valid result?
   unless (defined $blastresult && $blastresult =~ /\w+/) {
      print "Error in Blast analysis using database ${$args{auto_ini_ref}}{intergenicblastdb}\n";
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in Blast analysis using database ${$args{auto_ini_ref}}{intergenicblastdb}".
                  "\nCould not parse result: $blastresult\n\n";
      close WRITE;
      #print "\nError in Blast analysis using database ${$args{auto_ini_ref}}{intergenicblastdb} and Sequence: $nt_seq...\n...$aa_seq ...";
      $blastresult = 0;
      return (\$blastresult);
   };
   if ($blastresult =~ m/\*\*\*\*\* No hits found \*\*\*\*\*/) {
      $blastresult = 0;
      return (\$blastresult, $$aa_ref);
   }

   return (\$blastresult, $$aa_ref);
}

sub blast_ig_clustered {
   my %args = @_;
   my ($blastresult, $aa_seq, $aa_ref, $header, $nt_seq);

   #parse header and sequence, EXPECTED TO BE IN FASTA
   $args{sequence_nt} =~ m/>([^\n]*?)\n(.+)/s;
   $header = $1;
   $nt_seq = $2;

   #translate to aa sequence
   ($aa_seq) = &nt2aa_unchecked(
                                sequence     => $nt_seq
                               );

   #format into fasta
   ($aa_ref) = &seq2fasta(header   => $header,
                          sequence => $aa_seq
                         );
   #get current working directory
   my $curdir = "";
   $curdir = getcwd();

   #change to Blast directory
   if (${$args{auto_ini_ref}}{legacy_blast} == 1) {
      chdir ${$args{ini_ref}}{blast_executables};
      $blastresult = "";
      if (open RESULT, "-|") { # original process
         local $/;
         $blastresult = <RESULT>;
      } else { # child
         if (open STDIN, "-|") { # child
            exec "./blastpgp -a ${$args{auto_ini_ref}}{IG_CPU} -d ${$args{auto_ini_ref}}{intergenicblastdb} -b 5"; #child
            die "Cannot exec: $!";
         } else { # grandchild
            print $$aa_ref;
            CORE::exit;
         }
      }
   } elsif (${$args{auto_ini_ref}}{blast_plus} == 1) {
      chdir ${$args{ini_ref}}{blast_plus_executables};
      $blastresult = "";
      if (open RESULT, "-|") { # original process
         local $/;
         $blastresult = <RESULT>;
      } else { # child
         if (open STDIN, "-|") { # child
            exec "./blastp -num_threads ${$args{auto_ini_ref}}{IG_CPU} -db ${$args{auto_ini_ref}}{intergenicblastdb} -num_descriptions 5 -num_alignments 5"; #child
            die "Cannot exec: $!";
         } else { # grandchild
            print $$aa_ref;
            CORE::exit;
         }
      }
   }
   #change back to previous working dir
   chdir $curdir;

   #valid result?
   unless (defined $blastresult && $blastresult =~ /\w+/) {
      print "Error in Blast analysis using database ${$args{auto_ini_ref}}{intergenicblastdb}\n";
      open WRITE, ">>${$args{auto_ini_ref}}{work_dir}\/Error.log";
      print WRITE "Error in Blast analysis using database ${$args{auto_ini_ref}}{intergenicblastdb}".
                  "\nCould not parse result: $blastresult\n\n";
      close WRITE;
      #print "\nError in Blast analysis using database ${$args{auto_ini_ref}}{intergenicblastdb} and Sequence: $nt_seq...\n...$aa_seq ...";
      $blastresult = 0;
      return (\$blastresult);
   };
   if ($blastresult =~ m/\*\*\*\*\* No hits found \*\*\*\*\*/) {
      $blastresult = 0;
      return (\$blastresult, $$aa_ref);
   }

   return (\$blastresult, $$aa_ref);
}
1;
                