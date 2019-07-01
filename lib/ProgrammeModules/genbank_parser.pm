#!/opt/ActivePerl-5.8/bin/perl

#Genbank parser
#input arguments: main_window, progress_bar, gb_file
#returns : $header, $source, $core, $DNAseq, \@start, \@feature_list, \%genbank, $counter

#compile new GB file
#input arguments: main_window, progress_bar, ini_ref, auto_ini_ref, $header, $source, $DNAseq, \@start, \@feature_list, \%genbank
#writes new_gb_file to results folder;

package ProgrammeModules::genbank_parser;
use strict;
use Basics::Merge         qw( merge );
use Basics::progress_bar  qw(:DEFAULT);
use initialise::read_me   qw(:DEFAULT);
use initialise::gb_header qw(:DEFAULT);

use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&gb_parser &gb_parser_cmd &compile_gb $header $source $core $DNAseq %genbank @start @feature_list $counter );
use vars qw($header $source $core $DNAseq %genbank @start @feature_list $counter);

#local variables
my (%args, %seen, @features, @tmp, $GB_file, $temp, $i);

sub gb_parser {
   my %args = (gb_file           => '',
               process_gb        => 'extract',
               quiet             => 0,
               DNA_cleanup       => 1,
               @_
              );

   #slurp Genbank file into memory
   $GB_file = "";
   $header = "";
   $source = "";
   $core = "";
   $DNAseq = "";
   @features = ();
   @feature_list = ();
   @tmp = ();
   @start = ();
   %genbank = ();

   {
      local( $/, *GB_ENTRY ) ;
      open( GB_ENTRY, "<$args{gb_file}" ) or do {
         if ($args{quiet} == 0) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Cannot open submitted file $args{gb_file}",
                                                          -bitmap  => 'error',
                                                          -buttons => ['ok']);
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing Genbank files".
                           "\nCannot open submitted file $args{gb_file}\n\n";
            close ERRORLOG;
         }
         return (0, "Submitted Genbank file $args{gb_file} does not exist");
      };
      $GB_file = <GB_ENTRY>;
      close GB_ENTRY;
   }
   #remove empty space from beginning of file and re-write file
   if ($GB_file =~ m/^\s*/) {
      $GB_file =~ s/^\s*//;
      open  WRITE, "+>$args{gb_file}";
      print WRITE $GB_file;
      close WRITE;
   }

   #is it a genbank file?
   if ($GB_file !~ /^locus/i) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "The accessed file $args{gb_file} is not a valid Genbank file",
                                                    -bitmap  => 'error',
                                                    -buttons => ['ok']);
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error parsing Genbank files".
                     "\nThe accessed file $args{gb_file} is not a valid Genbank file\n\n";
      close ERRORLOG;
      return (0, "The accessed file $args{gb_file} is not a valid Genbank file");
   }

   ${$args{progress_bar}}->configure(-label=>"Reading Genbank file into memory");
   ${$args{main_window}}->update;
   #catch header
   'reset' =~ /reset/;
   $GB_file =~ m/^(.*? Location\/Qualifiers\s*\n)/s;
   $header = $1;

   #catch core
   'reset' =~ /reset/;
   $GB_file =~ m/Location\/Qualifiers\s*\n(     [a-zA-Z_]+\s+.*?\n(BASE COUNT|ORIGIN))/s;
   $core = $1;

   #catch source
   'reset' =~ /reset/;
   $core =~ m/((     source.*?\n)+)(     \S+|BASE COUNT|ORIGIN)/s;
   $source = $1;
   $core =~ s/\Q$source\E//s;
   $core =~ s/(BASE COUNT|ORIGIN)//gs;

   #catch DNA seq
   'reset' =~ /reset/;
   $GB_file =~ m/\n((BASE COUNT|ORIGIN).*)/s;
   $DNAseq = $1;
   undef $GB_file; #give back memory
   unless (defined $source) {$source = ''};

   #clean up DNA seq
   if ($args{DNA_cleanup} == 1) {
      $DNAseq =~ s/ORIGIN//s;
      $DNAseq =~ s/BASE COUNT[^\n]*?\n//is;
      $DNAseq =~ s/\d//gs;
      $DNAseq =~ s/\s//gs;
      $DNAseq =~ s/\\//g;
      $DNAseq =~ s/\///g;
      $DNAseq =~ s/[^a-z]+//igs;
   }

   ${$args{progress_bar}}->configure(-label=>"Reading Genbank file into memory...OK");
   ${$args{main_window}}->update;

   #if whole sequence is requested, exit here
   if ($args{process_gb} eq 'sequence') {
      #return references
      return ($header, $source, $core, $DNAseq, undef, undef, undef, undef);
   }

   #process into single features
   ${$args{progress_bar}}->configure(-label=>"Processing into single keys");
   ${$args{main_window}}->update;

   @tmp = split(/\n     (\S)/s, $core);
   #cleaning up;
   $tmp[0] =~ s/^\s+//;
   push (@features, $tmp[0]);

   #reset counter
   $counter = 0;
   foreach my $entry (@tmp) {
      if ($entry =~ /^\S$/) {
         my $temp = '';
         $temp = $tmp[$counter].$tmp[$counter+1];
         #generating valid entries list
         push (@features, $temp);
      }
      $counter++;
   }

   #release memory from temp array
   undef @tmp;

   #generate hash for entries; sort them later by start position
   my $gene_counter    = 100000;
   my $CDS_counter     = 100000;
   my %feature_counter = ();
   my %duplicate_gene  = ();
   my %duplicate_CDS   = ();
   foreach my $entry (@features) {
      next if ($entry !~ /\S+/);
      'reset' =~ m/reset/;
      my $feature = "";
      my $left_bd = "";
      my $right_bd = "";
      my $orientation;

      $entry =~ m/^\s*(\S+)\s*(complement)?\D*?(\d+)\D*?(\d+)/;
      $feature = $1;
      $orientation = $2;
      $left_bd = $3;
      $right_bd = $4;

      #skip source if escaped until here
      next if ($feature eq 'source');

      #modify right boundary for special keys with only left boundary entry
      unless (defined ($right_bd)) {
         #retry matching
         'reset' =~ m/reset/;
         $entry =~ m/^\s*(\S+)\s*(complement)?\D*?(\d+)/;
         $feature = $1;
         $orientation = $2;
         $left_bd = $3;
         #set boundaries equal
         $right_bd = $left_bd;
      }

      if (defined $orientation) {
         $orientation = 'antisense';
      } else {
         $orientation = 'sense';
      }
      unless (defined($feature)) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not grab values from $_ in file \n$args{gb_file}\.",
                                                       -bitmap  => 'error',
                                                       -buttons => ['ok']);
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing Genbank file $args{gb_file}".
                        "\nCould not grab values from entry $_\n\n";
         close ERRORLOG;
         next;
      }
      unless (defined($left_bd)) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not grab values from $_ in file \n$args{gb_file}\.",
                                                       -bitmap  => 'error',
                                                       -buttons => ['ok']);
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing Genbank file $args{gb_file}".
                        "\nCould not grab values from entry $_\n\n";
         close ERRORLOG;
         next;
      }

      #modify feature if necessary to make compatible with key-parsing
      if ($feature =~ m/_/) {
         $feature =~ s/_/~/g;
      }

      #do we have an ORF counter already? if other feature, increase ID count
      if ($feature =~ /^gene$/) {
         my $ORF_number = '';
         'reset' =~ m/reset/;
         $entry =~ m/\s+\/gene\=\".*?_\s*([\d\s]+)\s*\"/s;
         $ORF_number = $1;

         unless (defined ($ORF_number) && $ORF_number =~ /^\d+$/) {
            'reset' =~ m/reset/;
            $entry =~ m/\s+\/locus_tag\=.*?(\d+)\s*\"/s;
            $ORF_number = $1;
         }
         #no ORF number?
         unless (defined ($ORF_number) && $ORF_number =~ /^\d+$/) {
            $gene_counter++;
            $ORF_number = $gene_counter;
         }
         $ORF_number =~ s/\s+//gs;

         #remove leading '0's
         $ORF_number =~ s/^0+//;

         #nothing left? There have been cases where counts start with '0'
         unless ($ORF_number =~ m/\d+/) {
            $ORF_number = '0';
         }
         #does ORF number already exists?
         while (exists $duplicate_gene{$ORF_number}) {
            $ORF_number = $gene_counter;
            $gene_counter++;
         }
         $duplicate_gene{$ORF_number} = 1;

         #set gene counter
         $feature_counter{$feature} = $ORF_number;

      } elsif ($feature =~ /^CDS$/) {
         my $ORF_number = '';
         'reset' =~ m/reset/;
         $entry =~ m/\s+\/gene\=\".*?_\s*([\d\s]+)\s*\"/s;
         $ORF_number = $1;

         unless (defined ($ORF_number) && $ORF_number =~ /^\d+$/) {
            'reset' =~ m/reset/;
            $entry =~ m/\s+\/locus_tag\=.*?(\d+)\s*\"/s;
            $ORF_number = $1;
         }
         #no ORF number?
         unless (defined ($ORF_number) && $ORF_number =~ /^\d+$/) {
            $CDS_counter++;
            $ORF_number = $CDS_counter;
         }
         $ORF_number =~ s/\s+//gs;

         #remove leading '0's
         $ORF_number =~ s/^0+//;

         #nothing left? There have been cases where counts start with '0'
         unless ($ORF_number =~ m/\d+/) {
            $ORF_number = '0';
         }

         #does ORF number already exists?
         while (exists $duplicate_CDS{$ORF_number}) {
            $ORF_number = $CDS_counter;
            $CDS_counter++;
         }
         $duplicate_CDS{$ORF_number} = 1;

         #set gene/CDS counter
         $feature_counter{$feature} = $ORF_number;
      } else {
         $feature_counter{$feature}++;
      }

      $genbank{$left_bd.'_'.$right_bd.'_'.$feature.'_'.$feature_counter{$feature}} = $entry;
      push (@feature_list, $left_bd.'_'.$right_bd.'_'.$feature.'_'.$feature_counter{$feature});
   }

   #release memory from feature array
   @features = ();
   undef %duplicate_gene;
   undef %duplicate_CDS;

   #return references
   return ($header, $source, $core, $DNAseq, \@start, \@feature_list, \%genbank, $counter);

}


sub gb_parser_cmd { #duplicate to above, without graphical feedback
   my %args = (gb_file => "",
               process_gb => 'extract',
               @_);

   #slurp Genbank file into memory
   $GB_file = "";
   $header = "";
   $source = "";
   $core = "";
   $DNAseq = "";
   @features = ();
   @feature_list = ();
   @tmp = ();
   $counter = 0;
   @start = ();
   %genbank = ();

   {
      local( $/, *GB_ENTRY ) ;
      open( GB_ENTRY, "<$args{gb_file}" ) or do {
         print "Cannot open submitted file $args{gb_file}";
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error opening Genbank file $args{gb_file}".
                        "\nFile does not exist_\n\n";
         close ERRORLOG;
         return (0, "Submitted Genbank file $args{gb_file} does not exist");
      };
      $GB_file = <GB_ENTRY>;
      close GB_ENTRY;
   }

   #remove empty space from beginning of file and re-write file
   if ($GB_file =~ m/^\s*/) {
      $GB_file =~ s/^\s*//;
      open  WRITE, "+>$args{gb_file}";
      print WRITE $GB_file;
      close WRITE;
   }

   #is it a genbank file?
   if ($GB_file !~ /^locus/i) {
      print "The accessed file $args{gb_file} is no valid Genbank file";
      open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error parsing Genbank file $args{gb_file}".
                     "\nFile appears to be not in Genbank format\n\n";
      close ERRORLOG;
      return (0, "The accessed file $args{gb_file} is no valid Genbank file");
   }

   #catch header
   'reset' =~ /reset/;
   $GB_file =~ m/^(.*? Location\/Qualifiers\s*\n)/s;
   $header = $1;

   #catch core
   'reset' =~ /reset/;
   $GB_file =~ m/Location\/Qualifiers\s*\n(     [a-zA-Z_]+\s+.*?\n(BASE COUNT|ORIGIN))/s;
   $core = $1;

   #catch source
   'reset' =~ /reset/;
   $core =~ m/(     source.*?\n)(     [a-zA-Z]+|BASE COUNT|ORIGIN)/s;
   $source = $1;
   $core =~ s/$source//s;
   $core =~ s/(BASE COUNT|ORIGIN)//s;

   #catch DNA seq
   'reset' =~ /reset/;
   $GB_file =~ m/\n((BASE COUNT|ORIGIN).*)/s;
   $DNAseq = $1;
   undef $GB_file; #give back memory
   unless (defined $source) {$source = ''};

   #clean up DNA seq
   $DNAseq =~ s/ORIGIN//s;
   $DNAseq =~ s/BASE COUNT[^\n]*?\n//is;
   $DNAseq =~ s/\d//gs;
   $DNAseq =~ s/\s//gs;
   $DNAseq =~ s/\\//g;
   $DNAseq =~ s/\///g;
   $DNAseq =~ s/[^a-z]+//igs;

   #if whole sequence is requested, exit here
   if ($args{process_gb} eq 'sequence') {
      #return references
      return ($header, $source, $core, $DNAseq, undef, undef, undef, undef);
   }


   #process into single features
   @tmp = split(/\n     (\S)/s, $core);
   #cleaning up;
   $tmp[0] =~ s/^\s+//;
   push (@features, $tmp[0]);

   #reset counter
   $counter = 0;
   foreach (@tmp) {
      if ($_ =~ /^\S$/) {
         my ($temp);
         $temp = $tmp[$counter].$tmp[$counter+1];
         #generating valid entries list
         push (@features, $temp);
      }
      $counter++;
   }

   #release memory from temp array
   undef @tmp;

   #generate hash for entries; sort them later by start position
   my $gene_counter    = 100000;
   my $CDS_counter     = 100000;
   my %feature_counter = ();
   my %duplicate_gene  = ();
   my %duplicate_CDS   = ();
   foreach my $entry (@features) {
      next if ($entry !~ /\S+/);
      'reset' =~ m/reset/;
      my $feature = "";
      my $left_bd = "";
      my $right_bd = "";
      my $orientation;

      $entry =~ m/^\s*(\S+)\s*(complement)?\D*?(\d+)\D*?(\d+)/;
      $feature = $1;
      $orientation = $2;
      $left_bd = $3;
      $right_bd = $4;

      #skip source if escaped until here
      next if ($feature eq 'source');

      #modify right boundary for special keys with only left boundary entry
      unless (defined ($4)) {
         #retry matching
         'reset' =~ m/reset/;
         $entry =~ m/\s+\/locus_tag\=\"?.*?(\d+)\s*\"?$/s;
         $feature = $1;
         $orientation = $2;
         $left_bd = $3;
         #set boundaries equal
         $right_bd = $left_bd;
      }

      if (defined $orientation) {
         $orientation = 'antisense';
      } else {
         $orientation = 'sense';
      }
      unless (defined($feature)) {
         print "Could not grab values from $_...";
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing Genbank file $args{gb_file}".
                        "\nCould not grab values from entry $_\n\n";
         close ERRORLOG;
         next;
      }
      unless (defined($left_bd)) {
         print "Could not grab values from $_...";
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing Genbank file $args{gb_file}".
                        "\nCould not grab values from entry $_\n\n";
         close ERRORLOG;
         next;
      }

      #modify feature if necessary to make compatible with key-parsing
      if ($feature =~ m/_/) {
         $feature =~ s/_/~/g;
      }

      #do we have an ORF counter already? if other feature, increase ID count
      if ($feature =~ /^gene$/) {
         my $ORF_number;
         'reset' =~ m/reset/;
         $entry =~ m/\s+\/gene\=\".*?_\s*([\d\s]+)\s*\"/s;
         $ORF_number = $1;

         unless (defined ($ORF_number) && $ORF_number =~ /^\d+$/) {
            'reset' =~ m/reset/;
            $entry =~ m/\s+\/locus_tag\=.*?(\d+)\s*\"/s;
            $ORF_number = $1;
         }
         #no ORF number?
         unless (defined ($ORF_number) && $ORF_number =~ /^\d+$/) {
            $gene_counter++;
            $ORF_number = $gene_counter;
         }
         $ORF_number =~ s/\s+//gs;

         #remove leading '0's
         $ORF_number =~ s/^0+//;

         #nothing left? There have been cases where counts start with '0'
         unless ($ORF_number =~ m/\d+/) {
            $ORF_number = '0';
         }

         #does ORF number already exists?
         while (exists $duplicate_gene{$ORF_number}) {
            $ORF_number = $gene_counter;
            $gene_counter++;
         }
         $duplicate_gene{$ORF_number} = 1;

         #set gene counter
         $feature_counter{$feature} = $ORF_number;

      } elsif ($feature =~ /^CDS$/) {
         my $ORF_number;
         'reset' =~ m/reset/;
         $entry =~ m/\s+\/gene\=\".*?_\s*([\d\s]+)\s*\"/s;
         $ORF_number = $1;
         unless (defined ($ORF_number) && $ORF_number =~ /^\d+$/) {
            'reset' =~ m/reset/;
            $entry =~ m/\s+\/locus_tag\=.*?(\d+)\s*\"/s;
            $ORF_number = $1;
         }
         #no ORF number?
         unless (defined ($ORF_number) && $ORF_number =~ /^\d+$/) {
            $CDS_counter++;
            $ORF_number = $CDS_counter;
         }
         $ORF_number =~ s/\s+//gs;

         #remove leading '0's
         $ORF_number =~ s/^0+//;

         #nothing left? There have been cases where counts start with '0'
         unless ($ORF_number =~ m/\d+/) {
            $ORF_number = '0';
         }

         #does ORF number already exists?
         while (exists $duplicate_CDS{$ORF_number}) {
            $ORF_number = $CDS_counter;
            $CDS_counter++;
         }
         $duplicate_CDS{$ORF_number} = 1;

         #set gene/CDS counter
         $feature_counter{$feature} = $ORF_number;
      } else {
         $feature_counter{$feature}++;
      }

      $genbank{$left_bd.'_'.$right_bd.'_'.$feature.'_'.$feature_counter{$feature}} = $entry;
      push (@feature_list, $left_bd.'_'.$right_bd.'_'.$feature.'_'.$feature_counter{$feature});
   }

   #release memory from feature array
   @features = ();
   undef %duplicate_gene;
   undef %duplicate_CDS;

   #return references
   return ($header, $source, $core, $DNAseq, \@start, \@feature_list, \%genbank, $counter);

}

#generate new Genbank file
sub compile_gb {
   my %args = (update_allowed => 1,
               destination    => '',
               from_rotate    => 0,
               @_
              ); #update_allowed for updating existing genbank files during annotation.
                 #should be set to 0 if used from misc modules outside annotation
                 #otherwise will overwrite new feature changes
                 #from_rotate leaves writing the GB file to the module and only returns a GB reference
   my (%genbank, @gb_order, $gene_annotation, $CDS_annotation, $gene_locustag, $CDS_locustag,
       $CDS_product, $CDS_proteinid, $new_genbank, $counter, $gb_filename);
   $counter = 0;

   #define destination
   unless ($args{destination} =~ /\S+/) {
      $args{destination} = ${$args{ini_ref}}{results};
   }
   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling new Genbank file: $args{input_file}",
                   label        => 'Compiling features'
                  );
   &show_pbar_2;


   my $max_count = $#{$args{feature_list_ref}};
   ${$args{progress_bar}}->configure(-label=>"Compiling new Genbank file");
   ${$args{main_window}}->update;

   #parsing input genbank file if update option is selected
   if (${$args{auto_ini_ref}}{update_genbank} == 1 && $args{update_allowed} == 1) {
      ($gene_annotation, $CDS_annotation,
       $gene_locustag, $CDS_locustag,
       $CDS_product, $CDS_proteinid)  = &update_parsing (main_window      => $args{main_window},
                                                         progress_bar     => $args{progress_bar},
                                                         auto_ini_ref     => $args{auto_ini_ref},
                                                         ini_ref          => $args{ini_ref},
                                                         input_file       => $args{input_file},
                                                         genbank_ref      => $args{genbank_ref},
                                                         feature_list_ref => $args{feature_list_ref},
                                                        );
   }

   #re-create header
   &update_pbar_2(title        => "Compiling new Genbank file: $args{input_file}",
                  label        => "Compiling new Genbank file...generating header",
                  progress     => 1,
                 );
   chomp($args{genbank_source});
   chomp($args{genbank_header});

   #update date in gb header?
   if (${$args{auto_ini_ref}}{update_gb_date} == 1) {
      my $date = &set_date;
      $args{genbank_header} =~ s/\S+\s*\n/$date\n/;
   }

   #remove DBlink entry in gb header?
   if (${$args{auto_ini_ref}}{remove_dblink} == 1) {
      $args{genbank_header} =~ s/\nDBLINK[^\n]+?\n/\n/s;
   }

   #right format for genbank source?
   unless ($args{genbank_source} =~ /^     /) {
      $args{genbank_source} = '     '.$args{genbank_source};
   }

   #do we have feature/qualifier comment?
   unless ($args{genbank_header} =~ /FEATURES             Location\/Qualifiers/ ||
           $args{genbank_source} =~ /FEATURES             Location\/Qualifiers/) {
      $args{genbank_header} .= "\n".'FEATURES             Location/Qualifiers';
   }

   $new_genbank = $args{genbank_header}."\n".$args{genbank_source};

   #iterate over all entries
   #------------------------

   #check if we have exons
   if ($#{$args{auto_ini_ref}{exons}} >= 0) {
      foreach my $exon (@{$args{auto_ini_ref}{exons}}) {
         next if ($exon !~ /\w+/);
         my ($feature_range, @feature_range, $left_bd, $right_bd, $mRNA);

         #get intron boundaries
         'reset' =~ m/reset/;
         $exon =~ m/^(\d+)\.\.(\d+)/;
         $left_bd = $1;
         $right_bd = $2;

         #determine orientation
         @feature_range = grep { m/$left_bd\_$right_bd\_CDS_/ } @{$args{feature_list_ref}};

         if (${$args{genbank_ref}}{$feature_range[0]} =~ /CDS\s+complement/) {
            $mRNA = "\n     mRNA            complement\(join\($exon\)\)";
         } else {
            $mRNA = "\n     mRNA            join\($exon\)";
         }

         #add to genbank file
         $new_genbank .= $mRNA;
      }
   }

   #split Genbank key order into array
   @gb_order = split /;/,${$args{auto_ini_ref}}{genbank_key_order};

   &update_pbar_2(title        => "Compiling new Genbank file: $args{input_file}",
                  label        => "Compiling new Genbank file...generating entries",
                  progress     => 1,
                 );

   #sort feature_list_ref via start position and genbank_key_order
   $counter = 1;
   foreach my $key (@gb_order) {
      map { $_ =~ s/$key/$counter\-$key/ } @{$args{feature_list_ref}};
      $counter++;
   }
   #add generic number to those entries not listed in key order
   foreach my $entry (@{$args{feature_list_ref}}) {
      $entry =~ s/^(\d+_\d+_)(\D+\_.+)/$1$counter\-$2/;
      $counter++;
   }
   my @sorted =
         map  $_->[0] =>
         sort { $a->[1] <=> $b->[1] ||
                $a->[2] <=> $b->[2] }
         map  [ $_, m/^(\d+)\_/, m/^\d+_\d+_(\d+)\-.*?_\d+/ ]
         => @{$args{feature_list_ref}};

   $counter = 1;
   $max_count = $#sorted + 1;
   foreach my $entry (@sorted) {
      #remove genbank order key again
      $entry =~ s/_\d+\-(.+)$/\_$1/;

      #create access key for maintaining previous annotation
      my $short = $entry;
      $short =~ s/_\d+$//;

      #now add to GB file
      chomp(${$args{genbank_ref}}{$entry});
      ${$args{genbank_ref}}{$entry} =~ s/^\n//;

      #check if previous annotation exists
      if (${$args{auto_ini_ref}}{update_genbank} == 1 && $args{update_allowed} == 1 && ${$args{auto_ini_ref}}{keep_annotation} == 1) {
         #modify gene annotation
         if (defined $gene_annotation->{$short}            &&
             ${$args{genbank_ref}}{$entry} =~ m/^\s*gene/  &&
             ${$args{auto_ini_ref}}{maintain_gene}   == 1) {
            if (${$args{genbank_ref}}{$entry} =~ m/\/gene\=\".*?_\s*\d+\s*\"/s) {
               ${$args{genbank_ref}}{$entry} =~ s/\/gene\=\".*?(_\s*\d+)\s*\"/\/gene\=\"$gene_annotation->{$short}$1\"/s;
            } else {
               ${$args{genbank_ref}}{$entry} =~ s/\/gene\=\"[^\"]*?\"/\/gene\=\"$gene_annotation->{$short}\"/s;
            }
         }

         if (defined $gene_locustag->{$short}                &&
             ${$args{genbank_ref}}{$entry} =~ m/^\s*gene/    &&
             ${$args{auto_ini_ref}}{maintain_locus_tag} == 1) {
            ${$args{genbank_ref}}{$entry} =~ s/\/locus_tag\=\"?\S+/\/locus_tag\=\"$gene_locustag->{$short}\"/s;
         }

         #if it was not defined AND add_qualifier is NOT set, delete respective gene / locus_tag qualifiers
         if (!defined $gene_annotation->{$short} && ${$args{auto_ini_ref}}{add_qualifiers} == 0) {
            ${$args{genbank_ref}}{$entry} =~ s/                     \/gene\=\".+?\"\s*\n//s;
         }
         if (!defined $gene_locustag->{$short} && ${$args{auto_ini_ref}}{add_qualifiers} == 0) {
            ${$args{genbank_ref}}{$entry} =~ s/                     \/locus_tag\=\"?.+?\"?\s*\n//s;
         }

         #add gene qualifier if not present and fill-in is selected
         if (${$args{auto_ini_ref}}{add_qualifiers} == 1) {
            #create artificial CDS key assuming same start-stop positions for gene and CDS
            my $CDS_key = $short;
            $CDS_key =~ s/_gene/_CDS/;
            #easy case first, gene feature does not contain gene feature
            if (${$args{genbank_ref}}{$entry} =~ m/^\s*gene/ && ${$args{genbank_ref}}{$entry} !~ m/\/gene\=/) {
               if (defined $gene_annotation->{$short}) {
                  ${$args{genbank_ref}}{$entry} =~ s/(^\s*gene[^\n]+\n)/$1                     \/gene=\"$gene_annotation->{$short}\"\n/;#"
               } elsif (!defined $gene_annotation->{$short} && defined $CDS_annotation->{$CDS_key}) {
                  ${$args{genbank_ref}}{$entry} =~ s/(^\s*gene[^\n]+\n)/$1                     \/gene=\"$CDS_annotation->{$CDS_key}\"\n/;#"
               } elsif (!defined $gene_annotation->{$short} && defined $CDS_product->{$CDS_key}) {
                  ${$args{genbank_ref}}{$entry} =~ s/(^\s*gene[^\n]+\n)/$1                     \/gene=\"$CDS_product->{$CDS_key}\"\n/;#"
               }
            }
            #if there is gene qualifier already but no previous gene qualifier defined, replace with previous CDS gene or product qualifier
            if (${$args{genbank_ref}}{$entry} =~ m/^\s*gene/ && ${$args{genbank_ref}}{$entry} =~ m/\/gene\=/) {
               if (!defined $gene_annotation->{$short} && defined $CDS_annotation->{$CDS_key}) {
                  if (${$args{genbank_ref}}{$entry} =~ m/\/gene\=\".*?_\s*\d+\s*\"/s) {
                     ${$args{genbank_ref}}{$entry} =~ s/\/gene\=\".*?(_\s*\d+)\s*\"/\/gene\=\"$CDS_annotation->{$CDS_key}$1\"/s;
                  } else {
                     ${$args{genbank_ref}}{$entry} =~ s/\/gene\=\"[^\"]*?\"/\/gene\=\"$CDS_annotation->{$CDS_key}\"/s;
                  }
               } elsif (!defined $gene_annotation->{$short} && !defined $CDS_annotation->{$CDS_key} && defined $CDS_product->{$CDS_key}) {
                  if (${$args{genbank_ref}}{$entry} =~ m/\/gene\=\".*?_\s*\d+\s*\"/s) {
                     ${$args{genbank_ref}}{$entry} =~ s/\/gene\=\".*?(_\s*\d+)\s*\"/\/gene\=\"$CDS_product->{$CDS_key}$1\"/s;
                  } else {
                     ${$args{genbank_ref}}{$entry} =~ s/\/gene\=\"[^\"]*?\"/\/gene\=\"$CDS_product->{$CDS_key}\"/s;
                  }
               }
            }
         }
         #modify CDS annotation
         if (${$args{genbank_ref}}{$entry} =~ m/^\s*CDS/) {
            #gene annotation
            if (defined $CDS_annotation->{$short} && ${$args{genbank_ref}}{$entry} =~ m/\/gene\=/s && ${$args{auto_ini_ref}}{maintain_gene} == 1 ) {
               if (${$args{genbank_ref}}{$entry} =~ m/\/gene\=\".*?_\s*\d+\s*\"/s) {
                  ${$args{genbank_ref}}{$entry} =~ s/\/gene\=\".*?(_\s*\d+)\s*\"/\/gene\=\"$CDS_annotation->{$short}$1\"/s;
               } else {
                  ${$args{genbank_ref}}{$entry} =~ s/\/gene\=\"[^\"]*?\"/\/gene\=\"$CDS_annotation->{$short}\"/s;
               }
            } elsif (defined $CDS_annotation->{$short} && ${$args{genbank_ref}}{$entry} !~ m/\/gene\=/s && ${$args{auto_ini_ref}}{maintain_gene} == 1 ) {
               ${$args{genbank_ref}}{$entry} =~ s/\n/\n                     \/gene=\"$CDS_annotation->{$short}\"\n/s;
            }

            #locus_tag annotation
            if (defined $CDS_locustag->{$short} && ${$args{genbank_ref}}{$entry} =~ m/\/locus_tag\=/s && ${$args{auto_ini_ref}}{maintain_locus_tag} == 1 ) {
               ${$args{genbank_ref}}{$entry} =~ s/\/locus_tag\=\"?\S+/\/locus_tag\=\"$CDS_locustag->{$short}\"/s;
            } elsif (defined $CDS_locustag->{$short} && ${$args{genbank_ref}}{$entry} !~ m/\/locus_tag\=/s && ${$args{auto_ini_ref}}{maintain_locus_tag} == 1 ) {
               ${$args{genbank_ref}}{$entry} =~ s/\n/\n                     \/locus_tag\=\"$CDS_locustag->{$short}\"\n/s;
            }

            #product annotation
            if (defined $CDS_product->{$short} && ${$args{genbank_ref}}{$entry} =~ m/\/product\=/ && ${$args{auto_ini_ref}}{maintain_product} == 1 ) {
               ${$args{genbank_ref}}{$entry} =~ s/\/product\=\"[^\"]*?\"/\/product\=\"$CDS_product->{$short}\"/s;
            } elsif (defined $CDS_product->{$short} && ${$args{genbank_ref}}{$entry} !~ m/\/product\=/ && ${$args{auto_ini_ref}}{maintain_product} == 1 ) {
               ${$args{genbank_ref}}{$entry} =~ s/\n/\n                     \/product\=\"$CDS_product->{$short}\"\n/s;
            }

            #protein annotation
            if (defined $CDS_proteinid->{$short} && ${$args{genbank_ref}}{$entry} =~ m/\/protein_id\=/ && ${$args{auto_ini_ref}}{maintain_protein_id} == 1 ) {
               ${$args{genbank_ref}}{$entry} =~ s/\/protein_id\=\"[^\"]*?\"/\/protein_id\=\"$CDS_proteinid->{$short}\"/s;
            } elsif (defined $CDS_proteinid->{$short} && ${$args{genbank_ref}}{$entry} !~ m/\/protein_id/ && ${$args{auto_ini_ref}}{maintain_protein_id} == 1 ) {
               ${$args{genbank_ref}}{$entry} =~ s/\n/\n                     \/protein_id\=\"$CDS_proteinid->{$short}\"\n/s;
            }; #"
         }
         #add gene qualifier if not present and fill-in is selected
         if (${$args{auto_ini_ref}}{add_qualifiers} == 1) {
            #create artificial gene key assuming same start-stop positions for gene and CDS
            my $gene_key = $short;
            $gene_key =~ s/_CDS/_gene/;
            #easy case first, gene feature does not contain gene feature
            if (${$args{genbank_ref}}{$entry} =~ m/^\s*CDS/ && ${$args{genbank_ref}}{$entry} !~ m/\/gene\=/) {
               if (defined $CDS_annotation->{$short}) {
                  ${$args{genbank_ref}}{$entry} =~ s/(^\s*CDS[^\n]+\n)/$1                     \/gene=\"$CDS_annotation->{$short}\"\n/;#"
               } elsif (!defined $CDS_annotation->{$short} && defined $gene_annotation->{$gene_key}) {
                  ${$args{genbank_ref}}{$entry} =~ s/(^\s*gene[^\n]+\n)/$1                     \/gene=\"$gene_annotation->{$gene_key}\"\n/;#"
               } elsif (!defined $CDS_annotation->{$short} && defined $CDS_product->{$short}) {
                  ${$args{genbank_ref}}{$entry} =~ s/(^\s*gene[^\n]+\n)/$1                     \/gene=\"$CDS_product->{$short}\"\n/;#"
               }
            }
            #if there is gene qualifier already but no previous gene qualifier defined, replace with previous CDS gene or product qualifier
            if (${$args{genbank_ref}}{$entry} =~ m/^\s*CDS/ && ${$args{genbank_ref}}{$entry} =~ m/\/gene\=/) {
               if (!defined $CDS_annotation->{$short} && defined $gene_annotation->{$gene_key}) {
                  if (${$args{genbank_ref}}{$entry} =~ m/\/gene\=\".*?_\s*\d+\s*\"/s) {
                     ${$args{genbank_ref}}{$entry} =~ s/\/gene\=\".*?(_\s*\d+)\s*\"/\/gene\=\"$gene_annotation->{$gene_key}$1\"/s;
                  } else {
                     ${$args{genbank_ref}}{$entry} =~ s/\/gene\=\"[^\"]*?\"/\/gene\=\"$gene_annotation->{$gene_key}\"/s;
                  }
               } elsif (!defined $CDS_annotation->{$short} && !defined $gene_annotation->{$gene_key} && defined $CDS_product->{$short}) {
                  if (${$args{genbank_ref}}{$entry} =~ m/\/gene\=\".*?_\s*\d+\s*\"/s) {
                     ${$args{genbank_ref}}{$entry} =~ s/\/gene\=\".*?(_\s*\d+)\s*\"/\/gene\=\"$CDS_product->{$short}$1\"/s;
                  } else {
                     ${$args{genbank_ref}}{$entry} =~ s/\/gene\=\"[^\"]*?\"/\/gene\=\"$CDS_product->{$short}\"/s;
                  }
               }
            }
         }
      }

      #test if feature needs to be modified
      if ($entry =~ /~/) {
         ${$args{genbank_ref}}{$entry} =~ s/^(\s*\w+?)~/$1_/;
      }

      if (${$args{genbank_ref}}{$entry} !~ /^     /) {
         ${$args{genbank_ref}}{$entry} = "     ".${$args{genbank_ref}}{$entry};
      }

      $new_genbank .= "\n".${$args{genbank_ref}}{$entry};
      #delete key
      #delete ${$args{genbank_ref}}{$entry};
      if (($counter % 100) == 0) {
         &update_pbar_2(title        => "Compiling new Genbank file: $args{input_file}",
                        label        => "Compiling new Genbank file...generating entries",
                        progress     => ($counter / $max_count) * 100,
                       );
      }
      $counter++;
   }

   #format and add nucleotide sequence
   &update_pbar_2(title        => "Compiling new Genbank file: $args{input_file}",
                  label        => "Compiling new Genbank file...formatting nt sequence",
                  progress     => 1,
                 );
   ${$args{main_window}}->update;

   #format
   if (defined $args{DNAseq} && $args{DNAseq} =~ /\w+/) {
      #remove any empty space
      $args{DNAseq} =~ s/\s+//gs;
      $max_count = length ($args{DNAseq});
   } else {
      &update_pbar_2(title        => "Compiling new Genbank file: $args{input_file}",
                     label        => "Compiling new Genbank file...reading source nt sequence",
                     progress     => 1,
                    );
      #read input sequence file
      my ($seqref) = &slurp(main_window  => $args{main_window},
                            progress_bar => $args{progress_bar},
                            auto_ini_ref => $args{auto_ini_ref},
                            ini_ref      => $args{ini_ref},
                            directory    => ${$args{ini_ref}}{input_files},
                            filename     => $args{input_file}
                           );
      if ($seqref eq '0') {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error reading sequence file $args{input_file}.\nSkipping",
                                                       -bitmap  => 'error',
                                                       -buttons => ['ok']);
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error compiling Genbank files".
                        "\nError reading sequence file $args{input_file}.\nSkipping\n\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }
      $$seqref =~ s/^>[^\n\r]*?\n\r?//s;
      $args{DNAseq} = $$seqref;

      #remove any empty space
      $args{DNAseq} =~ s/\s+//gs;
      $max_count = length ($args{DNAseq});
      &update_pbar_2(title        => "Compiling new Genbank file: $args{input_file}",
                     label        => "Compiling new Genbank file...formatting nt sequence",
                     progress     => 1,
                    );
   }
   $args{DNAseq} = lc $args{DNAseq};
   my $num_A = $args{DNAseq} =~ s/(a)/$1/gi;
   my $num_C = $args{DNAseq} =~ s/(c)/$1/gi;
   my $num_G = $args{DNAseq} =~ s/(g)/$1/gi;
   my $num_T = $args{DNAseq} =~ s/(t)/$1/gi;
   my $num_N = $args{DNAseq} =~ s/([^acgt])/$1/gi;
   unless ($num_A =~ /\d+/) {$num_A = 0};
   unless ($num_C =~ /\d+/) {$num_C = 0};
   unless ($num_G =~ /\d+/) {$num_G = 0};
   unless ($num_T =~ /\d+/) {$num_T = 0};
   unless ($num_N =~ /\d+/) {$num_N = 0};

   $new_genbank .= "\nBASE COUNT    $num_A a   $num_C c   $num_G g   $num_T t  $num_N other\n".
                   "ORIGIN\n";
   $i=1;
   my $test="false";
   my $dna_len = 0;
   $dna_len = length ($args{DNAseq});
   my $dnasequence = "";
   #define format
   format DNA_SEQ =
@######### ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<<
$i,        $args{DNAseq}, $args{DNAseq}, $args{DNAseq}, $args{DNAseq}, $args{DNAseq}, $args{DNAseq}
.
   open DNA_SEQ, '>', \$dnasequence;

   while ($i <= $dna_len) {
      if ((($i-1) % 30000) == 0) {
         &update_pbar_2(title        => "Compiling new Genbank file: $args{input_file}",
                        label        => "Compiling new Genbank file...formatting nt sequence",
                        progress     => ($i / $max_count) * 100,
                       );
         ${$args{main_window}}->update;
      }
      write DNA_SEQ;
      $i = $i +60;
   }
   close DNA_SEQ;

   $new_genbank .= $dnasequence."\/\/\n";

   #write new Genbank file to result folder
   &update_pbar_2(title        => "Compiling new Genbank file: $args{input_file}",
                  label        => "Compiling new Genbank file...saving new Ganbank file",
                  progress     => 100,
                 );
   $args{input_file} =~ s/_GAMOLAdna$//;
   if ($args{input_file} =~ /\.(gbk|gb)$/) {
      $gb_filename = $args{destination}.'/'.$args{input_file};
   } else {
      $gb_filename = $args{destination}.'/'.$args{input_file}.'.gb';
   }

   #if parsing is for rotating GB files, return reference and leave writing the GB file to the module,
   #otherwise write new GB file
   if ($args{from_rotate} == 1) {
      &hide_pbar_2;
      return (\$new_genbank);
   } else {
      open WRITE, "+>$gb_filename" or die "\nError in saving new Genbank file $gb_filename";
      print WRITE $new_genbank;
      close WRITE;
      undef $new_genbank;
      &hide_pbar_2;
      return (1);
   }
}


sub update_parsing {
   my %args = @_;
   my (@keep, @update, $gene_annotation, $CDS_annotation, $gene_locustag, $CDS_locustag,
       $CDS_product, $CDS_proteinid, $gb_file, $counter, $max_count, $match);

   #first, parse genbank file
   $gb_file = $args{input_file};
   $gb_file =~ s/_GAMOLAdna$//;

   #return if input file was fasta
   return unless (-e ${$args{ini_ref}}{move_gb}.'/'.$gb_file);

   my ($header, $source, $core, $DNAseq, $start_ref,
       $feature_list_ref, $genbank_ref, $feature_counter) = &gb_parser (main_window   => $args{main_window},
                                                                        progress_bar  => $args{progress_bar},
                                                                        auto_ini_ref  => $args{auto_ini_ref},
                                                                        ini_ref       => $args{ini_ref},
                                                                        gb_file       => ${$args{ini_ref}}{move_gb}.'/'.$gb_file

   );
   #now remove all feature that have to be replaced/updated

   #get gene and CDS annotations if selected
   #if (${$args{auto_ini_ref}}{blast_selector} == 1 && ${$args{auto_ini_ref}}{maintain_annotation} =~ /yes/i) {
   if (${$args{auto_ini_ref}}{update_genbank} == 1 && ${$args{auto_ini_ref}}{keep_annotation} == 1) {
      #get gene annotation
      my @temp = grep (/_gene_/,@{$feature_list_ref});
      foreach my $entry (@temp) {
         #remove ID from entry to make it compatible
         my $short = $entry;
         $short =~ s/_\d+$//;

         #catch gene annotation
         if (${$args{auto_ini_ref}}{maintain_gene} == 1) {
            'reset' =~ m/reset/;
            ${$genbank_ref}{$entry} =~ m/\/gene\=\"([^\"]*?)\"/s;
            my $annotation = $1;
            if (defined $annotation) {
               $annotation =~ s/_\s*\d+\s*$//;
               $gene_annotation->{$short} = $annotation;
            }
         }

         #catch locus_tag if present
         if (${$args{auto_ini_ref}}{maintain_locus_tag} == 1) {
            'reset' =~ m/reset/;
            ${$genbank_ref}{$entry} =~ m/\/locus_tag\=\"?\s*(\S+)\s*\"?/s;
            my $locus_tag = $1;
            if (defined $locus_tag) {
               $locus_tag =~ s/\n.+//s;
               $locus_tag =~ s/\s+/ /gs;
               $locus_tag =~ s/\"//g;
               $gene_locustag->{$short} = $locus_tag;
            }
         }
      }
      #get CDS annotation
      undef @temp;
      @temp = grep (/_CDS_/,@{$feature_list_ref});
      foreach my $entry (@temp) {
         #remove ID from entry to make it compatible
         my $short = $entry;
         $short =~ s/_\d+$//;

         #catch gene annotation
         if (${$args{auto_ini_ref}}{maintain_gene} == 1) {
            'reset' =~ m/reset/;
            ${$genbank_ref}{$entry} =~ m/\/gene\=\"([^\"]*?)\"/s;
            my $annotation = $1;
            if (defined $annotation) {
               $annotation =~ s/_\s*\d+\s*$//;
               $CDS_annotation->{$short} = $annotation;
            }
         }

         #catch locus_tag if present
         if (${$args{auto_ini_ref}}{maintain_locus_tag} == 1) {
            'reset' =~ m/reset/;
            ${$genbank_ref}{$entry} =~ m/\/locus_tag\=\"?\s*(\S+)\s*\"?/s;
            my $locus_tag = $1;
            if (defined $locus_tag) {
               $locus_tag =~ s/\n.+//s;
               $locus_tag =~ s/\s+/ /gs;
               $locus_tag =~ s/\"//g;
               $CDS_locustag->{$short} = $locus_tag;
            }
         }

         #catch product if present
         if (${$args{auto_ini_ref}}{maintain_product} == 1) {
            'reset' =~ m/reset/;
            ${$genbank_ref}{$entry} =~ m/\/product\=\"([^\"]+?)\"/s;
            my $product = $1;
            if (defined $product) {
               $product =~ s/\s+/ /gs;
               $product =~ s/\"//g;
               $CDS_product->{$short} = $product;
            }
         }

         #catch protein_id if present
         if (${$args{auto_ini_ref}}{maintain_protein_id} == 1) {
            'reset' =~ m/reset/;
            ${$genbank_ref}{$entry} =~ m/\/protein_id\=\"([^\"]+?)\"/s;
            my $prot_id = $1;
            if (defined $prot_id) {
               $prot_id =~ s/\s+/ /gs;
               $prot_id =~ s/\"//g;
               $CDS_proteinid->{$short} = $prot_id;
            }
         }
      }
   }

   #define features to take out for updating
   $match = '(';
   if (${$args{auto_ini_ref}}{blast_selector} == 1) {
      $match .= '_gene_|';
      $match .= '_CDS_|';
   }
   if (${$args{auto_ini_ref}}{COG_selector} == 1) {
      $match .= '_COG_match_|_COG~match_|';
   }
   if (${$args{auto_ini_ref}}{Pfam_selector} == 1) {
      $match .= '_PFAM_match_|_PFAM~match_|';
   }
   if (${$args{auto_ini_ref}}{TIGRfam_selector} == 1) {
      $match .= '_TIGR_match_|_TIGR~match_|';
   }
   if (${$args{auto_ini_ref}}{trna_selector} == 1) {
      $match .= '_tRNA_|';
   }
   if (${$args{auto_ini_ref}}{rrna_selector} == 1) {
      $match .= '_rRNA_|';
   }
   if (${$args{auto_ini_ref}}{ncrna_selector} == 1) {
      $match .= '_ncRNA_|_misc~RNA_|';
   }
   if (${$args{auto_ini_ref}}{signalp_selector} == 1) {
      $match .= '_signalP_|';
   }
   if (${$args{auto_ini_ref}}{tmhmm_selector} == 1) {
      $match .= '_TMM_|';
   }
   if (${$args{auto_ini_ref}}{terminator_selector} == 1) {
      $match .= '_terminator_|';
   }
   if (${$args{auto_ini_ref}}{vector_screen_selector} == 1) {
      $match .= '_Vector_match_|_Vector~match_|';
   }
   if (${$args{auto_ini_ref}}{CRISPR_selector} == 1) {
      $match .= '_CRISPR_|';
   }
   $match =~ s/\|$//;
   $match .= ')';
   $match =~ s/\(\)//; #clear if empty

   #grep kept and updated features
   @keep   = grep (!/$match/,@{$feature_list_ref});
   @update = grep ( /$match/,@{$feature_list_ref});

   #delete features to be updated from original Genbank feature hash
   foreach my $entry (@update) {
      delete ${$genbank_ref}{$entry};
   }

   #merge feature list and remaining genbank hash
   #$args{feature_list_ref} contains new entries generated by this run
   #first work on feature_list_ref: remove all entries that are updated from original GB file, then add the ones to keep
   #remove gene model from args{feature_list_ref} -> replace with gene model from original Genbank file (feature_list_ref).
   #This will take care of all the tRNA and rRNA gene features

   #map tRNA/rRNA start stop positions to a temp hash
   #will be used to selectively keep tRNA/rRNA gene features
   my %rtRNA_hash = ();
   foreach my $entry (@{$args{feature_list_ref}}) {
      if ($entry =~ m/(_tRNA_|_rRNA_)/) {
         'reset' =~ m/reset/;
         $entry =~ m/^(\d+)_(\d+)_/;
         $rtRNA_hash{$1} = $2;
      }
   }

   #remove OLD/existing features that will be updated
   for my $i ( 0 .. $#{$feature_list_ref} ) {
      # undefine element
      my ($feature_key, $rna_left, $rna_right) = '';
      'reset' =~ m/reset/;
      ${$feature_list_ref}[$i] =~ m/^(\d+)_(\d+)_(\D+)_\d+/;
      ($rna_left, $rna_right, $feature_key) = ($1, $2, $3);
      ##remove features that will be updated but keep original gene model
      #undef ${$feature_list_ref}[$i] if ($match =~ m/$feature_key/ && $feature_key !~ m/(gene|CDS)/) ;

      #remove features that will be updated
      if ($feature_key =~ m/(gene|CDS)/ && $rna_right == $rtRNA_hash{$rna_left}) {
         #keep gene model in old results if coincides with RNA gene features
         next;
      }
      undef ${$feature_list_ref}[$i] if ($match =~ m/$feature_key/);
   }

   #remove gene/CDS from NEW results to keep tRNA/rRNA gene features if they exist
   for my $i ( 0 .. $#{$args{feature_list_ref}} ) {
      #undef ${$args{feature_list_ref}}[$i] if (${$args{feature_list_ref}}[$i] =~ m/(_gene_|_CDS_)/) ; #remove gene model to be replaced with original one
      'reset' =~ m/reset/;
      ${$args{feature_list_ref}}[$i] =~ m/^(\d+)_(\d+)_/;
      my ($rna_left, $rna_right) = ($1, $2);
      if (${$args{feature_list_ref}}[$i] =~ m/(_gene_|_CDS_)/ &&
          exists $rtRNA_hash{$rna_left}                       &&
          $rna_right == $rtRNA_hash{$rna_left}) {
         undef ${$args{feature_list_ref}}[$i]; #remove only NEW gene/CDS entries if they relate to a tRNA/rRNA feature in original genbank file
      }
      #remove "gene" feature if not in match string AND isn't a tRNA/rRNA
      if ($match !~ m/_gene_/ &&
          ${$args{feature_list_ref}}[$i] =~ m/_gene_/ &&
          $rna_right != $rtRNA_hash{$rna_left}) {
         undef ${$args{feature_list_ref}}[$i]
      }
   }
   # now remove undefined elements
   @{$feature_list_ref}       = grep{ defined } @{$feature_list_ref};
   @{$args{feature_list_ref}} = grep{ defined } @{$args{feature_list_ref}};

   #now add kept features
   #@{$args{feature_list_ref}} = (@{$args{feature_list_ref}}, @keep);
   push(@{$args{feature_list_ref}}, @{$feature_list_ref});

   #now update features: default behaviour for 'merge' is LEFT precedence: only new values from right hash will be added
   %{$args{genbank_ref}}      = %{ merge( \%{$args{genbank_ref}}, \%{$genbank_ref} ) };

   #remove potential duplicates
   %seen = ();
   @{$args{feature_list_ref}} = grep { ! $seen{$_} ++ } @{$args{feature_list_ref}};

   undef @keep;
   undef @update;
   return ($gene_annotation, $CDS_annotation, $gene_locustag, $CDS_locustag,
           $CDS_product, $CDS_proteinid);
}



1;