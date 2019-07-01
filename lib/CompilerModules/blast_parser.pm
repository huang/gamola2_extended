#!/opt/ActivePerl-5.8/bin/perl

#Blast result parser
#input arguments: main_window, progress_bar, gb_file, gene_model



package CompilerModules::blast_parser;
use strict;
use initialise::read_me        qw(:DEFAULT);
use initialise::gb_header      qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);

use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&blast_parser &gene_model_parser &blast_summary);
use vars qw();


#local variables
my (%args, %seen, $value);

format ENTRY =
~~                   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$value
.

sub blast_parser {
   my %args = @_;

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling Blast results",
                   label        => ''
                  );
   &show_pbar_2;

   #parse through Blast results
   foreach (my $i = 1; $i < $args{counter}; $i++) {
      my ($filename, $id, $left_bd, $right_bd, $orientation, $aa_seq, $all_hits, $best_hit, $evalue,
          $boundary_gene, $boundary_CDS, $locus_tag, $product,
          $gene, $note, $translation, $codon_table,
          $key_gene, $key_CDS, $feature_gene, $feature_CDS);

      #get filename, left_bd, right_bd, orientation
      'reset' =~ m/reset/;
      $args{orf_to_id}->{$i} =~ m/^(.*?)___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
      $filename    = $1;
      $id          = $2;
      $left_bd     = $3;
      $right_bd    = $4;
      $orientation = $5;
      unless (defined $orientation) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse from entry ${$args{orf_to_id}}->{$i}\.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing global Blast result hash".
                        "\nCould not parse from entry ${$args{orf_to_id}}->{$i}\.\n\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }

      #update status bar
      if (($i % 20) == 0) {
         &update_pbar_2(title        => "Compiling Blast results",
                        label        => "Compiling Blast results for $filename",
                        progress     => ($i / $args{counter}) * 100,
                       );
      }

      #read input file
      (my $file_ref) = slurp(main_window => $args{main_window},
                             directory   => ${$args{ini_ref}}{blast_results},
                             filename    => $filename.'_'.$id
                            );
      if ($file_ref eq '0') {&hide_pbar_2; return (0)};
      #catch aa sequence and best Blast hit
      if (${$file_ref} =~ /No Blast Hits found/s) {
         'reset' =~ m/reset/;
         ${$file_ref} =~ m/Deduced aminoacid sequence:\s*(\S*?)\s+/si;
         $aa_seq = $1;
         unless (defined $aa_seq) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse aa sequence from Blast entry $filename.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing Blast results for file $filename".
                           "\nCould not parse aa sequence from Blast entry $filename.\n\n";
            close ERRORLOG;
            &hide_pbar_2;
            return (0);
         }
         $best_hit = 'unknown';
         $evalue   = 1000000;
      } else {
         #get all blast hits
         'reset' =~ m/reset/;
         ${$file_ref} =~ m/Deduced aminoacid sequence:\s*(\S*?)\s+.*?Blast overview.*?(\>.*?)\nComplete list of Blast results/si;
         $aa_seq = $1;
         $all_hits = $2;
         unless (defined $all_hits) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not get all hits $all_hits \nfrom Blast entry $filename.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing Blast results for file $filename".
                           "\nCould not get all hits $all_hits \nfrom Blast entry $filename.\n\n";
            close ERRORLOG;
            next;
         }

         #create lookup hash
         my %seen;
         my @tmp = split/,/, ${$args{auto_ini_ref}}{filter_blast};
         foreach (@tmp) {$seen{$_}++};
         undef @tmp;

         #initialise filter
         $best_hit = '';
         if (${$args{auto_ini_ref}}{filter_blast} =~ /\w+/) {
            while ($best_hit eq '' && $all_hits =~ /\>/) {
               'reset'   =~ m/reset/;
               $all_hits =~ m/\>([^\n\r]*?)[\n\r]/;
               $best_hit = $1;
               unless (defined $best_hit) {
                  my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                                -text    => "Could not parse best hit from Blast entry $filename.",
                                                                -buttons => ['OK'],
                                                                -bitmap  => 'info');
                  $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
                  $error_msg-> Show();
                  open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
                  print ERRORLOG "Error parsing Blast results for file $filename".
                                 "\nCould not parse best hit from Blast entry $filename.\n\n";
                  close ERRORLOG;
                  $all_hits = '';
                  $best_hit = '';
                  last;
               }
               $all_hits =~ s/\>\Q$best_hit\E[\n\r]+//;

               #test if seen
               foreach my $filter (keys %seen) {
                  if ($best_hit =~ m/$filter/i) {
                     $best_hit = '';
                     last;
                  }
               }

               #skip if COG result; if selected
               if (${$args{auto_ini_ref}}{allow_COG_anno} == 1 &&
                  $best_hit =~ m/^\>.+? COG\d+\: /) {
                  $best_hit = '';
                  last;
               }

            }
         } else {
            'reset'   =~ m/reset/;
            $all_hits =~ m/\>([^\n\r]*?)[\n\r]/;
            $best_hit = $1;
            unless (defined $best_hit) {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Could not parse best hit from Blast entry $filename.",
                                                             -buttons => ['OK'],
                                                             -bitmap  => 'info');
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error parsing Blast results for file $filename".
                              "\nCould not parse best hit from Blast entry $filename.\n\n";
               close ERRORLOG;
               $all_hits = '';
               $best_hit = '';
               last;
            }
         }
         undef %seen;

         #no hits? set to unknown
         if ($best_hit eq '') {
            $best_hit = 'unknown';
            $evalue   = 1000000;
         }
         #otherwise parse for evalue
         else {
            #catch evalue
            'reset'   =~ m/reset/;
            $best_hit =~ m/Expect\=(.+)$/;
            $evalue   = $1;
            unless (defined $evalue) {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Could not parse evalue from entry $best_hit.",
                                                             -buttons => ['OK'],
                                                             -bitmap  => 'info');
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error parsing Blast results for file $filename".
                              "\nCould not parse evalue from entry $best_hit.\n\n";
               close ERRORLOG;
               next;
            }
         }
         if ($evalue =~ /^e/) {$evalue = '1'.$evalue};
      }

      #create boundary entry
      if ($orientation eq 'sense') {
         $boundary_gene = '     gene            '.$left_bd.'..'.$right_bd."\n";
         $boundary_CDS  = '     CDS             '.$left_bd.'..'.$right_bd."\n";
      } elsif ($orientation eq 'antisense') {
         $boundary_gene = '     gene            complement('.$left_bd.'..'.$right_bd."\)\n";
         $boundary_CDS  = '     CDS             complement('.$left_bd.'..'.$right_bd."\)\n";
      }

      #create product entry
      if ($best_hit eq 'unknown') {
         $product = '                     /product="unknown'."_$id".'"'."\n";
      } else {
         'reset' =~ m/reset/;
         if ( $best_hit =~ m/[^\|]*?\|[^\|]*?\|\s*/) {
            $best_hit =~ m/[^\|]*?\|[^\|]*?\|\s*(.+)/;
            $value = $1;
         } else {
            $best_hit =~ m/^\S+\s+(.+)/;
            $value = $1;
         }
         unless (defined $value) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Error while parsing for product annotation in $best_hit.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing Blast results for file $filename".
                           "\nError while parsing for product annotation in $best_hit.\n\n";
            close ERRORLOG;
            $value = $best_hit;
         }
         $value =~ s/\~+//g;
         $value = '/product="'.$value.'"'."\n";
         open  ENTRY, '>', \$product;
         write ENTRY;
         close ENTRY;
      }

      #create locus tag
      $locus_tag = '                     /locus_tag="'.$id.'"'."\n";

      #create gene tag
      if ($best_hit eq 'unknown') {
         $gene = '                     /gene="unknown'."_$id".'"'."\n";
      } elsif ($evalue > ${$args{ini_ref}}{Blast_cut} && $best_hit ne 'unknown') {
         $gene = '                     /gene="conserved hypothetical'."_$id".'"'."\n";
      } else {
         'reset' =~ /reset/;
         ($value) = &parse_entry(main_window  => $args{main_window},
                                 progress_bar => $args{progress_bar},
                                 auto_ini_ref => $args{auto_ini_ref},
                                 ini_ref      => $args{ini_ref},
                                 entry        => $best_hit
                                );
         $value =~ s/\~+//g;
         $value = '/gene="'.$value."_$id".'"'."\n";
         open  ENTRY, '>', \$gene;
         write ENTRY;
         close ENTRY;
      }

      #create codon start and translation table tags
      'reset' =~ m/reset/;
      ${$args{auto_ini_ref}}{translation_table} =~ m/\s*(\d+)\:/;
      my $translation_code = $1;
      unless(defined $translation_code) {$translation_code = 1};
      $codon_table = '                     /codon_start=1'."\n".
                     '                     /transl_table='.$translation_code."\n";

      #create note tags
      if ($best_hit eq 'unknown') {
         $note = '                     /note="unknown"'."\n";
      } else {
         'reset' =~ /reset/;
         $best_hit =~ m/(.*?)Length\=\d+/;
         $value = $1;
         unless (defined $value) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Error while parsing for note annotation in $best_hit.",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error parsing Blast results for file $filename".
                           "\nError while parsing for note annotation in $best_hit.\n\n";
            close ERRORLOG;
            $value = $best_hit;
         }
         $value =~ s/[\~\s]+$//g;
         $value =~ s/\s+/ /g;
         $value = '/note="'.$value.'"'."\n";
         open  ENTRY, '>', \$note;
         write ENTRY;
         close ENTRY;
      }

      #create translation tag
      $value = '/translation="'.$aa_seq.'"'."\n";
      open  ENTRY, '>', \$translation;
      write ENTRY;
      close ENTRY;

      #combine feature
      $feature_CDS =  $boundary_CDS.
                      $product.
                      $gene.
                      $locus_tag.
                      $codon_table.
                      $note.
                      $translation;
      $feature_gene = $boundary_gene.
                      $gene.
                      $locus_tag;
      $key_gene     = $left_bd.'_'.$right_bd.'_gene_'.$id; #maintain original uniqe ID instead of ORFnumber
      $key_CDS      = $left_bd.'_'.$right_bd.'_CDS_'.$id;

      ${$args{genbank_ref}}{$key_gene} = $feature_gene;
      ${$args{genbank_ref}}{$key_CDS}  = $feature_CDS;
      push (@{$args{feature_list_ref}}, $key_gene);
      push (@{$args{feature_list_ref}}, $key_CDS);

   }
   &hide_pbar_2;
   return;
}

sub gene_model_parser {
   my %args = @_;
   my (@list, $seq_ref);

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'Compiling Gene model results',
                   label        => ''
                  );
   &show_pbar_2;

   #slurp up sequence
   $seq_ref = &slurp(main_window  => $args{main_window},
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

   foreach (my $i = 1; $i < $args{counter}; $i++) {
      my ($filename, $id, $left_bd, $right_bd, $orientation, $boundary_gene,
          $boundary_CDS, $locus_tag, $product, $gene, $codon_table, $translation,
          $feature_CDS, $feature_gene, $key_gene, $key_CDS,
          $query_seq_nt,  $query_seq_aa);

      #get filename, left_bd, right_bd, orientation
      'reset' =~ m/reset/;
      $args{orf_to_id}->{$i} =~ m/^(.*?)___(\d+)___(\d+)___(\d+)___(sense|antisense)/;
      $filename    = $1;
      $id          = $2;
      $left_bd     = $3;
      $right_bd    = $4;
      $orientation = $5;
      unless (defined $orientation) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not parse from entry ${$args{orf_to_id}}->{$i}",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing global Blast results hash".
                        "\nCould not parse from entry ${$args{orf_to_id}}->{$i}\n\n";
         close ERRORLOG;
         &hide_pbar_2;
         return (0);
      }

      #update status bar
      if (($i % 20) == 0) {
         &update_pbar_2(title        => 'Compiling Gene model results',
                        label        => "Compiling Gene model results for $filename",
                        progress     => ($i / $args{counter}) * 100,
                       );
      }

      #create boundary entry
      if ($orientation eq 'sense') {
         $boundary_gene = '     gene            '.$left_bd.'..'.$right_bd."\n";
         $boundary_CDS  = '     CDS             '.$left_bd.'..'.$right_bd."\n";
      } elsif ($orientation eq 'antisense') {
         $boundary_gene = '     gene            complement('.$left_bd.'..'.$right_bd."\)\n";
         $boundary_CDS  = '     CDS             complement('.$left_bd.'..'.$right_bd."\)\n";
      }

      #create product entry
      $product = '                     /product="unknown'."_$i".'"'."\n";

      #extract nt and aa sequence from ORF boundaries
      $query_seq_nt = substr(${$seq_ref}, ($left_bd - 1), ($right_bd - $left_bd + 1));
      $query_seq_aa = nt2aa_unchecked(main_window  => $args{main_window},
                                      progress_bar => $args{progress_bar},
                                      auto_ini_ref => $args{auto_ini_ref},
                                      ini_ref      => $args{ini_ref},
                                      module       => 'Blast Parser',
                                      orientation  => $orientation,
                                      sequence     => $query_seq_nt,
                                      filename     => $args{input_file},
                                      left_bd      => $left_bd,
                                      right_bd     => $right_bd
                                     );
      if ($query_seq_aa eq '0') {
         &hide_pbar_2;
         $query_seq_aa = '';
      }

      #create translation tag
      $value = '/translation="'.$query_seq_aa.'"'."\n";
      open  ENTRY, '>', \$translation;
      write ENTRY;
      close ENTRY;

      #create locus tag
      $locus_tag = '                     /locus_tag="'.$i.'"'."\n";

      #create gene tag
      $gene = '                     /gene="unknown'."_$i".'"'."\n";

      #create codon start and translation table tags
      'reset' =~ m/reset/;
      ${$args{auto_ini_ref}}{translation_table} =~ m/\s*(\d+)\:/;
      my $translation_code = $1;
      unless(defined $translation_code) {$translation_code = 1};
      $codon_table = '                     /codon_start=1'."\n".
                     '                     /transl_table='.$translation_code."\n";

      #combine feature
      $feature_CDS =  $boundary_CDS.
                      $product.
                      $gene.
                      $locus_tag.
                      $codon_table.
                      $translation;
      $feature_gene = $boundary_gene.
                      $gene.
                      $locus_tag;
      $key_gene     = $left_bd.'_'.$right_bd.'_gene_'.$id; #maintain original uniqe ID instead of ORFnumber
      $key_CDS      = $left_bd.'_'.$right_bd.'_CDS_'.$id;

      ${$args{genbank_ref}}{$key_gene} = $feature_gene;
      ${$args{genbank_ref}}{$key_CDS}  = $feature_CDS;
      push (@{$args{feature_list_ref}}, $key_gene);
      push (@{$args{feature_list_ref}}, $key_CDS);

   }
   &hide_pbar_2;
   return;
}

sub parse_entry {
   my %args = @_;
   my ($entry);

   if ($args{entry} eq 'unknown') {
      return ($args{entry});
   }

   #support swisspprot and regular genbank
   if ($args{entry} =~ /\|[PQOA]\d/ || $args{entry} =~ /^sp\|/) { #swissprot
      'reset' =~ /reset/;
      if ($args{entry} =~ m/[^\|]*?\|[^\|]*?\|\s*\S+\s+/) {
         $args{entry} =~ m/[^\|]*?\|[^\|]*?\|\s*\S+\s+(.*?)Length\=\d+/;
         $entry = $1;
      } else {
         $args{entry} =~ m/^\S+\s+(.*?)Length\=\d+/;
         $entry = $1;
      }
      unless (defined $entry) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while parsing for gene annotation in $args{entry}.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing Blast entry".
                        "\nError while parsing for gene annotation in $args{entry}.\n\n";
         close ERRORLOG;
      }
      #clean up
      $entry =~ s/\s+$//g;
      $entry =~ s/~+//g;
      $entry =~ s/ - .+$//;
   } else { #assume Genbank
      'reset' =~ /reset/;
      if ($args{entry} =~ m/[^\|]*?\|[^\|]*?\|\s*/) {
         $args{entry} =~ m/[^\|]*?\|[^\|]*?\|\s*(.*?)Length\=\d+/;
         $entry = $1;
      } else {
         $args{entry} =~ m/^\S+\s+(.*?)Length\=\d+/;
         $entry = $1;
      }
      unless (defined $entry) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Error while parsing for gene annotation in $args{entry}.",
                                                       -buttons => ['OK'],
                                                       -bitmap  => 'info');
         $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
         $error_msg-> Show();
         open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
         print ERRORLOG "Error parsing Blast entry".
                        "\nError while parsing for gene annotation in $args{entry}.\n\n";
         close ERRORLOG;
         return ($args{entry});
      }
      #clean up
      $entry =~ s/\s+$//g;
      $entry =~ s/~+//g;
      $entry =~ s/\s\[[^\]]*?\]$//;
   }
   return ($entry);
}

sub blast_summary {
   my %args = @_;
   my ($max_count, $counter);
   my $date = &set_date;
   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Creating Blast summary",
                   label        => ''
                  );
   &show_pbar_2;

   #open Blast Summary File
   open WRITE, "+>${$args{ini_ref}}{blast_results}".'/Blast_Summary_'.$date.'.txt' or do {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Error creating Blast summary file.",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
      print ERRORLOG "Error processing Blast entries".
                     "\nError creating Blast summary file.\n\n";
      close ERRORLOG;
      &hide_pbar_2;
      return;
   };

   #get all Blast results in folder for each input file and process
   $max_count = $#{$args{input_list_ref}} + 1;
   $counter = 0;
   foreach my $input_file (@{$args{input_list_ref}}) {
      my (@inputfile, @sorted);
      ${$args{progress_bar}}->configure(-label=>" Processing results for $input_file");
      ${$args{main_window}}->update;
      $counter++;
      &update_pbar_2(title        => "Creating Blast summary",
                     label        => "Processing results for $input_file",
                     progress     => ($counter / $max_count) * 100,
                    );

      #grab all blast result files
      opendir RESULTS, ${$args{ini_ref}}{blast_results} or do {return};
      @inputfile = grep /^$input_file.*?_\d+$/, readdir(RESULTS);
      closedir RESULTS;

      #sort list by ORF number
      @sorted =
         map  $_->[0] =>
         sort { $a->[1] <=> $b->[1] }
         map  [ $_, m/_(\d+)$/ ]
         => @inputfile;
      undef @inputfile;

      #open each file and process
      foreach my $file (@sorted) {
         my ($Blastresult, $DNAseq, $aaseq, $genemodel, $shortlist);
         {
            local( $/, *GB_ENTRY ) ;
            open READ, "<${$args{ini_ref}}{blast_results}".'/'.$file or do {
               my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                             -text    => "Error reading Blast results file $file.",
                                                             -buttons => ['OK'],
                                                             -bitmap  => 'info');
               $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
               $error_msg-> Show();
               open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
               print ERRORLOG "Error processing Blast entries".
                              "\nError reading Blast results file $file.\n\n";
               close ERRORLOG;
               &hide_pbar_2;
               close WRITE;
               return;
            };
            $Blastresult = <READ>;
            close READ;
         }
         $Blastresult =~ m/(DNA sequence:\s+\S+)\s/;
         $DNAseq = $1;
         $Blastresult =~ m/(Deduced aminoacid sequence:\s+\S+)\s/;
         $aaseq = $1;
         $Blastresult =~ m/(Gene model summary:\s+.+sense)/;
         $genemodel = $1;
         if ($Blastresult =~ m/No Blast Hits found/) {
            $shortlist = 'No Blast Hits found';
         } else {
            $Blastresult =~ s/Complete list of Blast results.*//s;
            $Blastresult =~ s/.*Blast overview\s+//s;
            my @temp = split /\>/, $Blastresult;
            @temp = splice(@temp,0,${$args{auto_ini_ref}}{blast_list} + 1);
            $shortlist = "@temp";
            undef @temp;
         }
         print WRITE "Input File:\t$input_file\n".
                      $genemodel."\n".
                      $DNAseq."\n".
                      $aaseq."\n".
                      $shortlist."\n----------------------------------------\n";
         undef $genemodel;
         undef $DNAseq;
         undef $aaseq;
         undef $shortlist;
      }
   }
   close WRITE;
   &hide_pbar_2;
   return;
}

1;

