#!/opt/ActivePerl-5.8/bin/perl

#Vector sequence result parser
#input arguments: main_window, progress_bar, gb_file, gene_model
#list vector_hit_number results



package CompilerModules::vector_parser;
use strict;
use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);

use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&vector_parser);
use vars qw();


#local variables
my (%args, %seen, $value);

format ENTRY =
~~                   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$value
.

sub vector_parser {
   my %args = @_;
   my $max_number        = 0;
   my $counter           = 1;
   my @hit_summary       = ();
   my @result_files      = ();
   my @vector_summary    = ();
   my %vector_features   = ();
   my %vector_keys       = ();



   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling Vector Screen results",
                   label        => ''
                  );
   &show_pbar_2;

   #grab all results files for input file
   opendir INPUT, ${$args{ini_ref}}{vector_results};
   @result_files = grep /^$args{filename}/, readdir (INPUT);
   closedir INPUT;
   $max_number = @result_files;

   #return if no results
   if ($#result_files < 0) {
      &hide_pbar_2;
      return (1);
   }

   #iterate through each result file;
   foreach my $entry (@result_files) {
      my $contig_left_bd   = '';
      my $right_bd         = '';
      my $contig_right_bd  = '';
      my $score            = '';
      my $descriptor       = '';
      my $feature_boundary = '';
      my $feature_product  = '';
      my $feature_note     = '';
      my $feature_vector   = '';
      my $key_vector       = '';
      my @summary          = ();
      my $colour           = '';
      my $vector_score     = 0;
      my $hit_qual         = '';
      #slurp up sequence
      my $seq_ref = &slurp(main_window  => $args{main_window},
                           progress_bar => $args{progress_bar},
                           auto_ini_ref => $args{auto_ini_ref},
                           ini_ref      => $args{ini_ref},
                           directory    => ${$args{ini_ref}}{vector_results},
                           filename     => $entry
                         );
      if ($seq_ref eq '0') {
         &hide_pbar_2;
         return (0);
      }

      #skip entry if no hits
      next if ($$seq_ref =~ m/\*\*\*\*\* No hits found \*\*\*\*\*/);

      #get left boundary
      'reset'   =~ m/reset/;
      $$seq_ref =~ m/Query\=.+?\s+Left_bd\:\s+(\d+)\s+Right_bd\:\s+(\d+)/s;
      ($contig_left_bd, $contig_right_bd) = ($1, $2);

      #split result file into individual results blocks
      my @hits = split /\n>/, $$seq_ref;
      #remove header
      shift @hits;

      #iterate through all hits
      foreach my $vectorhit (@hits) {
         my $seq_align_length = 0;

         #grab descriptor
         'reset'     =~ m/reset/;
         $vectorhit  =~ m/\S+?\s+(.+?)Length\s*=/s;
         $descriptor = $1;
         $descriptor =~ s/\s+/ /gs;

         #grab score
         'reset'    =~ m/reset/;
         $vectorhit =~ m/Score\s*=\s+(\S+)\s+bits/s;
         $score     = $1;

         unless (defined $score) {
            my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                          -text    => "Could not parse Vector score from $vectorhit\nSkipping entry",
                                                          -buttons => ['OK'],
                                                          -bitmap  => 'info');
            $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
            $error_msg-> Show();
            open  ERRORLOG, ">>${$args{auto_ini_ref}}{work_dir}/Error.log";
            print ERRORLOG "Error in Vector Screen:".
                           "\nCould not parse Vector score from $vectorhit in file $entry\n\n";
            close ERRORLOG;
            next;
         }

         #skip if score is below threshold for weak hits
         next if ($score < 16);

         #get individual alignment blocks
         my @blocks = split/\n\s+Score =\s+/, $vectorhit;

         #parse through blocks and build final boundaries
         foreach my $block (@blocks) {
            next unless ($block =~ m/Query/);
            'reset' =~ m/reset/;
            $block  =~ m/^(\S+)\s.*?\nQuery\:?\s+(\d+)\s/s;
            my ($local_score, $local_left_bd) = ($1, $2);
            $block =~ m/.+\nQuery\:?\s+\d+\s.+?(\d+)/gs;
            my ($alignment_length) = $1;
            $alignment_length = int($alignment_length - $local_left_bd + 1);
            push (@summary, $local_score.'_'.$local_left_bd.'_'.$alignment_length);

         }
         undef @blocks;

         #multiple blocks? build joined feature 'Vector_match'
         if ($#summary > 0) {
            $feature_boundary = 'join(';
            my %seen = ();
            my @blocks = ();
            foreach my $block (@summary) {
               my ($block_lb, $block_rb);

               'reset' =~ m/reset/;
               $block =~ m/^([\d\.]+)_(\d+)_(\d+)/;
               my ($local_score, $local_left_bd, $alignment_length) = ($1, $2, $3);

               $block_lb = ($contig_left_bd + $local_left_bd - 1);
               $block_rb = ($contig_left_bd + $local_left_bd + $alignment_length - 1);

               #skip repeated blocks
               next if (exists $seen{($contig_left_bd + $local_left_bd - 1).'_'.($contig_left_bd + $local_left_bd + $alignment_length - 1)});
               $seen{($contig_left_bd + $local_left_bd - 1).'_'.($contig_left_bd + $local_left_bd + $alignment_length - 1)} = 1;

               push (@blocks, $block_lb.'_'.$block_rb.'_'.$local_score.'_'.$alignment_length);
            }

            #sort blocks by left and rigth boundaries
            my @blocks_sorted =
               map  $_->[0] =>
               sort { $a->[1] <=> $b->[1] ||
                      $a->[2] <=> $b->[2]}
               map  [ $_, m/^(\d+)_/, m/^\d+_(\d+)_/ ]
               => @blocks;

            #check if total length of the hit exceeds the maximum length of 10,000 nt.
            #if yes, only use the first block, discard the rest.
            'reset' =~ m/reset/;
            my ($first_lb, $last_rb);
            $blocks_sorted[0] =~ m/^(\d+)_/;
            $first_lb = $1;
            $blocks_sorted[-1] =~ m/^\d+_(\d+)_/;
            $last_rb = $1;
            if (($last_rb - $first_lb) > 10000) {
              splice (@blocks_sorted, 1);
            }

            #build
            foreach my $entry (@blocks_sorted) {
               $entry =~ m/^(\d+)_(\d+)_([\d\.]+)_(.+)/;
               my ($local_left_bd, $local_right_bd, $local_score, $alignment_length) = ($1, $2, $3, $4);

               $feature_boundary .= $local_left_bd.'..'.$local_right_bd.', ';
               $feature_note     .= 'Section Score: '.$local_score.', ';
               $right_bd          = $contig_left_bd + $local_left_bd + $alignment_length - 1;
               $vector_score     += $local_score;
               $seq_align_length += $alignment_length;
            }

            $feature_boundary  =~ s/\,\s*$//;
            $feature_note      =~ s/\,\s*$//;
            $feature_boundary .= "\)\n";
            $feature_note      = "Total Score: $vector_score; ".$feature_note;
            $feature_note     .= "\n";
            $vector_score      = int ($vector_score);
            undef @blocks;
            undef @blocks_sorted;
            undef %seen;

         }
         #if not, create single feature
         else {
            'reset' =~ m/reset/;
            $summary[0] =~ m/^([\d\.]+)_(\d+)_(\d+)/;
            my ($local_score, $local_left_bd, $alignment_length) = ($1, $2, $3);
            $feature_boundary = ($contig_left_bd + $local_left_bd - 1).'..'.($contig_left_bd + $local_left_bd + $alignment_length - 1)."\n";
            $feature_note     = 'Section_Score: '.$local_score;
            $right_bd         = $contig_left_bd + $local_left_bd + $alignment_length - 1;
            $vector_score     = int ($local_score);
            $seq_align_length = $alignment_length;
         }

         #test for internal/terminal location and strength of hit
         'reset' =~ m/reset/;
         $feature_boundary =~ m/^.*?(\d+).*\.\.(\d+)\)?\n/;
         my ($feature_lb, $feature_rb) = ($1, $2);

         #calculate ratio, minimum 1
         my $ratio = (($contig_right_bd - $contig_left_bd) / 350000);
         if ($ratio < 1) {$ratio = 1};

         if (($feature_lb - $contig_left_bd) <= 25 || ($contig_right_bd - $feature_rb) <= 25) { #terminal hit
            if (($vector_score / $ratio)  >= 24) {
               $hit_qual = 'Strong Match';
            } elsif (24 > ($vector_score / $ratio) && ($vector_score / $ratio) >= 19) {
               $hit_qual = 'Moderate Match';
            } elsif (19 > ($vector_score / $ratio) && ($vector_score / $ratio) >= 16) {
               $hit_qual = 'Weak Match';
            } else {
               $vector_score     = 0;
               $feature_product  = '';
               $feature_boundary = '';
               $feature_note     = '';
               $feature_vector   = '';
               $key_vector       = '';
               $value            = '';
               @summary          = ();
               next;
            }
         } else { #internal match
            if (($vector_score / $ratio) >= 30) {
               $hit_qual = 'Strong Match';
            } elsif (30 > ($vector_score / $ratio) && ($vector_score / $ratio) >= 25) {
               $hit_qual = 'Moderate Match';
            } elsif (25 > ($vector_score / $ratio) && ($vector_score / $ratio) >= 23) {
               $hit_qual = 'Weak Match';
            } else {
               $vector_score     = 0;
               $feature_product  = '';
               $feature_boundary = '';
               $feature_note     = '';
               $feature_vector   = '';
               $key_vector       = '';
               $value            = '';
               @summary          = ();
               next;
            }
         }

         #create note
         $value = '/note="'.$feature_note.'; Hit quality: '.$hit_qual.'"'."\n";

         open  ENTRY, '>', \$feature_note;
         write ENTRY;
         close ENTRY;
         #create product
         $value = '/product="'.$descriptor.'"'."\n";
         open  ENTRY, '>', \$feature_product;
         write ENTRY;
         close ENTRY;
         #create colour
         $colour = "                     \/colour\=17\n";
         #combine feature
         $feature_vector =  '     Vector_match    '.$feature_boundary.
                            $feature_product.
                            $feature_note.
                            $colour;
         $key_vector      = $feature_lb.'_'.$feature_rb.'_Vector_match_'.$counter;

         #increase internal counter, will later be replaed with global feature counter $args{counter}
         $counter++;

         #put everyting safely away in array and hash to sort and retrieve later
         push (@vector_summary, $vector_score.'___'.$seq_align_length.'___'.$key_vector);
         $vector_keys{$vector_score.'___'.$seq_align_length.'___'.$key_vector}     = $key_vector;
         $vector_features{$vector_score.'___'.$seq_align_length.'___'.$key_vector} = $feature_vector;

         $vector_score     = 0;
         $feature_product  = '';
         $feature_boundary = '';
         $feature_note     = '';
         $feature_vector   = '';
         $key_vector       = '';
         $value            = '';
         @summary          = ();

      }
   }

   #sort through all results: vector_score -> alignment_length -> left_bd
   my @sorted =
      map  $_->[0] =>
      sort { $b->[1] <=> $a->[1] ||
             $b->[2] <=> $a->[2] ||
             $a->[3] <=> $b->[3]}
      map  [ $_, m/^(\d+)_/, m/^\d+___(\d+)_/, m/^\d+___\d+___(\d+)_/ ]
      => @vector_summary;

   #remove duplicate hits with same start-stop boundary or same lb or rb (only keep best hit)
   my %seen_lb_rb   = ();
   my %seen_lb      = ();
   my %seen_rb      = ();
   my @nonredundant = ();
   foreach my $entry (@sorted) {
      'reset' =~ m/reset/;
      $entry =~ m/^\d+___\d+___(\d+)_(\d+)_/;
      next if (exists $seen_lb_rb{$1.'_'.$2});
      next if (exists $seen_lb   {$1});
      next if (exists $seen_rb   {$2});
      $seen_lb_rb{$1.'_'.$2} = 1;
      $seen_lb   {$1}        = 1;
      $seen_rb   {$2}        = 1;
      push (@nonredundant, $entry);
   }

   #shrink sorted nonredundant array to vector_hit_number length
   splice (@nonredundant,${$args{auto_ini_ref}}{vector_hit_number});

   #create keys and hashes for features
   foreach my $entry (@nonredundant) {
      ${$args{genbank_ref}}{$vector_keys{$entry}} = $vector_features{$entry};
      push (@{$args{feature_list_ref}}, $vector_keys{$entry});
   }

   undef %seen;
   undef %vector_keys;
   undef %vector_features;
   undef @vector_summary;
   undef @nonredundant;

   &hide_pbar_2;
   return (1);

}

1;

