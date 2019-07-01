#!/opt/ActivePerl-5.8/bin/perl
#transfer annotation from one GAMOLA genbank file
#to new one
#changes are logged in separate file

package Miscellaneous::transfer_annotation;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&transfer_annotation);
use vars qw();

use Tk;
use Tk::Pane;
use initialise::read_me              qw(:DEFAULT);
use ProgrammeModules::sequence       qw(:DEFAULT);
use ProgrammeModules::genbank_parser qw(:DEFAULT);
use Basics::progress_bar             qw(:DEFAULT);
use Cwd;

#local variables
my (%args, $tl, $new_width, $old_width, $current_dir, $local_blastdb_dir,
    @new_selected_files, %new_all_files, @new_short_files,
    @old_selected_files, %old_all_files, @old_short_files,
    $transfer_gene_anno, $transfer_product_anno, $transfer_locustag_anno,
    @ORFnumbers
    );


sub transfer_annotation {
   my %args = @_;

   #remember last directory for adding files
   $current_dir = ${$args{auto_ini_ref}}{work_dir};
   $local_blastdb_dir = ${$args{ini_ref}}{blast_db_path};

   my ($progress_bar, $frame_left, $frame_right, $frame_middle, $frame_bottom,
       $old_anno_lb, $old_width, $old_add_files, $old_clear,
       $new_anno_lb, $new_width, $new_add_files, $new_clear,
       $identity_cutoff);

   #set defaults
   $old_width          = 20;
   $new_width          = 20;
   $identity_cutoff    = 90;
   $transfer_gene_anno = 1;
   $transfer_product_anno  = 0;

   if (defined $tl && Tk::Exists($tl)) {
      $tl->state('normal');
      $progress_bar = $tl->Frame(-borderwidth => 2, -relief => 'groove');
      $frame_bottom = $tl->Frame(-relief => 'groove')                   ;
      $frame_left   = $tl->Frame(-relief => 'raised')                   ;
      $frame_middle = $tl->Frame()                                      ;
      $frame_right  = $tl->Frame(-relief => 'raised')                   ;
      $progress_bar->configure(-label => "Current Blast DB path: $local_blastdb_dir"); #create dummy for space
      $tl->update();
   } else {
      $tl = ${$args{main_window}}->DialogBox(-title   => 'Transfer Annotation between Genbank files',
                                             -buttons => [ 'Exit' ]
                                            );
      $tl->minsize(20, 10);
      #$tl->maxsize(120, 20);
      $progress_bar = $tl->Frame(-borderwidth => 2, -relief => 'groove') -> pack(-side => 'bottom', -fill => 'x');
      $frame_bottom = $tl->Frame(-relief => 'groove')                    -> pack(-side => 'bottom', -fill => 'x',    -expand => 1);
      $frame_left   = $tl->Frame(-relief => 'raised')                    -> pack(-side => 'left',   -fill => 'both', -expand => 1);
      $frame_middle = $tl->Frame()                                       -> pack(-side => 'left');
      $frame_right  = $tl->Frame(-relief => 'raised')                    -> pack(-side => 'left',   -fill => 'both', -expand => 1);
   }

   $progress_bar->configure(-label => "Current Blast DB path: $local_blastdb_dir"); #create dummy for space
   $tl->update();

   #left frame -> old annotation
   {
      $frame_left                       -> Label (-text => "Existing Genbank annotation\n")
         ->grid(-row => 0, -column => 0, -columnspan => 2, -sticky => 'nsew');
      $old_add_files = $frame_left      ->Button(-text    => 'Add files',
                                                 -command => sub {&populate_old(main_window  => \$tl,
                                                                                auto_ini_ref => $args{auto_ini_ref},
                                                                                textbox      => $old_anno_lb,
                                                                               )
                                                                  }
                                                 )
         ->grid(-row => 1, -column => 0);

      $frame_left                       ->Label(-text => 'Double left-click entry to remove.')
         ->grid(-row => 1, -column => 1, -sticky => 'nsew');
      $old_anno_lb         = $frame_left->Scrolled("Listbox",
                                                    -setgrid    => 1,
                                                    -selectmode => 'single',
                                                    -width      => $old_width
                                                    )
         ->grid(-row => 2, -column => 0, -columnspan => 2, -sticky => 'nsew');
      $old_clear         = $frame_left ->Button(-text => 'Clear current selection',
                                               -command => sub {
                                                                 %old_all_files = ();
                                                                 @old_short_files = ();
                                                                 $old_anno_lb->delete(0, 'end');
                                                                 $tl->update;
                                                                }
                                              )
      ->grid(-row => 3, -column => 0, -columnspan => 2, -sticky => 'nsew');
      $old_anno_lb                      ->configure(-width   => $old_width,
                                                    -setgrid => 1 );
   }

   #middle frame: optics
   $frame_middle -> Label (-text => "Transfer annotation to\n    ===>")->grid(-row => 0, -column => 1,  -sticky => 'nsew');

   #right frame -> new annotation
   {
      $frame_right                       -> Label (-text => "New Genbank annotation\n")
         ->grid(-row => 0, -column => 0, -columnspan => 2, -sticky => 'nsew');
      $new_add_files = $frame_right      ->Button(-text => 'Add files',
                                                  -command => sub {&populate_new(main_window  => \$tl,
                                                                                 auto_ini_ref => $args{auto_ini_ref},
                                                                                 textbox      => $new_anno_lb,
                                                                                )
                                                                  }
                                                 )
         ->grid(-row => 1, -column => 0);

      $frame_right                       ->Label(-text => 'Double right-click entry to remove.')
         ->grid(-row => 1, -column => 1, -sticky => 'nsew');
      $new_anno_lb         = $frame_right->Scrolled("Listbox",
                                                    -setgrid    => 1,
                                                    -selectmode => 'single',
                                                    -width      => $new_width
                                                    )
         ->grid(-row => 2, -column => 0, -columnspan => 2, -sticky => 'nsew');
      $new_clear         = $frame_right  ->Button(-text    => 'Clear current selection',
                                                  -command => sub {
                                                                    %new_all_files = ();
                                                                    @new_short_files = ();
                                                                    $new_anno_lb->delete(0, 'end');
                                                                    $tl->update;
                                                                   }
                                                 )
      ->grid(-row => 3, -column => 0, -columnspan => 2, -sticky => 'nsew');
      $new_anno_lb                       ->configure(-width   => $new_width,
                                                     -setgrid => 1 );
   }

   #bottom frame: start transfer
   {
      $frame_bottom -> BrowseEntry(-label   =>'Identity cutoff [%]: ',
                                    -choices  =>[1..100],
                                    -width    => 5,
                                    -variable =>\$identity_cutoff
                                   )->grid(-row => 0, -column => 0, -sticky => 'w');;

      $frame_bottom ->Button(-text => 'Change Blast DB path',
                             -command => sub {$local_blastdb_dir= &dir_select(main_window  => \$tl,
                                                                              progress_bar => \$progress_bar,
                                                                              auto_ini_ref => $args{auto_ini_ref},
                                                                              ini_ref      => $args{ini_ref},
                                                                              directory    => $local_blastdb_dir
                                                                             );
                                              $progress_bar->configure(-label => "Current Blast DB path: $local_blastdb_dir");
                                              $tl->update;
                                             }
                            )->grid(-row => 0, -column => 1, -sticky => 'nsew');

      $frame_bottom ->Button(-text => 'Start Annotation Transfer',
                             -command => sub {&transfer(main_window            => \$tl,
                                                        progress_bar           => \$progress_bar,
                                                        ini_ref                => $args{ini_ref},
                                                        auto_ini_ref           => $args{auto_ini_ref},
                                                        old_anno_lb            => \$old_anno_lb,
                                                        new_anno_lb            => \$new_anno_lb,
                                                        identity_cutoff        => $identity_cutoff,
                                                        transfer_gene_anno     => $transfer_gene_anno,
                                                        transfer_product_anno  => $transfer_product_anno,
                                                        transfer_locustag_anno => $transfer_locustag_anno
                                                       );
                                             }
                            )->grid(-row => 0, -column => 2, -sticky => 'nsew');
      $frame_bottom -> Checkbutton(-text     => 'Transfer gene annotation',
                                   -variable => \$transfer_gene_anno
                                  )
                             ->grid(-row => 1, -column => 0, -sticky => 'nsew');
      $frame_bottom -> Checkbutton(-text     => 'Transfer product annotation',
                                   -variable => \$transfer_product_anno
                                  )
                             ->grid(-row => 1, -column => 1, -sticky => 'nsew');
      $frame_bottom -> Checkbutton(-text     => 'Transfer locus_tag annotation',
                                   -variable => \$transfer_locustag_anno
                                  )
                             ->grid(-row => 1, -column => 2, -sticky => 'nsew');

   }
   #bind left and right mouse button to selection and removal of entries
   $old_anno_lb->bind('<Double-ButtonRelease-1>' => sub {&remove_entry_old(textbox      => \$old_anno_lb,
                                                                           main_window  => \$tl,
                                                                           progress_bar => \$progress_bar,
                                                                           ini_ref      => $args{ini_ref},
                                                                           auto_ini_ref => $args{auto_ini_ref}
                                                                          )
                                                        });
   $new_anno_lb->bind('<Double-ButtonRelease-3>' => sub {&remove_entry_new(textbox      => \$new_anno_lb,
                                                                           main_window  => \$tl,
                                                                           progress_bar => \$progress_bar,
                                                                           ini_ref      => $args{ini_ref},
                                                                           auto_ini_ref => $args{auto_ini_ref}
                                                                          )
                                                        });
   my $wait = $tl->Show();
   $new_anno_lb->delete(0, 'end');
   $old_anno_lb->delete(0, 'end');
   $progress_bar->configure(-label => "Current Blast DB path: $local_blastdb_dir");
   $tl->update;
   undef @new_selected_files;
   undef @old_selected_files;
   undef %new_all_files;
   undef %old_all_files;
   undef @new_short_files;
   undef @old_short_files;
   undef $current_dir;
   if (-d ${$args{ini_ref}}{annotation_transfer_dir}.'/Blast_Results') {
      cleanup(${$args{ini_ref}}{annotation_transfer_dir}.'/Blast_Results');
   }
   if ($wait eq 'OK') {
      $tl->state('withdrawn');
   }

}

sub cleanup {
   my $dir = shift;
   local *DIR;

   opendir DIR, $dir or die "opendir $dir: $!";
   my $found = 0;
   while ($_ = readdir DIR) {
      next if /^\.{1,2}$/;
      my $path = "$dir/$_";
      unlink $path if -f $path;
      cleanup($path) if -d $path;
   }
   closedir DIR;
   rmdir $dir or print "error - $!";
}

sub transfer {
   my %args = @_;
   my (@old_list, @old_hits, @old_nohits, %old_tags, %old_header, %old_source, %old_core, %old_DNAseq, %old_start_ref, %old_feature_list_ref, %old_genbank_ref, %old_feature_counter,
       @new_list, @new_hits, @new_nohits, %new_tags, %new_header, %new_source, %new_core, %new_DNAseq, %new_start_ref, %new_feature_list_ref, %new_genbank_ref, %new_feature_counter); #hash references for genbank files
   my (%seen, @multi_hits, @old_multi_hits, @new_multi_hits, $feature_counter, $max_count, $progress,
       $status);


   #return if no files selected
   if (scalar (keys %old_all_files) < 1 || scalar (keys %new_all_files) < 1) {

      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Insufficient file selection",
                                        -icon    => 'error',
                                        -type    => 'ok');
      return;
   };

   #initialise progress bar
   &progress_bar_1(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => 'Transfer annotation',
                   label        => 'Transfer'
                  );
   &show_pbar_1;

   #cleanup Annotation transfer directory
   if (-d ${$args{ini_ref}}{annotation_transfer_dir}) {
      &cleanup(${$args{ini_ref}}{annotation_transfer_dir});
   }
   #re-create Annotation_transfer directories
   mkdir (${$args{ini_ref}}{annotation_transfer_dir}, 0777) or do {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Cannot create Annotation transfer directory ${$args{ini_ref}}{annotation_transfer_dir}",
                                        -icon    => 'error',
                                        -type    => 'ok');
      &hide_pbar_1;
      return;
   };
   mkdir (${$args{ini_ref}}{annotation_transfer_dir}.'/Blast_Results', 0777) or do {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Cannot create Blast directory for Annotation transfer",
                                        -icon    => 'error',
                                        -type    => 'ok');
      &hide_pbar_1;
      return;
   };

   #create databases and generate feature hashes for each selected genbank file
   #hashes wil hold entries for all files
   $status = &old_database(main_window          => $args{main_window},
                           progress_bar         => $args{progress_bar},
                           auto_ini_ref         => $args{auto_ini_ref},
                           ini_ref              => $args{ini_ref},
                           old_header           => \%old_header,
                           old_source           => \%old_source,
                           old_core             => \%old_core,
                           old_DNAseq           => \%old_DNAseq,
                           old_start_ref        => \%old_start_ref,
                           old_feature_list_ref => \%old_feature_list_ref,
                           old_genbank_ref      => \%old_genbank_ref,
                           old_feature_counter  => \%old_feature_counter,
                           old_list             => \@old_list,
                           old_tags             => \%old_tags
                         );
   #error in generating databases?
   if ($status == 0) {
      &hide_pbar_1;
      return (0);
   }

   $status = &new_database(main_window          => $args{main_window},
                           progress_bar         => $args{progress_bar},
                           auto_ini_ref         => $args{auto_ini_ref},
                           ini_ref              => $args{ini_ref},
                           new_header           => \%new_header,
                           new_source           => \%new_source,
                           new_core             => \%new_core,
                           new_DNAseq           => \%new_DNAseq,
                           new_start_ref        => \%new_start_ref,
                           new_feature_list_ref => \%new_feature_list_ref,
                           new_genbank_ref      => \%new_genbank_ref,
                           new_feature_counter  => \%new_feature_counter,
                           new_list             => \@new_list,
                           new_tags             => \%new_tags,
                           );

   #error in generating databases?
   if ($status == 0) {
      &hide_pbar_1;
      return (0);
   }

   #iterate through old annotation and try finding corresponding ORFs in new genbank files
   $max_count = $#old_list + 1;
   $progress = 1;
   foreach my $old_entry (@old_list) {
      my ($old_aa, $old_file, @new_subset, $blastresult,
          $multi_gene_anno, $multi_product_anno, $multi_locus_tag);
      my ($old_gene_annotation, $old_product_annotation, $old_gene_feature, $old_locus_tag, $old_ORF_number);
      'reset'    =~ m/reset/;
      $old_entry =~ m/(.*?)_\d+_(.+)/;
      $old_aa    =  $1;
      $old_file  =  $2;
      unless (defined $old_file) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not parse from 'OLD' entry $old_entry\.",
                                           -icon    => 'error',
                                           -type    => 'OK');
         &hide_pbar_1;
         return (0);
      }
      $old_aa =~ s/\*//gs; #get rid of stop codon
      $old_aa =~ s/^.//;   #get rid of start codon

      &update_pbar_1(title       => 'Transfer annotation',
                     label       => "Transfer existing annotation from $old_file",
                     progress    => ($progress / $max_count) * 100
                    );
      $progress++;

      #skip hits with no gene or product qualifiers
      next if ($old_genbank_ref{$old_file}->{$old_tags{$old_entry}} !~ /\/(gene|product)\=/s);

      #try and find unchanged ORFs in new Genabnk files
      @new_subset = grep (/$old_aa/, @new_list);

      #no hits? try Blasting
      if ($#new_subset < 0) {
         #change to Blast directory
         my $curdir = getcwd();
         chdir ${$args{ini_ref}}{blast_executables};

         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               exec "./blastall -a ${$args{auto_ini_ref}}{CPU} -p blastp -F F -d $local_blastdb_dir/new.aa"; #child
               die "Cannot exec: $!";
            } else { # grandchild
               print $old_aa;
               CORE::exit;
            }
         }

         #change back to previous working dir
         chdir $curdir;

         #grab best hit
         'reset' =~ m/reset/;
         $blastresult =~ m/\>gb\|\d+\|\s*([^_]*?)_(\d+).*?Identities =.*?\((\d+)\%\).*?Query/s;
         my $best_seq = $1;
         my $query_id = $2;
         my $identity = $3;
         $best_seq    =~ s/\*//;

         if (defined $best_seq && defined $query_id) {
            open WRITE, "+>${$args{ini_ref}}{annotation_transfer_dir}\/Blast_Results\/$best_seq\_$query_id\.txt" or die;
            print WRITE $blastresult;
            close WRITE;
         }

         if (!defined $identity && $blastresult !~ /\*\*\*\*\* No hits found \*\*\*\*\*\*/) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Error in parsing Blast results for 'old' entry $old_tags{$old_entry}\n\n$blastresult\n",
                                              -icon    => 'error',
                                              -type    => 'ok');
            &hide_pbar_1;
            return;
         }

         #test if identity is bigger then identity cutoff
         if (defined $identity && $identity >= $args{identity_cutoff}) {
            @new_subset = grep (/$best_seq/, @new_list);
         }
      }

      #still no hits? enter into old negative
      if ($#new_subset < 0) {
         push (@old_nohits, $old_file .'___'.$old_tags{$old_entry});
         next;
      }
      #if NEW hits are found, record for later - will be eliminated for reciprocal search
      else {
         @new_hits = (@new_hits, @new_subset);
      }

      #capture old gene qualifier annotation if selected
      'reset' =~ m/reset/;
      $old_genbank_ref{$old_file}->{$old_tags{$old_entry}} =~ (/(gene|CDS).*?\//s);
      $old_gene_feature    = $1;

      #capture old ORF number
      'reset' =~ m/reset/;
      $old_genbank_ref{$old_file}->{$old_tags{$old_entry}} =~ (/$old_gene_feature.*?\/gene\=\".*_(\d+)\s*\"/s);
      $old_ORF_number      = $1;


      if ($args{transfer_gene_anno} == 1 && $old_genbank_ref{$old_file}->{$old_tags{$old_entry}} =~ /(gene|CDS).*?\/gene\=\"([^\"]*?)\"/s) {
         'reset' =~ /reset/;
         $old_genbank_ref{$old_file}->{$old_tags{$old_entry}} =~ /(gene|CDS).*?\/gene\=\"([^\"]*?)\"/s;
         $old_gene_annotation = $2;

         unless (defined $old_gene_annotation) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Error in parsing gene annotation from entry ->\n$old_file\n$old_genbank_ref{$old_file}->{$old_tags{$old_entry}}",
                                              -icon    => 'error',
                                              -type    => 'ok');
            &hide_pbar_1;
            return;
         }
         $multi_gene_anno     = $old_gene_annotation;
         $old_gene_annotation =~ s/\s+/ /gs;
         $old_gene_annotation =~ s/_\d+\s*$//;
         $old_gene_annotation =~ s/\>\>//;
         $old_gene_annotation =~ s/\<\<//;
      }

      #capture old product qualifier annotation if selected
      if ($args{transfer_product_anno} == 1 && $old_genbank_ref{$old_file}->{$old_tags{$old_entry}} =~ /\/product\=\"([^\"]*?)\"/s) {
         'reset' =~ /reset/;
         $old_genbank_ref{$old_file}->{$old_tags{$old_entry}} =~ /(gene|CDS).*?\/product\=\"([^\"]*?)\"/s;
         $old_product_annotation = $2;

         unless (defined $old_product_annotation) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Error in parsing product annotation from entry ->\n$old_file\n$old_genbank_ref{$old_file}->{$old_tags{$old_entry}}",
                                              -icon    => 'error',
                                              -type    => 'ok');
            &hide_pbar_1;
            return;
         }
         $multi_product_anno = $old_product_annotation;
         $old_product_annotation =~ s/\s+/ /gs;
         $old_product_annotation =~ s/_\d+\s*$//;
         $old_product_annotation =~ s/\>\>//;
         $old_product_annotation =~ s/\<\<//;
      }

      #capture old locus tag  if selected
      if ($args{transfer_locustag_anno} == 1 && $old_genbank_ref{$old_file}->{$old_tags{$old_entry}} =~ /\/locus_tag\=\"([^\"]*?)\"/s) {
         'reset' =~ /reset/;
         $old_genbank_ref{$old_file}->{$old_tags{$old_entry}} =~ /(gene|CDS).*?\/locus_tag\=\"([^\"]*?)\"/s;
         $old_locus_tag = $2;

         unless (defined $old_locus_tag) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Error in parsing locus tag from entry ->\n$old_file\n$old_genbank_ref{$old_file}->{$old_tags{$old_entry}}",
                                              -icon    => 'error',
                                              -type    => 'ok');
            &hide_pbar_1;
            return;
         }
         $multi_locus_tag = $old_locus_tag;
         $old_locus_tag =~ s/\s+/ /gs;
      }

      #transfer annotation from old to new
      $feature_counter    = 0;
      $multi_gene_anno    =~ s/\s+/ /gs;
      $multi_product_anno =~ s/\s+/ /gs;
      $multi_locus_tag =~ s/\s+/ /gs;

      foreach my $new_entry (@new_subset) {
         my ($new_file);
         #defined respective file
         'reset'    =~ m/reset/;
         $new_entry =~ m/.*?_\d+_(.+)/;
         $new_file  =  $1;
         unless (defined $new_file) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Error in parsing annotation from new entry ->\n$new_file\n$new_entry",
                                              -icon    => 'error',
                                              -type    => 'ok');
            &hide_pbar_1;
            return;
         }

         #skip if wrong feature
         next if ($new_genbank_ref{$new_file}->{$new_tags{$new_entry}} !~ /^\s*$old_gene_feature\s+/);

         #add in case of multiple hits
         $feature_counter++;
         #push (@multi_hits, $old_file.'___'.$old_tags{$old_entry}.'___'.$multi_anno);
         if (defined $multi_gene_anno && $multi_gene_anno =~ /\w+/) {
            push (@multi_hits, $new_file.'___'.$new_tags{$new_entry}.'___'.$multi_gene_anno);
         }
         if (defined $multi_product_anno && $multi_product_anno =~ /\w+/) {
            push (@multi_hits, $new_file.'___'.$new_tags{$new_entry}.'___'.$multi_product_anno);
         }
         if (defined $multi_locus_tag && $multi_locus_tag =~ /\w+/) {
            push (@multi_hits, $new_file.'___'.$new_tags{$new_entry}.'___'.$multi_locus_tag);
         }

         #transfer annotation
         'reset' =~ m/reset/;
         if (defined $old_gene_annotation && $old_gene_annotation =~ /\w+/)  {
            $new_genbank_ref{$new_file}->{$new_tags{$new_entry}}  =~ m/\/gene\=\".*_(\d+)\s*\"/s;
            my $new_ORFnumber = $1;
            #add locus tag where possible
            my $locus_tag;
            if ($new_genbank_ref{$new_file}->{$new_tags{$new_entry}} =~ m/\/locus_tag\=/) {
               $new_genbank_ref{$new_file}->{$new_tags{$new_entry}} =~ m/\/locus_tag\=\"([^\"]+)\"/s;
               $locus_tag = $1;
            } else {
               $locus_tag = '';
            }

            if (defined $new_ORFnumber) {
               push (@ORFnumbers, 'Old ORF number: '.$old_ORF_number.' New ORF number: '.$new_ORFnumber.' Locus_tag: '.$locus_tag."\n");
            } else {
               push (@ORFnumbers, 'Old ORF number: '.$old_ORF_number.' New ORF number: '.$new_ORFnumber.' Locus_tag: '.$locus_tag."\t".$new_tags{$new_entry}."\t".$new_genbank_ref{$new_file}->{$new_tags{$new_entry}}."\n");
            }
         }
         if ($args{transfer_gene_anno} == 1 && defined $old_gene_annotation && $old_gene_annotation =~ /\w+/)  {
            $new_genbank_ref{$new_file}->{$new_tags{$new_entry}}  =~ s/\/gene\=\"[^\"]*?(_\s*\d+)?\s*\"/\/gene\=\"$old_gene_annotation$1\"/s;
         }
         if ($args{transfer_product_anno} == 1 && defined $old_product_annotation && $old_product_annotation =~ /\w+/)  {
            $new_genbank_ref{$new_file}->{$new_tags{$new_entry}}  =~ s/\/product\=\"[^\"]*?(_\s*\d+)?\s*\"/\/product\=\"$old_product_annotation$1\"/s;
         }
         if ($args{transfer_locustag_anno} == 1 && defined $old_locus_tag && $old_locus_tag =~ /\w+/)  {
            if ($new_genbank_ref{$new_file}->{$new_tags{$new_entry}}  =~ /\/locus_tag/) {
               $new_genbank_ref{$new_file}->{$new_tags{$new_entry}}  =~ s/(locus_tag[^\n]+?\n)/$1                     \/note\=\"$old_locus_tag\"\n/s; #"
            } else {
               $new_genbank_ref{$new_file}->{$new_tags{$new_entry}}  =~ s/($old_gene_feature[^\n]+?\n)/$1                     \/locus_tag\=\"$old_locus_tag\"\n/s;#"
            }
         }
      }

      #multiple hits?
      if ($feature_counter > 1) {
         @old_multi_hits = (@old_multi_hits, @multi_hits);
      }
      undef @multi_hits;
   }

   #find ORFs in new Genbank files that have not been assigned and try finding them in the old annotation, in case sequence has changed slightly
   foreach my $new_entry (@new_list, @new_hits) { $seen{$new_entry}++ };
   @new_hits = ();
   foreach my $new_entry (keys %seen) {
      if ($seen{$new_entry} == 1) {push(@new_hits, $new_entry)}; #new_hits now contains unmapped orfs from new genbank files
   }

   #iterate over unmapped new ORFs
   $max_count = $#new_hits + 1;
   $progress = 1;
   foreach my $new_entry (@new_hits) {
      my ($new_aa, $new_file, $new_feature, @old_subset, $blastresult, $old_annotation);

      #parse new entry
      'reset'   =~ m/reset/;
      $new_entry    =~ m/(.*?)_\d+_(.+)/;
      $new_aa   =  $1;
      $new_file =  $2;
      unless (defined $new_file) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not parse from 'OLD' entry $new_entry\.",
                                           -icon    => 'error',
                                           -type    => 'OK');
         &hide_pbar_1;
         return (0);
      }

      $new_aa =~ s/\*//gs; #get rid of stop codon
      $new_aa =~ s/^.//;   #get rid of start codon

      &update_pbar_1(title       => 'Transfer annotation',
                     label       => "Reverse Transfer of existing annotation from $new_file",
                     progress    => ($progress / $max_count) * 100
                    );
      $progress++;

      #try and find unchanged ORFs in new Genabnk files
      @old_subset = grep (/$new_aa/, @old_list);

      #no hits? try Blasting
      if ($#old_subset < 0) {
         #change to Blast directory
         my $curdir = getcwd();
         chdir ${$args{ini_ref}}{blast_executables};

         if (open RESULT, "-|") { # original process
            local $/;
            $blastresult = <RESULT>;
         } else { # child
            if (open STDIN, "-|") { # child
               print "\nBlasting: ";
               exec "./blastall -a ${$args{auto_ini_ref}}{CPU} -p blastp -F F -d $local_blastdb_dir/old.aa"; #child
               die "Cannot exec: $!";
            } else { # grandchild
               print $new_aa;
               CORE::exit;
            }
         }

         #change back to previous working dir
         chdir $curdir;

         #grab best hit
         'reset' =~ m/reset/;
         $blastresult =~ m/\>gb\|\d+\|\s*([^_]*?)_(\d+).*?Identities =.*?\((\d+)\%\).*?Query/s;
         my $best_seq = $1;
         my $query_id = $2;
         my $identity = $3;

         if (defined $best_seq && defined $query_id) {
            open WRITE, "+>${$args{ini_ref}}{annotation_transfer_dir}\/Blast_Results\/$best_seq\_$query_id\.txt" or die;
            print WRITE $blastresult;
            close WRITE;
         }

         if (!defined $identity && $blastresult !~ /\*\*\*\*\* No hits found \*\*\*\*\*\*/) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Error in parsing Blast results for 'new' entry $new_tags{$new_entry}\n\n$blastresult\n",
                                              -icon    => 'error',
                                              -type    => 'ok');
            &hide_pbar_1;
            return;
         }

         #test if identity is bigger then identity cutoff
         if (defined $identity && $identity >= $args{identity_cutoff}) {
            @old_subset = grep (/$best_seq/, @old_list);
         }
      }

      #still no hits? enter into new negative
      if ($#old_subset < 0) {
         push (@new_nohits, $new_file.'___'.$new_tags{$new_entry});
         next;
      }

      #get new feature
      'reset' =~ m/reset/;
      $new_genbank_ref{$new_file}->{$new_tags{$new_entry}}  =~ m/\s*(gene|CDS)\s+/;
      $new_feature = $1;
      unless (defined $new_feature) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not parse from 'NEW' entry ->\n$new_genbank_ref{$new_file}->{$new_tags{$new_entry}}\.",
                                           -icon    => 'error',
                                           -type    => 'OK');
         &hide_pbar_1;
         return (0);
      }
      $new_feature =~ s/\s+/ /gs;

      #transfer annotation from old to new
      $feature_counter = 0;
      foreach my $old_entry (@old_subset) {
         my ($old_aa, $old_file);
         my ($old_gene_annotation, $old_product_annotation, $old_gene_feature, $old_product_feature, $old_locus_tag,
             $multi_gene_anno, $multi_product_anno, $multi_locus_tag);
         #get old parameters
         'reset'        =~ m/reset/;
         $old_entry =~ m/(.*?)_\d+_(.+)/;
         $old_aa        =  $1;
         $old_file      =  $2;
         unless (defined $old_file) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not parse from 'OLD' entry $old_entry\.",
                                              -icon    => 'error',
                                              -type    => 'OK');
            &hide_pbar_1;
            return (0);
         }

         #ignore if not gene feature
         #next if ($old_genbank_ref{$old_file}->{$old_tags{$old_entry}} !~ /\s*(gene|CDS).*?\/gene\=/s);

         #capture old gene annotation
         if ($args{transfer_gene_anno} == 1 && $old_genbank_ref{$old_file}->{$old_tags{$old_entry}} =~ /(gene|CDS).*?\/gene\=\"([^\"]*?)\"/s) {
            'reset' =~ /reset/;
            $old_genbank_ref{$old_file}->{$old_tags{$old_entry}} =~ /\s*(gene|CDS).*?\/gene\=\"([^\"]*?)\"/s;
            $old_gene_feature    = $1;
            $old_gene_annotation = $2;
            unless (defined $old_gene_annotation) {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Error in parsing gene annotation from entry ->\n$old_genbank_ref{$old_file}->{$old_tags{$old_subset[0]}}",
                                                 -icon    => 'error',
                                                 -type    => 'ok');
               &hide_pbar_1;
               return;
            }
            $old_gene_annotation =~ s/\s+/ /gs;
            $multi_gene_anno     = $old_gene_annotation;
            $old_gene_annotation =~ s/_\d+\s*$//;
            $old_gene_annotation =~ s/\>\>//;
            $old_gene_annotation =~ s/\<\<//;
         }

         #capture old product qualifier annotation if selected
         if ($args{transfer_product_anno} == 1 && $old_genbank_ref{$old_file}->{$old_tags{$old_entry}} =~ /(gene|CDS).*?\/product\=\"([^\"]*?)\"/s) {

            'reset' =~ /reset/;
            $old_genbank_ref{$old_file}->{$old_tags{$old_entry}} =~ /(gene|CDS).*?\/product\=\"([^\"]*?)\"/s;
            $old_product_annotation = $2;

            unless (defined $old_product_annotation) {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Error in parsing product annotation from entry ->\n$old_file\n$old_genbank_ref{$old_file}->{$old_tags{$old_entry}}",
                                                 -icon    => 'error',
                                                 -type    => 'ok');
               &hide_pbar_1;
               return;
            }
            $multi_product_anno = $old_product_annotation;
            $old_product_annotation =~ s/\s+/ /gs;
            $old_product_annotation =~ s/_\d+\s*$//;
            $old_product_annotation =~ s/\>\>//;
            $old_product_annotation =~ s/\<\<//;
         }

         #capture old locus tag if selected
         if ($args{transfer_locustag_anno} == 1 && $old_genbank_ref{$old_file}->{$old_tags{$old_entry}} =~ /(gene|CDS).*?\/locus_tag\=\"([^\"]*?)\"/s) {

            'reset' =~ /reset/;
            $old_genbank_ref{$old_file}->{$old_tags{$old_entry}} =~ /(gene|CDS).*?\/locus_tag\=\"([^\"]*?)\"/s;
            $old_locus_tag = $2;

            unless (defined $old_locus_tag) {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Error in parsing locus tag from entry ->\n$old_file\n$old_genbank_ref{$old_file}->{$old_tags{$old_entry}}",
                                                 -icon    => 'error',
                                                 -type    => 'ok');
               &hide_pbar_1;
               return;
            }
            $multi_locus_tag = $old_locus_tag;
            $old_locus_tag =~ s/\s+/ /gs;
         }

         #skip if wrong feature
         next if ($old_genbank_ref{$old_file}->{$old_tags{$old_entry}} !~ /^\s*$new_feature\s+/);

         #add in case of multiple hits
         $feature_counter++;
         #push (@multi_hits, $old_file.'___'.$old_tags{$old_entry}.'___'.$multi_anno);
         if (defined $multi_gene_anno && $multi_gene_anno =~ /\w+/) {
            push (@multi_hits, $new_file.'___'.$new_tags{$new_entry}.'___'.$multi_gene_anno);
         }
         if (defined $multi_product_anno && $multi_product_anno =~ /\w+/) {
            push (@multi_hits, $new_file.'___'.$new_tags{$new_entry}.'___'.$multi_product_anno);
         }
         if (defined $multi_locus_tag && $multi_locus_tag =~ /\w+/) {
            push (@multi_hits, $new_file.'___'.$new_tags{$new_entry}.'___'.$multi_locus_tag);
         }

         #transfer annotation
         'reset' =~ m/reset/;
         #$new_genbank_ref{$new_file}->{$new_tags{$new_entry}}  =~ s/\/gene\=\"[^\"]*?(_\s*\d+)\s*\"/\/gene\=\"$old_annotation$1\"/s;
         if ($args{transfer_gene_anno} == 1 && defined $old_gene_annotation && $old_gene_annotation =~ /\w+/)  {

            $new_genbank_ref{$new_file}->{$new_tags{$new_entry}}  =~ s/\/gene\=\"[^\"]*?(_\s*\d+)?\s*\"/\/gene\=\"$old_gene_annotation$1\"/s;
         }
         if ($args{transfer_product_anno} == 1 && defined $old_product_annotation && $old_product_annotation =~ /\w+/)  {
            $new_genbank_ref{$new_file}->{$new_tags{$new_entry}}  =~ s/\/product\=\"[^\"]*?(_\s*\d+)?\s*\"/\/product\=\"$old_product_annotation$1\"/s;
         }
         if ($args{transfer_locustag_anno} == 1 && defined $old_locus_tag && $old_locus_tag =~ /\w+/)  {
            if ($new_genbank_ref{$new_file}->{$new_tags{$new_entry}}  =~ /\/locus_tag/) {
               $new_genbank_ref{$new_file}->{$new_tags{$new_entry}}  =~ s/(locus_tag[^\n]+?\n)/$1                     \/note\=\"$old_locus_tag\"\n/s;#"
            } else {
               $new_genbank_ref{$new_file}->{$new_tags{$new_entry}}  =~ s/($old_gene_feature[^\n]+?\n)/$1                     \/locus_tag\=\"$old_locus_tag\"\n/s;#"
            }
         }
      }

      #multiple hits?
      if ($feature_counter > 1) {
         @new_multi_hits = (@new_multi_hits, @multi_hits);
      }
      undef @multi_hits;

   }

   # generate new Genbank file
   $max_count = keys %new_all_files;
   $progress = 1;
   while ( my ($key, $new_file) = each %new_all_files) {
      my (%genbank, @feature_list, $short_file);
      &update_pbar_1(title       => 'Transfer annotation',
                     label       => "Building new Genbank file for $new_file",
                     progress    => ($progress / $max_count) * 100
                    );
      $progress++;

      #re-constitute merged arrays and hashes
      while (my ($key,$value) = each %{ $new_genbank_ref{$new_file} }) {
         $genbank{$key} = $value;
      }
      foreach my $entry (@{ $new_feature_list_ref{$new_file} }) {
         push (@feature_list, $entry);
      }

      #determine short file name
      'reset' =~ m/reset/;
      $new_file =~ m/.+[\/\\](.+)$/;
      $short_file = $1;
      unless (defined $short_file) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Error in parsing file name from ->\n$new_file",
                                           -icon    => 'error',
                                           -type    => 'ok');
         &hide_pbar_1;
         return;
      }

      &compile_gb(main_window      => $args{main_window},
                  progress_bar     => $args{progress_bar},
                  auto_ini_ref     => $args{auto_ini_ref},
                  ini_ref          => $args{ini_ref},
                  input_file       => $short_file,
                  DNAseq           => $new_DNAseq{$new_file},
                  genbank_ref      => \%genbank,
                  feature_list_ref => \@feature_list,
                  genbank_header   => $new_header{$new_file},
                  genbank_source   => $new_source{$new_file},
                  destination      => ${$args{ini_ref}}{annotation_transfer_dir}
                 );
      undef %genbank;
      undef @feature_list;

   }

   #write unassigned and multiple hits to files
   &write_ambiguous(main_window      => $args{main_window},
                    progress_bar     => $args{progress_bar},
                    auto_ini_ref     => $args{auto_ini_ref},
                    ini_ref          => $args{ini_ref},
                    old_nohits       => \@old_nohits,
                    new_nohits       => \@new_nohits,
                    ORFnumbers       => \@ORFnumbers,
                    old_multi_hits   => \@old_multi_hits,
                    new_multi_hits   => \@new_multi_hits,
                    new_genbank_ref  => \%new_genbank_ref,
                    old_genbank_ref  => \%old_genbank_ref,
                   );
   &hide_pbar_1;
   undef @old_nohits;
   undef @new_nohits;
   undef @new_multi_hits;
   undef @old_multi_hits;

   #unlink temporary databases
   my @suffix = qw(phr pin pnd pni psd psi psq pal);
   foreach my $suffix (@suffix) {
      unlink $local_blastdb_dir.'/new.aa.'.$suffix;
      unlink $local_blastdb_dir.'/old.aa.'.$suffix;
   }

   #done
   ${$args{main_window}}->messageBox(-title   => 'Info',
                                     -message => "Annotation transfer completed",
                                     -icon    => 'info',
                                     -type    => 'ok');
   return;

}

sub write_ambiguous {
   my %args = @_;

   #write old ORfnumbers to new ORF number
   open WRITE, "+>${$args{ini_ref}}{annotation_transfer_dir}\/Old_to_new_ORFnumbers.txt";
   foreach my $entry (sort @{ $args{ORFnumbers} }) {
      print WRITE $entry;
   }
   close WRITE;

   if ($#{$args{old_nohits}} >= 0) {
      open WRITE, "+>${$args{ini_ref}}{annotation_transfer_dir}\/Unassigned_ORFs_from_existing_annotation.txt";
      foreach my $entry (sort @{ $args{old_nohits} }) {
         my ($file, $ID, $annotation);
         #parse entry
         'reset' =~ m/reset/;
         $entry  =~ m/(.+?)___(.+)/;
         $file   =  $1;
         $ID     =  $2;
         unless (defined $ID) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Error in parsing entry from OLD_nohits array ->\n$entry",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return;
         }

         #get annotation
         'reset' =~ m/reset/;
         ${$args{old_genbank_ref}}{$file}->{$ID}  =~ m/\/gene\=\"([^\"]*?)\"/;
         $annotation = $1;

         #no hit? try product
         unless (defined $annotation) {
            'reset' =~ m/reset/;
            ${$args{old_genbank_ref}}{$file}->{$ID}  =~ m/\/product\=\"([^\"]*?)\"/;
            $annotation = $1;
         }
         #no hit? try locus_tag
         unless (defined $annotation) {
            'reset' =~ m/reset/;
            ${$args{old_genbank_ref}}{$file}->{$ID}  =~ m/\/locus_tag\=\"([^\"]*?)\"/;
            $annotation = $1;
         }
         unless (defined $annotation) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Error in parsing annotation entry from OLD_nohits array ->\n${$args{old_genbank_ref}}{$file}->{$ID}",
                                              -icon    => 'error',
                                              -type    => 'ok');
         }
         $annotation =~ s/\s+/ /gs;

         #make short filename
         $file =~ s/^.*\/(.+)$/$1/;

         print WRITE $file.'-->'.$ID.'-->'.$annotation."\n";
      }
      close WRITE;
   }

   if ($#{$args{new_nohits}} >= 0) {
      open WRITE, "+>${$args{ini_ref}}{annotation_transfer_dir}\/Unassigned_ORFs_in_new_Genbank_files.txt";
      foreach my $entry (sort @{ $args{new_nohits} }) {
         my ($file, $ID, $annotation);
         #parse entry
         'reset' =~ m/reset/;
         $entry  =~ m/(.+?)___(.+)/;
         $file   =  $1;
         $ID     =  $2;

         unless (defined $ID) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Error in parsing entry from NEW_nohits array ->\n$entry",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return;
         }

         #get annotation
         'reset' =~ m/reset/;
         if (${$args{new_genbank_ref}}{$file}->{$ID}  =~ m/\/(gene|product)\=\"/) {
            ${$args{new_genbank_ref}}{$file}->{$ID}  =~ m/\/(gene|product)\=\"([^\"]*?)\"/;
            $annotation = $2;

            unless (defined $annotation) {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Error in parsing annotation entry from NEW_nohits array ->\n${$args{new_genbank_ref}}{$file}->{$ID}",
                                                 -icon    => 'error',
                                                 -type    => 'ok');
            }
         } else {
            $annotation = 'N/A';
         }
         $annotation =~ s/\s+/ /gs;

         #make short filename
         $file =~ s/^.*\/(.+)$/$1/;

         print WRITE $file.'-->'.$ID.'-->'.$annotation."\n";
      }
      close WRITE;
   }

   if ($#{$args{new_multi_hits}} >= 0) {
      open WRITE, "+>${$args{ini_ref}}{annotation_transfer_dir}\/Multiple_hits_in_new_Genbank_files.txt";
      foreach my $entry (sort @{ $args{new_multi_hits} }) {
         my ($file, $ID, $annotation);
         'reset'     =~ m/reset/;
         $entry      =~ m/(.*?)___(.*?)___(.+)/;
         $file       =  $1;
         $ID         =  $2;
         $annotation =  $3;

         unless (defined $annotation) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Error in parsing entry from new_multi_hits array ->\n$$entry",
                                              -icon    => 'error',
                                              -type    => 'ok');
         }
         $annotation =~ s/\s+/ /gs;

         #make short filename
         $file =~ s/^.*\/(.+)$/$1/;

         print WRITE $file.'-->'.$ID.'-->'.$annotation."\n";
      }
      close WRITE;
   }

   if ($#{$args{old_multi_hits}} >= 0) {
      open WRITE, "+>${$args{ini_ref}}{annotation_transfer_dir}\/Multiple_hits_in_existing_annotation.txt";
      foreach my $entry (sort @{ $args{old_multi_hits} }) {
         my ($file, $ID, $annotation);
         'reset'     =~ m/reset/;
         $entry      =~ m/(.*?)___(.*?)___(.+)/;
         $file       =  $1;
         $ID         =  $2;
         $annotation =  $3;

         unless (defined $annotation) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Error in parsing entry from old_multi_hits array ->\n$$entry",
                                              -icon    => 'error',
                                              -type    => 'ok');
         }
         $annotation =~ s/\s+/ /gs;

         #make short filename
         $file =~ s/^.*\/(.+)$/$1/;

         print WRITE $file.'-->'.$ID.'-->'.$annotation."\n";
      }
      close WRITE;
   }

}

sub old_database {
   my %args = @_;
   my (@old_coding, $old_id, %old_seen);

   #write temporary msfasta file to create Blast database
   open WRITE, "+>${$args{ini_ref}}{annotation_transfer_dir}\/old_seq.fasta";

   #build Blast database from old Genbank files
   $old_id = 1;
   while (my ($keys, $file) = each %old_all_files) {
      my ($old_header, $old_source, $old_core, $old_DNAseq,
          $old_start_ref, $old_feature_list_ref, $old_genbank_ref,
          $old_feature_counter) = &gb_parser (main_window   => $args{main_window},
                                              progress_bar  => $args{progress_bar},
                                              auto_ini_ref  => $args{auto_ini_ref},
                                              ini_ref       => $args{ini_ref},
                                              gb_file       => $file,
                                             );

      #dereference and COPY to new targets
      $args{old_header}{$file}           = $old_header;
      $args{old_source}{$file}           = $old_source;
      $args{old_core}{$file}             = $old_core;
      $args{old_DNAseq}{$file}           = $old_DNAseq;
      $args{old_start_ref}{$file}        = [ @{$old_start_ref}] ;
      $args{old_feature_list_ref}{$file} = [ @{$old_feature_list_ref} ];
      while (my ($key, $value) = each %{$old_genbank_ref}) {
          $args{old_genbank_ref}->{$file}->{$key} = $value;
       }
      $args{old_feature_counter}{$file}  = $old_feature_counter;

      #get gene/CDS features
      @old_coding = grep {/_(gene|CDS)_/} @{$old_feature_list_ref};
      next if ($#old_coding < 0); #skip if no coding regions found in Genebank file

      #iterate through coding features
      foreach my $entry (@old_coding) {
         my ($ID, $left_bd, $right_bd, $orientation, $query_seq_nt, $query_seq_aa);
         'reset' =~ m/reset/;
         $entry  =~ m/^(\d+)_(\d+)_/;
         $left_bd         = $1;
         $right_bd        = $2;

         unless (defined $right_bd) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Error in parsing boundaries for 'old' entry $entry",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return;
         }

         #get orientation
         'reset' =~ m/reset/;
         if (${$old_genbank_ref}{$entry} =~ m/(gene|CDS)[^\n\r]*?complement/) {
            $orientation = 'antisense';
         } else {
            $orientation = 'sense';
         }

         #extract nt sequence from ORF
         $query_seq_nt = substr($old_DNAseq, ($left_bd - 1), ($right_bd - $left_bd + 1));

         #translate nt sequence to aa if selected
         $query_seq_aa = nt2aa_unchecked(main_window  => $args{main_window},
                                         progress_bar => $args{progress_bar},
                                         auto_ini_ref => $args{auto_ini_ref},
                                         ini_ref      => $args{ini_ref},
                                         orientation  => $orientation,
                                         sequence     => $query_seq_nt,
                                         filename     => $file,
                                         left_bd      => $left_bd,
                                         right_bd     => $right_bd
                                        );
         next if ($query_seq_aa eq '0');

         #entry already exists?
         if (exists $old_seen{$query_seq_aa}) {
            push (@{$args{old_list}}, $query_seq_aa.'_'.$old_id.'_'.$file);
            ${$args{old_tags}}{$query_seq_aa.'_'.$old_id.'_'.$file} = $entry; #add feature to hash to capure it for later
            $old_id++;
            next;
         }
         $old_seen{$query_seq_aa}++;

         #add line to fasta file
         #print WRITE '>gi|'.$old_id.'|gb|'.$old_id.'| '.substr($query_seq_aa,(length($query_seq_aa) - 30)).'_'.$old_id."\n".$query_seq_aa."\n";
         print WRITE '>gi|'.$old_id.'|gb|'.$old_id.'| '.substr($query_seq_aa,0,50).'_'.$old_id."\n".$query_seq_aa."\n";

         #add to tag array and hash to speed up annotation transfer
         push (@{$args{old_list}}, $query_seq_aa.'_'.$old_id.'_'.$file);
         ${$args{old_tags}}{$query_seq_aa.'_'.$old_id.'_'.$file} = $entry;

         $old_id++;
      }
   }
   undef %old_seen;
   close WRITE;

   #no entries for database? no gene model->abort
   if (-s ${$args{ini_ref}}{annotation_transfer_dir}.'/old_seq.fasta' < 1) {
       ${$args{main_window}}->messageBox(-title   => 'Error',
                                         -message => "No entries for Blast databases.\nProbably no gene models present in input files.",
                                         -icon    => 'error',
                                         -type    => 'OK');
       return (0);
   }

   #create Blast database
   `${$args{ini_ref}}{blast_executables}/formatdb -i ${$args{ini_ref}}{annotation_transfer_dir}/old_seq.fasta -p T -o T -n $local_blastdb_dir/old.aa`;

   #no database created? abort
   unless (-e $local_blastdb_dir.'/old.aa.pin') {
       ${$args{main_window}}->messageBox(-title   => 'Error',
                                         -message => "No Blast database created for source files.\n",
                                         -icon    => 'error',
                                         -type    => 'OK');
       return (0);
   }

   unlink "${$args{ini_ref}}{annotation_transfer_dir}\/old_seq.fasta";

   return (1);
}

sub new_database {
   my %args = @_;
   my (@new_coding, $new_id, %new_seen);

   #write temporary msfasta file to create Blast database
   open WRITE, "+>${$args{ini_ref}}{annotation_transfer_dir}\/new_seq.fasta";

   #build Blast database from new Genbank files
   $new_id = 1;
   while (my ($key, $file) =  each %new_all_files) {
      my ($new_header, $new_source, $new_core, $new_DNAseq,
          $new_start_ref, $new_feature_list_ref, $new_genbank_ref,
          $new_feature_counter) = &gb_parser (main_window   => $args{main_window},
                                              progress_bar  => $args{progress_bar},
                                              auto_ini_ref  => $args{auto_ini_ref},
                                              ini_ref       => $args{ini_ref},
                                              gb_file       => $file,
                                             );

       #dereference and COPY to new targets
       $args{new_header}{$file}           = $new_header;
       $args{new_source}{$file}           = $new_source;
       $args{new_core}{$file}             = $new_core;
       $args{new_DNAseq}{$file}           = $new_DNAseq;
       $args{new_start_ref}{$file}        = [ @{$new_start_ref} ];
       $args{new_feature_list_ref}{$file} = [ @{$new_feature_list_ref} ];
       while (my ($key, $value) = each %{$new_genbank_ref}) {
          $args{new_genbank_ref}->{$file}->{$key} = $value;
       }
       $args{new_feature_counter}{$file}  = $new_feature_counter;


      #get gene/CDS features
      @new_coding = grep {/_(gene|CDS)_/} @{$new_feature_list_ref};
      next if ($#new_coding < 0); #skip if no coding regions found in Genbank file

      #iterate through coding features
      foreach my $entry (@new_coding) {
         my ($ID, $left_bd, $right_bd, $orientation, $query_seq_nt, $query_seq_aa);
         'reset' =~ m/reset/;
         $entry  =~ m/^(\d+)_(\d+)_/;
         $left_bd         = $1;
         $right_bd        = $2;

         unless (defined $right_bd) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Error in parsing boundaries for 'new' entry $entry",
                                              -icon    => 'error',
                                              -type    => 'ok');
            return;
         }

         #get orientation
         'reset' =~ m/reset/;
         if (${$new_genbank_ref}{$entry} =~ m/(gene|CDS)[^\n\r]*?complement/) {
            $orientation = 'antisense';
         } else {
            $orientation = 'sense';
         }

         #extract nt sequence from ORF
         $query_seq_nt = substr($new_DNAseq, ($left_bd - 1), ($right_bd - $left_bd + 1));

         #translate nt sequence to aa if selected
         $query_seq_aa = nt2aa_unchecked(main_window  => $args{main_window},
                                         progress_bar => $args{progress_bar},
                                         auto_ini_ref => $args{auto_ini_ref},
                                         ini_ref      => $args{ini_ref},
                                         orientation  => $orientation,
                                         sequence     => $query_seq_nt,
                                         filename     => $file,
                                         left_bd      => $left_bd,
                                         right_bd     => $right_bd
                                        );
         next if ($query_seq_aa eq '0');

         #entry already exists?
         if (exists $new_seen{$query_seq_aa}) {
            push (@{$args{new_list}}, $query_seq_aa.'_'.$new_id.'_'.$file);
            $args{new_tags}{$query_seq_aa.'_'.$new_id.'_'.$file} = $entry; #add feature to hash to capure it for later
            $new_id++;
            next;
         }
         $new_seen{$query_seq_aa}++;

         #add line to fasta file
         #print WRITE '>gi|'.$new_id.'|gb|'.$new_id.'| '.substr($query_seq_aa,(length($query_seq_aa) - 30)).'_'.$new_id."\n".$query_seq_aa."\n";
         print WRITE '>gi|'.$new_id.'|gb|'.$new_id.'| '.substr($query_seq_aa,0, 50).'_'.$new_id."\n".$query_seq_aa."\n";

         #add to tag array and hash to speed up annotation transfer
         push (@{$args{new_list}}, $query_seq_aa.'_'.$new_id.'_'.$file);
         $args{new_tags}{$query_seq_aa.'_'.$new_id.'_'.$file} = $entry;

         $new_id++;
      }
   }
   undef %new_seen;
   close WRITE;

   #no entries for database? no gene model->abort
   if (-s ${$args{ini_ref}}{annotation_transfer_dir}.'/new_seq.fasta' < 1) {
       ${$args{main_window}}->messageBox(-title   => 'Error',
                                         -message => "No entries for Blast databases.\nProbably no gene models present in input files.",
                                         -icon    => 'error',
                                         -type    => 'OK');
       return (0);
   }

   #create Blast database
   `${$args{ini_ref}}{blast_executables}/formatdb -i ${$args{ini_ref}}{annotation_transfer_dir}/new_seq.fasta -p T -o T -n $local_blastdb_dir/new.aa`;

   #no database created? abort
   unless (-e $local_blastdb_dir.'/new.aa.pin') {
       ${$args{main_window}}->messageBox(-title   => 'Error',
                                         -message => "No Blast database created for target files.\n",
                                         -icon    => 'error',
                                         -type    => 'OK');
       return (0);
   }

   unlink "${$args{ini_ref}}{annotation_transfer_dir}\/new_seq.fasta";

   return (1);
}

sub populate_new {
   my %args = @_;

   (my $new_selected_files) = ${$args{main_window}}->getOpenFile(-initialdir => $current_dir,
                                                                 -title      => 'Select files for custom database',
                                                                 -multiple   => 1);

   if (defined(@{$new_selected_files})) {
      #generate short_list
      foreach my $file (@{$new_selected_files}) {
         'reset' =~ m/reset/;
         $file =~ m/(.+)\/([^\/]+)$/; #display only file_name
         $current_dir = $1;
         my $short = $2;
         push (@new_short_files, $short);
         #configure new width for textbox
         if (length($short) > $new_width) {
            $new_width = length($short);
         }
      }

      #remove duplicates
      &delete_duplicates(array => \@new_short_files);

      #populate textbox with short files
      $args{textbox}->delete(0, 'end');
      foreach my $entry (@new_short_files) {
         $args{textbox}->insert('end', $entry);
      }

      #merge selected files with overall selection array
      foreach (@{$new_selected_files}) {
         'reset' =~ m/reset/;
         m/\/([^\/]+)$/; #display only file_name
         my $short = $1;
         $new_all_files{$short} = $_;
      }

      #update display
      #$args{textbox}->configure(-width => $new_width);
      $tl->update;
   }
}

sub populate_old {
   my %args = @_;

   (my $old_selected_files) = ${$args{main_window}}->getOpenFile(-initialdir => $current_dir,
                                                                 -title      => 'Select files for custom database',
                                                                 -multiple   => 1);

   if (defined(@{$old_selected_files})) {
      #generate short_list
      foreach my $file (@{$old_selected_files}) {
         'reset' =~ m/reset/;
         $file =~ m/(.+)\/([^\/]+)$/; #display only file_name
         $current_dir = $1;
         my $short = $2;
         push (@old_short_files, $short);
         #configure new width for textbox
         if (length($short) > $old_width) {
            $old_width = length($short);
         }
      }

      #remove duplicates
      &delete_duplicates(array => \@old_short_files);

      #populate textbox with short files
      $args{textbox}->delete(0, 'end');
      foreach my $entry (@old_short_files) {
         $args{textbox}->insert('end', $entry);
      }

      #merge selected files with overall selection array
      foreach (@{$old_selected_files}) {
         'reset' =~ m/reset/;
         m/\/([^\/]+)$/; #display only file_name
         my $short = $1;
         $old_all_files{$short} = $_;
      }

      #update display
      #$args{textbox}->configure(-width => $old_width);
      $tl->update;
   }
}

sub remove_entry_new {
   my %args = @_;
   return if ($#new_short_files < 0);
   my ($current_index) = ${$args{textbox}}->curselection();
   return unless ($current_index =~ /\S+/);
   #remove from shortlist array
   splice(@new_short_files, $current_index, 1);
   #remove from all files hash
   my $sel_key = ${$args{textbox}}->get($current_index);
   foreach my $key (keys %new_all_files) {
      if ($key =~ m/$sel_key$/) {
         delete $new_all_files{$key};
         ${$args{progress_bar}}->Label(-text => "Deleted entry $key")-> grid (-row => 0, -column => 0, -sticky => 'esnw');
         ${$args{main_window}}->update;
         last;
      }
   }
   #remove from display
   ${$args{textbox}}->delete($current_index);
   ${$args{main_window}}->update;

}

sub remove_entry_old {
   my %args = @_;
   return if ($#old_short_files < 0);
   my ($current_index) = ${$args{textbox}}->curselection();
   return unless ($current_index =~ /\S+/);
   #remove from shortlist array
   splice(@old_short_files, $current_index, 1);
   #remove from all files hash
   my $sel_key = ${$args{textbox}}->get($current_index);

   foreach my $key (keys %old_all_files) {
      if ($key =~ m/$sel_key$/) {
         delete $old_all_files{$key};
         ${$args{progress_bar}}->Label(-text => "Deleted entry $key")-> grid (-row => 0, -column => 0, -sticky => 'esnw');
         ${$args{main_window}}->update;
         last;
      }
   }
   #remove from display
   ${$args{textbox}}->delete($current_index);
   ${$args{main_window}}->update;

}

sub delete_duplicates {
   my %args = @_;
   my %seen = ();
   @{$args{array}} = grep { ! $seen{$_} ++ } @{$args{array}};
}




1;




