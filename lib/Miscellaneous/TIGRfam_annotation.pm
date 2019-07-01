#!/opt/ActivePerl-5.8/bin/perl

#rotate_gb: rotates a selected Genbank fileno clockwise or counterclockwise
#input arguments: main_window, ini_ref, auto_ini_ref,


package Miscellaneous::TIGRfam_annotation;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&TIGRfam_annotation);
use vars qw();

use Tk;
use initialise::read_me              qw(:DEFAULT);
use ProgrammeModules::sequence       qw(:DEFAULT);
use ProgrammeModules::genbank_parser qw(:DEFAULT);
use CompilerModules::tigrfam_parser  qw(:DEFAULT);
use Basics::progress_bar             qw(:DEFAULT);
use Cwd;

#local variables
my (%args, $tl, %all_files, @short_files,  @selected_files, $width, $current_dir,
    $TIGR_EConly, $TIGR_ECCDS, $select_files, $file_lb, $clear_selection);

#initialise

sub TIGRfam_annotation {
   my %args = @_;
   $width = 40;
   my ($file_lb, $EConly, $ECCDS, $TIGR_ECdiscard, $compile,
       $frame_left, $frame_right, $progress_bar);


   if (defined $tl && Tk::Exists($tl)) {
      $tl->state('normal');
      $progress_bar = $tl->Frame(-borderwidth => 2, -relief => 'groove');
      $frame_left   = $tl->Frame()                                      ;
      $frame_right  = $tl->Frame(-relief => 'groove')                   ;
   } else {
      $tl = ${$args{main_window}}->DialogBox(-title   => 'Transfer TIGRfam annotation',
                                             -buttons => [ 'OK' ]
                                            );
      $tl->minsize(20, 10);
      $tl->maxsize(120, 20);
      $progress_bar = $tl->Frame(-borderwidth => 2, -relief => 'groove') -> pack(-side => 'bottom',-fill => 'x');
      $frame_left   = $tl->Frame()                                       -> pack(-side => 'left',  -fill => 'y');
      $frame_right  = $tl->Frame(-relief => 'groove')                    -> pack(-side => 'left',  -fill => 'both', -expand => 1);
   }

   $progress_bar->configure(-label => " "); #create dummy for space

   #set some defaults
   $TIGR_EConly = 0;
   $TIGR_ECCDS  = 1;

   #entries into frame
   {
      $frame_left                	->Label(-text => 'Transfer TIGRfam_annotation')-> grid (-row => 0, -column => 0, -columnspan => 2, -sticky => 'esnw',-pady => 25); #dummy

      $EConly			= $frame_left	->Checkbutton(-text   	=> 'Transfer EC numbers only',
                                                 -variable	=> \$TIGR_EConly,
                                                 -command 	=> sub {
                                                               if ($TIGR_EConly == 1) {
                                                                  $TIGR_ECCDS = 0;
                                                               } elsif ($TIGR_EConly == 0) {
                                                                  $TIGR_ECCDS = 1;
                                                               }
                                                             }
                                       )
      -> grid(-row => 1, -column => 0, -sticky => 'w');

      $ECCDS      	= $frame_left->Checkbutton(-text     => 'Transfer EC numbers and gene annotation',
                                               -variable => \$TIGR_ECCDS,
                                               -command  => sub {
                                                                  if ($TIGR_ECCDS == 1) {
                                                                     $TIGR_EConly = 0;
                                                                  } elsif ($TIGR_ECCDS == 0) {
                                                                     $TIGR_EConly = 1;
                                                                  }
                                                                }
                                              )
      -> grid (-row => 2, -column => 0, -sticky => 'w');

      $frame_left->Checkbutton(-text     => 'Discard existing EC number entries',
                               -variable => \$TIGR_ECdiscard,
                              )
      -> grid (-row => 3, -column => 0, -sticky => 'w');

      $compile   = $frame_left ->Button(-text => 'Transfer annotation',
                                        -command => sub {
                                                          &compile(main_window       => \$tl,
                                                                   progress_bar      => \$progress_bar,
                                                                   ini_ref           => $args{ini_ref},
                                                                   auto_ini_ref      => $args{auto_ini_ref},
                                                                   TIGR_EConly       => $TIGR_EConly,
                                                                   TIGR_ECdiscard    => $TIGR_ECdiscard
                                                                  );
                                                        }
                                          )->grid(-row => 9, -column => 0, -sticky => 'sew');
      $frame_left->gridRowconfigure   (9, -weight => 1);
   }


   #Textbox and button in right frame
   {
      $select_files = $frame_right       ->Button(-text => 'Add files',
                                                  -command => sub {&populate(main_window  => \$tl,
                                                                             auto_ini_ref => $args{auto_ini_ref},
                                                                             textbox      => $file_lb)}
                                                 )->grid(-row => 0, -column => 0);
      $frame_right                       ->Label(-text => 'Double click entry to remove.') -> grid(-row => 0, -column => 1, -sticky => 'nsew');
      $file_lb         = $frame_right    ->Scrolled("Listbox",
                                                    -setgrid    => 1,
                                                    -selectmode => 'single',
                                                    -width      => $width
                                                    )->grid(-row => 1, -column => 0, -columnspan => 2, -sticky => 'nsew');
      $file_lb                           ->configure(-width   => $width,
                                                     -setgrid => 1 );
      $clear_selection = $frame_right ->Button(-text => 'Clear current selection',
                                               -command => sub {
                                                                 %all_files = ();
                                                                 @short_files = ();
                                                                 $file_lb->delete(0, 'end');
                                                                 $tl->update;
                                                                }
                                              )->grid(-row => 2, -column => 0, -columnspan => 2, -sticky => 'nsew');

   }
   #bind right mouse button to selection and removal of entry
   $file_lb->bind('<Double-ButtonRelease-1>' => sub {&remove_entry(textbox      => \$file_lb,
                                                                   main_window  => \$tl,
                                                                   progress_bar => \$progress_bar,
                                                                   ini_ref      => $args{ini_ref},
                                                                   auto_ini_ref => $args{auto_ini_ref}
                                                                  )
                                                    });
   my $wait = $tl->Show();
   if ($wait eq 'OK') {
      $file_lb->delete(0, 'end');
      $tl->update;
      undef %all_files;
      undef @short_files;
      undef @selected_files;
      $tl->state('withdrawn');
   }
}

sub compile {
   my %args = @_;
   my ($max_count, $counter);

   #create status box
   &progress_bar_1(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Processing Genbank files",
                   label        => 'Processing features'
                  );
   &show_pbar_1;

   $max_count = $#short_files + 1;
   $counter = 1;
   while (my ($key, $value) = each (%all_files)) {
      my (@tigrs, @CDS, $CDScounter, $max_cds);

      #parse gb file
      &update_pbar_1(title        => "Processing Genbank files",
                     label        => "Reading entries for: $key",
                     progress     => ($counter / $max_count) * 100,
                    );
      $counter++;

      #get Genbank features
      my ($header, $source, $core, $DNAseq, $start_ref,
          $feature_list_ref, $genbank_ref, $gb_counter) = &gb_parser(main_window   => $args{main_window},
                                                                     progress_bar  => $args{progress_bar},
                                                                     ini_ref       => $args{ini_ref},
                                                                     auto_ini_ref  => $args{auto_ini_ref},
                                                                     gb_file       => $value,
                                                                    );

      ${$args{progress_bar}}->configure(-label=>"Grabbing required features");
      ${$args{main_window}}->update;

      #first, grab TIGRfams from feature list ref
      @tigrs = grep {/TIGR_match/} @{$feature_list_ref};

      #no TIGRfams? next;
      next if ($#tigrs < 0);

      #second, grab CDS features
      @CDS = grep {/CDS/} @{$feature_list_ref};

      #create status box
      &progress_bar_3(main_window  => $args{main_window},
                      progress_bar => $args{progress_bar},
                      auto_ini_ref => $args{auto_ini_ref},
                      ini_ref      => $args{ini_ref},
                      title        => "Processing features",
                      label        => 'Processing features'
                     );
      &show_pbar_3;

      #iterate through individual CDS features ->take care of joined features?
      ${$args{progress_bar}}->configure(-label=>"Iterating through features");
      ${$args{main_window}}->update;
      $max_cds = $#CDS + 1;
      $CDScounter = 1;
      foreach my $CDS (@CDS) {
         my ($boundary, @pairs, @gene, $best_evalue, $best_tigr, $CDS_left, $CDS_right);

         #parse CDSs
         if ($CDScounter % 100 == 0) {
            &update_pbar_3(title        => "Processing features",
                           label        => "Processing features",
                           progress     => ($CDScounter / $max_cds) * 100,
                          );
         }
         $CDScounter++;

         #get CDS feature boundaries to fetch corresponding gene feature
         'reset'    =~ m/reset/;
         $CDS       =~ m/(\d+)\_(\d+)\_/;
         $CDS_left  = $1;
         $CDS_right = $2;

         #fetch gene feature
         @gene = grep (/$CDS_left\_$CDS_right\_gene\_/, @{$feature_list_ref});

         #get boundary line from entry
         'reset' =~ m/reset/;
         ${$genbank_ref}{$CDS} =~ m/CDS(.*?)\//s;
         $boundary = $1;

         #split into individual boundaries if joined
         $boundary =~ s/\s+//gs;
         $boundary =~ s/[\(\)]//g;
         $boundary =~ s/[\>\<]//g;
         $boundary =~ s/(join|complement)//g;
         @pairs = split/\,/,$boundary;

         #iterate through individual boundaries
         $best_evalue = 100;
         foreach my $set (@pairs) {
            my ($left_bd, $right_bd, @local);
            'reset' =~ m/reset/;
            $set =~ m/(\d+)\.\.(\d+)/;
            $left_bd = $1;
            $right_bd = $2;
            #get corresponding TIGRfams
            @local = grep { m/^(\d+)\_/;
                            $1 >= $left_bd;} @tigrs;
            @local = grep { m/^\d+\_(\d+)\_/;
                            $1 <= $right_bd;} @local;

            #find best TIGRfam hit
            foreach my $tigr (@local) {
               my ($evalue, $product);
               'reset'                   =~ m/reset/;
               ${$genbank_ref}{$tigr} =~ m/\/product\=\"([^\"]*?)\"/s;
               $product = $1;
               $product =~ s/\s+//gs;
               'reset'  =~ m/reset/;
               $product =~ m/e-value:\s*(.+)/;
               $evalue  = $1;
               $evalue  =~ s/\s+//g;
               unless (defined $evalue) {$evalue = 1}; #if not defiend, set to an arbitrary number to catch anyway
               if ($evalue < $best_evalue) {
                  $best_evalue = $evalue;
                  $best_tigr   = ${$genbank_ref}{$tigr};
               }
            }
         }
         #submit to TIGRparsing
         my ($label, $EC, $descriptor) = &tigr_annotation(main_window  => $args{main_window},
                                                          progress_bar => $args{progress_bar},
                                                          auto_ini_ref => $args{auto_ini_ref},
                                                          ini_ref      => $args{ini_ref},
                                                          annotation   => $best_tigr,
                                                         );

         #no results? next
         next if ($label !~ /\w+/ && $EC !~ /\w+/ && $descriptor !~ /\w+/);

         #update annotation
         if ($args{TIGR_EConly} == 0) { #apply gene/cds TIGR annotation
            #first gene
            if (defined $gene[0] && ${$genbank_ref}{$gene[0]} =~ /\/gene/ && $label =~ /\w+/) {
               ${$genbank_ref}{$gene[0]} =~ s/(\/gene\=\").*?(_?\s*\d*\s*\")/$1$label$2/s;
            } elsif (defined $gene[0] && ${$genbank_ref}{$gene[0]} =~ /\/note/ && $label =~ /\w+/) {
               ${$genbank_ref}{$gene[0]} =~ s/(\/note\=\").*?(_?\s*\d*\s*\")/$1$label$2/s;
            }
            #then CDS
            if (${$genbank_ref}{$CDS} =~ /\/gene/  && $descriptor =~ /\w+/) {
               if ($label =~ /\w+/) {$label = ', '.$label};
               ${$genbank_ref}{$CDS} =~ s/(\/gene\=\").*?(_\s*\d+\s*\")/$1$descriptor$label$2/s;
            } elsif (${$genbank_ref}{$CDS} =~ /\/product/  && $descriptor =~ /\w+/) {
               if ($label =~ /\w+/) {$label = ', '.$label};
               ${$genbank_ref}{$CDS} =~ s/(\/product\=\").*?(_\s*\d+\s*\")/$1$descriptor$label$2/s;
            } elsif (${$genbank_ref}{$CDS} =~ /\/note/  && $descriptor =~ /\w+/) {
               if ($label =~ /\w+/) {$label = ', '.$label};
               ${$genbank_ref}{$CDS} =~ s/(\/note\=\").*?(_\s*\d+\s*\")/$1$descriptor$label$2/s;
            }
         }

         #generate ECnumber entries
         my @EC = split /\s+/, $EC; #if there are multiple EC numbers, put each in a separate entry
         $EC    = '';
         foreach my $entry (@EC) {
            next unless ($entry =~ /\S+/);
            $EC .= '                     /EC_number='.$entry."\n";
         }
         $EC =~ s/^                     \/EC_number\=//;
         $EC =~ s/\n$//s;

         #delete existing EC number entries if selected
         if (${$genbank_ref}{$CDS} =~ /\/EC_number/ && $EC =~ /\w+/ && $args{TIGR_ECdiscard} == 1) {
            ${$genbank_ref}{$CDS} =~ s/\/EC_number.*?(\s+?\/)/$1/gs;
         }

         #add EC number
         if (${$genbank_ref}{$CDS} =~ /\/(gene|product|note)/ && ${$genbank_ref}{$CDS} !~ /\/EC_number/ && $EC =~ /\w+/) {
            if (${$genbank_ref}{$CDS} =~ /\/gene/) {
               ${$genbank_ref}{$CDS} =~ s/(\/gene\=\"[^\"]*?\")/$1\n                     \/EC_number\=$EC/s;
            } elsif (${$genbank_ref}{$CDS} =~ /\/product/) {
               ${$genbank_ref}{$CDS} =~ s/(\/product\=\"[^\"]*?\")/$1\n                     \/EC_number\=$EC/s;
            } elsif (${$genbank_ref}{$CDS} =~ /\/note/) {
               ${$genbank_ref}{$CDS} =~ s/(\/note\=\"[^\"]*?\")/$1\n                     \/EC_number\=$EC/s;
            } else {
               ${$genbank_ref}{$CDS} =~ s/(                     \/)/                     \/EC_number\=$EC\n$1/s;
            }
         } elsif (${$genbank_ref}{$CDS} =~ /\/EC_number/ && $EC =~ /\w+/ && $args{TIGR_ECdiscard} == 1) {
            ${$genbank_ref}{$CDS} =~ s/(\/EC_number\=\"?.*?\"?)\n/\/EC_number\=$EC/s;
         } elsif ($EC =~ /\w+/) {
            ${$genbank_ref}{$CDS} =~ s/(                     \/)/                     \/EC_number\=$EC\n$1/s;
         }
      }
      &hide_pbar_3;

      &update_pbar_1(title        => "Processing Genbank files",
                     label        => "Writing new genbank file",
                     progress     => ($counter / $max_count) * 100,
                    );

      #write new genbank file
      ${$args{progress_bar}}->configure(-label=>"Writing new Genbank file");
      ${$args{main_window}}->update;
      &compile_gb(main_window      => $args{main_window},
                  progress_bar     => $args{progress_bar},
                  auto_ini_ref     => $args{auto_ini_ref},
                  ini_ref          => $args{ini_ref},
                  input_file       => $key,
                  genbank_ref      => $genbank_ref,
                  feature_list_ref => $feature_list_ref,
                  genbank_header   => $header,
                  genbank_source   => $source,
                  DNAseq           => $DNAseq,
                  update_allowed   => 0
                 );
   }
   &hide_pbar_1;
   ${$args{progress_bar}}->configure(-label=>"Finished");
   ${$args{main_window}}->update;
   ${$args{main_window}}->messageBox(-title   => 'Info',
                                     -message => "Successfully finished annotation transfer",
                                     -icon    => 'info',
                                     -type    => 'ok');
   return;
}



sub populate {
   my %args = @_;

   (my $selected_files) = ${$args{main_window}}->getOpenFile(-initialdir => $current_dir,
                                                             -title      => 'Select files for custom database',
                                                             -multiple   => 1);

   if (defined(@{$selected_files})) {
      #generate short_list
      foreach my $file (@{$selected_files}) {
         'reset' =~ m/reset/;
         $file =~ m/(.+)\/([^\/]+)$/; #display only file_name
         $current_dir = $1;
         my $short = $2;
         push (@short_files, $short);
         #configure new width for textbox
         if (length($short) > $width) {
            $width = length($short);
         }
      }

      #remove duplicates
      &delete_duplicates(array => \@short_files);

      #populate textbox with short files
      $args{textbox}->delete(0, 'end');
      foreach my $entry (@short_files) {
         $args{textbox}->insert('end', $entry);
      }

      #merge selected files with overall selection array
      foreach (@{$selected_files}) {
         'reset' =~ m/reset/;
         m/\/([^\/]+)$/; #display only file_name
         my $short = $1;
         $all_files{$short} = $_;
      }

      #update display
      $args{textbox}->configure(-width => $width);
      $tl->update;
   }
}

sub remove_entry {
   my %args = @_;
   return if ($#short_files < 0);
   my ($current_index) = ${$args{textbox}}->curselection();
   #remove from shortlist array
   splice(@short_files, $current_index, 1);
   #remove from all files hash
   my $sel_key = ${$args{textbox}}->get($current_index);
   foreach my $key (keys %all_files) {
      if ($key =~ /\/$sel_key$/) {
         delete $all_files{$key};
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
