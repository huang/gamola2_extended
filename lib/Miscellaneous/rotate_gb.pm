#!/opt/ActivePerl-5.8/bin/perl

#rotate_gb: rotates a selected Genbank fileno clockwise or counterclockwise
#input arguments: main_window, ini_ref, auto_ini_ref,


package Miscellaneous::rotate_gb;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&rotate_gb);
use vars qw();

use Tk;
use Tk::ProgressBar;
use Tk::LabEntry;
use initialise::read_me              qw(:DEFAULT);
use ProgrammeModules::sequence       qw(:DEFAULT);
use ProgrammeModules::genbank_parser qw(:DEFAULT);
use Cwd;

#local variables
my (%args, $rotation, $c_rot, $cc_rot, $width, $tl, $new_start, $selected_file, $short_form,
    $header, $source, $core, $DNAseq, $start_ref, $feature_list_ref, $genbank_ref, $counter,
    $bar_progress);

#initialise

sub rotate_gb {
   my %args = @_;
   $width = 40;
   my ($file_lb, $clockwise, $counterclockwise, $progress_bar, $frame_bottom, $frame_left,
       $progress) = '';
   $bar_progress = 0;

   if (defined $tl && Tk::Exists($tl)) {
      $tl->state('normal');
      $selected_file = '';
      $progress_bar  = $tl	->Frame(-borderwidth => 2, -relief => 'groove');
      $frame_bottom  = $tl	->Frame(-borderwidth => 2, -relief => 'groove');
      $progress      = $progress_bar->ProgressBar(-width    => 20,
                                                  -length   => 300,
                                                  -from     => 0,
                                                  -to       => 100,
                                                  -blocks   => 100,
                                                  -variable => \$bar_progress
                                                  );

      $frame_left    = $tl->Frame();
   } else {
      $selected_file = '';
      $tl = ${$args{main_window}}->Toplevel(-title 				=> 'Rotate Genbank file',
                                           );
      $progress_bar  = $tl	->Frame(-borderwidth => 2, -relief => 'groove')->pack(-side => 'bottom',-fill => 'x');
      $frame_bottom  = $tl	->Frame(-borderwidth => 2, -relief => 'groove')->pack (-side => 'bottom', -fill => 'x');
      $progress      = $progress_bar->ProgressBar(-width    => 20,
                                                  -length   => 300,
                                                  -from     => 0,
                                                  -to       => 100,
                                                  -blocks   => 100,
                                                  -variable => \$bar_progress
                                                  )
                             -> grid (-row => 0, -column => 0, -sticky => 'esnw');
      $frame_bottom	 			->Label(-text => "")	-> grid (-row => 0, -column => 0, -sticky => 'esnw');
      $frame_left    = $tl->Frame()->pack (-side => 'left', -fill => 'y');
   }

   #entries into frame
   $frame_left                		->Label(-text => 'Parameters for Genbank rotation')-> grid (-row => 0, -column => 0, -columnspan => 2, -sticky => 'esnw',-pady => 25); #dummy

   $clockwise			= $frame_left	->Checkbutton(-text   	=> 'Shift forward',
                                                  -variable	=> \$c_rot,
                                                  -command 	=> sub {
                                                                if ($c_rot == 1) {
                                                                   $cc_rot = 0;
                                                                   $rotation = 'clockwise';
                                                                } elsif ($c_rot == 0) {
                                                                   $cc_rot = 1;
                                                                   $rotation = 'counterclockwise';
                                                                }
                                                              }
                                    )
   -> grid ($counterclockwise	= $frame_left->Checkbutton(-text     => 'Shift backwards',
                                                         -variable => \$cc_rot,
                                                         -command  => sub {
                                                                            if ($cc_rot == 1) {
                                                                               $c_rot = 0;
                                                                               $rotation = 'counterclockwise';
                                                                            } elsif ($cc_rot == 0) {
                                                                               $c_rot = 1;
                                                                               $rotation = 'clockwise';
                                                                            }
                                                                          }
                                                        )
        )
   -> grid (-row => 1, -sticky => 'nsew');
   #$frame_left	->Label(-text	=> 'Enter new startposition')
   #                                 -> grid (-row => 2, -column => 0);
   my $new_pos			= $frame_left   ->LabEntry(-label        => "Enter new start-position",
                                               -textvariable 	=> \$new_start,
                                               -labelPack    => [ -side => 'left' ]
                                    )-> grid (-row => 2, -column => 0, -columnspan => 2, -sticky => 'nsew');

   my $select_gb		= $frame_left ->Button(-text => 'Select Genbank file',
                                        -command => sub {($selected_file) = ${$args{main_window}}->getOpenFile(-initialdir => ${$args{auto_ini_ref}}{work_dir},
                                                                                                                      -title      => 'Select files for custom database',
                                                                                                                      -multiple   => 0);
                                                          $selected_file =~ m/.+\/(.+)$/;
                                                          $short_form = $1;
                                                         }
                                    )#-> grid (-row => 3, -column => 0);
   ->grid($frame_left ->Label(-textvariable	=> \$short_form,
                     -relief				=> 'sunken',
                     -width 				=> $width),
          -sticky => 'nsew', -row => 3);

   $frame_left->Button(-text => "Rotate selected Genbank file",
                       -command => sub {  #first, check all variables
                                          $frame_bottom->Label(-text=>"Verifying parameters") ->grid (-row => 0, -column => 0, -sticky => 'ew');
                                          ${$args{main_window}}->update;
                                          my ($status) = &check_variables(main_window  => $args{main_window},
                                                                          ini_ref      => $args{ini_ref},
                                                                          auto_ini_ref => $args{auto_ini_ref},
                                                                          progress_bar => \$frame_bottom,
                                                                          rotation     => $rotation,
                                                                          new_start    => $new_start,
                                                                          gb_file      => $selected_file);
                                          if ($status == 0) {
                                             $tl->destroy;
                                             $selected_file = '';
                                             &rotate_gb(main_window  => $args{main_window},
                                                        ini_ref      => $args{ini_ref},
                                                        auto_ini_ref => $args{auto_ini_ref},
                                                        progress_bar => \$frame_bottom,
                                                        );
                                             return;
                                          }
                                          #then parse GB file
                                          $frame_bottom->Label(-text=>"Parsing Genbank file") ->grid (-row => 0, -column => 0, -sticky => 'ew');
                                          ${$args{main_window}}->update;
                                          ($header, $source, $core, $DNAseq, $start_ref, $feature_list_ref, $genbank_ref, $counter) = gb_parser(main_window  => $args{main_window},
                                                                                                                                                ini_ref      => $args{ini_ref},
                                                                                                                                                auto_ini_ref => $args{auto_ini_ref},
                                                                                                                                                progress_bar => \$frame_bottom,
                                                                                                                                                gb_file		 => $selected_file
                                                                                                                                                );
                                          $frame_bottom->Label(-text=>"Rotating Genbank features") ->grid (-row => 0, -column => 0, -sticky => 'ew');
                                          ${$args{main_window}}->update;
                                          #finally, rotate -> check if features are cut, create joined features
                                          my ($rotated_genbank_ref, $rotated_feature_list_ref, $status_1) = &rotate(main_window      => $args{main_window},
                                                                                                                    ini_ref          => $args{ini_ref},
                                                                                                                    auto_ini_ref     => $args{auto_ini_ref},
                                                                                                                    progress_bar     => \$frame_bottom,
                                                                                                                    start_ref        => $start_ref,
                                                                                                                    feature_list_ref => $feature_list_ref,
                                                                                                                    genbank_ref      => $genbank_ref,
                                                                                                                    DNAlength        => length($DNAseq),
                                                                                                                    DNAseq_ref       => \$DNAseq,
                                                                                                                    rotation         => $rotation,
                                                                                                                    new_start        => $new_start,
                                                                                                                   );
                                          if ($status_1 == 0) {
                                             $tl->destroy;
                                             $selected_file = '';
                                             &rotate_gb(main_window  => $args{main_window},
                                                        ini_ref      => $args{ini_ref},
                                                        auto_ini_ref => $args{auto_ini_ref},
                                                        );
                                             return;
                                          }
                                          #compile new GB file
                                          $frame_bottom->Label(-text=>"Compiling new Genbank file") ->grid (-row => 0, -column => 0, -sticky => 'ew');
                                          ${$args{main_window}}->update;
                                          my ($new_gb_ref) = &compile_gb(main_window      => \$tl,
                                                                         ini_ref          => $args{ini_ref},
                                                                         auto_ini_ref     => $args{auto_ini_ref},
                                                                         progress_bar     => \$frame_bottom,
                                                                         p                => \$progress,
                                                                         update_allowed   => 0,
                                                                         genbank_source   => $source,
                                                                         genbank_header   => $header,
                                                                         feature_list_ref => $rotated_feature_list_ref,
                                                                         genbank_ref      => $rotated_genbank_ref,
                                                                         DNAseq           => $DNAseq,
                                                                         from_rotate      => 1
                                                                        );

                                          #write new Genbank file
                                          $frame_bottom->Label(-text=>"Compiling new Genbank file") ->grid (-row => 0, -column => 0, -sticky => 'ew');
                                          ${$args{main_window}}->update;
                                          #modify filename
                                          if ($selected_file =~ /\.gb[k]?$/) {
                                             $selected_file =~ s/\.gb[k]?$/_rotated\.gb/;
                                          } else {
                                             $selected_file .= '_rotated';
                                          }
                                          open WRITE, "+>$selected_file" or do {
                                             ${$args{main_window}}->messageBox(-title   => 'Error',
                                                                               -message => "Could not open file $selected_file for writing.",
                                                                               -type    => 'OK',
                                                                               -icon    => 'error');
                                             $tl->destroy;
                                             $selected_file = '';
                                             &rotate_gb(main_window  => $args{main_window},
                                                        ini_ref      => $args{ini_ref},
                                                        auto_ini_ref => $args{auto_ini_ref},
                                                        );
                                             return;
                                          };
                                          print WRITE $$new_gb_ref;
                                          close WRITE;
                                          undef $new_gb_ref;

                                          #success
                                          ${$args{main_window}}->messageBox(-title   => 'Success',
                                                                            -message => "Successfully rotated Genbank file.",
                                                                            -type    => 'OK',
                                                                            -icon    => 'info');
                                          #$tl->state('withdrawn');
                                          $tl->destroy;
                                          $selected_file = '';
                                          return (1);
                                       }
                       )

   ->grid($frame_left->Button(-text => "Cancel",
                              -command => sub {
                                               #$tl->state('withdrawn');
                                               $tl->destroy;
                                               return (1);
                                              }
                             ),
         -row => 4, -sticky => 'nsew', -pady => 25);
}

sub rotate {
   my %args = @_;
   my ($new_left_bd, $new_right_bd, %rotated_genbank, @rotated_feature_list, $counter);
   $counter = 0;
   my $percent_done = $#{$args{feature_list_ref}};
   ${$args{progress_bar}}->Label(-text=>'Rotating Genbank features') ->grid (-row => 0, -column => 0, -sticky => 'ew');
   ${$args{main_window}}->update;

   #set up new start position depending on rotation orientation
   if ($args{rotation} eq 'counterclockwise') {
      $args{new_start} = $args{DNAlength} - $args{new_start};
   }

   #rotate all features->check for borderbreach
   foreach my $key (@{$args{feature_list_ref}}) {
      #update progress bar
      $counter++;
      if ($counter % 20 == 0) {
         $bar_progress = ($counter / $percent_done) * 100;
         ${$args{main_window}}->update;
      }

      my (@boundaries, @rotated_boundary, $left_bd, $right_bd, $orientation);
      'reset' =~ m/reset/;

      #catch boundary line from entry
      ${$args{genbank_ref}}{$key} =~ m/^([^\/]*?)\//s;
      my $boundaries = $1;
      #no hit? maybe only a one liner anyway, try working on whole entry
      unless (defined $boundaries && $boundaries =~ /\w+/) {
         $boundaries = ${$args{genbank_ref}}{$key};
      }
      #remove leading characters
      $boundaries =~ s/^\s*\S+\s+\D+//;
      #remove trailing characters
      $boundaries =~ s/\D+$//;
      #remove whitespaces
      $boundaries =~ s/\s+//gs;

      @boundaries = ($boundaries =~ m/(\d+\D*?\d+)/g);
      if ($#boundaries < 0) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not parse feature boundaries from key -> $key .\n...${$args{genbank_ref}}{$key} ...",
                                           -type    => 'OK',
                                           -icon    => 'error');
         return (undef,undef,0);
      }

      foreach my $boundary (@boundaries) {
         my ($left_bd, $right_bd,$new_left_proximal, $new_right_proximal, $new_left_distal, $new_right_distal);
         'reset'=~m/reset/;
         $boundary =~ m/(\d+)\D*?(\d+)/;
         $left_bd  = $1;
         $right_bd = $2;
         unless (defined $right_bd && defined $left_bd) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not parse feature boundaries from key -> $key .\nBoundary: $boundary...",
                                              -type    => 'OK',
                                              -icon    => 'error');
            return (undef,undef,0);
         }
         if ($left_bd <= $args{new_start} && $right_bd >= $args{new_start}) { #new start within feature?
            $new_left_proximal  = 1;
            $new_right_proximal = $right_bd - $args{new_start};
            $new_left_distal    = $args{DNAlength} - ($args{new_start} - $left_bd);
            $new_right_distal   = $args{DNAlength};
         } else { #regular rotation
            if ($left_bd < $args{new_start} && $right_bd <= $args{new_start}) {
               $new_left_proximal  = ($args{DNAlength} - $args{new_start}) + $left_bd;
               $new_right_proximal = ($args{DNAlength} - $args{new_start}) + $right_bd;
            } elsif ($left_bd >= $args{new_start} && $right_bd > $args{new_start}) {
               $new_left_proximal  = ($left_bd  - $args{new_start});
               $new_right_proximal = ($right_bd - $args{new_start});
            }
            $new_left_distal  = '-';
            $new_right_distal = '-';
         }

         push (@rotated_boundary, $new_left_proximal.'..'.$new_right_proximal);
         if ($new_left_distal =~ /\d+/) {
            push (@rotated_boundary, $new_left_distal.'..'.$new_right_distal);
         }
      }
      #recreate boundaries
      if ($#rotated_boundary == 0) { #one entry only
         ${$args{genbank_ref}}{$key} =~ s/^(\s*\S+\s+\D+)\d+\D*?\d+/$1$rotated_boundary[0]/;
      } else {
         my $new_boundary = join(',', @rotated_boundary);
         #one liner?
         if (${$args{genbank_ref}}{$key} !~ /\s+\//s) {
            ${$args{genbank_ref}}{$key} =~ s/^(\s*\S+\s+\D+).*/$1$new_boundary\)\n/s;
         } else {
            #${$args{genbank_ref}}{$key} =~ s/^(\s*\D+)(join\D*?)?\(?[^\/]*?\//$1 join\($new_boundary\)\n                     \//s;
            ${$args{genbank_ref}}{$key} =~ s/^(\s*\S+\s+\D+)[^\/]*?\//$1 join\($new_boundary\)\n                     \//s;
            ${$args{genbank_ref}}{$key} =~ s/complement\(join\( join/complement\(join/;
            ${$args{genbank_ref}}{$key} =~ s/join\( join/join/;
            if (${$args{genbank_ref}}{$key} =~ m/complement\(join/) {
                ${$args{genbank_ref}}{$key} =~ s/\)+/\)\)/;
            }
         }
      }
      #change key
      my $transfer = ${$args{genbank_ref}}{$key};
      'reset' =~ m/reset/;
      $rotated_boundary[0] =~ m/(\d+)\D*?(\d+)/;
      $left_bd = $1;
      $right_bd = $2;
      unless (defined $right_bd && defined $left_bd) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Could not parse feature boundaries from key -> $key .\nRotated boundary: $rotated_boundary[0]...",
                                           -type    => 'OK',
                                           -icon    => 'error');
         return (undef,undef,0);
      }
      $key =~ s/^\d+_/$left_bd\_/;
      $key =~ s/^(\d+_)\d+_/$1$right_bd\_/;
      $rotated_genbank{$key} = $transfer;
      push (@rotated_feature_list, $key);
      undef @boundaries;
      undef @rotated_boundary;

      # {
         # $key =~ m/^(\d+)_(\d+)_/;
         # $left_bd = $1;
         # $right_bd = $2;
         # unless (defined $right_bd && defined $left_bd) {
            # ${$args{main_window}}->messageBox(-title   => 'Error',
                                              # -message => "Could not parse feature boundaries from key -> $key .",
                                              # -type    => 'OK',
                                              # -icon    => 'error');
            # return (undef,undef,0);
         # }
         # #${$args{progress_bar}}->Label(-text=>"Rotating Genbank features at position $left_bd") ->grid (-row => 0, -column => 0, -sticky => 'ew');
         # #${$args{main_window}}->update;
   #
         # #rotate features
         # if ($left_bd < $args{new_start} && $right_bd > $args{new_start}) { #new start within feature?
            # my $new_left_proximal = 1;
            # my $new_right_proximal = $right_bd - $args{new_start};
            # my $new_left_distal = $args{DNAlength} - ($args{new_start} - $left_bd);
            # my $new_right_distal = $args{DNAlength};
   #
            # #modify entry and check for sense and antisense direction for bounday ordering
   #
            # ${$args{genbank_ref}}{$key} =~ s/$left_bd\.\.$right_bd/join\($new_left_proximal\.\.$new_right_proximal\,$new_left_distal\.\.$new_right_distal\)/;
            # my $transfer = ${$args{genbank_ref}}{$key};
            # $key =~ s/$left_bd/$new_left_proximal/;
            # $key =~ s/$right_bd/$new_right_proximal/;
            # #transfer into new datasets
            # $rotated_genbank{$key} = $transfer;
            # push (@rotated_feature_list, $key);
            # #push (@rotated_start, $new_left_proximal);
         # } else { #regular rotation
            # if ($left_bd < $args{new_start} && $right_bd < $args{new_start}) {
               # $new_left_bd  = ($args{DNAlength} - $args{new_start}) + $left_bd;
               # $new_right_bd = ($args{DNAlength} - $args{new_start}) + $right_bd;
            # } elsif ($left_bd > $args{new_start} && $right_bd > $args{new_start}) {
               # $new_left_bd  = ($left_bd  - $args{new_start});
               # $new_right_bd = ($right_bd - $args{new_start});
            # }
            # #modify entry
            # ${$args{genbank_ref}}{$key} =~ s/$left_bd/$new_left_bd/;
            # ${$args{genbank_ref}}{$key} =~ s/$right_bd/$new_right_bd/;
            # my $transfer = ${$args{genbank_ref}}{$key};
            # $key =~ s/$left_bd/$new_left_bd/;
            # $key =~ s/$right_bd/$new_right_bd/;
            # #transfer into new datasets
            # $rotated_genbank{$key} = $transfer;
            # push (@rotated_feature_list, $key);
            # #push (@rotated_start, $new_left_bd);
         # }
      # }
   }

   #rotate DNA sequence
   my $new_distal_DNA = substr(${$args{DNAseq_ref}},0,$args{new_start});
   my $new_proximal_DNA = substr(${$args{DNAseq_ref}},$args{new_start});
   ${$args{DNAseq_ref}} = $new_proximal_DNA.$new_distal_DNA;

   return (\%rotated_genbank, \@rotated_feature_list, '1');

}

sub check_variables {
   my %args = @_;
   my ($content, $DNAseq);

   #check status of rotation
   ${$args{progress_bar}}->Label(-text=>"Verifying orientation") ->grid (-row => 0, -column => 0, -sticky => 'ew');
   ${$args{main_window}}->update;
   unless (defined $rotation && ($rotation eq 'clockwise' || $rotation eq 'counterclockwise')) {
           ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "No orientation selected",
                                              -type    => 'OK',
                                              -icon    => 'error');
           return (0);
   }

   #check if selected file exists and is Genbank file -> also extract DNA sequence
   ${$args{progress_bar}}->Label(-text=>"Verifying Genbank file") ->grid (-row => 0, -column => 0, -sticky => 'ew');
   ${$args{main_window}}->update;
   unless (defined $args{gb_file} && -e $args{gb_file}) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Selected file $args{gb_file} does not exist",
                                        -type    => 'OK',
                                        -icon    => 'error');
      return (0);
   }
   {
      local ($/, *READ);
      open READ, "<$args{gb_file}" or return (0);
      $content = <READ>;
      close READ;
   }

   #gb format - quick check
   if ($content !~ /^locus/i) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Selected file $args{gb_file} is not a Genbank file",
                                        -type    => 'OK',
                                        -icon    => 'error');
      return (0);
   }

   #catch DNA sequence
   ${$args{progress_bar}}->Label(-text=>"Verifying new start position") ->grid (-row => 0, -column => 0, -sticky => 'ew');
   ${$args{main_window}}->update;
   'reset' =~ /reset/;
   $content =~ m/\n((BASE COUNT|ORIGIN).*)/is;
   $DNAseq = $1;
   unless (defined $DNAseq) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Could not identify DNA sequence in file $args{gb_file}",
                                        -type    => 'OK',
                                        -icon    => 'error');
      return (0);
   }
   #clean up DNA seq
   $DNAseq =~ s/ORIGIN//s;
   $DNAseq =~ s/BASE COUNT[^\n]*?\n//is;
   $DNAseq =~ s/\d//gs;
   $DNAseq =~ s/\s//gs;
   $DNAseq =~ s/\\//g;
   undef $content;

   #check if new start is larger than length of DNAseq
   if ($args{new_start} =~ /\D/igs) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Selected new start position is not a number: $args{new_start}",
                                        -type    => 'OK',
                                        -icon    => 'error');
      return (0);
   }
   if ($args{new_start} > length($DNAseq) || $args{new_start} < 0) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Selected new start position is outside Genbank boundaries: $args{new_start}",
                                        -type    => 'OK',
                                        -icon    => 'error');
      return (0);
   }

   if ($args{new_start} == 0) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "No rotation possible with a value of '0'",
                                        -type    => 'OK',
                                        -icon    => 'error');
      return (0);
   }
   ${$args{progress_bar}}->Label(-text=>"Verification OK") ->grid (-row => 0, -column => 0, -sticky => 'ew');
   ${$args{main_window}}->update;
   return (1);
}

1;
