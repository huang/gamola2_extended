#!/opt/ActivePerl-5.8/bin/perl

#progress bar
#input arguments: main_window, progress_bar, timeout

package Basics::progress_bar;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;

$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&progress_bar_1 &progress_bar_2 &progress_bar_3 &update_pbar_1 &update_pbar_2 &update_pbar_3
             &hide_pbar_1 &hide_pbar_2 &hide_pbar_3 show_pbar_1 show_pbar_2 show_pbar_3
             &kill_all_pbar) ;
use vars qw();

use Tk;
use Basics::StatusBar;

use strict;
my ($p_tl_1, $p_progress_1, $p_label_1, $sb_1, $p_1, $frame_1, $label_1, $width_1, $previous_width_1,
    $p_tl_2, $p_progress_2, $p_label_2, $sb_2, $p_2, $frame_2, $label_2, $width_2, $previous_width_2,
    $p_tl_3, $p_progress_3, $p_label_3, $sb_3, $p_3, $frame_3, $label_3, $width_3, $previous_width_3
    );

my @colors = (   0, '#ff002a',  1, '#ff0014',  2, '#ff000a',  3, '#ff0500',  4, '#ff1000',
                 5, '#ff1b00',  6, '#ff3000',  7, '#ff3b00',  8, '#ff4600',  9, '#ff5100',
                10, '#ff6100', 11, '#ff7600', 12, '#ff8100', 13, '#ff8c00', 14, '#ff9700',
                15, '#ffa100', 16, '#ffbc00', 17, '#ffc700', 18, '#ffd200', 19, '#ffdd00',
                20, '#ffe700', 21, '#fffd00', 22, '#f0ff00', 23, '#e5ff00', 24, '#dbff00',
                25, '#d0ff00', 26, '#baff00', 27, '#afff00', 28, '#9fff00', 29, '#95ff00',
                30, '#8aff00', 31, '#74ff00', 32, '#6aff00', 33, '#5fff00', 34, '#54ff00',
                35, '#44ff00', 36, '#2eff00', 37, '#24ff00', 38, '#19ff00', 39, '#0eff00',
                40, '#03ff00', 41, '#00ff17', 42, '#00ff21', 43, '#00ff2c', 44, '#00ff37',
                45, '#00ff42', 46, '#00ff57', 47, '#00ff67', 48, '#00ff72', 49, '#00ff7d',
                50, '#00ff87', 51, '#00ff9d', 52, '#00ffa8', 53, '#00ffb8', 54, '#00ffc3',
                55, '#00ffcd', 56, '#00ffe3', 57, '#00ffee', 58, '#00fff8', 59, '#00faff',
                60, '#00eaff', 61, '#00d4ff', 62, '#00c9ff', 63, '#00bfff', 64, '#00b4ff',
                65, '#00a9ff', 66, '#008eff', 67, '#0083ff', 68, '#0079ff', 69, '#006eff',
                70, '#0063ff', 71, '#004eff', 72, '#003eff', 73, '#0033ff', 74, '#0028ff',
                75, '#001dff', 76, '#0008ff', 77, '#0200ff', 78, '#1200ff', 79, '#1d00ff',
                80, '#2800ff', 81, '#3d00ff', 82, '#4800ff', 83, '#5300ff', 84, '#5d00ff',
                85, '#6e00ff', 86, '#8300ff', 87, '#8e00ff', 88, '#9900ff', 89, '#a300ff',
                90, '#ae00ff', 91, '#c900ff', 92, '#d400ff', 93, '#df00ff', 94, '#e900ff',
                95, '#f400ff', 96, '#ff00f3', 97, '#ff00e3', 98, '#ff00d9', 99, '#ff00ce'
            );

sub progress_bar_1 {
   my %args = (title     => '',
               label     => '',
               progress  => '0',
               @_
              );

   $p_progress_1 = $args{progress};
   $p_label_1    = $args{label};

   $p_tl_1       = ${$args{main_window}}->Toplevel (title => $args{title});
   $p_tl_1       ->state('withdrawn');
   $sb_1         = $p_tl_1->StatusBar();
   $frame_1      = $p_tl_1->Frame(-borderwidth => 2, -relief => 'groove')->pack(-expand =>1, -side => 'top', -fill =>'x');

   #define length of bar
   unless (defined $width_1) {
      $width_1 = length($p_label_1) + 2 ;
      if ($width_1 < (length($args{title}) + 2) ) {
         $width_1 = (length($args{title}) +2);
      }
      if ($width_1 < 20) {$width_1 = 20};
   }
   $label_1 = $frame_1->Label(-width          => $width_1,
                              -relief         => 'flat',
                              -textvariable   => \$p_label_1,
   )->pack();
   $sb_1->addLabel(-width          => length($args{title}),
                   -relief         => 'flat',
                   -textvariable   => \$args{title},
   );

   $p_1 = $sb_1->addProgressBar(-from     => 0,
                                -to       => 100,
                                -blocks   => 100,
                                -colors   => \@colors,
                                -variable => \$p_progress_1,
   );

   return (1);
}

sub progress_bar_2 {
   my %args = (title     => '',
               label     => '',
               progress  => '0',
               @_
              );

   $p_progress_2 = $args{progress};
   $p_label_2    = $args{label};

   $p_tl_2  = ${$args{main_window}}->Toplevel (title => $args{title});
   $p_tl_2  ->state('withdrawn');
   $sb_2    = $p_tl_2->StatusBar();
   $frame_2 = $p_tl_2->Frame(-borderwidth => 2, -relief => 'groove')->pack(-expand =>1, -side => 'top', -fill =>'x');

   #define length of bar
   unless (defined $width_2) {
      $width_2 = length($p_label_2) + 2 ;
      if ($width_2 < (length($args{title}) + 2) ) {
         $width_2 = (length($args{title}) +2);
      }
      if ($width_2 < 20) {$width_2 = 20};
   }
   $label_2 = $frame_2->Label(-width          => $width_2,
                              -relief         => 'flat',
                              -textvariable   => \$p_label_2,
   )->pack();
   $sb_2->addLabel(-width          => length($args{title}),
                   -relief         => 'flat',
                   -textvariable   => \$args{title},
   );

   $p_2 = $sb_2->addProgressBar(-from     => 0,
                                -to       => 100,
                                -blocks   => 100,
                                -colors   => \@colors,
                                -variable => \$p_progress_2,
   );

   return (1);
}

sub progress_bar_3 {
   my %args = (title     => '',
               label     => '',
               progress  => '0',
               @_
              );

   $p_progress_3 = $args{progress};
   $p_label_3    = $args{label};

   $p_tl_3  = ${$args{main_window}}->Toplevel (title => $args{title});
   $p_tl_3  ->state('withdrawn');
   $sb_3    = $p_tl_3->StatusBar();
   $frame_3 = $p_tl_3->Frame(-borderwidth => 2, -relief => 'groove')->pack(-expand =>1, -side => 'top', -fill =>'x');


   #define length of bar
   unless (defined $width_3) {
      $width_3 = length($p_label_3) + 2 ;
      if ($width_3 < (length($args{title}) + 2) ) {
         $width_3 = (length($args{title}) +2);
      }
      if ($width_3 < 20) {$width_3 = 20};
   }
   $label_3 = $frame_3->Label(-width          => $width_3,
                              -relief         => 'flat',
                              -textvariable   => \$p_label_3,
   )->pack();
   $sb_3->addLabel(-width          => length($args{title}),
                   -relief         => 'flat',
                   -textvariable   => \$args{title},
   );

   $p_3 = $sb_3->addProgressBar(-from     => 0,
                                -to       => 100,
                                -blocks   => 100,
                                -colors   => \@colors,
                                -variable => \$p_progress_3,
   );

   return (1);
}

sub progress_bar_auto {
   my %args = @_;
   my $timeout;
   unless (defined ($timeout)) {$timeout = 0};
   my @colors = (   0, '#ff002a',  1, '#ff0014',  2, '#ff000a',  3, '#ff0500',  4, '#ff1000',
                    5, '#ff1b00',  6, '#ff3000',  7, '#ff3b00',  8, '#ff4600',  9, '#ff5100',
                   10, '#ff6100', 11, '#ff7600', 12, '#ff8100', 13, '#ff8c00', 14, '#ff9700',
                   15, '#ffa100', 16, '#ffbc00', 17, '#ffc700', 18, '#ffd200', 19, '#ffdd00',
                   20, '#ffe700', 21, '#fffd00', 22, '#f0ff00', 23, '#e5ff00', 24, '#dbff00',
                   25, '#d0ff00', 26, '#baff00', 27, '#afff00', 28, '#9fff00', 29, '#95ff00',
                   30, '#8aff00', 31, '#74ff00', 32, '#6aff00', 33, '#5fff00', 34, '#54ff00',
                   35, '#44ff00', 36, '#2eff00', 37, '#24ff00', 38, '#19ff00', 39, '#0eff00',
                   40, '#03ff00', 41, '#00ff17', 42, '#00ff21', 43, '#00ff2c', 44, '#00ff37',
                   45, '#00ff42', 46, '#00ff57', 47, '#00ff67', 48, '#00ff72', 49, '#00ff7d',
                   50, '#00ff87', 51, '#00ff9d', 52, '#00ffa8', 53, '#00ffb8', 54, '#00ffc3',
                   55, '#00ffcd', 56, '#00ffe3', 57, '#00ffee', 58, '#00fff8', 59, '#00faff',
                   60, '#00eaff', 61, '#00d4ff', 62, '#00c9ff', 63, '#00bfff', 64, '#00b4ff',
                   65, '#00a9ff', 66, '#008eff', 67, '#0083ff', 68, '#0079ff', 69, '#006eff',
                   70, '#0063ff', 71, '#004eff', 72, '#003eff', 73, '#0033ff', 74, '#0028ff',
                   75, '#001dff', 76, '#0008ff', 77, '#0200ff', 78, '#1200ff', 79, '#1d00ff',
                   80, '#2800ff', 81, '#3d00ff', 82, '#4800ff', 83, '#5300ff', 84, '#5d00ff',
                   85, '#6e00ff', 86, '#8300ff', 87, '#8e00ff', 88, '#9900ff', 89, '#a300ff',
                   90, '#ae00ff', 91, '#c900ff', 92, '#d400ff', 93, '#df00ff', 94, '#e900ff',
                   95, '#f400ff', 96, '#ff00f3', 97, '#ff00e3', 98, '#ff00d9', 99, '#ff00ce' );

   $args{p_bar} = ${$args{progress_bar}}->ProgressBar(-width => 10,
                                                -length =>350,
                                                -from => 0,
                                                -to => 100,
                                                -blocks => 100,
                                                -colors => \@colors,
                                                -variable => \$timeout) ->grid (-row => 10, -column => 0, -sticky => 'ew');
   #${args{main_window}} -> update;

   for ($timeout = 0; $timeout <= 100; $timeout ++) {
      if ($timeout > 99) {$timeout = 0};
      ${$args{progress_bar}}->update;
      ${$args{progress_bar}}->after (100);
   }

}

sub update_pbar_1 {
   my %args =  (title    => '',
                label    => '',
                progress => 0,
                @_);

   $p_tl_1->configure(-title => $args{title});

   unless (defined ($width_1) && $width_1 =~ /\d+/) {
      $width_1 = 20;
   }
   unless (defined $previous_width_1) {
      $previous_width_1 = 0;
   }
   $p_progress_1 = $args{progress};
   $p_label_1    = $args{label};
   $width_1      = length($p_label_1) + 2 ;
   if ($width_1 != $previous_width_1) {
      if ($width_1 < (length($args{title}) + 2) ) {
         $width_1 = (length($args{title}) +2);
      }
      if ($width_1 < 20) {$width_1 = 20};
      $label_1->configure(-width => $width_1);
      $previous_width_1 = $width_1;
   }
   $label_1->configure(-textvariable   => \$p_label_1);
   $p_tl_1->update;
}

sub update_pbar_2 {
   my %args =  (title    => '',
                label    => '',
                progress => 0,
                @_);

   $p_tl_2->configure(-title => $args{title});

   unless (defined ($width_2) && $width_2 =~ /\d+/) {
      $width_2 = 20;
   }
   unless (defined $previous_width_2) {
      $previous_width_2 = 0;
   }
   $p_progress_2 = $args{progress};
   $p_label_2    = $args{label};
   $width_2      = length($p_label_2) + 2 ;
   if ($width_2 != $previous_width_2) {
      if ($width_2 < (length($args{title}) + 2) ) {
         $width_2 = (length($args{title}) +2);
      }
      if ($width_2 < 20) {$width_2 = 20};
      $label_2->configure(-width => $width_2);
      $previous_width_2 = $width_2;
   }
   $label_2->configure(-textvariable   => \$p_label_2);
   $p_tl_2->update;
}

sub update_pbar_3 {
   my %args =  (title    => '',
                label    => '',
                progress => 0,
                @_);

   $p_tl_3->configure(-title => $args{title});

   unless (defined ($width_3) && $width_3 =~ /\d+/) {
      $width_3 = 20;
   }
   unless (defined $previous_width_3) {
      $previous_width_3 = 0;
   }
   $p_progress_3 = $args{progress};
   $p_label_3    = $args{label};
   $width_3      = length($p_label_3) + 2 ;
   if ($width_3 != $previous_width_3) {
      if ($width_3 < (length($args{title}) + 2) ) {
         $width_3 = (length($args{title}) +2);
      }
      if ($width_3 < 20) {$width_3 = 20};
      $label_3->configure(-width => $width_3);
      $previous_width_3 = $width_3;
   }
   $label_3->configure(-textvariable   => \$p_label_3);
   $p_tl_3->update;
}

sub hide_pbar_1 {
   my $state = $p_tl_1->state();
   if ($state eq 'withdrawn') {
      return(1);
   } else {
      $p_tl_1->state('withdrawn');
      $p_tl_1->update;
      return(1);
   }
}

sub hide_pbar_2 {
   my $state = $p_tl_2->state();
   if ($state eq 'withdrawn') {
      return(1);
   } else {
      $p_tl_2->state('withdrawn');
      $p_tl_2->update;
      return(1);
   }
}

sub hide_pbar_3 {
   my $state = $p_tl_3->state();
   if ($state eq 'withdrawn') {
      return(1);
   } else {
      $p_tl_3->state('withdrawn');
      $p_tl_3->update;
      return(1);
   }
}

sub show_pbar_1 {
   my $state = $p_tl_1->state();
   if ($state eq 'normal') {
      return(1);
   } else {
      $p_tl_1->state('normal');
      $p_tl_1->update;
      return(1);
   }
}

sub show_pbar_2 {
   my $state = $p_tl_2->state();
   if ($state eq 'normal') {
      return(1);
   } else {
      $p_tl_2->state('normal');
      $p_tl_2->update;
      return(1);
   }
}

sub show_pbar_3 {
   my $state = $p_tl_3->state();
   if ($state eq 'normal') {
      return(1);
   } else {
      $p_tl_3->state('normal');
      $p_tl_3->update;
      return(1);
   }
}

sub kill_all_pbar {
   if (Tk::Exists $p_tl_1) {
      $p_tl_1->destroy();
   }
   if (Tk::Exists $p_tl_2) {
      $p_tl_2->destroy();
   }
   if (Tk::Exists $p_tl_3) {
      $p_tl_3->destroy();
   }
}

1;