#!/opt/ActivePerl-5.8/bin/perl

#progress bar
#input arguments: title, optional text, optional time scale

use FindBin qw($Bin);
use lib "$Bin";
use Tk;
use StatusBar;

use strict;
use vars qw($Progress @colors $mw $sb $p $frame);

#my %args = @ARGV;
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

#define empty string for label if not provided
unless (defined ($ARGV[1])) {$ARGV[1] = ""};
#define default time scale if not provided
unless (defined ($ARGV[2])) {$ARGV[2] = '80'};

my $mw = new MainWindow (title => $ARGV[0]);
$sb    = $mw->StatusBar();
$frame = $mw->Frame(-borderwidth => 2, -relief => 'groove')->pack(-expand =>1, -side => 'top', -fill =>'x');

$frame->Label(
           -relief         => 'flat',
           -textvariable   => \$ARGV[1],
)->pack();

my $width = length ($ARGV[0]);
$sb->addLabel(
              -width          => $width,
              -relief         => 'flat',
              -textvariable   => \$ARGV[0],
);

my $length = length{$ARGV[1]};
if ($length <100) {$length = 100};
$p = $sb->addProgressBar(
              -from           => 0,
              -to             => 100,
              -blocks         => 100,
              -colors         => \@colors,
              -variable       => \$Progress,
);

$mw->repeat($ARGV[2], sub {
                             $Progress = 0 if (++$Progress > 100);
                          });

MainLoop();
