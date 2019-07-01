#!/opt/ActivePerl-5.8/bin/perl5.8.8

package initialise::blast_type;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&blast_select &file_select $dirselect);
use vars qw($dirselect);

use Tk;
use Basics::Dirlist;
use initialise::error_popup;
use Cwd;


#local variables
use vars qw($title $mw $fr $ok $newdir $dl $top $reselect);

sub blast_select {
   my $title = shift;
   my $mw;
   my $reselect = 0;
   close (STDOUT);
   open (STDOUT, ">/dev/null");


   if (defined $mw && Tk::Exists($mw)) {
       $mw->state('normal');
   } else {
      $mw = MainWindow->new;
   }

   $mw->title($title);
   $mw->focus;
   my $fr = $mw->Frame->pack(-fill => "x", -side => "bottom");

   undef $dirselect;
   my $ok = 0;
   my $newdir = cwd();

   $dl = $mw->Scrolled
       ('Dirlist',
        -scrollbars => 'osoe',
        -width => 25,
        -height => 25,

        -browsecmd => sub {
                            $newdir = shift;
                            if (-d $newdir) {
                                  $dl->configure(-directory => $newdir);
                                  $dirselect = $newdir;
                            }
                          });

   $dl->pack(-fill => "both", -expand => 1);
   $fr->Button(-text => 'Ok',
         -command => sub { $ok =  1 })->pack(-side => 'left');
   $fr->Button(-text => 'Cancel',
         -command => sub { $ok = -1 })->pack(-side => 'left');

   # wait until a file is selected
   $fr->waitVariable(\$ok);
   unless (defined ($dirselect) || $ok == -1 ) {
      if (defined $top && Tk::Exists($top)) {
          $top->state('normal');
      } else {
         $top = $mw->Toplevel(-title=>'ERROR');
      }
      $top->Label (-text    => "No directory selected\ntry again")->pack();
      $top->Button(-text    => "Close",
                   -command => sub {$mw->state('withdrawn');
                                    $top->state('withdrawn');
                                    &dir_select ($title);
                                  })->pack();
      $top->focus;
   } else {
      $mw->state('withdrawn');
      return ($dirselect);
   }
}

