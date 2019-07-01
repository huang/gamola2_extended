#!/opt/ActivePerl-5.8/bin/perl

package Basics::error_popup;
use Tk;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&error);


my ($em, $message);


sub error {
   $message = shift;
   $em = MainWindow->new();
   $em->title('Error Message');
   $em->Label( -text => $message)->pack();
   $em->Button(-text =>'OK',
               -command => sub {
                                 $em->destroy;
                                 return (1);
                               }) -> pack();

   MainLoop;
}
1;