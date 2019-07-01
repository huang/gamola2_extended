#!/opt/ActivePerl-5.8/bin/perl5.8.8

#read directories without and return selection
#required input parameters for dir_select: main_window, label, ini_ref, auto_ini_ref,
#and directory as plain text (hash key from ref hashes; i.e. 'pfam_db_path')

#read file content
#required input parameters for slurp: main_window, directory and filename

package initialise::read_me;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&dir_select &slurp &slurp_cmd &select_blast_filter $dirselect);
use vars qw($dirselect);

use Tk;
use Tk::Toplevel;
use Tk::DirTree;
use Basics::Dirlistnofile;
use Cwd;


#local variables
use vars qw(%args $tl $fr $ok $newdir $dl $top $reselect $file);


sub dir_select {
   my %args = (directory => '',
               label     => 'Choose directory',
               @_);
   my ($curr_dir, $OK, $t, $d, $f);
   if (!Tk::Exists($t)) {
      $t = ${$args{main_window}}->Toplevel();
   } else {
      $t->state('normal');
   }
   $t->title($args{label});

   $ok = 0; # flag: "1" means OK, "-1" means cancelled

   # Create Frame widget before the DirTree widget, so it's always visible
   # if the window gets resized.
   $f = $t->Frame->pack(-fill => "x", -side => "bottom");

   if ($args{directory} =~ /\S+/ && defined (${$args{ini_ref}}{$args{directory}})) {
      $curr_dir = ${$args{ini_ref}}{$args{directory}};
   } elsif ($args{directory} =~ /\S+/ && defined (${$args{auto_ini_ref}}{$args{directory}})) {
      $curr_dir = ${$args{auto_ini_ref}}{$args{directory}};
   } else {
      $curr_dir = Cwd::cwd();
   }

   #-width => 35,
   $d = $t->Scrolled('DirTree',
           -scrollbars => 'osoe',
           -height => 20,
           -selectmode => 'browse',
           -exportselection => 1,
           -browsecmd => sub { $curr_dir = shift },

           # With this version of -command a double-click will
           # select the directory
           #-command   => sub { $ok = 1 },

           # With this version of -command a double-click will
           # open a directory. Selection is only possible with
           # the Ok button.
           -command   => sub { $d->opencmd($_[0]) },
          )->pack(-fill => "both", -expand => 1);
   # Set the initial directory
   $d->chdir($curr_dir);
   $d->configure(-width => (length($args{label}) * 2));

   $f->Button(-text => 'Ok',
         -command => sub { $ok =  1 })->pack(-side => 'left');
   $f->Button(-text => 'Cancel',
         -command => sub { $ok = -1 })->pack(-side => 'left');

   # You probably want to set a grab. See the Tk::FBox source code for
   # more information (search for grabCurrent, waitVariable and
   # grabRelease).
   $f->waitVariable(\$ok);
   if (defined (${$args{ini_ref}}{$args{directory}})) {
      ${$args{ini_ref}}{$args{directory}} = $curr_dir;
   } elsif (defined (${$args{auto_ini_ref}}{$args{directory}})) {
      ${$args{auto_ini_ref}}{$args{directory}} = $curr_dir;
   }

   $t->state('withdrawn');
   return ($curr_dir);
}



sub slurp {
   my %args = (quiet => 0,
               @_
              );
   my $timeout = 7000;
   if (defined ${$args{auto_ini_ref}}{msg_timeout}) {
      $timeout = ${$args{auto_ini_ref}}{msg_timeout};
   }

   local( $/, *ENTRY ) ;
   open( ENTRY, $args{directory}.'/'.$args{filename} ) or do
   {
      if ($args{quiet} == 0) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "Could not open file $args{filename} in directory $args{directory}.",
                                                       -bitmap  => 'info',
                                                       -buttons => ['OK']
                                                      );
         $error_msg->after($timeout, sub {$error_msg->Exit()});
         $error_msg-> Show();
      }
      return (0);
   };
   my $file = <ENTRY>;
   close ENTRY;
   unless ($file =~ /\w+/) {
      if ($args{quiet} == 0) {
         my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                       -text    => "File $args{filename} appears to be empty.",
                                                       -bitmap  => 'info',
                                                       -buttons => ['OK']
                                                       );
         $error_msg->after($timeout, sub {$error_msg->Exit()});
         $error_msg-> Show();
      }
      return (0);
   };
   return (\$file);
}

sub slurp_cmd { #as above, without graphical feedback
   my %args = (quiet => 0,
               @_
              );

   local( $/, *ENTRY ) ;
   open( ENTRY, $args{directory}.'/'.$args{filename} ) or do
   {
      if ($args{quiet} == 0) {print "Error: Could not open file $args{filename} in directory $args{directory}."};
      return (0);
   };
   my $file = <ENTRY>;
   close ENTRY;
   unless ($file =~ /\w+/) {
      if ($args{quiet} == 0) {print "Error: File $args{filename} appears to be empty."};
      return (0);
   };
   return (\$file);
}

sub selection {
   my %args = @_;
   my $top;

   if (! defined ($args{dirselect}) ) {
      $top = ${$args{window}}->Dialog(-title          => 'ERROR',
                                      -text           => "No directory selected\ntry again",
                                      -default_button => 'Close',
                                      -buttons        => ['Close'],
                                      -bitmap         => 'question') -> Show();

      if ($top eq 'Close') {
         if (Tk::Exists(${$args{window}})) {${$args{window}}->withdraw()};
         return (0);
      }
   } elsif  ( defined ($args{dirselect})) {
      if (Tk::Exists(${$args{window}})) {${$args{window}}->destroy()};
      return (1);
   } else {
      $top = ${$args{window}}->Dialog(-title          => 'ERROR',
                                      -text           => "No directory selected\ntry again",
                                      -default_button => 'Close',
                                      -buttons        => ['Close'],
                                      -bitmap         => 'question') -> Show();

      if ($top eq 'Close') {
         if (Tk::Exists(${$args{window}})) {${$args{window}}->withdraw()};
         return (0);
      }
   }
}

sub select_blast_filter {
   my %args = @_;

   (my $selected_file) = ${$args{main_window}}->getOpenFile(-initialdir => ${$args{auto_ini_ref}}{work_dir},
                                                             -title      => "Select file for Blast filter. White space parsing",
                                                             -multiple   => 0);
   unless (defined ($selected_file) && (-e $selected_file)) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "No file selected or non-accessible",
                                                    -bitmap  => 'error',
                                                    -buttons => ['ok']);
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      return;
   }

   #pull filepath apart
   'reset' =~ m/reset/;
   $selected_file =~ m/^(.+)\/([^\/]+)$/;
   my $local_path = $1;
   my $filename = $2;
   #read file to memory
   my $file_ref = &slurp(main_window   => $args{main_window},
                         progress_bar  => $args{progress_bar},
                         ini_ref       => $args{ini_ref},
                         auto_ini_ref  => $args{auto_ini_ref},
                         directory     => $local_path,
                         filename      => $filename
                        );

   #reformat textfile into a space separated single string
   $$file_ref =~ s/\s+/ /gs;
   if ($$file_ref =~ /[^\w ]/) {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'Error',
                                                    -text    => "Invalid characters found",
                                                    -bitmap  => 'error',
                                                    -buttons => ['ok']);
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      return;
   }
   $args{auto_ini_ref}{filter_blast} = $$file_ref;
   return(1);
}
1;
