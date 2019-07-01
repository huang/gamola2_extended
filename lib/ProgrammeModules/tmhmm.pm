#!/opt/ActivePerl-5.8/bin/perl

# This is version 2.0c of tmhmm


# Give ONE fasta file on cmdline OR use stdin
# A single sequence can be given WITHOUT the ID line (">ID")
# Such a sequence will be called "WEBSEQUENCE"

# OPTION PARSING ##########################################
package ProgrammeModules::tmhmm;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&new);
use vars qw();
use ProgrammeModules::tmhmmformat;

sub new {
   # $opt_basedir: basis directory for TMHMM package
   # $opt_scrdir: Script directory (defaults basedir/bin)
   # $opt_bindir: Bin directory (defaults basedir/bin)
   # $opt_libdir: Library directory (defaults basedir/lib)


   # Process options

   my ($opt_basedir, $opt_d, $opt_workdir, $opt_wwwdir, $opt_serverhome, $opt_html, $opt_short, $opt_plot, $opt_v1, $hash_key, $wd) = @_;
   my ($result, $format) = "";

   $opt_basedir .= "/" if ($opt_basedir !~ /\/$/ );
   $opt_scrdir = $opt_basedir."bin" if (!defined($opt_scrdir));
   $opt_bindir = $opt_basedir."bin" if (!defined($opt_bindir));
   $opt_libdir = $opt_basedir."lib" if (!defined($opt_libdir));

   ###########################################

   # Choose old model if requested
   if ( $opt_v1) { $modelfile="$opt_libdir/TMHMM1.0.model"; }
   else {          $modelfile="$opt_libdir/TMHMM2.0.model"; }

   $optfile="$opt_libdir/TMHMM2.0.options";

   # Debugging?
   if ($opt_d) { $err = ""; }
   else { $err = "2>/dev/null"; }

   # Programs to run
   $tmhmmformat = "$opt_scrdir/tmhmmformat.pl -workdir $wd -wwwdir $wd -serverhome $opt_serverhome";
   $decode = "$opt_bindir/decodeanhmm $modelfile $wd/$hash_key -f $optfile -plp";

   # options to $tmhmmformat
   $tmhmmformat .= " -html" if ($opt_html);
   $tmhmmformat .= " -v1" if ($opt_v1);
   if ( $opt_short ) {
       $tmhmmformat .= " -short -noplot";
   }
   elsif (!$opt_plot) {
       $tmhmmformat .= " -noplot";
   }

   #run decode
   $result = `$decode`;

   #run format
   $opt_html = 0;       # Make HTML output
   $opt_short = 0;      # Make short output format
   $opt_plot = 1;       # Make plots
   $opt_v1 = 0;         # Use old model (version 1)
   $opt_signal = 10;    # Cut-off used for producing a warning of possible
                        # Signal sequence
   $opt_Nterm = 60;     # Number of bases to consider for signal peptides
                        # in the Nterm


   $opt_workdir = $wd;  # Working dir.
   $opt_wwwdir = $wd;   # The place where the www server looks for files
                        # (The www name for the working dir)
   #$opt_serverhome = "$opt_serverhome"; # This is the location of various stationary files
                          # (Help files etc)

   $format= &tmhmmformat($opt_html, $opt_short, $opt_plot, $opt_v1, $opt_signal, $opt_Nterm, $opt_workdir, $opt_wwwdir, $opt_wwwdir, $opt_serverhome, $result);

   return ($format);

   ##################
   ##################

   # You MUST give a file handle as argument
   sub do_file {
      my $header = shift;
      $nodd=0;

      while ( @entry = &read_fasta($header) ) {
         if (!$entry[1]) { next; }
         $id = $entry[0];
         $id =~ s/>//;
         $id =~ s/\s.*$//;

         # Translate to upper case.
         $entry[1] =~ tr/a-z/A-Z/;
         $nodd = ( $entry[1] =~ s/[^ACDEFGHIKLMNPQRSTWVYBZX]/X/g ) ;

         print_entry(\@entry);
      }
      return ($id);
   }




   # print_entry ##########################################

   # Print an entry
   # First arg is the entry, second (optional) is the filehandle
   # Eg.
   #         print_entry(\@entry,\*OFILE);
   sub print_entry {
      #if (defined($_[1])) { $fh = $_[1]; }
      #else {$fh = \*STDOUT;}
      $fh = \*STDOUT;
      $linelen = 70;
      @entry = @{$_[0]};
      print $fh "$entry[0]\n";
      $movein = "";
      if (defined($entry[2])) {
         $movein = '  ';
         $labpref = '# '
      }
      $l = length($entry[1]);
      for ($i=0; $i<$l; $i += $linelen) {
         print $fh $movein . substr($entry[1],$i,$linelen) . "\n";
         if (defined($entry[2])) {
            print $fh $labpref . substr($entry[2],$i,$linelen) . "\n";
         }
      }
   }



   # read_fasta ##########################################

   # Read a FASTA entry
   sub read_fasta {
       my $seq = shift;
       @entry = ();

       # Individual characters in the sequence fitting this reg. exp.
       # are deleted
       $ignore = '[^a-zA-Z]';

       #catch header line
       $seq =~ m/(\>[^\n]*?)\n/;
       $entry[0] = $1;

       #catch seq
       $seq =~ m/\>[^\n]*?\n(.+)/s;
       $entry[1] = $1;

       #cleanup seq
       $entry[1] =~ s/\s//gs;
       $entry[1] =~ s/$ignore//g;
       $entry[1] =~ s/\W//g;

       return @entry;
   }
}

1;
