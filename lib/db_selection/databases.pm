#!/opt/ActivePerl-5.8/bin/perl

#select databases for Blast

package db_selection::databases;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
use Cwd;
use Tk;

$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&select_glimmer_model &select_blast_db &select_cog_db &select_pfam_db &select_rfam_db &test_pfam_db_selection &select_tigrfam_db &clear_db_selection $newfile $short_file &select_meta_db) ;
use vars qw();

#local modules used
use initialise::read_me qw(:DEFAULT);
use ProgrammeModules::ncrna qw(:DEFAULT);

#local variables
my ($status, $tl, $set);
$status = 0;

#select Glimmer model
sub select_glimmer_model {
   my %args = @_;
   my ($default_model, $newfile, $types, $short_file);

   if (${$args{auto_ini_ref}}{use_glimmer2} == 1) {
      $types = [
               ["Model2 Files", ['.model2']],
               ["Model3 Files", ['.model3']],
               ["All Files"   , ["*"]]
               ];
      $default_model = '.model2';
   } elsif (${$args{auto_ini_ref}}{use_glimmer3} == 1) {
      $types = [
               ["Model3 Files", ['.model3']],
               ["Model2 Files", ['.model2']],
               ["All Files"   , ["*"]]
               ];
      $default_model = '.model3';
   }

   $newfile = ${$args{main_window}}->getOpenFile(-filetypes        => $types,
                                                 -initialdir       => ${$args{ini_ref}}{glimmer_model},
                                                 -initialfile      => ${$args{auto_ini_ref}}{gl_short_file},
                                                 -defaultextension => \$default_model);
   unless (defined $newfile && -e $newfile) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => 'No model file selected, maintaining original selection',
                                        -type    => 'OK',
                                        -icon    => 'error');
      return;
   };

   $newfile =~ m/([^\/]*?)$/;
   $short_file = $1;
   ${$args{auto_ini_ref}}{selected_gl_model} = $newfile;
   ${$args{auto_ini_ref}}{gl_short_file} = $short_file;

   return;

}

#determine selected databses for Blast, Pfam
sub select_blast_db {
   my %args = (critica_select => 0,
               vector_db      => 0,
               @_);
   my ($default_type, $initial_file, $types, $newfile, $corename, $full_corename,
       $short_file);
   #REDO:

   #special case for PSI blast
   #if (${$args{auto_ini_ref}}{blast_type} eq 'PSI-Blast' && $args{vector_db} == 0) {
   #   &select_psi_blast_dbs(main_window  => $args{main_window},
   #                         ini_ref      => $args{ini_ref},
   #                         auto_ini_ref => $args{auto_ini_ref}
   #                        );
   #   return(1);
   #}

   #define default type
   if ((${$args{auto_ini_ref}}{blast_type} eq 'BlastP' || ${$args{auto_ini_ref}}{blast_type} eq 'gappedBlastP'
                                                       || ${$args{auto_ini_ref}}{blast_type} eq 'PSI-Blast')
       && ${$args{auto_ini_ref}}{critica_select} == 0
       && $args{vector_db} == 0) {
      $default_type = '.phr';
      $initial_file = ${$args{auto_ini_ref}}{blast_db};
      $types = [
               ["Amino Acid db", ['.phr', '.pal']],
               ["Nucleotide db", ['.nhr', '.nal']],
               ["All Files"   , ["*"]]
               ];
   } elsif ((${$args{auto_ini_ref}}{blast_type} eq 'BlastN' || ${$args{auto_ini_ref}}{blast_type} eq 'tBlastX')
            && ${$args{auto_ini_ref}}{critica_select} == 0
            && $args{vector_db} == 0) {
      $default_type = '.nhr';
      $initial_file = ${$args{auto_ini_ref}}{blast_db};
      $types = [
               ["Nucleotide db", ['.nhr', '.nal']],
               ["Amino Acid db", ['.phr', '.pal']],
               ["All Files"   , ["*"]]
               ];
   }
   if ($args{critica_select} == 1) {
      $default_type = '.nhr';
      $initial_file = ${$args{auto_ini_ref}}{selected_critica_db};
      $types = [
               ["Nucleotide db", ['.nhr', '.nal']],
               ["Amino Acid db", ['.phr', '.pal']],
               ["All Files"   , ["*"]]
               ];
   }
   if ($args{vector_db} == 1) {
      $default_type = '.nhr';
      $initial_file = ${$args{auto_ini_ref}}{selected_vector_db};
      $types = [
               ["Nucleotide db", ['.nhr', '.nal']],
               ["All Files"   , ["*"]]
               ];
   }

   $newfile = ${$args{main_window}}->getOpenFile(-filetypes        => $types,
                                                 -initialdir       => ${$args{ini_ref}}{blast_db_path},
                                                 -initialfile      => $initial_file,
                                                 -defaultextension => $default_type,
                                                 -multiple         => 1);
   #was a file selected?
   unless (defined $newfile && -e $$newfile[0]) {
      ${$args{main_window}}->messageBox(-title   => 'No selection',
                                        -message => 'No db file selected, maintaining original selection',
                                        -type    => 'OK',
                                        -icon    => 'info');
      ${$args{main_window}}->update;
      return ();
   }
   #multiple files selected?
   if ($#{$newfile} > 0) {
      ${$args{main_window}}->messageBox(-title   => 'Multiple selection',
                                        -message => "More than one db file was selected.\nUsing first selected one, ignoring all others",
                                        -type    => 'OK',
                                        -icon    => 'info');
   }

   $$newfile[0] =~ m/([^\/]*?)$/;
   $short_file = $1;
   unless (defined $short_file) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Error parsing filename for selected database $$newfile[0]\.",
                                        -type    => 'OK',
                                        -icon    => 'info');
   }

   #check for amino acid format db for BlastP, gapped BlastP
   if ((${$args{auto_ini_ref}}{blast_type} eq 'BlastP' || ${$args{auto_ini_ref}}{blast_type} eq 'gappedBlastP'
                                                       || ${$args{auto_ini_ref}}{blast_type} eq 'PSI-Blast')
       && $args{critica_select} == 0
       && $args{vector_db}      == 0) {
      ($corename, $$newfile[0]) = &verify_aa_db(main_window  => $args{main_window},
                                                ini_ref      => $args{ini_ref},
                                                auto_ini_ref => $args{auto_ini_ref},
                                                newfile      => $$newfile[0],
                                                short_file   => $short_file
                                               );
      if ($corename eq '0') {
         &select_blast_db(main_window  => $args{main_window},
                          ini_ref      => $args{ini_ref},
                          auto_ini_ref => $args{auto_ini_ref},
                          vector_db    => 0
                         );
      }
   } elsif ((${$args{auto_ini_ref}}{blast_type} eq 'BlastN' || ${$args{auto_ini_ref}}{blast_type} eq 'tBlastX')
            && $args{critica_select} == 0
            && $args{vector_db}      == 0) {
      #check for nt format db for BlastP, gapped BlastP
      ($corename, $$newfile[0]) = &verify_nt_db(main_window  => $args{main_window},
                                                ini_ref      => $args{ini_ref},
                                                auto_ini_ref => $args{auto_ini_ref},
                                                newfile      => $$newfile[0],
                                                short_file   => $short_file
                                               );
      if ($corename eq '0') {
         &select_blast_db(main_window  => $args{main_window},
                          ini_ref      => $args{ini_ref},
                          auto_ini_ref => $args{auto_ini_ref},
                          vector_db    => 0
                         );
      }
   }
   #check for nt format db for Critica and vector_screen
   if ($args{critica_select} == 1 || $args{vector_db} == 1) {
      ($corename, $$newfile[0]) = &verify_nt_db(main_window  => ${$args{main_window}},
                                                ini_ref      => $args{ini_ref},
                                                auto_ini_ref => $args{auto_ini_ref},
                                                newfile      => $$newfile[0],
                                                short_file   => $short_file
                                               );
      if ($corename eq '0') {#goto REDO};
         &select_blast_db(main_window    => $args{main_window},
                          ini_ref        => $args{ini_ref},
                          auto_ini_ref   => $args{auto_ini_ref},
                          vector_db      => $args{vector_db},
                          critica_select => $args{critica_select},
                          vector_db      => 1
                         );
      }
   }

   #get path to Blast DB
   'reset'      =~ m/reset/;
   $$newfile[0] =~ m/^(.+)\/(.+)$/;
   my $path     = $1;
   #get rid of default extensions
   $$newfile[0] =~ s/\.(nhr|nal|phr|pal)$//;
   $corename    =~ s/\.(nhr|nal|phr|pal)$//;

   if (${$args{auto_ini_ref}}{blast_selector} == 1
       && $args{critica_select}               == 0
       && $args{vector_db}                    == 0) { #modify Blast database
      ${$args{auto_ini_ref}}{blast_db}      = $corename;
      ${$args{auto_ini_ref}}{full_blast_db} = $$newfile[0];
      ${$args{ini_ref}}{blast_db_path}      = $path;
      ${$args{main_window}}->update;
      return (1);
   }
   if ($args{critica_select} == 1) { #modify Critica Blast db
      ${$args{auto_ini_ref}}{selected_critica_db} = $corename;
      ${$args{auto_ini_ref}}{full_critica_db}     = $$newfile[0];
      ${$args{main_window}}->update;
      return (1);
   }
   if ($args{vector_db} == 1) { #modify vectorBlast database
      ${$args{auto_ini_ref}}{selected_vector_db}      = $corename;
      ${$args{auto_ini_ref}}{full_vector_db}          = $$newfile[0];
      ${$args{ini_ref}}{vector_db_path}               = $path;
   }
   ${$args{main_window}}->update;
   return (1);
}

sub select_cog_db {
   my %args = @_;
   my ($newfile, $short_file, $corename);
   my $types = [
               ["COG db", ['.phr', '.pal']],
               ["All Files"   , ["*"]]
               ];

   #Cog databases are expected to be clustered in on root directory
   #${$args{ini_ref}}{COG_db_path} is the path to root directory
   $newfile = ${$args{main_window}}->getOpenFile(-title            => 'Select COG database',
                                                 -filetypes        => $types,
                                                 -initialdir       => ${$args{ini_ref}}{COG_db_path},
                                                 -initialfile      => ${$args{auto_ini_ref}}{COG_db},
                                                 -defaultextension => '.phr',
                                                 -multiple         => 1
                                                 );

   #was a file selected?
   unless (defined $newfile && -e $$newfile[0]) {
      ${$args{main_window}}->messageBox(-title   => 'No selection',
                                        -message => 'No COG db file selected, maintaining original selection',
                                        -type    => 'OK',
                                        -icon    => 'info');
      return ();
   }
   #multiple files selected?
   if ($#{$newfile} > 0) {
      ${$args{main_window}}->messageBox(-title   => 'Multiple selection',
                                        -message => "More than one db file was selected.\nUsing first selected one, ignoring all others",
                                        -type    => 'OK',
                                        -icon    => 'info');
   }

   $$newfile[0] =~ m/([^\/]*?)$/;
   $short_file = $1;

   #check for amino acid format db for BlastP, gapped BlastP
   ($corename, $newfile) = &verify_COG_db(main_window  => $args{main_window},
                                          auto_ini_ref => $args{auto_ini_ref},
                                          ini_ref      => $args{ini_ref},
                                          newfile      => $$newfile[0],
                                          short_file   => $short_file
                                        );
   if ($corename eq '0') {
      ${$args{main_window}}->messageBox(-title   => 'No selection',
                                        -message => 'No valid COG db file selected, maintaining original selection',
                                        -type    => 'OK',
                                        -icon    => 'info');
      return ();
   }

   #$corename =~ s/.*?\/(.+)$/$1/;
   'reset' =~ m/reset/;
   $newfile =~ m/(.*)\/.+?\/.+$/;
   my $path = $1;
   ${$args{auto_ini_ref}}{COG_db}      = $corename;
   ${$args{auto_ini_ref}}{full_COG_db} = $newfile;
   ${$args{ini_ref}}{COG_db_path}      = $path;

   #warn if parsing files are not present in COG dir
   if (${$args{auto_ini_ref}}{standardCOG} == 1) {
      unless (-e ${$args{ini_ref}}{COG_db_path}.'/fun.txt' && -e ${$args{ini_ref}}{COG_db_path}.'/myva=gb' &&
              -e ${$args{ini_ref}}{COG_db_path}.'/org.txt' && -e ${$args{ini_ref}}{COG_db_path}.'/whog') {
         ${$args{main_window}}->messageBox(-title   => 'No COG parsing files',
                                           -message => 'COG parsing files missing in COG directory. '.
                                                       'The following files are required: fun.txt, '.
                                                       'myva=gb, org.txt and whog',
                                           -type    => 'OK',
                                           -icon    => 'warning');
      }
   } elsif (${$args{auto_ini_ref}}{COG2008} == 1) {
      ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
      my $COG2008_dir = $1;
      unless (-e $COG2008_dir.'/cogdef.txt' &&
              -e $COG2008_dir.'/domdat.txt' &&
              -e $COG2008_dir.'/genomes.txt') {
         ${$args{main_window}}->messageBox(-title   => 'No COG2008 parsing files',
                                           -message => 'COG2008 parsing files missing in COG directory '. $COG2008_dir . '.'.
                                                       'The following files are required: cogdef.txt, '.
                                                       'domdat.txt, org.txt and genomes.txt',
                                           -type    => 'OK',
                                           -icon    => 'warning');
      }
   } elsif (${$args{auto_ini_ref}}{arCOG} == 1) {
      ${$args{auto_ini_ref}}{full_COG_db} =~ m/^(.+)\/\w+$/;
      my $arCOG_dir = $1;
      unless (-e $arCOG_dir.'/arCOG2.txt' &&
              -e $arCOG_dir.'/arCOGdef.tab') {
         ${$args{main_window}}->messageBox(-title   => 'No arCOG parsing files',
                                           -message => 'arCOG parsing files missing in COG directory '. $arCOG_dir .'. '.
                                                       'The following files are required: '.
                                                       'arCOG2.txt and arCOGdef.tab',
                                           -type    => 'OK',
                                           -icon    => 'warning');
      }
   }

   return (1);
}

sub select_pfam_db {
   my %args = @_;
   my ($newfile, $types);

   #set default filetypes for selected Pfam version
   if (${$args{auto_ini_ref}}{use_Pfam2} == 1) {
      $types = [
               ["Pfam2 db", ['Pfam*']],
               ["All Files"   , ["*"]]
               ];
   } elsif (${$args{auto_ini_ref}}{use_Pfam3} == 1) {
     $types = [
               ["Pfam3 db", ['*.h3f']],
               ["All Files"   , ["*"]]
               ];
   }

   undef $newfile;
   $newfile = ${$args{main_window}}->getOpenFile(-title            => 'Select Pfam databases',
                                                 -initialdir       => ${$args{ini_ref}}{pfam_db_path}.'/',
                                                 -filetypes        => $types,
                                                 -defaultextension => '',
                                                 -multiple         => 1);


   #was a file selected?
   unless (defined $newfile ) {
      ${$args{main_window}}->messageBox(-title   => 'No selection',
                                        -message => 'No PFam db file selected, maintaining original selection',
                                        -type    => 'OK',
                                        -icon    => 'info');
      return ();
   }

   #strip pathways from selected dbs
   foreach my $file (@{$newfile}) {
      my ($db_path, $corename);
      'reset' =~ /reset/;
      $file =~ m/(.*\/)(.+)$/;
      $db_path   = $1;
      $corename  = $2;

      #clean up Pfam3 database name where necessary
      $corename =~ s/\.(h3f|h3i|h3m|h3p)$//;
      $file     =~ s/\.(h3f|h3i|h3m|h3p)$//;
      $db_path =~ s/[\\\/]$//g;

      ${$args{auto_ini_ref}}{Pfam_db}      .= ' ; '.$corename;
      ${$args{auto_ini_ref}}{full_Pfam_db} .= ' ; '.$file;
      ${$args{ini_ref}}{pfam_db_path}       = $db_path;
   }

   ${$args{auto_ini_ref}}{full_Pfam_db} =~ s/^ ; //;
   ${$args{auto_ini_ref}}{Pfam_db}      =~ s/^ ; //;

   return (1);
}

sub select_rfam_db {
   my %args = @_;
   my ($newfile, $types, $infernal_version);

   #check with version of Infernal is installed
   {
      chdir ${$args{ini_ref}}{infernal_dir};
      system "./cmsearch -h > version.txt";
      open READ, "<${$args{ini_ref}}{infernal_dir}\/version.txt";
      while (<READ>) {
         next unless ($_ =~ m/^\#?\s*infernal\s+(\S+)\s/i);
         $_ =~ m/^\#?\s*infernal\s+(\S+)\s/i;
         $infernal_version = $1;
         last;
      }
      close READ;
      unlink "${$args{ini_ref}}{infernal_dir}\/version.txt";
      chdir ${$args{auto_ini_ref}}{work_dir};
   }

   $types = [
            ["Rfam db", ['Rfam*']],
            ["All Files"   , ["*"]]
            ];

   undef $newfile;
   $newfile = ${$args{main_window}}->getOpenFile(-title            => 'Select Rfam databases',
                                                 -initialdir       => ${$args{ini_ref}}{rfam_db_path}.'/',
                                                 -filetypes        => $types,
                                                 -defaultextension => '',
                                                 -multiple         => 1);


   #was a file selected?
   unless (defined $newfile ) {
      ${$args{main_window}}->messageBox(-title   => 'No selection',
                                        -message => 'No RFam db file selected, maintaining original selection',
                                        -type    => 'OK',
                                        -icon    => 'info');
      return ();
   }

   #clear variables:
   ${$args{auto_ini_ref}}{Rfam_db}      = '';
   ${$args{auto_ini_ref}}{full_Rfam_db} = '';

   #strip pathways from selected dbs
   foreach my $file (@{$newfile}) {
      my ($db_path, $corename);
      'reset' =~ /reset/;
      $file =~ m/(.*\/)(.+)$/;
      $db_path   = $1;
      $corename  = $2;

      #clean up Rfam database name where necessary
      my $db_file = $file;
      $corename =~ s/\.(cm)$//;
      #$db_file  =~ s/\.(nhr|nsq|nin|nsd|nsi)$//;
      $db_path  =~ s/[\\\/]$//g;

      ${$args{auto_ini_ref}}{Rfam_db}      .= ' ; '.$corename;
      ${$args{auto_ini_ref}}{full_Rfam_db} .= ' ; '.$db_file;
      ${$args{ini_ref}}{rfam_db_path}       = $db_path;

      #was compressed or unformatted file selected?
      if ($file =~ m/(^.*\/Rfam\.fasta$|^.*\/Rfam\.fasta\.gz$)/) {
         my ($status) = &recompile_Rfam_fasta(main_window  => $args{main_window},
                                              progress_bar => $args{progress_bar},
                                              auto_ini_ref => $args{auto_ini_ref},
                                              ini_ref      => $args{ini_ref},
                                              Rfam_file    => $file
                                             );
         if ($status == 0) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Could not compile Rfam.fasta nt database in directory ${$args{ini_ref}}{rfam_db_path}",
                                              -type    => 'OK',
                                              -icon    => 'info');
            return (0); #exit to main
         }
         #call database selection again
         &select_rfam_db(main_window  => $args{main_window},
                         auto_ini_ref => $args{auto_ini_ref},
                         ini_ref      => $args{ini_ref}
                        );
      }
   }

   ${$args{auto_ini_ref}}{full_Rfam_db} =~ s/^ ; //;
   ${$args{auto_ini_ref}}{Rfam_db}      =~ s/^ ; //;

   return (1);
}

sub test_pfam_db_selection {
   my %args = @_;
   my @dbs;

   #iterate through each selected database and return if possible problem is found
   @dbs = split (/ ; /, ${$args{auto_ini_ref}}{full_Pfam_db});

   #test if a potential Pfam3 db has been seelcted when running Hmmer2
   foreach my $entry (@dbs) {
      if (${$args{auto_ini_ref}}{use_Pfam2} == 1 &&
          ($entry =~ /\.(h3f|h3i|h3m|h3p)$/ ||
           -e $entry.'.h3f' ||
           -e $entry.'.h3i' ||
           -e $entry.'.h3m' ||
           -e $entry.'.h3p'
          )
         ) {
         return ('Pfam2error');
      }
      if (${$args{auto_ini_ref}}{use_Pfam3} == 1) {
         #clean up filename for testing.
         $entry =~ s/\.(h3f|h3i|h3m|h3p)$//;
         unless (-e $entry.'.h3f' &&
                 -e $entry.'.h3i' &&
                 -e $entry.'.h3m' &&
                 -e $entry.'.h3p'
                ) {
            return ('Pfam3error');
         }
         #make sure the Pfam descriptor file is present, and all are in the same directory.
         #Adjust pfam_descriptor pathway

         if (-e ${$args{ini_ref}}{pfam_db_path}.'/pfamA.txt') {
            ${$args{ini_ref}}{Pfam_descriptor} = ${$args{ini_ref}}{pfam_db_path};
         } else {
            my $newfile = '';
            undef $newfile;
            while (!defined $newfile && !-e $newfile) {
               my $types = [
                            ["PfamA.txt", ['.txt']],
                            ["All Files"   , ["*"]]
                           ];
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Pfam descriptor file is missing. Select file",
                                                 -type    => 'OK',
                                                 -icon    => 'info');
               $newfile = ${$args{main_window}}->getOpenFile(-title            => 'Select PfamA.txt file',
                                                             -initialdir       => ${$args{ini_ref}}{Pfam_descriptor}.'/',
                                                             -defaultextension => '.txt',
                                                             -filetypes        => $types,
                                                             -multiple         => 0);
            }
            'reset'  =~ /reset/;
            $newfile =~ m/(.*\/).+$/;
            ${$args{ini_ref}}{Pfam_descriptor}   = $1;
         }

      }
   }
   return (1);
}

sub select_tigrfam_db {
   my %args = @_;
   my ($newfile);
   my $types = [
               ["TIGRfam db", ['.LIB']],
               ["All Files"   , ["*"]]
               ];
   $newfile = ${$args{main_window}}->getOpenFile(-title            => 'Select TIGRfam databases',
                                                 -initialdir       => ${$args{ini_ref}}{TIGRfam_db_path}.'/',
                                                 -filetypes        => $types,
                                                 -defaultextension => '.LIB',
                                                 -multiple         => 1);

   #was a file selected?
   unless (defined $newfile ) {
      ${$args{main_window}}->messageBox(-title   => 'No selection',
                                        -message => 'No TGIRfam db file selected, maintaining original selection',
                                        -type    => 'OK',
                                        -icon    => 'info');
      return (0);
   }

   #strip pathways from selected dbs
   foreach my $file (@{$newfile}) {
      my ($db_path, $corename);
      'reset' =~ /reset/;
      $file =~ m/(.*\/)(.+)$/;
      $db_path   = $1;
      $corename  = $2;

      $db_path =~ s/[\\\/]$//g;

      ${$args{auto_ini_ref}}{TIGRfam_db}      .= ' ; '.$corename;
      ${$args{auto_ini_ref}}{full_TIGRfam_db} .= ' ; '.$file;
      ${$args{ini_ref}}{TIGRfam_db_path}       = $db_path;
   }

   ${$args{auto_ini_ref}}{full_TIGRfam_db} =~ s/^ ; //;
   ${$args{auto_ini_ref}}{TIGRfam_db}      =~ s/^ ; //;

   return (1);
}

sub select_psi_blast_dbs {
   my %args = @_;
   my ($newfile, $short_file, $db_path, $corename, $psi_aa, $full_psi_aa,
       $psi_nt, $full_psi_nt, $full_corename);

   #select amino acid db first
   PSI_AA_REDO:
   my $types = [
               ["Amino Acid db", ['.phr', '.pal']]
               ];

   $newfile = ${$args{main_window}}->getOpenFile(-title       => 'Select aa database for PSI-Blast',
                                                 -filetypes   => $types,
                                                 -initialdir  => ${$args{ini_ref}}{blast_db_path},
                                                 -initialfile => ${$args{auto_ini_ref}}{blast_db},
                                                 -defaultextension => '.phr');
   #was a file selected?
   unless (defined $newfile && -e $newfile) {
      ${$args{main_window}}->messageBox(-title   => 'No selection',
                                        -message => 'No db file selected, maintaining original selection',
                                        -type    => 'OK',
                                        -icon    => 'info');
      return (0);
   }

   $newfile    =~ m/^(.+)\/([^\/]*?)$/;
   ($db_path, $short_file) = ($1, $2);

   #set pathways; Both aa and nt databases MUST be in the same directory
   ${$args{ini_ref}}{blast_db_path} = $db_path;

   #verify aa selection
   ($corename, $newfile) = &verify_aa_db(main_window  => $args{main_window},
                                         ini_ref      => $args{ini_ref},
                                         auto_ini_ref => $args{auto_ini_ref},
                                         newfile      => $newfile,
                                         short_file   => $short_file
                                        );

   if ($corename eq '0') {goto PSI_AA_REDO};
   $corename    =~ s/.*\/(.+)$/$1/;
   $psi_aa      = $corename;
   $full_psi_aa = $newfile;
   $full_psi_aa =~ s/\.pal$//;
   $corename    = '';

   #select nt db second
   PSI_NT_REDO:
   undef $newfile;
   $types = [
            ["Nucleotide db", ['.nhr', '.nal']]
            ];

   $newfile = ${$args{main_window}}->getOpenFile(-title            => 'Select nt database for PSI-Blast',
                                                 -filetypes        => $types,
                                                 -initialdir       => ${$args{ini_ref}}{blast_db_path},
                                                 -initialfile      => '',
                                                 -defaultextension => '.nhr');
   #was a file selected?
   unless (defined $newfile && -e $newfile) {
      ${$args{main_window}}->messageBox(-title   => 'No selection',
                                        -message => 'No db file selected, maintaining original selection',
                                        -type    => 'OK',
                                        -icon    => 'info');
      return (0);
   }

   $newfile    =~ m/([^\/]*?)$/;
   $short_file = $1;

   #verify nt selection
   ($corename, $newfile) = &verify_nt_db(main_window  => $args{main_window},
                                         ini_ref      => $args{ini_ref},
                                         auto_ini_ref => $args{auto_ini_ref},
                                         newfile      => $newfile,
                                         short_file   => $short_file
                                        );
   if ($corename eq '0') {goto PSI_NT_REDO};
   $corename      =~ s/.*\/(.+)$/$1/;
   $psi_nt        = $corename;
   $full_psi_nt   = $newfile;
   $full_psi_nt   =~ s/\.nal$//;
   $corename      = '';
   $corename      = $psi_aa.' ; '.$psi_nt;
   $full_corename = $full_psi_aa.' ; '.$full_psi_nt;

   ${$args{auto_ini_ref}}{blast_db} = $corename;
   ${$args{auto_ini_ref}}{full_blast_db} = $full_corename;

   return (1);

}

sub verify_aa_db {
   my %args = @_;
   my ($path, $corename, $extension, $multi_db, $dbparts);

   #check if correct extension
   if ($args{short_file} !~ /\.phr$/ && $args{short_file} !~ /\.pin$/ && $args{short_file} !~ /\.pnd$/ &&
       $args{short_file} !~ /\.pni$/ && $args{short_file} !~ /\.psd$/ && $args{short_file} !~ /\.psi$/ &&
       $args{short_file} !~ /\.psq$/ && $args{short_file} !~ /\.pal$/) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Selected database file $args{short_file} does not appear to represent an amino-acid database.\nSelect a different database",
                                        -type    => 'OK',
                                        -icon    => 'error');
      return (0);
   }

   #check if multipart db is present for selected file
   $corename      = '';
   $args{newfile} =~ m/^(.+)\/(.+)(\..+)$/;
   ($path, $corename, $extension) = ($1, $2, $3);
   $corename      =~ s/\.\d+$//;

   if (-e $path.'/'.$corename.'.pal' && $extension !~ /\.pal/) {
      ${$args{main_window}}->messageBox(-title   => 'Info',
                                        -message => "Selected db file $args{short_file} appears to be part of a multi-part database.\nAutomatically switching to correct file $corename".'.pal',
                                        -type    => 'OK',
                                        -icon    => 'info');
      $args{short_file} = $corename.'.pal';
      $args{newfile} =~ s/\.\d+$extension$/\.pal/;
   }

   #check if all files exist
   unless (-e $path.'/'.$corename.'.phr' && -e $path.'/'.$corename.'.pin' && -e $path.'/'.$corename.'.pnd' &&
           -e $path.'/'.$corename.'.pni' && -e $path.'/'.$corename.'.psd' && -e $path.'/'.$corename.'.psi' &&
           -e $path.'/'.$corename.'.psq' ) {
      #take into account a different database format
      unless (-e $path.'/'.$corename.'.psq' && -e $path.'/'.$corename.'.psi' && -e $path.'/'.$corename.'.psd' &&
              -e $path.'/'.$corename.'.pog' && -e $path.'/'.$corename.'.pin' && -e $path.'/'.$corename.'.phr') {
         if ($args{short_file} !~ /\.pal$/) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => 'Selected amino acid database appears to be incomplete. 7 files should be present',
                                              -type    => 'OK',
                                              -icon    => 'error');
            return (0);
         }
      }
   }

   #check if database is broken up into several parts
   if ($args{newfile} =~ /\.pal$/) {
      local( $/, *ENTRY ) ;
      open( ENTRY, $args{newfile} );
      $multi_db = <ENTRY>;
      close ENTRY;

      $multi_db =~ m/\nDBLIST\s(.*?)\n/s;
      $dbparts  = $1;
      $dbparts  =~ s/\"//gs;

      while ($dbparts =~ m/([^\s]+)/g) {
         $set = $path.'/'.$&;

         unless (-e $set.'.phr' && -e $set.'.pin' && -e $set.'.pnd' &&
                 -e $set.'.pni' && -e $set.'.psd' && -e $set.'.psi' &&
                 -e $set.'.psq' ) {
            #take into account a different database format
            unless (-e $set.'.psq' && -e $set.'.psi' && -e $set.'.psd' &&
                    -e $set.'.pog' && -e $set.'.pin' && -e $set.'.phr') {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Selected multi-part db $args{short_file} appears to be broken.\nReselect a different database",
                                                 -type    => 'OK',
                                                 -icon    => 'error');
               return (0);
            }
         }
      }
   }

   return ($corename, $args{newfile});
}

sub verify_nt_db {
   my %args = @_;
   my ($path, $corename, $extension, $multi_db, $dbparts);

   #check if correct extension
   if ($args{short_file} !~ /\.nhr$/ && $args{short_file} !~ /\.nin$/ && $args{short_file} !~ /\.nnd$/ &&
       $args{short_file} !~ /\.nni$/ && $args{short_file} !~ /\.nsd$/ && $args{short_file} !~ /\.nsi$/ &&
       $args{short_file} !~ /\.nsq$/ && $args{short_file} !~ /\.nal$/) {
      ${$args{main_window}}->messageBox(-title   => 'Error',
                                        -message => "Selected database file $args{short_file} does not appear to represent an nucleotide database.\nSelect a different database",
                                        -type    => 'OK',
                                        -icon    => 'error');
      return (0);
   }

   #check if multipart db is present for selected file
   $corename      = '';
   $args{newfile} =~ m/^(.+)\/(.+)(\..+)$/;
   ($path, $corename, $extension) = ($1, $2, $3);
   $corename      =~ s/\.\d+$//;

   if (-e $path.'/'.$corename.'.nal' && $extension !~ /\.nal/) {
      ${$args{main_window}}->messageBox(-title   => 'Info',
                                        -message => "Selected db file $args{newfile} appears to be part of a multi-part database.\nAutomatically switching to correct file $corename".'.nal',
                                        -type    => 'OK',
                                        -icon    => 'info');
      $args{short_file} = $corename.'.nal';
      $args{newfile}    =~ s/\.\d+$extension$/\.nal/;
   }
   #check if all files exist
   unless (-e $path.'/'.$corename.'.nhr' && -e $path.'/'.$corename.'.nin' && (-e $path.'/'.$corename.'.nnd' || -e $path.'/'.$corename.'.nhd') &&
           -e $path.'/'.$corename.'.nsd' && -e $path.'/'.$corename.'.nsi' && (-e $path.'/'.$corename.'.nni' || -e $path.'/'.$corename.'.nhi') &&
           (-e $path.'/'.$corename.'.nsq'|| -e $path.'/'.$corename.'.nog' )) {
      #take into account a different database format
      unless (-e $path.'/'.$corename.'.nsq' && -e $path.'/'.$corename.'.nsi' && -e $path.'/'.$corename.'.nsd' &&
              -e $path.'/'.$corename.'.nin' && -e $path.'/'.$corename.'.nhr') {
         if ($args{short_file} !~ /\.nal$/) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Selected nucleotide database $args{short_file} appears to be incomplete.\n7 files should be present",
                                              -type    => 'OK',
                                              -icon    => 'error');
            return (0);
         }
      }
   }

   #check if database is broken up into several parts
   if ($args{newfile} =~ /\.nal$/) {
      local( $/, *ENTRY ) ;
      open( ENTRY, $args{newfile} );
      $multi_db = <ENTRY>;
      close ENTRY;

      $multi_db =~ m/\nDBLIST\s(.*?)\n/s;
      $dbparts  = $1;
      $dbparts  =~ s/\"//gs;

      while ($dbparts =~ m/([^\s]+)/g) {
         $set = $path.'/'.$&;

         unless (-e $set.'.nhr' && -e $set.'.nin' && -e $set.'.nnd' &&
                 -e $set.'.nni' && -e $set.'.nsd' && -e $set.'.nsi' &&
                 -e $set.'.nsq' ) {
            #take into account a different database format
            unless (-e $set.'.nsq' && -e $set.'.nsi' && -e $set.'.nsd' &&
                    -e $set.'.nin' && -e $set.'.nhr') {
               ${$args{main_window}}->messageBox(-title   => 'Error',
                                                 -message => "Selected multi-part db $args{short_file} appears to be broken.\nReselect a different database",
                                                 -type    => 'OK',
                                                 -icon    => 'error');
               return (0);
            }
         }
      }
   }
   return ($corename, $args{newfile});
}

sub verify_COG_db {
   my %args = @_;
   my ($filepath, $corename, $extension, $multi_db, $dbparts);

   #check if correct extension
   if ($args{short_file} !~ /\.(phr|pin|psd|psi|psq|pal)$/)  {
       ${$args{main_window}}->messageBox(-title   => 'Error',
                                         -message => "Selected COG database file $args{short_file} does not appear to represent a valid COG database. Select a different database",
                                         -type    => 'OK',
                                         -icon    => 'error');
      return (0);
   }

   #check if multipart db is present for selected file
   $args{newfile} =~ m/(.+)\/(.+)(\..+)$/;

   $filepath  = $1;
   $corename  = $2;
   $extension = $3;
   $corename  =~ s/\.\d+$//;

   if (-e $filepath.'/'.$corename.'.pal' && $extension !~ /\.pal/) {
      ${$args{main_window}}->messageBox(-title   => 'Info',
                                        -message => "Selected COG db file $args{short_file} appears to be part of a multi-part database. Automatically switching to correct file $corename\.pal",
                                        -type    => 'OK',
                                        -icon    => 'info');
      $args{short_file} = $corename.'.pal';
      $args{newfile}    =~ s/\d+\.$extension$/\.pal/;
   }

   #check if all files exist
   unless (-e $filepath.'/'.$corename.'.phr' && -e $filepath.'/'.$corename.'.pin' &&
           -e $filepath.'/'.$corename.'.psd' && -e $filepath.'/'.$corename.'.psi' &&
           -e $filepath.'/'.$corename.'.psq' ) {
      if ($args{short_file} !~ /\.pal$/) {
         ${$args{main_window}}->messageBox(-title   => 'Error',
                                           -message => "Selected COG database $corename appears to be incomplete. 5 files should be present",
                                           -type    => 'OK',
                                           -icon    => 'error');
         return (0);
      }
   }

   #check if database is broken up into several parts
   if ($args{newfile} =~ /\.pal$/) {
      local( $/, *ENTRY ) ;
      open( ENTRY, $args{newfile} );
      $multi_db = <ENTRY>;
      close ENTRY;

      $multi_db =~ m/\nDBLIST\s(.*?)\n/s;
      $dbparts = $1;

      while ($dbparts =~ m/([^\s]+)/g) {
         $set = ${$args{ini_ref}}{COG_dir}.'/'.$&;

         unless (-e $filepath.'/'.$set.'.phr' && -e $filepath.'/'.$set.'.pin' &&
                 -e $filepath.'/'.$set.'.psd' && -e $filepath.'/'.$set.'.psi' &&
                 -e $filepath.'/'.$set.'.psq' ) {
            ${$args{main_window}}->messageBox(-title   => 'Error',
                                              -message => "Selected multi-part COG db $args{short_file} appears to be broken. Reselect a different database",
                                              -type    => 'OK',
                                              -icon    => 'error');
            return (0);
         }
      }
   }

   return ($corename, $filepath.'/'.$corename );
}

sub clear_db_selection  {
   my %args = @_;

   ${$args{auto_ini_ref}}{$args{short_db}} = '';
   ${$args{auto_ini_ref}}{$args{full_db}}  = '';
   ${$args{main_window}}->update;
   return;
}

sub select_meta_db {
   my %args = @_;
   #my ($critica_sel, $nw, $ini_ref, $auto_ini_ref) = @_;
   my ($default_type, $initial_file, $types, $newfile, $corename, $full_corename,
       $short_file);
   REDO:

   #define default type
   $default_type = '.nhr';
   $initial_file = ${$args{auto_ini_ref}}{blast_db};
   $types = [
            ["Nucleotide db", ['.nhr', '.nal']],
            ["Amino Acid db", ['.phr', '.pal']],
            ["All Files"   , ["*"]]
            ];

   $newfile = ${$args{main_window}}->getOpenFile(-filetypes        => $types,
                                                 -initialdir       => ${$args{ini_ref}}{blast_db_path},
                                                 -initialfile      => $initial_file,
                                                 -defaultextension => $default_type,
                                                 -multiple         => 1);
   #was a file selected?
   unless (defined $newfile && -e $$newfile[0]) {
      ${$args{main_window}}->messageBox(-title   => 'No selection',
                                        -message => 'No db file selected, maintaining original selection',
                                        -type    => 'OK',
                                        -icon    => 'info');
      return ();
   }
   #multiple files selected?
   if ($#{$newfile} > 0) {
      ${$args{main_window}}->messageBox(-title   => 'Multiple selection',
                                        -message => "More than one db file was selected.\nUsing first selected one, ignoring all others",
                                        -type    => 'OK',
                                        -icon    => 'info');
   }

   #strip corename from path
   'reset' =~ m/reset/;
   $$newfile[0] =~ m/([^\/]*?)$/;
   $corename = $1;

   #get rid of default extensions
   $$newfile[0] =~ s/\.(nhr|nal|phr|pal)$//;
   $corename =~ s/\.(nhr|nal|phr|pal)$//;

   return ($$newfile[0], $corename);

}

1;