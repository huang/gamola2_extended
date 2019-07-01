#!/opt/ActivePerl-5.8/bin/perl

#tRNScan result parser
#input arguments: main_window, progress_bar, gb_file, gene_model



package CompilerModules::trnascan_parser;
use strict;
use initialise::read_me        qw(:DEFAULT);
use ProgrammeModules::sequence qw(:DEFAULT);
use Basics::progress_bar       qw(:DEFAULT);

use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&trnascan_parser);
use vars qw();


#local variables
my (%args, %seen, $value);

format ENTRY =
~~                   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$value
.

sub trnascan_parser {
   my %args = @_;
   my ($max_number, $counter);

   #test if result file is present
   unless (-e ${$args{ini_ref}}{trnascan_results}.'/'.$args{filename}.'_tRNAscan') {
      my $error_msg = ${$args{main_window}}->Dialog(-title   => 'No file found',
                                                    -text    => "tRNAscan result file is not present in folder ${$args{ini_ref}}{trnascan_results}",
                                                    -buttons => ['OK'],
                                                    -bitmap  => 'info');
      $error_msg->after(${$args{auto_ini_ref}}{msg_timeout}, sub {$error_msg->Exit()});
      $error_msg-> Show();
      return;
   }

   #return if empty file
   return (1) if (-s ${$args{ini_ref}}{trnascan_results}.'/'.$args{filename}.'_tRNAscan' < 1);

   #create status box
   &progress_bar_2(main_window  => $args{main_window},
                   progress_bar => $args{progress_bar},
                   auto_ini_ref => $args{auto_ini_ref},
                   ini_ref      => $args{ini_ref},
                   title        => "Compiling tRNAscan results",
                   label        => ''
                  );
   &show_pbar_2;

   #determine length of Summary file
   $max_number = `wc -l ${$args{ini_ref}}{trnascan_results}\/$args{filename}\_tRNAscan`;
   if ($max_number < 1) {$max_number = 1};

   #open summary file
   $counter = 0;
   open READ, ${$args{ini_ref}}{trnascan_results}.'/'.$args{filename}.'_tRNAscan';
   while (<READ>) {
      my (@temp, $boundary_trna, $gene, $product,
          $feature_trna, $key_trna);
      if ($. < 4) {
         next;
      }

      #update progress bar
      $counter++;
      if (($counter % 20) == 0) {
         &update_pbar_2(title        => "Compiling tRNAscan results",
                        label        => "Reading entry from $args{input_file}",
                        progress     => ($counter / $max_number) * 100,
                       );
      }

      #increase ID
      $args{counter}++;

      #parse entry
      @temp = split /\s+/;
      if ($temp[4] =~ /Pseudo/) {
         next;
      }

      #create boundary entry
      if ($temp[2] > $temp[3]) {
         $boundary_trna = '     tRNA            complement('.$temp[3].'..'.$temp[2]."\)\n";
      } else {
         $boundary_trna = '     tRNA            '.$temp[2].'..'.$temp[3]."\n";
      }

      #create gene structure
      $value = '/gene="'.$temp[4].'"'."\n";
      open  ENTRY, '>', \$gene;
      write ENTRY;
      close ENTRY;

      #create product
      $value = '/product="tRNA: '.$temp[4].'; Anticodon: '.$temp[5].'"'."\n";
      open  ENTRY, '>', \$product;
      write ENTRY;
      close ENTRY;

      #combine feature
      $feature_trna =  $boundary_trna.
                       $gene.
                       $product;
      $key_trna     =  $temp[2].'_'.$temp[3].'_tRNA_'.$args{counter}; #create  uniqe ID instead of ORFnumber

      ${$args{genbank_ref}}{$key_trna} = $feature_trna;
      push (@{$args{feature_list_ref}}, $key_trna);
   }
   close READ;
   &hide_pbar_2;
   return (1);

}



1;

