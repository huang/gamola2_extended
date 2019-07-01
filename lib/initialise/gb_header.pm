#!/opt/ActivePerl-5.8/bin/perl5.8.8

#setup for genbank header


package initialise::gb_header;
use strict;
use vars qw($VERSION @ISA @EXPORT);
use Exporter;
$VERSION = '0.01';
@ISA = ('Exporter');
#exported items
@EXPORT = qw(&gb_header &default_header &set_date) ;
use vars qw();

#local modules used
use Tk;
use Tk::ROText;
use initialise::read_me qw(:DEFAULT);

#local variables
my (%args, $value, $gbkey);

format ENTRY =
@<<<<<<<<<  ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$gbkey,     $value
~~          ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            $value
.


#build genbank headers
sub gb_header {
   my %args = (filename => 'NONE',
               seq_length => '0',
               definition => 'NONE',
               @_);
   my ($gbheader, $top_frame, $bottom_frame, $text_frame, $entries, $entries_fields, $entries_label,
       $defined_fields, $selected_field, $lb_frame,
       $max_width, $counter, $rm_field, $line_add, $field, $subfield, $fieldhash
      );
   my (@browse, @header);

   #create and show toplevel is necessary
   if (defined $gbheader && Tk::Exists($gbheader)) {
       $gbheader->state('normal');
   } else {
      $gbheader = ${$args{main_window}}->Toplevel(-title          => "Genbank header",
                                                 );
   }

   #predefine set of default fields
   foreach my $field (@{${$args{auto_ini_ref}}{'gb_order'}}) {
      my ($type, $value);
      'reset' =~ m/reset/;
      $field =~ m/([FS])__(.+)/;
      $type = $1;
      $value = $2;

      $fieldhash->{$value} = $type;
      push (@browse, $value);
      push (@header, $value);
   }

   #create top frame
   $selected_field = '';
   $max_width = 1;
   $counter = 1;
   $top_frame = $gbheader->Frame(-borderwidth => 2, -relief => 'groove')->pack(-side => 'top', -fill => 'x');

   #create Bottom frame with exit button
   $bottom_frame = $gbheader->Frame(-borderwidth => 2, -relief => 'groove')->pack(-side => 'bottom', -fill => 'x');

   #create right handed text frame
   $text_frame = $gbheader->Frame(-borderwidth => 2, -relief => 'groove')->pack(-side => 'right', -fill => 'x');
   $lb_frame   = $gbheader->Frame(-borderwidth => 2, -relief => 'groove')->pack(-side => 'left', -fill => 'x');
   #create Pane widget for header fields
   $entries    = $lb_frame->Scrolled("ROText",
                                    -scrollbars => 'osw',
                                    -width => 100,
                                   )->pack (-expand => 1, -fill => 'y', -side=>'bottom', -anchor => 'nw');
   $entries->menu(undef);

   #parse exsting header if defined
   if (defined $args{predefined_header}) {
      &parse_header(main_window       => $args{main_window},
                    progress_bar      => $args{progress_bar},
                    auto_ini_ref      => $args{auto_ini_ref},
                    ini_ref           => $args{ini_ref},
                    predefined_header => $args{predefined_header},
                    defined_fields    => \$defined_fields,
                    fieldhash         => \$fieldhash,
                    browse            => \@browse,
                    header            => \@header,
                    filename          => $args{filename}
                    );
   }
   #create a set of standard fields
   else {
      foreach my $entry (@header) {
         if ($max_width < length($entry)) {$max_width = length($entry)};
         $defined_fields->{$entry} = ${$args{auto_ini_ref}}{$fieldhash->{$entry}.'__'.$entry};
      }
         #modify variable default paramters
         $defined_fields->{'date'}        = &set_date;
         $defined_fields->{'designation'} = $args{filename};
         $defined_fields->{'accession'}   = $args{filename};
         $defined_fields->{'definition'}  = $args{definition};
         $defined_fields->{'length'}      = $args{seq_length};


      $max_width += 10;
   }

   #populate bottom frame
   {
      $bottom_frame->Label(-text => "Add field: Left double click\nRemove field: Middle double click\nSelect line: Left click")
                                                   ->grid (-row => 0, -column => 0, -sticky => 'e', -columnspan => 3);
      $bottom_frame->Button(-text => "OK",
                            -command => sub {
                                             my $final_header =  &compile_header(defined_fields    => \$defined_fields,
                                                                                 fieldhash         => \$fieldhash,
                                                                                 header            => \@header,
                                                                                 main_window       => \$gbheader,
                                                                                 filename          => $args{filename}
                                                                                );
                                             $gbheader->state('withdrawn');
                                             return ($final_header);
                                            })     ->grid (-row => 0, -column => 3);
      $bottom_frame->Button(-text => "Cancel",
                            -command => sub {
                                             $gbheader->state('withdrawn');
                                             return ();
                                            })     ->grid (-row => 0, -column => 4);

      $bottom_frame->Button(-text    => 'Set as default header',
                            -command => sub {
                                              #reset default array
                                              @{${$args{auto_ini_ref}}{'gb_order'}} = ();
                                              #iterate over current array
                                              foreach my $entry (@header) {
                                                 ${$args{auto_ini_ref}}{$fieldhash->{$entry}.'__'.$entry} = $defined_fields->{$entry};
                                                 push (@{${$args{auto_ini_ref}}{'gb_order'}}, $fieldhash->{$entry}.'__'.$entry);
                                              }
                                              $gbheader->messageBox(-title   => 'Default Genbank header',
                                                                    -message => 'Set current setup as default header',
                                                                    -type    => 'OK'
                                                                    );

                                              ${$args{main_window}}->update;
                                            }
                           )
                                                   ->grid(-row => 0, -column => 5);
   }

   #add text widget
   $text_frame->Label(-text => 'Enter field value')->pack();
   my $text = $text_frame->Text(-state => 'normal',
                                -width => 40,
                                )->pack(-fill => 'y');

   #add tools to top frame
   $top_frame->Button(-text    => 'Add field',
                      -command => sub {
                                       &add_entry(entries           => \$entries,
                                                  text              => \$text,
                                                  ini_ref           => $args{ini_ref},
                                                  auto_ini_ref      => $args{auto_ini_ref},
                                                  entries_fields    => \$entries_fields,
                                                  entries_label     => \$entries_label,
                                                  defined_fields    => \$defined_fields,
                                                  header            => \@header,
                                                  browse            => \@browse,
                                                  field             => $field,
                                                  fieldhash         => $fieldhash,
                                                  line_add          => $line_add,
                                                  selected_field    => $selected_field,
                                                  counter           => $counter,
                                                  max_width         => $max_width
                                                 );

                                       ${$args{main_window}}->update;
                                      }
                     )
      ->grid(-row => 0, -column => 0, -sticky => 'e');



   my $main_field = $top_frame-> BrowseEntry(-choices   => \@browse,
                                             -width     => 20,
                                             -variable  => \$selected_field,
                                             -command   => sub {
                                                                if (defined $fieldhash->{$selected_field}) {
                                                                   if ($fieldhash->{$selected_field} eq 'F') {
                                                                      $field    = 1;
                                                                      $subfield = 0;
                                                                   } elsif ($fieldhash->{$selected_field} eq 'S') {
                                                                      $field    = 0;
                                                                      $subfield = 1;
                                                                   }
                                                                }
                                                               }
                                            )
      ->grid (-row => 0, -column => 1, -sticky => 'e');

   $top_frame->Checkbutton(-text     => 'Field',
                           -variable => \$field,
                           -command  => sub {
                                              if ($field == 1) {
                                                 $subfield = 0;
                                              } elsif ($field == 0) {
                                                 $subfield = 1;
                                              }
                                            }
                           )
         ->grid (-row => 0, -column => 2, -sticky => 'e');

   $top_frame->Checkbutton(-text     => 'Subfield',
                           -variable => \$subfield,
                           -command  => sub {
                                              if ($subfield == 1) {
                                                 $field = 0;
                                              } elsif ($subfield == 0) {
                                                 $field = 1;
                                              }
                                            }
                           )
      ->grid (-row => 0, -column => 3, -sticky => 'e');





   #add fields to Pane:
   foreach my $key (@header) {
      #indent for subfields
      if ($fieldhash->{$key} eq 'S') {
         $entries->insert('end', "   ");
      }

      #define height
      my $height = ($defined_fields->{$key} =~ tr/\\z//s);

      #create text widget
      $entries_fields->{$key} = $entries->Text (-width  => 40,
                                                -height => ($height + 1)
                                               );
      #insert text
      $entries_fields->{$key}->insert('end', $defined_fields->{$key});

      $entries->windowCreate('end', -window => $entries_fields->{$key});
      $entries->insert('end', $key);
      $entries->insert('end', "\n");

   }

   #add button bindings
   $gbheader->bind('<Double-ButtonPress-2>' => sub {
                                                       &remove_entry(entries        => \$entries,
                                                                     ini_ref        => $args{ini_ref},
                                                                     auto_ini_ref   => $args{auto_ini_ref},
                                                                     entries_fields => \$entries_fields,
                                                                     defined_fields => \$defined_fields,
                                                                     header         => \@header,
                                                                     );
                                                       ${$args{main_window}}->update;
                                                    });
   $gbheader->bind('<Double-ButtonPress-1>' => sub {
                                                       &add_entry(entries           => \$entries,
                                                                  text              => \$text,
                                                                  ini_ref           => $args{ini_ref},
                                                                  auto_ini_ref      => $args{auto_ini_ref},
                                                                  entries_fields    => \$entries_fields,
                                                                  defined_fields    => \$defined_fields,
                                                                  header            => \@header,
                                                                  browse            => \@browse,
                                                                  field             => $field,
                                                                  fieldhash         => $fieldhash,
                                                                  line_add          => $line_add,
                                                                  selected_field    => $selected_field,
                                                                  counter           => $counter,
                                                                  max_width         => $max_width
                                                                  );
                                                       ${$args{main_window}}->update;
                                                    });
   $gbheader->bind('<ButtonRelease-1>' => sub {
                                                       $line_add = $entries->index('current');
                                                       $line_add =~ s/\.\d+//;
                                                       ${$args{main_window}}->update;
                                                    });
}


sub remove_entry {
   my %args = @_;
   my $remove_entry;

   #determine current index
   my ($current_index) = ${$args{entries}}->index('current');
   $current_index =~ s/\.\d+//;

   #remove from header array
   $remove_entry = splice (@{$args{header}}, ($current_index + 1), 1);

   #remove from field lists
   delete ${$args{entries_fields}}->{$remove_entry};
   delete ${$args{defined_fields}}->{$remove_entry};

   #remove from display
   ${$args{entries}}->delete($current_index.'.0 linestart', $current_index.'.0 lineend' );
   ${$args{entries}}->delete($current_index.'.0');
   ${$args{entries}}->update;


}

sub add_entry {
   my %args = @_;
   my ($height, $field_text);

   #return if header field is empty
   if (!defined $args{selected_field}) {return};

   #define variables if necessary
   unless (defined $args{selected_field}) {$args{selected_field}    = ''};

   #adjust for multiple instances
   while (exists ${$args{defined_fields}}->{$args{selected_field}}) {
      $args{counter}++;
      $args{selected_field} =~ s/_\d+$//;
      $args{selected_field} = $args{selected_field}.'_'.$args{counter};
   }

   #define field or subfield
   if ($args{field} == 1) {
      $args{fieldhash}->{$args{selected_field}} = 'F';
   } elsif ($args{field} == 0) {
      $args{fieldhash}->{$args{selected_field}} = 'S';
   } else {
      return;
   }
   #add to header array
   if ($args{selected_field} =~ /\S+/) {
      #splice (@{$args{header}}, ($args{line_add} - 1), 0, $args{selected_field});
      splice (@{$args{header}}, $args{line_add}, 0, $args{selected_field});
      #add to browse list
      push (@{$args{browse}}, $args{selected_field});
   }

   #retrieve possible text from text widget
   $field_text = ${$args{text}}->get("1.0",'end');

   #add entry to list
   if ($args{selected_field} =~ /\w+/) {
      if ($field_text =~ /\S+/) {
         ${$args{defined_fields}}->{$args{selected_field}} = $field_text;
      } else {
         ${$args{defined_fields}}->{$args{selected_field}} = '';
      }
      $height = ($field_text =~ tr/\n//s);

      #create line for field or subfield
      if ($args{field} == 1) {
         ${$args{entries}}->insert(($args{line_add} + 1).'.0 linestart', "\n");
      } elsif ($args{field} == 0) {
         ${$args{entries}}->insert(($args{line_add} + 1).'.0 linestart', "   \n");
      }
   }

   #create text widget
   ${$args{entries_fields}}->{$args{selected_field}} = ${$args{entries}}->ROText (-width  => 40,
                                                                                  -height => ($height + 1),
                                                                                  );
   #insert text
   ${$args{entries_fields}}->{$args{selected_field}}->insert('end', ${$args{defined_fields}}->{$args{selected_field}});

   ${$args{entries}}->windowCreate(($args{line_add} + 1).'.0 lineend', -window => ${$args{entries_fields}}->{$args{selected_field}});
   ${$args{entries}}->insert(($args{line_add} + 1).'.0 lineend', $args{selected_field});


   #reset reference fields
   if ($args{selected_field} =~ /reference/i) {
      $args{selected_field} = 'reference';
   }

   return;
}


sub parse_header {
   my %args = @_;
   my (@local_header, @tmp);

   $args{predefined_header} = "\n".$args{predefined_header};
   @tmp = split/\n\r?/, $args{predefined_header};
   my $section = shift(@tmp);

   #define single header entries
   foreach my $entry (@tmp) {
      if ($entry =~ /^\s*[A-Z\_\/]+\s+/) {
         push (@local_header, $section);
         $section = $entry;
      } else {
         $section .= $entry;
      }
   }
   push (@local_header, $section);

   #push header entries into fields
   foreach my $entry (@local_header) {
      my ($qualifier, $value);
      $entry =~ s/[\n\r]//gs;
      #$entry =~ s/\s+/ /gs;

      if ($entry =~ /^locus/i) {
         'reset' =~ m/reset/;
         $entry =~ m/locus\s+(\S+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
         my $acc    = $1;
         my $length = $2;
         my $type   = $3;
         my $format = $4;
         my $domain = $5;
         my $date   = $6;
         unless (defined $date) {
            $acc    = '';
            $length = '';
            $type   = '';
            $format = '';
            $domain = '';
            $date   = '';
         }

         ${$args{defined_fields}}->{'designation'} = $acc;
         ${$args{defined_fields}}->{'length'}      = $length;
         ${$args{defined_fields}}->{'type'}        = $type;
         ${$args{defined_fields}}->{'format'}      = $format;
         ${$args{defined_fields}}->{'domain'}      = $domain;
         ${$args{defined_fields}}->{'date'}        = $date;
         ${$args{fieldhash}}     ->{'designation'} = 'S';
         ${$args{fieldhash}}     ->{'length'} = 'S';
         ${$args{fieldhash}}     ->{'type'} = 'S';
         ${$args{fieldhash}}     ->{'format'} = 'S';
         ${$args{fieldhash}}     ->{'domain'} = 'S';
         ${$args{fieldhash}}     ->{'date'} = 'S';
         ${$args{fieldhash}}     ->{'locus'} = 'F';

         my @temp = qw[locus designation length type format domain date];
         @{$args{browse}} = (@{$args{browse}}, @temp);
         @{$args{header}} = (@{$args{header}}, @temp);

         next;
      }
      'reset' =~ m/reset/;
      $entry =~ m/\s*([A-Z\_\.]+)\s+(.+)/;
      $qualifier = $1;
      $value = $2;

      unless (defined $value) {next};
      ${$args{defined_fields}}->{(lc $qualifier)} = $value;
      if ($entry =~ /^\s+/) {
         ${$args{fieldhash}}->{(lc $qualifier)} = 'S';
      } else {
         ${$args{fieldhash}}->{(lc $qualifier)} = 'F';
      }
      push (@{$args{browse}}, (lc $qualifier));
      push (@{$args{header}}, (lc $qualifier));
   }
   return (1);
}

sub set_date {
   my (@months, $day, $month, $year, $date);
   # determine the actual date
   ($day, $month, $year) = (localtime) [3,4,5];
   @months = ('JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC');
   $date =  $day."-".$months[$month]."-".($year+1900);

   return ($date)
}

sub compile_header {
   my %args = @_;
   my ($full_header, %seen, $formatted, $new_locus);
   my ($designation, $length, $type, $format, $domain, $date);

   foreach my $entry (@{$args{header}}) {
      $gbkey = $entry;
      $value = ${$args{defined_fields}}->{$entry};
      if (${$args{fieldhash}}->{$entry} eq 'S') {
         $gbkey = '  '.$gbkey;
      }
      $gbkey = uc($gbkey);
      $gbkey =~ s/_\d+$//;

      open  ENTRY, '>',\$formatted;
      write ENTRY;
      close ENTRY;

      $full_header .= $formatted;
   }

   #re-constitute locus line
   'reset' =~ m/reset/;
   $full_header =~ m/(LOCUS.*?\n\r?)[A-Z\_\.]+\s+/s;
   my $locus = $1;
   unless (defined $locus) {
      ${$args{main_window}}->messageBox(-title   => 'Error parsing Genbank header',
                                        -message => 'Cannot find Locus features',
                                        -type    => 'OK'
                                       );
      return ($full_header);
   }
   $locus =~ m/designat\w*\s+([^\n\r]*?)\n\r?/is;
   $designation = $1;
   unless (defined $designation) {$designation = ''};
   unless ($designation =~ /\w+/) {$designation = ${$args{filename}}};
   $locus =~ m/length\s+([^\n\r]*?)\n\r?/is;
   $length = $1;
   unless (defined $length) {$length = ''};
   if ($length =~ /\d+/) {$length .= ' bp'};
   $locus =~ m/type\s+([^\n\r]*?)\n\r?/is;
   $type = $1;
   unless (defined $type) {$type = ''};
   $locus =~ m/format\s+([^\n\r]*?)\n\r?/is;
   $format = $1;
   unless (defined $format) {$format = ''};
   $locus =~ m/domain\s+([^\n\r]*?)\n\r?/is;
   $domain = $1;
   unless (defined $domain) {$domain = ''};
   $locus =~ m/date\s+([^\n\r]*?)\n\r?/is;
   $date = $1;
   unless (defined $date) {$date = ''};

   $new_locus = 'LOCUS       '.$designation.'   '.$length.'    '.$type.'   '.
                $format.'  '.$domain.'       '.$date."\n";

   $full_header =~ s/LOCUS.*?\n([A-Z\_\.]+\s+)/$new_locus$1/s;
   $full_header .= 'FEATURES             Location/Qualifiers';

   return($full_header);
}

sub default_header {
   my %args = (filename   => 'NONE',
               seq_length => '0',
               definition => 'NONE',
               domain     => 'BCT',
               @_
              );
   my ($defined_fields, $fieldhash, @header);


   #predefine set of default fields
   foreach my $field (@{${$args{auto_ini_ref}}{'gb_order'}}) {
      my ($type, $value);
      'reset' =~ m/reset/;
      $field =~ m/([FS])__(.+)/;
      $type = $1;
      $value = $2;
      $fieldhash->{$value}      = $type;
      $defined_fields->{$value} = ${$args{auto_ini_ref}}{$field};
      push (@header, $value);
   }

   #modify to individual specs
   $defined_fields->{'length'}      = $args{seq_length};
   $defined_fields->{'designation'} = $args{definition};
   $defined_fields->{'domain'}      = $args{domain};
   $defined_fields->{'definition'}  = $args{definition};
   $defined_fields->{'date'}        = &set_date;

   return (&compile_header(defined_fields    => \$defined_fields,
                           fieldhash         => \$fieldhash,
                           header            => \@header,
                           main_window       => $args{main_window},
                           filename          => $args{filename}
                          )
          );
}




1;
