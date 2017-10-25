package Bio::Gorap::Gorap;

use version; our $VERSION = "2.3.1";

1;

=head1 NAME

Bio::Gorap - A Perl distribution for genome wide ncRNA annotation based on various software 
wrappers, NCBI taxonomy functions and Stockholm alignment file manipulation.

head1 DESCRIPTION

Bio::Gorap is a distribution of Perl modules to build and extend a scripted main 
pipeline to efficiently annotate non-coding RNAs in genomic sequences. It screens for 
user defined or whole database provided queries and is highly configurable. It deals with 
different software, databases (NCBI, Rfam, Silva, mirBASe), offers an artificial evaluation 
study, HTML based visualization of results and can be easily extended by new annotation tools.

The following modules are part of the Bio::Gorap package.

=over 

=item L<Bio::Gorap::DB::Fasta>: Routines for storing, splitting with overlap and accessing genomic
sequences implemented via a Moose to L<Bio::DB::Fasta>.

=back

=head1 AUTHORS

=over

=item Konstantin Riege E<lt>konstantin.riege@uni-jena.deE<gt>

=back

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 Konstantin Riege E<lt>konstantin.riege@uni-jena.deE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.14.4 or,
at your option, any later version of Perl 5 you may have available.

This software is distributed without any warranty.

=cut
