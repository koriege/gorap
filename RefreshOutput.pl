#! /usr/bin/perl

use strict;
use warnings;

#TODO remove 
use lib 'lib';

use Cwd 'abs_path';
use Bio::Gorap::ThrListener;
use Bio::Gorap::Parameter;
use Bio::Gorap::DB::GFF;
use Bio::Gorap::DB::STK;
use Bio::Gorap::DB::Fasta;
use Bio::Gorap::DB::Taxonomy;
use Bio::Gorap::Evaluation::HTML;
use Bio::Gorap::Functions::ToolParser;
use File::Spec::Functions;
use File::Basename;
use Try::Tiny;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$year = $year + 1900;
$mon += 1;
print "$mday.$mon.$year-$hour:$min:$sec\n";
my $stamp = "$mday.$mon.$year-$hour:$min:$sec";

#set defaults and read CLI parameter or parameter file
my $parameter = Bio::Gorap::Parameter->new(
	#pwd => dirname(abs_path(__FILE__)),
	pwd => $ENV{PWD},
	pid => $$,
	commandline => 1
);

#gorap gff3 storage database initialization without loss of existing data in output directory 
my $gffdb = Bio::Gorap::DB::GFF->new(
	parameter => $parameter
); 

my $fastadb = Bio::Gorap::DB::Fasta->new(
	parameter => $parameter,
	do_chunks => 0
);

my $stkdb = Bio::Gorap::DB::STK->new(
	parameter => $parameter
);

sub get_overlaps {
	my ($f) = @_;
		
	my @tmp = split /\./, $f->seq_id;
	my $abbr = $tmp[0];
	pop @tmp;
	my $id = join '.' , @tmp;

	my @features;
	for ($gffdb->db->{$abbr}->features()){		
		next unless $_->display_name eq '!';
		@tmp = split /\./, $_->seq_id;
		pop @tmp;		
		next unless join('.' , @tmp) eq $id;
		my ($start, $stop, $strand) = $f->intersection($_);
		push @features , $_ if $start && ($stop - $start) > 10;
	}
	
	return @features;
}

for my $cfg (@{$parameter->queries}){
	$parameter->set_cfg($cfg);
	my $type = $parameter->cfg->rf_rna;
	my $hold;
	for my $f (@{$gffdb->get_all_features($type)}){
		my $seq;
		$seq = $stkdb->db->{$type}->get_seq_by_id($f->seq_id) if exists $stkdb->db->{$type};
		if ($seq){
			if($parameter->force && ($f->type=~/_Afu/ || $f->type=~/_SNORD/ || $f->type=~/_sn?o?s?n?o?[A-WYZ]+[a-z]?\d/)){
				my $higherscore;
				for (&get_overlaps($f)){
					if($_->type=~/_Afu/ || $_->type=~/_SNORD/ || $_->type=~/_sn?o?s?n?o?[A-WYZ]+[a-z]?\d/){
						$higherscore = 1 if $_->score > $f->score;
					}
				}
				if ($higherscore){
					$gffdb->update_filter($f->seq_id,$type,'O');
					$stkdb->db->{$type}->remove_seq($seq);
				} else {
					$hold=1;
					$gffdb->update_filter($f->seq_id,$type,'!');	
				}
			} else {
				$hold=1;
				$gffdb->update_filter($f->seq_id,$type,'!');
			}
		} else {
			$gffdb->update_filter($f->seq_id,$type,'X') if $f->display_name eq '!' && $f->primary_tag!~/SU_rRNA/;
		}		
	}	
	if ( ! $hold && exists $stkdb->db->{$type} && $parameter->force){
		unlink $stkdb->idToPath->{$type};
		delete $stkdb->db->{$type};
	}
}
print "Storing changes\n";
$stkdb->store;
$gffdb->store_overlaps;

Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$fastadb->oheaderToDBsize,$stkdb->idToPath,"$mday.$mon.$year-$hour:$min:$sec-amended");

unlink $_ for glob catfile($parameter->tmp,$parameter->pid.'.*');

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$year = $year + 1900;
$mon += 1;
print "$mday.$mon.$year-$hour:$min:$sec\n";

__END__

=head1 NAME

GORAP - Genomewide ncRNA Annotation Pipeline

=head1 SYNOPSIS

RefreshOutput.pl [OPTION]...
  
example: RefreshOutput.pl -o results -force

=head1 DESCRIPTION

Updates annotation GFF and FASTA files for manual amended Stockholm alignment files

B<-h>, B<--help>	

	this (help) message

B<-file>, B<--file>=F<FILE> 
	
	(optional)
	(default: E. coli example with parameter file at: $GORAP/parameter/parameter.txt) 
	run GORAP with a parameter file
	

-------------------------
	
		
B<-i>, B<--fastas>=F<FILE>,...

	(default: E. coli example genome) 
	(regex) path(s) of comma separated species FASTA file(s)
		
B<-a>, B<--abbreviations>=I<abbreviation,> ...	

	(default: build from FASTA file name(s)) 
	list of comma separated abbreviations as unique identifiers
	
B<-o>, B<--output>=F<PATH>	

	(default: pwd/gorap_out) 
	output directory

B<-t>, B<--tmp>=F<PATH>	

	(default: $TMPDIR or /tmp or $GORAP/tmp) 
	set the temporary directory
	
B<-c>, B<--cpu>=I<INT>	

	(default: 1) 
	count of cpu cores to use
	
B<-k>, B<--kingdom>=I<bac,arc,euk,fungi,virus>

	(default: all) 
	list of comma separated kingdoms to screen for kingdom specific ncRNAs only
	
B<-q>, B<--queries>=I<RF00001:RF00005,RF00008,> ...

	(default: all Rfam families)
	list of comma separated Rfam ids/numbers or ranges by ':'
	to enable additional features only, skip annotation with 0
	 
B<-r>, B<--rank>=[I<INT/STRING>]	

	(optional) 
	NCBI taxonomy matching scientific name or taxonomy id of rank/genus/... for given species
	please escape white spaces i.e. Bacillus\\ subtilis\\ group or write 'Bacillus subtilis group'
	
B<-s>, B<--species>=[I<INT/STRING>]

	(optional) 
	NCBI taxonomy matching scientific name or taxonomy id of given species
	please escape white spaces i.e. Bacillus\\ subtilis or write 'Bacillus subtilis'
								
B<-og>, B<--outgroup>=F<FILE>

	(optional) 
	path to an additional species FASTA file to be used as outgroup for
	SSU rRNA based phylogeny reconstruction including all given FASTA files of B<-i>

B<-oga>, B<--ogabbreviation>=F<FILE>
	
	(optional)
	(default: build from FASTA file name of B<-og> 
	abbreviation as unique identifiers
	
B<-b>, B<--bam>=F<FILE>,...

	(optional) 
	(regex) paths(s) of comma separated mapping results as SAM/BAM file(s)	

-force, --force
	
	(optional)
	remove alignment files of no query sequence hits remains

=head1 AUTHOR

Konstantin Riege, E<lt>konstantin.riege@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 by Konstantin Riege

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.14.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 SEE ALSO

The full description of all in-/outputs and parameters is maintained as PDF manual. 
You can access it on L<www.rna.uni-jena.de/software.php>.

=cut