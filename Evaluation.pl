#! /usr/bin/perl

use strict;
use warnings;

#TODO remove 
use lib 'lib';

use Bio::Gorap::Evaluation::Clustering;
use Bio::Gorap::Evaluation::Statistics;
use Bio::Gorap::Parameter;
use Bio::Gorap::DB::Taxonomy;
use File::Spec::Functions;
use POSIX qw(:sys_wait_h);
use Cwd qw(abs_path);

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time());
$year = $year + 1900;
$mon += 1;
print "$mday.$mon.$year-$hour:$min:$sec\n";

#push gorap tools to $PATH
my $PATHtools;
for(reverse glob(catdir($ENV{GORAP},'*','bin'))){
	$PATHtools .= $PATHtools ? ":$_" : $_;
}
local $ENV{PATH} = $ENV{PATH} ? "$PATHtools:$ENV{PATH}" : $PATHtools;
my ($trnalib) = reverse glob catdir($ENV{GORAP},'tRNAscan-SE*');
local $ENV{PERL5LIB} = $ENV{PERL5LIB} ? "$trnalib:$ENV{PERL5LIB}" : $trnalib if $trnalib;

#set defaults
my $parameter = Bio::Gorap::Parameter->new(	
	pwd => $ENV{PWD},
	pid => $$,	
	output => 'gorap_eval',
	commandline => 0
);

#$parameter->set_queries(2544);
&clustering_rfam(1);
&artificial_genome(1);

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time());
$year = $year + 1900;
$mon += 1;
print "$mday.$mon.$year-$hour:$min:$sec\n";

sub clustering_rfam {
	my ($threads) = @_;

	$parameter->threads($threads);

	my $clustering = Bio::Gorap::Evaluation::Clustering->new(
		parameter => $parameter
	);

	my $c=0;
	my $thrs={};
	for my $cfg (@{$parameter->queries}){
		$c++;		
		$parameter->set_cfg($cfg);	
		
		if (scalar(keys %{$thrs}) >= $parameter->threads){
			my $pid = wait();
			delete $thrs->{$pid};		
		}
		print $c.' of ',$#{$parameter->queries}+1,' - '.$parameter->cfg->rf."\n";

		if (my $pid = fork()) {
			$thrs->{$pid}++;
		} else {
			$clustering->write_iid_60_sets();
			exit;
		}
	}
	for (keys %{$thrs} ) {
		my $pid = wait();
		delete $thrs->{$pid};
	}

	for my $cfg (@{$parameter->queries}){		
		$parameter->set_cfg($cfg);	
		$clustering->build_cm_from_sets();		
	}
}

sub artificial_genome {
	my ($threads) = @_;

	$parameter->threads($threads);

	my $taxdb = Bio::Gorap::DB::Taxonomy->new(
		parameter => $parameter
	);
	
	my $statistic_functions = Bio::Gorap::Evaluation::Statistics->new(
		parameter => $parameter,
		taxdb => $taxdb
	);	

	$parameter->set_genomes($statistic_functions->create_genomes());

	for (0..$#{$parameter->genomes}){		
		my $genome = ${$parameter->genomes}[$_];
		my $abbr = ${$parameter->abbreviations}[$_];
		my $kingdom = $abbr;
		print "perl Gorap_noTaxonomy.pl -i $genome -a $abbr -k $kingdom -o $parameter->output -c $threads\n";
	}
}