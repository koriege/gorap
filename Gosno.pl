#!/usr/bin/perl
use strict;
use warnings;
use v5.10;

use Bio::Gorap::Functions::STK;
use Bio::SeqFeature::Generic;
use Bio::Gorap::CFG;
use Bio::AlignIO;

use Getopt::Long;
use File::Spec::Functions;

(Getopt::Long::Parser->new)->getoptions (
	's|scores=s' => \my $scorefile,
	'stk|stk=s' => \my $stkfile,
	'r|regex=s'	=> \my $regex,
	'rf|rfid=s' => \my $rf,
	'v|verbose' => \my $verbose,
	'h|help' => \my $help
) or &usage();

&usage() if $help || ! $stkfile || ! $rf || ! $ENV{GORAP};

sub usage(){
	print <<EOF;
	GOSNO filters a GORAP based C/D-Box snoRNA Stockholm alignment
	-> install GORAP and export all variables

	-stk  | --stk        => (reqiered) <FILE> Rfam based snoRNA STK
	-s    | --scores     => (reqiered) <FILE> related GFF or cmalign scores file
	-rf   | --rfid       => (reqiered) <STRING> related Rfam id (accession number)
	-r    | --regex      => (optional) <'STRING'> regular expression for sequence ids
	-v    | --verbose    => print to STDOUT
	-h    | --help       => this message
EOF
	exit;
}

$regex = '.*' unless $regex;

my ($cfgfile) = glob catfile($ENV{GORAP},'config',$rf.'*.cfg');
say $cfgfile;
my $cfg = Bio::Gorap::CFG->new( cfg => $cfgfile);

say $stkfile if $verbose;
die 'No C/D-Box snoRNA STK '.$! if $#{$cfg->constrains} == -1;

my $stk = (Bio::AlignIO->new(-format => 'stockholm', -file => $stkfile, -verbose => -1))->next_aln or die "Corrupt STK file $stkfile ".$!;
$stk->set_displayname_flat;

my $scores;
open S , '<'.$scorefile or die $!;
while (<S>){	
	next if $_=~/^#/ || $_=~/^\s*$/;
	chomp $_;
	$_=~s/^\s+//;
	my @l = split /\s+/ , $_;
	$scores->{$l[1]} = $l[6];
	$scores->{$l[0]} = $l[5];
}
close S;

my @features;
for my $seq ($stk->each_seq){		
	next unless $seq->display_id=~/$regex/;	
	say $seq->display_id.' '.$scores->{$seq->display_id} if $verbose;		
	push @features , Bio::SeqFeature::Generic->new(-seq_id => $seq->display_id, -primary_tag => $cfg->rf_rna, -score => $scores->{$seq->display_id});
}
($stk) = Bio::Gorap::Functions::STK->user_filter($stk, \@features, $cfg->constrains, $cfg->cs, $cfg->stk);
if ($verbose){
	say 'filtered:';
	for my $seq ($stk->each_seq){			
		next unless $seq->display_id=~/$regex/;		
		say $seq->display_id;		
	}
	say '';
}

(Bio::AlignIO->new(-format => 'stockholm', -file => '>'.$stkfile, -verbose => -1))->write_aln($stk) or die $!;
