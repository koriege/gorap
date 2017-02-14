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
	'h|help' => \my $help,
	'o|out=s' => \my $outstk
) or &usage();

&usage() if $help || ! $stkfile || ! $outstk || ! $rf || ! $ENV{GORAP};

sub usage(){
	print <<EOF;
	GOSNO filters a GORAP based C/D-Box or HACA snoRNA Stockholm alignment
	-> install GORAP and export all variables

	-stk  | --stk        => (reqiered) <FILE> Rfam based snoRNA STK
	-o    | --out        => (reqiered) <FILE> output STK
	-s    | --scores     => (reqiered) <FILE> related GFF or cmalign scores file
	-rf   | --rfid       => (reqiered) <STRING> related Rfam id (accession number: RF012345)
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

my $features;
for my $seq ($stk->each_seq){		
	next unless $seq->display_id=~/$regex/;	
	say $seq->display_id.' '.$scores->{$seq->display_id} if $verbose && exists $scores->{$seq->display_id};
	push @{$features} , Bio::SeqFeature::Generic->new(-seq_id => $seq->display_id, -primary_tag => $cfg->rf_rna, -score => $scores->{$seq->display_id});
}

($stk, $features) = Bio::Gorap::Functions::STK->structure_filter($stk, $features);
$stk = &remove_gap_columns_and_write($stk,$outstk);

($stk, $features) = Bio::Gorap::Functions::STK->user_filter($stk, $features, $cfg->constrains, $cfg->cs, $cfg->stk);
$stk = &remove_gap_columns_and_write($stk,$outstk);


if ($verbose){
	say 'remaining:';
	for my $seq ($stk->each_seq){			
		next unless $seq->display_id=~/$regex/;		
		say $seq->display_id if exists $scores->{$seq->display_id};
	}
}

sub remove_gap_columns_and_write {
	my ($stk, $file) = @_;
	my $tmpstk = $stk;	
	
	my ($ss , $sc) = Bio::Gorap::Functions::STK->get_ss_cs_from_object($stk);
	my @ss = split // , $ss;
	my @sc = split // , $sc;
	
	my $del = 0 ;
	$tmpstk->map_chars('\.','-');
	my @gapcolm = @{$stk->gap_col_matrix};
	for (my $i=$#gapcolm; $i>=0; $i-- ){	
		my $absent = 0;
		for my $k (keys %{$gapcolm[$i]}){				
			$absent += $gapcolm[$i]->{$k};
		}		
		
		if($absent == $tmpstk->num_sequences){
			$del = 1;
			$tmpstk = $tmpstk->remove_columns([$i,$i]);
			splice(@sc, $i, 1);
			splice(@ss, $i, 1);		
		}	
	}
	
	$tmpstk->set_displayname_flat;
	
	if ($del){
		my $out;		
		open my $READER, '>', \$out;
		(Bio::AlignIO->new(-format  => 'stockholm', -fh => $READER, -verbose => -1 ))->write_aln($tmpstk);
		close $READER;
		my @lines = split /\n/ , $out;
		$lines[-1] = '#=GC RF '.join('',(' ') x ($tmpstk->maxdisplayname_length()-7)).' '.join('' , @sc);
		push @lines , '#=GC SS_cons '.join('',(' ') x ($tmpstk->maxdisplayname_length()-12)).' '.join('' , @ss);
		push @lines, '//';
		
		open STK , '>'.$file or die $!;
		print STK $_."\n" for @lines;
		close STK;

		my $retaln = (Bio::AlignIO->new(-file => $file, -format => 'stockholm'))->next_aln;
		$retaln->set_displayname_flat;
		return $retaln;
	} else {
		(Bio::AlignIO->new(-format => 'stockholm', -file => '>'.$file, -verbose => -1))->write_aln($stk);
		return $stk;
	}
		
}