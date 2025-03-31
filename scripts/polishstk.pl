#! /usr/bin/env perl
# (c) Konstantin Riege

use feature ":5.10";
use strict;
use warnings;
use version; our $version = qv('1.0.0');

use Getopt::Long;
use Bio::AlignIO;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use File::Basename;
use Switch;
use POSIX qw(ceil);

# test
# ./scripts/polishstk.pl -i test/polish_test.stk | ./scripts/polishstk.pl
# test ACA-box zero mismatches constrain
# ./scripts/polishstk.pl -i test/polish_test.stk -c test/polish_test.cfg --debug
# ./scripts/polishstk.pl -i test/polish_test.stk -c test/polish_test.cfg -dbnodes data/taxonomy/nodes.dmp -dbnames data/taxonomy/names.dmp --debug

sub usage {
	print <<~EOF;
		DESCRIPTION
		removes gap columns and optional sorts stockholm alignment by taxid

		VERSION
		$version

		OPTIONS
		-h       | --help                 : this message
		-i       | --input       [path]   : default: stdin. path to input stk (seqid: sciname.taxid.n)
		-o       | --output      [path]   : default: stdout. path to output stk
		-dbnodes | --dbncbinodes [path]   : optional. to sort stk by lineage provide path to nodes file of a local NCBI taxonomy
		-dbnames | --dbncbinames [path]   : optional. to sort stk by lineage provide path to names file of a local NCBI taxonomy
		-c       | --constrain   [path]   : optional. path to Rfam cfg file to filter sequences by mismatch constrains
		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
		EOF
	exit 0;
}

my $sort=0;
my $debug=0;
my ($stkfile, $stkout, $nodesfile, $namesfile, $cfgfile);
(Getopt::Long::Parser->new)->getoptions(
	'h|help' => sub{ &usage },
	'i|input=s' => \$stkfile,
	'o|output=s' => \$stkout,
	'dbnodes|dbncbinodes=s' => \$nodesfile,
	'dbnames|dbncbinames=s' => \$namesfile,
	'c|constrain=s' => \$cfgfile,
	'd|debug' => \$debug
) or print $! and exit 1;


my $aln;
if ($stkfile){
	$aln = (Bio::AlignIO->new(-format => 'stockholm', -file => $stkfile, -verbose => -1))->next_aln;
} else {
	$aln = (Bio::AlignIO->new(-format => 'stockholm', -fh => \*STDIN, -verbose => -1))->next_aln;
}
# set flat i,e. remove "/sta-sto" from id when printing and convert gap charter to default "-" note: "." cannot be used in remove_gaps function
$aln->set_displayname_flat;
$aln->map_chars('\.','-');
my $ameta = $aln->consensus_meta;


############# apply filter

if($cfgfile){
	############# contrain filter
	# 0: conserved_orig
	# 1: structure
	# 2: conserved
	# 3: constrain
	# ...
	my @l = ([split//,uc($ameta->named_meta("RF"))=~tr/./-/r]);
	my $len = $l[0]->$#*;
	my $type;
	open CFG,"<",$cfgfile or die $!;
	while(<CFG>){
		chomp;
		next unless /\S/;
		my @F=split/(=|\s+)/;
		$type=$F[-1] if $F[0] eq "type";
		push @l,[(split//,uc($F[-1]))] if $F[0]=~/^(structure|conserved|constrain)/;
	}
	my $i=0;
	for(0..$len){
		if(! defined $l[2]->[$i] || $l[0]->[$i] ne $l[2]->[$i]){
			splice $l[$_]->@*,$i,0,"." for 1..$#l;
		}
		$i++;
	}
	if ($debug) {
		say join("",$l[$_]->@*) for (0..$#l);
	}

	my (@ugaug, @cuga);
	for my $c (3..$#l) {
		my @constrains = grep { $l[$c]->[$_] eq "|" } 0..$len;
		for ($i=0; $i<=$#constrains; $i+=2){
			my $j=$i+1;
			my @ref = $l[0]->@[$constrains[$i]+1..$constrains[$j]-1];
			if($debug){
				say "";
				say $l[$c]->@[$constrains[$i]..$constrains[$j]];
				say @ref;
				say "";
			}

			(join"",$l[$c]->@[$constrains[$i]+1..$constrains[$j]-1])=~/(\d+)/;
			my $allowed_mismatches = $1;
			for my $s ($aln->each_seq()) {
				next unless $s->id =~ /^GORAP/ || $debug;

				my @s = split//,uc($s->subseq($constrains[$i]+2,$constrains[$j]));
				# my $mm = scalar grep { $ref[$_] ne "-" && $s[$_] ne "-" && $s[$_] ne $ref[$_] } 0..$#ref; # allows indels
				# my $mm = scalar grep { $s[$_] ne "-" && $s[$_] ne $ref[$_] } 0..$#ref; # allows deletions but no insertions
				my $mm = scalar grep { $s[$_] ne $ref[$_] } 0..$#ref; # count indels as mismatches
				$aln->remove_seq($s) if $mm > $allowed_mismatches;

				if ($debug){
					say $s->id;
					say @s;
					say @ref;
					say "mm: $mm - allowed: $allowed_mismatches";
				}
			}

			if($type =~ /CD-box/) {
				my $ref = uc((join"",@ref)=~s/\W//gr);
				if ($ref =~/^[AGU]GAU\w\w{0,1}$/){
					@ugaug = ($constrains[$i]+1,$constrains[$j]-1);
				} elsif ($ref =~/^GAU\w\w{0,1}$/){
					@ugaug = ($constrains[$i],$constrains[$j]-1);
				} elsif ($ref =~/^\wUGA$/) {
					@cuga = ($constrains[$i]+1,$constrains[$j]-1);
				} elsif ($ref =~/^UGA$/) {
					@cuga = ($constrains[$i],$constrains[$j]-1);
				}
			}
		}
	}

	############# refined constrain filter for CD-box snoRNAs
	# CD-box kink-turn filter
	# U      != C
	# G - A  0 mm
	# A - G  0 mm
	# U - U  <- check possible pairings
	# G - C  <- check possible pairings
	#(A - U) <- bp only present in 70% sequences
	if($type =~ /CD-box/ && $#ugaug!=-1 && $#cuga !=-1) {
		my @ref_c = $l[0]->@[$ugaug[0]..$ugaug[1]];
		my @ref_d = $l[0]->@[$cuga[0]..$cuga[1]];
		my @idx_c = grep { $ref_c[$_] ne "-" } 0..$#ref_c;
		my @idx_d = grep { $ref_d[$_] ne "-" } 0..$#ref_d;
		if($debug){
			say "";
			say "kink-turn";
			say @ref_c;
			say @ref_d;
			say @ref_c[@idx_c];
			say @ref_d[@idx_d];
			say "";
		}

		for my $s ($aln->each_seq()) {
			next unless $s->id =~ /^GORAP/ || $debug;
			my $mm = 0;
			my @c = (split//,uc($s->subseq($ugaug[0]+1,$ugaug[1]+1)))[@idx_c]; # this selection allows insertions
			my @d = (split//,uc($s->subseq($cuga[0]+1,$cuga[1]+1)))[@idx_d];
			$mm++ if $c[0] eq "C";
			# allows also deletions
			$mm++ unless "$c[1]$d[3]" =~ /[-G][-A]/;
			$mm++ unless "$c[2]$d[2]" =~ /[-A][-G]/;
			$mm++ unless "$c[3]$d[1]" =~/([-U][-U]|[-A][-U]|[-C][-G]|[-G][-C]|[-U][-A])/;
			$mm++ unless "$c[4]$d[0]" =~/([-A][-U]|[-C][-G]|[-G][-C]|[-U][-A])/;
			# $mm++ unless "$c[1]$d[3]" =~ /GA/; # (100%)
			# $mm++ unless "$c[2]$d[2]" =~ /AG/; # (100%)
			# $mm++ unless "$c[3]$d[1]" =~/(UU|AU|CG|GC|UA)/; # UU (95%) or canonical (5%)
			# $mm++ unless "$c[4]$d[0]" =~/(AU|CG|GC|UA)/;    # canonical (100%)
			$aln->remove_seq($s) if $mm > 0;

			if($debug){
				say $s->id;
				say @c;
				say " ".(reverse @d);
				say "mm: $mm";
			}
		}
	}

	############# conservation filter
	if($debug || ($type =~ /(miRNA|CRISPR)/)) {
		my $atmp = $aln;
		for my $s ($atmp->each_seq()) {
			$atmp->remove_seq($s) if $s->id =~ /^GORAP/;
		}
		my @consensus = split//,uc($atmp->consensus_string(90));
		if($debug){
			say "";
			say "consensus";
			say @consensus;
			say "";
		}
		my $allowed_mismatches = ceil(0.1 * scalar grep { $consensus[$_] ne "?" } 0..$#consensus);
		if ($allowed_mismatches >= 1) {
			for my $s ($aln->each_seq()) {
				next unless $s->id =~ /^GORAP/ || $debug;
				my @s = split//,uc($s->seq);
				# do not skip gaps via && $s[$_] ne "-" &&
				my $mm = scalar grep { $consensus[$_] ne "?" && $s[$_] && $consensus[$_] ne $s[$_] } 0..$#s;
				$aln->remove_seq($s) if $mm > $allowed_mismatches;
				if($debug){
					say $s->id;
					say "mm: $mm - allowed: $allowed_mismatches ";
					say "";
				}
			}
		}
	}

	############# structure presence filter
	#       ( ( ( (
	#           ) )
	#       ( ( ( ( print # print if context switch
	#       ( ( ( ( ( (
	#           ) ) ) )
	#     print ) ) ) )   # print if former printed reached + jump over
	# print ) ) ) ) ) )   # print if bottom reached
	#
	# "((.(()).(()).))"
	# 3-4 : 5-6
	# 8-9 : 10-11
	# 0-1 : 13-14
	if($debug || ($type =~ /./)){
		my $ref = uc($ameta->named_meta("RF"))=~tr/./-/r;
		my $s = $ameta->named_meta("SS_cons");
		$s=~tr/[{</(((/;
		$s=~tr/]}>/)))/;
		my @ss=split//,$s;

		my ($c, @open, @close, @regions);
		my $stac=-1;
		for my $i (0..$#ss) {
			if( $ss[$i] eq "(" ){
				if($stac!=-1){
					$c++;
					# say $open[$c]."-".$open[$stac]." : ".$close[$stac]."-".$close[$c];
					push @regions,[$open[$c],$open[$stac],$close[$stac],$close[$c]];
					$stac=-1;
				}
				push @open, $i;
				$c=$#open;
			}
			if( $ss[$i] eq ")" ){
				$stac=$c if $stac==-1;
				$close[$c--]=$i;
				if ($c==-1 || defined $close[$c]){
					$c++;
					# say $open[$c]."-".$open[$stac]." : ".$close[$stac]."-".$close[$c];
					push @regions,[$open[$c],$open[$stac],$close[$stac],$close[$c]];
					$stac=-1;
					$c-- while defined $close[$c] && $c>=0;
				}
			}
		}

		if($debug){
			say "structure";
			say @ss;
			say "open: ".$_->[0]." to ".$_->[1]." close: ".$_->[2]." to ".$_->[3] for @regions;
			say "";
		}

		my $ref_acaregion = uc(substr($ref, $regions[-1]->[3]-3, $aln->length)=~s/-//gr);
		for my $s ($aln->each_seq()) {
			next unless $s->id =~ /^GORAP/ || $debug;
			my @s = split//,$s->seq;
			my $covered = 0;
			for (@regions){
				if( $_==0 && (length($s->subseq($_->[0]+1,$_->[1]+1)=~s/-//gr)==0 || length($s->subseq($_->[2]+1,$_->[3]+1)=~s/-//gr)==0) ){
					# sequence can be truncated
					$covered++ if length((join"",@s[$_->[0]..$_->[1]])=~s/-//gr) >=1 || length((join"",@s[$_->[2]..$_->[3]])=~s/-//gr) >=1;
				} else {
					# inner stem needs minimum size of 2bp
					$covered++ if length((join"",@s[$_->[0]..$_->[1]])=~s/-//gr) >=2 || length((join"",@s[$_->[2]..$_->[3]])=~s/-//gr) >=2;
				}
			}
			$aln->remove_seq($s) if $covered < ($#regions+1)/2;

			if ($debug){
				say $s->id;
				say "covered: $covered - stems: ",($#regions+1);
				say "";
			}

			if ($type =~ /HACA-box/ && $#regions > 0 && $ref_acaregion=~/(GCA|A[CUA]A)/) {
				# 2 stems and HACA type can be used to filter for at least presence of ACA. ANANNA is too hard to find
				# ACA oftem AUA but can be AAA and GCA (the latter RF00045)
				my $mm=1;
				unless (uc($s->subseq($regions[-1]->[3]-2,$aln->length)=~s/-//gr)=~/(GCA|A[CUA]A)/){
					$aln->remove_seq($s);
					$mm=0;
				}
				if($debug){
					say "HACA found: $mm";
					say uc(substr($ref, $regions[-1]->[3]-2, $aln->length));
					say uc($s->subseq($regions[-1]->[3]-1,$aln->length));
					say "";
				}
			}
		}
	}
}


############# removing gap columns
# rfam family sscons contains    (),-.:<>[]_{}2aAbBcCdD
# map to aln accepted char set   abcdefghijklmnopqrstuv
# remove gap columns
# re-add cleaned sequences to keep integer stk with annotations and meta

# generate char maps from Bio::Seq::Meta
my $ss=$ameta->named_meta("SS_cons");
say $ss if $debug;
my @chars=('.','-','a'..'z','A'..'Z');
my (%s2c, %c2s);
my %idx=('.' => 0, '-' => 1);
for (split//,$ss){
	$idx{$_} = exists $idx{$_} ? $idx{$_} : scalar(keys %idx);
	$s2c{$_}=$chars[$idx{$_}];
	$c2s{$chars[$idx{$_}]}=$_;
}
$ss=~s/(.)/$s2c{$1}/g;
if ($debug){
	say $ss;
	my $tmp=$ss;
	$tmp=~s/(.)/$c2s{$1}/g;
	say $tmp;
	say "";
}

# add fake seq to aln and remove gap columns
$aln->add_seq(Bio::LocatableSeq->new(-seq => $ameta->named_meta("RF"), -id  => "RF", -start => 1, -end => length($ameta->named_meta("RF")=~s/[\.-]//gr)));
$aln->add_seq(Bio::LocatableSeq->new(-seq => $ss, -id  => "SS_cons", -start => 1, -end => length($ss=~s/[\.-]//gr)));
$aln->map_chars('\.','-');
my $atmp = $aln->remove_gaps('-',1); # returns new object if altered, else reference!
# update meta
$ameta->named_meta("RF",$atmp->get_seq_by_id("RF")->seq);
$ameta->named_meta("SS_cons",$atmp->get_seq_by_id("SS_cons")->seq=~s/(.)/$c2s{$1}/gr);
# because atmp can be a reference, remove after updating meta
$aln->remove_seq($aln->get_seq_by_id("RF"));
$aln->remove_seq($aln->get_seq_by_id("SS_cons"));
# update sequences
# since input within pipeline does not have acc/sta-sto ids, no get_nse must not be used here
$_->seq($atmp->get_seq_by_id($_->id)->seq) for $aln->each_seq();


############# sorting

my (%taxIDsToLineage, $taxdb);
sub getLineage {
	my ($taxid) = @_;
	my @nodes;
	# check for reuse of prior looked up lineages
	return $taxIDsToLineage{$taxid} if exists $taxIDsToLineage{$taxid};
	my $taxon = $taxdb->get_taxon(-taxonid => $taxid);
	my $tree_functions = Bio::Tree::Tree->new( -verbose => -1);
	push @nodes, $_->id for $tree_functions->get_lineage_nodes($taxon);
	push @nodes, $taxid;
	$taxIDsToLineage{$taxid} = join ",",@nodes;
	return join ",",@nodes;
}
if ($nodesfile){
	$taxdb = Bio::DB::Taxonomy->new(
		-source => 'flatfile',
		-nodesfile => $nodesfile,
		-namesfile => $namesfile,
		-directory => dirname($nodesfile),
		-force => 0,
		-verbose => -1
	);

	my @tosort;
	for ($aln->each_seq()){
		push @tosort,$_;
		$aln->remove_seq($_);
	}
	$aln->add_seq($_) for sort { getLineage((split/\./,$a->id)[-2]) cmp getLineage((split/\./,$b->id)[-2]) } @tosort;
}

############# printing

my $str="";
open my $STKFD,">",\$str;
(Bio::AlignIO->new(-format  => "stockholm", -fh => $STKFD, -verbose => -1))->write_aln($aln);
if ($stkout) {
	open OUT,">",$stkout or die $!;
	for (split/\n/,$str){say OUT unless /^#=GS/};
	close OUT;
} else {
	for (split/\n/,$str){say unless /^#=GS/};
}
