#! /usr/bin/env perl
# (c) Konstantin Riege

use feature ":5.10";
use strict;
use warnings;
use version; our $version = qv('1.0.0');

use Tree::Simple;
use Bio::Tree::Node;
use Bio::TreeIO;
# use Bio::DB::EUtilities;
# use Try::Tiny;
# use XML::Simple;
use Getopt::Long;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use Bio::DB::Fasta;
use Bio::SeqIO;
use File::Basename;
use List::Util qw(min max sum);
use Switch;

# example for taxonomy threshold
# ./scripts/getcutoff.pl -dbnodes data/taxonomy/nodes.dmp -dbnames data/taxonomy/names.dmp -silva data/taxonomy/silva.newick -taxid 5052 -family data/rfam/RF00001.tsv -fasta data/rfam/RF00001.full.fa
# example for different kingdom fungi
# ./scripts/getcutoff.pl -dbnodes data/taxonomy/nodes.dmp -dbnames data/taxonomy/names.dmp -silva data/taxonomy/silva.newick -taxid 1423 -family data/rfam/RF01846.tsv
# example for same kingdom but overrepresentation
# ./scripts/getcutoff.pl -dbnodes data/taxonomy/nodes.dmp -dbnames data/taxonomy/names.dmp -silva data/taxonomy/silva.newick -taxid 1423 -family data/rfam/RF01772.tsv

sub usage {
	print <<~EOF;
		DESCRIPTION
		Filters Rfam family members by a given ancestor taxid and thus returns a taxonomic trusted bitscore cutoff.
		Given a fasta file, according to the input taxid representative descendants sequences will be sampled.

		VERSION
		$version

		OPTIONS
		-h       | --help                    : this message
		-dbnodes | --dbncbinodes [path]      : required. set path to nodes file of a local NCBI taxonomy
		-dbnames | --dbncbinames [path]      : required. set path to names file of a local NCBI taxonomy
		-silva   | --silvanewick [path]      : optional. path to silva newick
		-family  | --familytsv   [path]      : required. path to Rfam tsv file (6 columns: acc taxid sciname rfamid bitscore rfamname)
		-fasta   | --familyfasta [path]      : optional. path to Rfam fasta file (seqid: >sciname.taxid acc)
		-faout   | --fastaoutput [path]      : optional. path to filtered fasta file
		-taxid   | --taxid       [taxid]     : required. input for taxon id. mutually exclusive to -name
		-name    | --name        [species]   : required. input for scientific name. mutually exclusive to -taxid.
		                                       INFO: use one of bac|arc|euk|fungi|plants|virus for simple kingdom filter

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
		EOF
	exit 0;
}

usage if $#ARGV == -1;

my ($nodesfile, $namesfile, $silvafile, $familyfile, $fastafile, $fastaout, $query, $taxid);
(Getopt::Long::Parser->new)->getoptions(
	'h|help' => sub{ &usage },
	'dbnodes|dbncbinodes=s' => \$nodesfile,
	'dbnames|dbncbinames=s' => \$namesfile,
	'silva|silvanewick=s' => \$silvafile,
	'family|familytsv=s' => \$familyfile,
	'fasta|familyfasta=s' => \$fastafile,
	'faout|fastaoutput=s' => \$fastaout,
	'name|name=s' => \$query,
	'taxid|taxid=s' => \$taxid
) or print $! and exit 1;

############## initialize databases

my ($fadb, $sio, %acc2faid);
sub parse {
	my ($header) = @_;
	$header=~/^>\s*(\S+).*\s+(\S+)$/;
	push $acc2faid{$2}->@*,$1;
	return $1;
}
if ($fastafile) {
	# need forcing re-indexing to trigger id parser
	$fadb = Bio::DB::Fasta->new($fastafile, -makeid => \&parse, -index_name => "$fastafile.idx", -reindex => 1, clean => 1);
	if ($fastaout){
		$sio = Bio::SeqIO->new(-file => ">$fastaout", -format => "fasta", -width => 999999);
	} else {
		$sio = Bio::SeqIO->new(-fh => \*STDOUT, -format => "fasta", -width => 999999);
	}
}

my $silva;
if ($silvafile){
	$silva = (Bio::TreeIO->new(-format => 'newick', -file => $silvafile, -verbose => -1))->next_tree;
} else {
	$silva = Bio::Tree::Tree->new( -verbose => -1);
}

my $taxdb = Bio::DB::Taxonomy->new(
	-source => 'flatfile',
	-nodesfile => $nodesfile,
	-namesfile => $namesfile,
	-directory => dirname($nodesfile),
	-force => 0,
	-verbose => -1
);

############## function definitions

sub flatname {
	my ($name) = @_;
	$name=~s/^[[({](\w+)[]})]/$1/;
	$name=~s/(\W).*\1$//;
	$name=~s/[\W._-]*[[({].*[]})]$//;
	$name=~s/[\W._-]+$//;
	$name=~s/^[\W._-]+$//;
	$name=~s/[\W_-]+/_/g;
	$name=~s/_+/_/g;
	$name=~s/\.+/./g;
	return $name;
}

sub flatquery {
	my ($name) = @_;
	$name=lc($name);
	$name=~s/(^[\W._-]+|[\W._-]+$)//g;
	$name=~s/_+/ /g;
	$name=~s/\.+/\./g;
	return $name;
}

my %taxIDsToLineage;
sub getLineageNodes {
	my ($taxon) = @_;
	my @nodes;
	if ($taxon) {
		my $taxid = $taxon->id;
		# check for reuse of prior looked up lineages
		return $taxIDsToLineage{$taxid} if exists $taxIDsToLineage{$taxid};
		my $tree_functions = Bio::Tree::Tree->new( -verbose => -1);
		push @nodes, $_->id for $tree_functions->get_lineage_nodes($taxon);
		push @nodes, $taxid;
		push $taxIDsToLineage{$taxid}->@*, $_ for @nodes;
	}
	return \@nodes;
}

my %sci2tax;
sub query2taxid {
	my ($name) = @_;

	return $sci2tax{$name} if exists $sci2tax{$name};

	my $taxon = $taxdb->get_taxon(-name => $name);
	my $taxid;

	if ($taxon){
		$taxid = $taxon->id;
		$name = flatname($taxon->scientific_name);
		$sci2tax{$name} = $taxid;
	}
	# too fragile for production
	# else {
	# 	my ($try, $error, $tries) = (0, 1, 10);
	# 	while ($error && (++$try) < $tries){
	# 		$error = 0;
	# 		try {
	# 			my $factory = Bio::DB::EUtilities->new(
	# 				-eutil => 'esearch',
	# 				-db => 'taxonomy',
	# 				-email => 'mymail@foo.bar',
	# 				-term => $name,
	# 				-verbose => -1
	# 			);
	# 			($taxid) = $factory->get_ids;
	# 		} catch {
	# 			$error = 1;
	# 		};
	# 	}
	# }
	return $taxid;
}

############## options sanity checks

unless ($taxid){
	switch ($query) {
		case "virus" { $taxid=10239; }
		case "euk" { $taxid=2759; }
		case "arc" { $taxid=2157; }
		case "bac" { $taxid=2; }
		case "fungi" { $taxid=4751; }
		case "plants" { $taxid=33090; }
		else {
			$query = flatquery($query);
			$taxid = query2taxid($query);
			unless ($taxid){
				say STDERR ":ERROR: no taxon id found for flattened query $query";
				exit 1;
			}
		}
	}
}

my $root = $taxdb->get_taxon(-taxonid => $taxid);
my $rootnodes = getLineageNodes($root);

############## parse info file and check for overrepresentation

my (%acc2sci, %acc2tax, %acc2bs, %acc2rf, %kingdoms, %overrepresented, $types, $family);
my $minbs="+inf";
my $checked=0;
my $familysize=0;
my $thfactor=0.8;

open MAP, '<' , $familyfile or die $!;
while(<MAP>){
	chomp;
	my ($acc,$tax,$sci,$rf,$bs,$t) = split /\t/;
	$types=$t;
	$family = $rf;
	$acc2sci{$acc} = $sci;
	$acc2tax{$acc} = $tax;
	$acc2bs{$acc} = min(exists $acc2bs{$acc} ? $acc2bs{$acc} : "+inf", $bs);
	$sci2tax{$sci} = $tax;
	$minbs = min($minbs,$bs);

	unless ($checked){
		$checked=1;

		switch($types) {
			case /plant/i {
				$kingdoms{"plants"} = "+inf";
			}
			case /(^afu_|yeast|fungi)/i {
				$kingdoms{"fungi"} = "+inf";
			}
			case /euka/i {
				$kingdoms{"euk"} = "+inf";
				$kingdoms{"fungi"} = "+inf";
				$kingdoms{"plants"} = "+inf";
			}
			case /(metazo|protoz)/i {
				$kingdoms{"euk"} = "-inf";
				# prevent adding plants and fungi below
			}
			case /arch/i {
				$kingdoms{"arc"} = "+inf";
			}
			case /bact/i {
				$kingdoms{"bac"} = "+inf";
			}
			case /;CRISPR/ {
				$kingdoms{"arc"} = "+inf";
				$kingdoms{"bac"} = "+inf";
			}
			case /(;splicing|^U\d(\d|atac))/ {
				$kingdoms{"euk"} = "+inf";
				$kingdoms{"fungi"} = "+inf";
				$kingdoms{"plants"} = "+inf";
			}
			case /;snRNA/ {
				$kingdoms{"arc"} = "+inf";
				$kingdoms{"euk"} = "+inf";
				$kingdoms{"fungi"} = "+inf";
				$kingdoms{"plants"} = "+inf";
			}
			case /;miRNA/ {
				$kingdoms{"euk"} = "+inf";
				$kingdoms{"plants"} = "+inf";
			}
		}
	}

	my $nodes = getLineageNodes($taxdb->get_taxon(-taxonid => $tax));
	if ($nodes->$#* != -1){
		unless ($nodes->[0] == 28384 || $nodes->[0] == 12908){
			$familysize++;
			$kingdoms{"virus"}++ if $nodes->[0] == 10239;
			$kingdoms{"euk"}++ if $nodes->[1] == 2759;
			$kingdoms{"arc"}++ if $nodes->[1] == 2157;
			$kingdoms{"bac"}++ if $nodes->[1] == 2;
			$kingdoms{"fungi"}++ if $nodes->[3] && $nodes->[3] == 4751;
			$kingdoms{"plants"}++ if $nodes->[2] && $nodes->[2] == 33090;

			$overrepresented{$nodes->[-2]}++;
			$overrepresented{$nodes->[-3]}++;
		}
	}
}
close MAP;

if ($familysize != 0 && exists $kingdoms{"virus"} && $kingdoms{"virus"} >= $familysize*0.95){
	delete $kingdoms{"euk"};
	delete $kingdoms{"arc"};
	delete $kingdoms{"bac"};
	delete $kingdoms{"fungi"};
	delete $kingdoms{"plants"};
}
if ($familysize != 0 && exists $kingdoms{"fungi"} && $kingdoms{"fungi"} >= $familysize*0.95){
	delete $kingdoms{"virus"};
	delete $kingdoms{"euk"};
	delete $kingdoms{"arc"};
	delete $kingdoms{"bac"};
	delete $kingdoms{"plants"};
}
if ($familysize != 0 && exists $kingdoms{"plants"} && $kingdoms{"plants"} >= $familysize*0.95){
	delete $kingdoms{"virus"};
	delete $kingdoms{"euk"};
	delete $kingdoms{"arc"};
	delete $kingdoms{"bac"};
	delete $kingdoms{"fungi"};
}
if ($familysize != 0 && exists $kingdoms{"arc"} && $kingdoms{"arc"} >= $familysize*0.95){
	delete $kingdoms{"virus"};
	delete $kingdoms{"euk"};
	delete $kingdoms{"bac"};
	delete $kingdoms{"fungi"};
	delete $kingdoms{"plants"};
}
if ($familysize != 0 && exists $kingdoms{"bac"} && $kingdoms{"bac"} >= $familysize*0.95){
	delete $kingdoms{"virus"};
	delete $kingdoms{"euk"};
	delete $kingdoms{"arc"};
	delete $kingdoms{"fungi"};
	delete $kingdoms{"plants"};
}
if ($familysize != 0 && exists $kingdoms{"euk"} && $kingdoms{"euk"} >= $familysize*0.95){
	delete $kingdoms{"virus"};
	delete $kingdoms{"arc"};
	delete $kingdoms{"bac"};
	$kingdoms{"fungi"} = 1;
	$kingdoms{"plants"} = 1;
}

if ($rootnodes->$#* != -1){
	my $abort = 0;
	if ($rootnodes->[0] == 10239){
		$abort = 1 unless exists $kingdoms{"virus"};
	} elsif	($rootnodes->[1] == 2157) {
		$abort = 1 unless exists $kingdoms{"arc"};
	} elsif	($rootnodes->[1] == 2) {
		$abort = 1 unless exists $kingdoms{"bac"};
	} elsif ($rootnodes->[3] && $rootnodes->[3] == 4751) {
		$abort = 1 unless exists $kingdoms{"fungi"};
	} elsif ($rootnodes->[2] && $rootnodes->[2] == 33090) {
		$abort = 1 unless exists $kingdoms{"plants"};
	} elsif ($rootnodes->[1] == 2759) {
		$abort = 1 unless exists $kingdoms{"euk"};
	}

	if ($abort) {
		say join " ",(
			'family="'.$family.'";',
			'kingdoms="'.(join" ",sort {$a cmp $b} keys %kingdoms).'";',
			'cutoff="'."NA".'";',
			'minimum="'.$minbs.'";',
			'info="'."no cutoff. family belongs to different kingdom".'";'
		);
		exit 0;
	}

	my @overrepresented = reverse sort { $overrepresented{$a} <=> $overrepresented{$b} } keys %overrepresented;
	if ($familysize != 0 && $overrepresented[0] && ($overrepresented{$overrepresented[0]} == $familysize || $overrepresented{$overrepresented[0]} > max(1,$familysize/3))){

		$abort = 1;
		for my $overtax (@overrepresented){
			my $overtaxon = $taxdb->get_taxon(-taxonid => $overtax);
			my $overnodes = getLineageNodes($overtaxon);
			# lineage of input taxonid does not match any of the overrepresented species
			my $i = min($rootnodes->$#*,$overnodes->$#*);
			if ($rootnodes->[$i] == $overnodes->[$i]){
				$abort = 0;
				last;
			}
		}
		if ($abort){
			say join " ",(
				'family="'.$family.'";',
				'kingdoms="'.(join" ",sort {$a cmp $b} keys %kingdoms).'";',
				'cutoff="'."NA".'";',
				'minimum="'.$minbs.'";',
				'info="'."no cutoff. family dominated by ".flatname($taxdb->get_taxon(-taxonid => $overrepresented[0])->scientific_name)." ($overrepresented[0])".'";'
			);
			exit 0;
		}
	}
}

if($rootnodes->$#* <= 3){
	say join " ",(
		'family="'.$family.'";',
		'kingdoms="'.(join" ",sort {$a cmp $b} keys %kingdoms).'";',
		'cutoff="'.$minbs.'";',
		'minimum="'.$minbs.'";',
		'info="'."trusted cutoff".'";'
	);
	exit 0;
}

############## get linages from ncbi

# if ancestors are linear, go up first
my $taxon = $root;
my $ancestor = $taxdb->ancestor($root);
while( scalar $taxdb->each_Descendent($ancestor) == 1){
	$taxon = $ancestor;
	$ancestor = $taxdb->ancestor($ancestor);
}

my %relatedIDsToLineage;
for my $t ($taxon,$taxdb->get_all_Descendents($taxon)){
	$relatedIDsToLineage{$t->id} = getLineageNodes($t);
}

############## silva based refinement

# faster access than repetetive $silva->find_node(-id => $id);
my %silvaIDsToNodes;
for ($silva->get_nodes){
	$silvaIDsToNodes{$_->id} = $_ if $_->id;
}

# add missing via BFT: find inner node in silva, iterate over silva subtree to add missing species
my @q;
push @q, $taxon;
while($#q > -1){
	my $v = shift @q;
	my $s = $silvaIDsToNodes{$v->id};
	if ($s) {
		for ($s,$s->get_all_Descendents()){
			next unless $_->id;
			unless (exists $relatedIDsToLineage{$_->id}){
				my $t = $taxdb->get_taxon(-taxonid => $_->id);
				$relatedIDsToLineage{$t->id} = getLineageNodes($t);
			}
		}
	} else {
		push @q, $taxdb->each_Descendent($v);
	}
}

# remove misplaced: find ancestor in silva, remove species if not in silva subtree but in whole tree
my $v = $root;
my $found = 0;
while($v){
	if (exists $silvaIDsToNodes{$v->id}){
		$found = 1;
		my %nodes;
		for ($silvaIDsToNodes{$v->id}->get_all_Descendents()){
			next unless $_->id;
			$nodes{$_->id}=1
		}
		for (keys %relatedIDsToLineage){
			if (exists $silvaIDsToNodes{$_} && ! exists $nodes{$_}){
				delete $relatedIDsToLineage{$_};
			}
		}
		last;
	}
	$v = $taxdb->ancestor($v);
}

# if no ancestor found in silva: compute lca from species in silva, remove species if not in silva subtree but in whole tree
unless ($found){
	my @nodes;
	for (keys %relatedIDsToLineage){
		push @nodes, $silvaIDsToNodes{$_} if exists $silvaIDsToNodes{$_};
	}
	if ($#nodes != -1){
		my $lca = $silva->get_lca(-nodes => \@nodes);
		my %tmp;
		for ($lca->get_all_Descendents()){
			next unless $_->id;
			$tmp{$_->id}=1;
		}
		for (keys %relatedIDsToLineage){
			if (exists $silvaIDsToNodes{$_} && ! exists $tmp{$_}){
				delete $relatedIDsToLineage{$_};
			}
		}
	}
}

############## finally

my $threshold = "+inf";
my %threshold2acc;
my @thresholds;

for my $acc (sort keys %acc2tax){
	if (exists $relatedIDsToLineage{$acc2tax{$acc}}){
		# attention a single accession can have multiple sequences with heavily differing scores. in acc2bs hash is only the minimum
		# $sio->write_seq($fadb->get_Seq_by_id($_)) for $acc2faid{$acc}->@*;
		push $threshold2acc{$acc2bs{$acc}}->@*,$_ for $acc2faid{$acc}->@*;
		$threshold=min($threshold,$acc2bs{$acc});
	}
}

if ($threshold == "+inf") {
	say join " ",(
		'family="'.$family.'";',
		'kingdoms="'.(join" ",sort {$a cmp $b} keys %kingdoms).'";',
		'cutoff="'.$minbs.'";',
		'minimum="'.$minbs.'";',
		'info="'."trusted cutoff. no related species found".'";'
	);
} else {
	@thresholds = sort {$a <=> $b} keys %threshold2acc;
	$sio->write_seq($fadb->get_Seq_by_id($_)) for $threshold2acc{$thresholds[0]}->@*; # contains minimum
	if ($#thresholds>0) {
		$sio->write_seq($fadb->get_Seq_by_id($_)) for $threshold2acc{$thresholds[-1]}->@*; # contains max
	}
	if ($#thresholds>1) {
		my ($q1,$q3) = (0,0);
		# if($#thresholds>2){
		# 	my $q1 = int(($#thresholds+1)*0.25);
		# 	my $q3 = int(($#thresholds+1)*0.75);
		# 	$sio->write_seq($fadb->get_Seq_by_id($_)) for $threshold2acc{$thresholds[$q1]}->@*; # contains q1
		# 	$sio->write_seq($fadb->get_Seq_by_id($_)) for $threshold2acc{$thresholds[$q3]}->@*; # contains q3
		# }
		my $medix = int(($#thresholds+1)*0.5);
		if ($medix != $q1 && $medix != $q3) {
			$sio->write_seq($fadb->get_Seq_by_id($_)) for $threshold2acc{$thresholds[$medix]}->@*; # contains median
		}
		my $mean = sum(@thresholds)/($#thresholds+1);
		pop @thresholds while $thresholds[-1] > $mean;
		if($#thresholds!=$medix && $#thresholds!=$q1 && $#thresholds!=$q3){
			$sio->write_seq($fadb->get_Seq_by_id($_)) for $threshold2acc{$thresholds[-1]}->@*; # contains mean
		}
	}

	say join " ",(
		'family="'.$family.'";',
		'kingdoms="'.(join" ",sort {$a cmp $b} keys %kingdoms).'";',
		'cutoff="'.$threshold.'";',
		'minimum="'.$minbs.'";',
		'info="'."taxonomy based cutoff. minimum score from descendants of ".flatname($root->scientific_name)." (".$root->id.")".'";'
	);
}
