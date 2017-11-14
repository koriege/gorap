package Bio::Gorap::DB::Taxonomy;

use Moose;
use File::Spec::Functions;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use Bio::TreeIO;
use Bio::DB::EUtilities;
use List::Util qw(min max);
use Try::Tiny;

#uses the gorap parameter object to initialize the
#database of a Bio::DB::Taxonomy object
has 'parameter' => (
	is => 'ro',
	isa => 'Bio::Gorap::Parameter',
	required => 1,
	trigger => \&_set_db
);

#NCBI Taxonomy based Bio::DB::Taxonomy databases
has 'ncbi' => (
	is => 'rw',
	isa => 'Bio::DB::Taxonomy',
);

has 'silva' => (
	is => 'rw',
	isa => 'Bio::Tree::Tree',
	lazy => 1,
	default => sub { (Bio::TreeIO->new(-format => 'newick', -file => catfile($ENV{GORAP},'gorap','data','taxonomy','silva.newick') , -verbose => -1))->next_tree }
);

has 'accToTaxid' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
);

has 'accToName' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
);

has 'nameToTaxid' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
);

has 'taxidToName' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
);

has 'speciesID' => (
	is => 'rw',
	isa => 'Int'
);

has 'rankID' => (
	is => 'rw',
	isa => 'Int'
);

has 'rankIDlineage' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub {[]}
);

has 'taxIDsToLineage' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
);

has 'relatedRankIDsToLineage' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
);

has 'relatedSpeciesIDsToLineage' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
);

#needs gorap data
sub _set_db {
	my ($self) = @_;

	#initialize taxonomy database
	#directory was parameter-tmp -> unique id, but now we will reuse created indexes (force 0)
	print "No index found. Indexing Taxonomy DB - may take a while\n" unless (-e catfile($ENV{GORAP},'gorap','data','taxonomy','parents'));
	$self->ncbi(Bio::DB::Taxonomy->new(-source => 'flatfile', -nodesfile => catfile($ENV{GORAP},'gorap','data','taxonomy','nodes.dmp'), -namesfile => catfile($ENV{GORAP},'gorap','data','taxonomy','names.dmp'), -directory => catdir($ENV{GORAP},'gorap','data','taxonomy') , -force => 0 , -verbose => -1 ));

	#read in silva phylogeny accession numbers, already mapped to ncbi taxonomy
	#in case of coming from previous gorap version or using the ambiguous silva updater
	if (-e catfile($ENV{GORAP},'gorap','data','silvaNcbi.txt')){
		open MAP, '<'.catfile($ENV{GORAP},'gorap','data','silvaNcbi.txt') or die $!;
		while(<MAP>){
			chomp $_;
			my ($acc,$tax) = split(/\s+/,$_);
			$self->accToTaxid->{$acc} = $tax;
		}
		close MAP;
	}

	#read in rfam and silva accession numbers, already mapped to ncbi taxonomy
	if (-e catfile($ENV{GORAP},'gorap','data','accSciTax.txt')){
		open MAP, '<'.catfile($ENV{GORAP},'gorap','data','accSciTax.txt') or die $!;
		while(<MAP>){
			chomp $_;
			my @tmp=split(/\s+/,$_);
			$self->accToName->{$tmp[0]}=$tmp[1];
			$self->accToTaxid->{$tmp[0]}=$tmp[2];
			$self->nameToTaxid->{$tmp[1]}=$tmp[2];
			$self->taxidToName->{$tmp[2]}=$tmp[1];
		}
		close MAP;
	}
}

sub reindex {
	my ($self) = @_;
	print "Indexing Taxonomy DB - may take a while\n";
	$self->ncbi(Bio::DB::Taxonomy->new(-source => 'flatfile', -nodesfile => catfile($ENV{GORAP},'gorap','data','taxonomy','nodes.dmp'), -namesfile => catfile($ENV{GORAP},'gorap','data','taxonomy','names.dmp'), -directory => catdir($ENV{GORAP},'gorap','data','taxonomy') , -force => 1 , -verbose => -1 ));
}

sub flatname {
	my ($self, $name) = @_;

	$name=~s/(^\s+|\s+$)//g;
	$name=~s/\s+/_/g;
	$name=~s/\.\.+/\./g;
	$name=~s/__+/_/g;
	$name=~s/[^a-zA-Z0-9_]*//g;
	$name=ucfirst($name);

	return $name;
}

sub flatquery {
	my ($self, $name) = @_;

	$name=lc($name);
	$name=~s/(^\s+|\s+$)//g;
	$name=~s/__*/ /g;
	$name=~s/\.\.*/\./g;

	return $name;
}

#from taxid or Bio::Taxon
sub getLineageNodes { #does not include node of it self
	my ($self , $taxid) = @_;

	my $taxon;
	if (ref(\$taxid) eq 'SCALAR'){
		$taxon = $self->ncbi->get_taxon(-taxonid => $taxid);
	} else {
		$taxon = $taxid;
		$taxid = $taxon->id;
	}
	# check for reuse of prior crated valid lineages
	return $self->taxIDsToLineage->{$taxid} if exists $self->taxIDsToLineage->{$taxid};

	my @nodes = ();
	if ($taxon){
		my $tree_functions = Bio::Tree::Tree->new( -verbose => -1);
		@nodes = $tree_functions->get_lineage_nodes($taxon);
		# push @nodes, $taxid;
	} else { # try to reach web api
		my ($error, $try, $tries) = (1, 0, 5);
		while($error && (++$try) < $tries){
			$error = 0;
			try{
				my $factory = Bio::DB::EUtilities->new(
					-eutil => 'efetch',
					-email => 'mymail@foo.bar',
					-db => 'taxonomy',
					-id => $taxid,
					-verbose => -1
				);
				my $res = $factory->get_Response->content;
				# parse xml return into hash
				my $data = XMLin($res);
				if (ref $data){
					for (@{$data->{Taxon}->{LineageEx}->{Taxon}}){
						push @nodes, $_->{TaxId};
						my $name = $self->flatname($_->{ScientificName});
						my $id = $_->{TaxId};
						$self->taxidToName->{$id} = $name;
						$self->nameToTaxid->{$name} = $id;
					}
				}
			} catch {
				$error = 1;
			};
		}
		# push @nodes, $taxid;
	}
	push @{$self->taxIDsToLineage->{$taxid}}, @nodes if $#nodes > -1;

	return \@nodes;
}

sub getIDfromName {
	my ($self,$query) = @_;

	return 0 unless $query;
	return $self->nameToTaxid->{$query} if exists $self->nameToTaxid->{$query};
	$query = $self->flatquery($query);

	my $taxon;
	if ($query=~/^\d+$/){
		$taxon = $self->ncbi->get_taxon(-taxonid => $query);
	} else {
		$taxon = $self->ncbi->get_taxon(-name => $query);
	}

	my $taxid;
	if ($taxon){
		$taxid = $taxon->id;
		my $name = $self->flatname($taxon->scientific_name);
		$self->taxidToName->{$taxid} = $name;
		$self->nameToTaxid->{$name} = $taxid;
	} else {
		my ($try, $error, $tries) = (0, 1, 5);
		while ($error && (++$try) < $tries){
			$error = 0;
			try { # try to reach web api
				my $factory = Bio::DB::EUtilities->new(
					-eutil => 'esearch',
					-db => 'taxonomy',
					-email => 'mymail@foo.bar',
					-term => $query,
					-verbose => -1
				);
				($taxid) = $factory->get_ids;
			} catch {
				$error = 1;
			};
		}
	}

	if ($taxid){
		$self->getNameFromID($taxid);
		return $taxid;
	} else {
		return 0
	}
}

sub getNameFromID {
	my ($self,$taxid) = @_;

	return $self->taxidToName->{$taxid} if exists $self->taxidToName->{$taxid};

	my $name;
	my $taxon = $self->ncbi->get_taxon(-taxonid => $taxid);
	if ($taxon) {
		$name = $taxon->scientific_name;
	} else {
		my ($try, $error, $tries) = (0, 1, 5);
		while($error && (++$try) < $tries){
			$error = 0;
			try { # try to reach web api
				my $factory = Bio::DB::EUtilities->new(
					-eutil => 'esummary',
					-email => 'mymail@foo.bar',
					-db => 'taxonomy',
					-id => $taxid,
					-verbose => -1
				);
				($name) = $factory->next_DocSum->get_contents_by_name('ScientificName');
			} catch {
				$error = 1;
			};
		}
	}

	if ($name){
		$name = $self->flatname($name);
		$self->taxidToName->{$taxid} = $name;
		$self->nameToTaxid->{$name} = $taxid;
		return $name;
	} else {
		return 0
	}
}

sub getIDfromAccession {
	my ($self, $acc) = @_;

	return $self->accToTaxid->{$acc} if exists $self->accToTaxid->{$acc};

	my ($taxid, $name);
	my ($try, $error, $tries) = (0, 1, 5);
	while ($error && (++$try) < $tries){
		$error = 0;
		try {
			my $factory = Bio::DB::EUtilities->new(
				-eutil => 'esearch',
				-email => 'mymail@foo.bar',
				-db => 'nuccore',
				-term => $acc,
				-verbose => -1
			);
			my @uids = $factory->get_ids;
			$factory->reset_parameters(
				-eutil => 'esummary',
				-email => 'mymail@foo.bar',
				-db => 'nuccore',
				-id => \@uids
			);
			($taxid) = $factory->next_DocSum->get_contents_by_name('TaxId');
		} catch {
			$error = 1;
		};
	}

	if ($taxid){
		$self->accToTaxid->{$acc} = $taxid;
		my $name = $self->getNameFromID($taxid);
		$self->accToName->{$acc} = $name if $name;
		return $taxid;
	} else {
		return 0
	}
}

sub findRelatedSpecies {
	my ($self) = @_;

	print "Searching taxonomy for related species\n" if $self->parameter->has_species || $self->parameter->has_rank && $self->parameter->verbose;

	if ($self->parameter->has_species){
		$self->speciesID($self->getIDfromName($self->parameter->species));
		if ($self->speciesID){
			$self->relatedSpeciesIDsToLineage($self->findRelatedIDs($self->speciesID));
		} else {
			print ":WARNING: Given species not found" if $self->parameter->verbose;
		}
	}
	if ($self->parameter->has_rank){
		$self->rankID($self->getIDfromName($self->parameter->rank));
		if ($self->rankID){
			push @{$self->rankIDlineage} , $_->id for @{$self->getLineageNodes($self->rankID)};
			push @{$self->rankIDlineage} , $self->rankID;
			$self->relatedRankIDsToLineage($self->findRelatedIDs($self->rankID));
		} else {
			print ":WARNING: Given rank not found" if $self->parameter->verbose;
		}
	}


}

sub findRelatedIDs {
	my ($self, $taxid) = @_;

	my $taxonorig = $self->ncbi->get_taxon(-taxonid => $taxid);

	my $taxon = $taxonorig;
	my $ancestor = $self->ncbi->ancestor($taxon);
	while( scalar $self->ncbi->each_Descendent($ancestor) == 1){
		$taxon = $ancestor;
		$ancestor = $self->ncbi->ancestor($ancestor);
	}

	my $relatedIDsToLineage;
	for my $t ($taxon,$self->ncbi->get_all_Descendents($taxon)){
		my $nodes = $self->getLineageNodes($t);
		$_ = $_->id for @{$nodes};
		$relatedIDsToLineage->{$t->id} = $nodes;
	}
	# print scalar keys %{$relatedIDsToLineage}; print "\n";
	# print join(":",(@{$relatedIDsToLineage->{$_}},$_))."\n" for keys %{$relatedIDsToLineage};
	# exit;

	### silva based refinement
	my $silvaIDsToNodes; # faster access than $self->silva->find_node(-id => $id);
	for ($self->silva->get_nodes){
		$silvaIDsToNodes->{$_->id} = $_ if $_->id;
	}

	#add missing
	my @q;
	push @q, $taxon;
	while($#q > -1){ #BFT: find inner node in silva, iterate over silva subtree to add missing species
		my $v = shift @q;
		my $s = $silvaIDsToNodes->{$v->id};
		if ($s) {
			for (($s,$s->get_all_Descendents())){
				next unless $_->id;
				unless (exists $relatedIDsToLineage->{$_->id}){
					my $t = $self->ncbi->get_taxon(-taxonid => $_->id);
					my $nodes = $self->getLineageNodes($t);
					$_ = $_->id for @{$nodes};
					$relatedIDsToLineage->{$t->id} = $nodes;
				}
			}
		} else {
			push @q, $self->ncbi->each_Descendent($v);
		}
	}

	#remove misplaced
	my $v = $taxonorig;
	my $found = 0;
	while($v){
		if (exists $silvaIDsToNodes->{$v->id}){ #find ancestor in silva, remove species if not in silva subtree but in whole tree
			$found = 1;
			my $nodes;
			for ($silvaIDsToNodes->{$v->id}->get_all_Descendents()){
				next unless $_->id;
				$nodes->{$_->id}=1
			}
			for (keys %{$relatedIDsToLineage}){
				if (exists $silvaIDsToNodes->{$_} && ! exists $nodes->{$_}){
					delete $relatedIDsToLineage->{$_};
				}
			}
			last;
		}
		$v = $self->ncbi->ancestor($v);
	}

	unless ($found){ #if no ancestor found in silva: compute lca from species in silva, remove species if not in silva subtree but in whole tree
		my @nodes;
		for (keys %{$relatedIDsToLineage}){
			push @nodes, $silvaIDsToNodes->{$_} if exists $silvaIDsToNodes->{$_};
		}
		if ($#nodes != -1){
			my $lca = $self->silva->get_lca(-nodes => \@nodes);
			my $tmp;
			for ($lca->get_all_Descendents()){
				next unless $_->id;
				$tmp->{$_->id}=1
			}
			for (keys %{$relatedIDsToLineage}){
				if (exists $silvaIDsToNodes->{$_} && ! exists $tmp->{$_}){
					delete $relatedIDsToLineage->{$_};
				}
			}
		}
	}

	return $relatedIDsToLineage;
}


sub findRelatedIDs_old {
	my ($self,$taxid) = @_;

	my $relatedIDsToLineage;
	#get query taxon from ncbi
	my $taxon = $self->ncbi->get_taxon(-taxonid => $taxid);
	my $ancestor = $self->ncbi->ancestor($taxon);
	#update query taxon to its parent if its the only child
	while( scalar $self->ncbi->each_Descendent($ancestor) == 1){
		$taxon = $ancestor;
		$ancestor = $self->ncbi->ancestor($ancestor);
	}

	#store all children of the query taxon as related taxons with its lineages
	for my $t ($taxon,$self->ncbi->get_all_Descendents($taxon)){
		my $nodes = $self->getLineageNodes($t);
		$_ = $_->id for @{$nodes};
		$relatedIDsToLineage->{$t->id} = $nodes;
	}

	#if query taxon belongs to bacteria, go on refining the related list by silva tree
	return $relatedIDsToLineage unless ${$relatedIDsToLineage->{$taxon->id}}[1] == 2;

	#store silva tree nodes in hashmap for faster access than depth first search
	my $silvaIDsToNodes;
	do {$silvaIDsToNodes->{$_->id} = $_ if $_->id} for $self->silva->get_nodes;

	#find taxon or its offsprings (#childs=1) in silva tree
	my ($node) = $self->silva->find_node(-id => $taxon->id);
	my $tmpTaxon = $taxon;
	while ( !$node && scalar $self->ncbi->each_Descendent($tmpTaxon) == 1){
		($tmpTaxon) = $self->ncbi->each_Descendent($tmpTaxon);
		$node = $silvaIDsToNodes->{$tmpTaxon->id};
		#($node) = $self->silva->find_node(-id => $tmpTaxon->id);
	}

	my @heightSortRelatedIDs = sort {$#{$relatedIDsToLineage->{$a}} <=> $#{$relatedIDsToLineage->{$b}}} keys %{$relatedIDsToLineage};

	#find all ncbi related ids in silva tree
	my $relatedIDinSilva;
	my $c=0;
	for (@heightSortRelatedIDs){
		next if exists $relatedIDinSilva->{$_};
		my $tmpNode = $silvaIDsToNodes->{$_};
		#my ($tmpNode) = $self->silva->find_node(-id => $_);
		if ($tmpNode){
			$relatedIDinSilva->{$_} = 1;
			do { $relatedIDinSilva->{$_->id} = 1 if $_->id && $_->id=~/^\d+$/} for $tmpNode->get_all_Descendents;
		}
	}
	#if taxon was not found in silva tree compute lca of all ncbi ids found in silva
	my $lca = 0;
	unless ($node){
		if (scalar keys %{$relatedIDinSilva} > 1){
			$lca = 1;
			my @nodes;
			push @nodes , $silvaIDsToNodes->{$_} for keys %{$relatedIDinSilva};
			$node = $self->silva->get_lca(-nodes => \@nodes);
		} elsif (scalar keys %{$relatedIDinSilva} == 1){
			$node = $silvaIDsToNodes->{(keys %{$relatedIDinSilva})[0]};
		}
	}

	return $relatedIDsToLineage unless $node;

	#destinguish between ncbi ids supported by silva tree as offsprings and those not
	my $silvaRelatedIDs;
	$silvaRelatedIDs->{$node->id} = 1 if $node->id && $node->id=~/^\d+$/;
	do { $silvaRelatedIDs->{$_->id} = 1 if $_->id && $_->id=~/^\d+$/ } for $node->get_all_Descendents;

	#remove all ids from final relatedIDsToLineage not supported as offsprings by silva tree
	for (@heightSortRelatedIDs)	{
		next unless exists $relatedIDsToLineage->{$_}; # already removed
		#if in $relatedIDinSilva but not in $silvaRelatedIDs -> remove itself and all children from $relatedIDsToLineage
		if (exists $relatedIDinSilva->{$_} && ! exists $silvaRelatedIDs->{$_}){
			delete $relatedIDsToLineage->{$_};
			#remove also all childs which may be present in $relatedIDsToLineage
			delete $relatedIDsToLineage->{$_->id} for $self->ncbi->get_all_Descendents($self->ncbi->get_taxon(-taxonid => $_));
		}
	}

	#add all missing ids to final relatedIDsToLineage supported as offsprings by silva tree, if rank id was found
	#in silva, otherwise there will be added too many false positives
	unless ($lca){
		for (keys %{$silvaRelatedIDs}){
			next if exists $relatedIDsToLineage->{$_};
			my $taxon = $self->ncbi->get_taxon(-taxonid => $_);
			next unless $taxon;
			for my $t ($taxon , $self->ncbi->get_all_Descendents($taxon)){
				next if exists $relatedIDsToLineage->{$t->id};
				my $nodes = $self->getLineageNodes($t);
				$_ = $_->id for @{$nodes};
				$relatedIDsToLineage->{$t->id} = $nodes;
			}
		}
	}

	return $relatedIDsToLineage;
}

sub sort_stk {
	my ($self, $stk) = @_;

	# return $stk unless $self->rankID || $self->speciesID;

	my @tosort;
	my @remaining;
	# my $taxID = $self->speciesID ? $self->speciesID : $self->rankID;
	for ( $stk->each_seq() ) {
		if (exists $self->nameToTaxid->{(split(/\./,$_->id))[0]}){
			push @tosort , $_;
		} else {
			push @remaining , $_;
		}
		$stk->remove_seq($_);
	}
	try{
		my @sort = sort { $self->lineageNodesToString($self->getLineageNodes($self->nameToTaxid->{(split(/\./,$a->id))[0]}))
			cmp $self->lineageNodesToString($self->getLineageNodes($self->nameToTaxid->{(split(/\./,$b->id))[0]}))
		} @tosort;
		$stk->add_seq($_) for @sort;
	} catch {
		$stk->add_seq($_) for @tosort;
	};
	$stk->add_seq($_) for @remaining;

	$stk->set_displayname_flat;
	return $stk;
}

sub lineageNodesToString {
	my ($self, $nodes) = @_;

	my @s;
	push @s, $_->id for @{$nodes};

	return join(',',@s);
}

1;