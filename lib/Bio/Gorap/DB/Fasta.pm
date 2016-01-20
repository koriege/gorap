package Bio::Gorap::DB::Fasta;

use Moose;
use Bio::Index::Fasta;
use File::Spec::Functions;
use Encode;

#uses the gorap parameter object to initialize the 
#database of a Bio::Index::Fasta object
has 'parameter' => (
	is => 'ro',
	isa => 'Bio::Gorap::Parameter',
	required => 1 ,
	trigger => \&_set_db 
);

#genome file based hashmap of Bio::Index::Fasta databases
has 'db' => (
	is => 'rw',
    isa => 'Bio::Index::Fasta'
);

has 'do_chunks' => (
	is => 'ro',
	isa => 'Bool',
	default => 1
);

has 'chunks' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub {[]}
);

has 'oheaderToDBsize' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub {[]}
);

#for fast access during ToolI deletions
has 'oheaders' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
);

has 'nheaders' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub {[]}
);

my $set_db_abbr;
my @oheader;
my $doubleEntry;
sub _set_db {
	my ($self) = @_;	

	$self->db(Bio::Index::Fasta->new(-filename => catfile($self->parameter->tmp,$self->parameter->pid.'.faidx'), -write_flag => 1 , -verbose => -1));				
	
	print "Reading FASTA files\n" if $self->parameter->verbose;
	
	for my $i (0..$#{$self->parameter->genomes}){
		my $genome = ${$self->parameter->genomes}[$i];
		print $genome."\n" if $self->parameter->verbose;
		$set_db_abbr = ${$self->parameter->abbreviations}[$i];
		#make headers uniq by adding an abbreviation/filname in front of \S+
		&add_fasta($self,$genome,\&_parse_id); 
		my $residues=0;		
		for (@oheader){	
			$self->oheaders->{$_}=1;	
			push @{$self->nheaders} , $set_db_abbr.'.'.$_;
			$residues += ($self->db->fetch($set_db_abbr.'.'.$_))->length;			
		}		
		push @{$self->oheaderToDBsize} , [$oheader[0] , $residues];				
		@oheader=();
	}	
	$self->chunks(&chunk($self,$self->nheaders,$self->parameter->threads)) if $self->do_chunks && ! $self->parameter->skip_comp;	
}

sub add_fasta {
	my ($self,$file,$idparser) = @_;
	
	$self->db->id_parser(\&{$idparser}) if $idparser;
	$self->db->make_index($file);	
}

#gorap specific id parser
sub _parse_id {
	my ($header) = @_;	

	$header=~/^>\s*(\S+)/;
	push @oheader , $1;
	return $set_db_abbr.'.'.$1;	
}


#gorap specific subseq extraction
sub get_gff3seq {
	my ($self,$s) = @_;

	my ($id,$source,$type,$start,$stop,$score,$strand,$phase,$attributes) = @{$s};
	
	my ($abbr,@orig) = split /\./ , $id;
	pop @orig;

	my $seq=$self->db->fetch(join('.',($abbr,@orig)));	
	if ($strand eq "+"){		
		my $s = $seq->subseq($start,$stop);
		$s =~s/[tT]/U/g;
		return $s;	
	} else {	
		my $s = ((Bio::Seq->new( -seq => $seq->subseq($start,$stop) , -verbose => -1))->revcom)->seq;
		$s =~s/[tT]/U/g;
		return $s;											
	}	
}

#subsequence extraction from this db
sub get_subseq {
	my ($self,$id,$start,$stop,$strand,$dna) = @_;

	my $seq=$self->db->fetch($id);	
	if ($strand eq "+"){		
		my $s = $seq->subseq($start,$stop);
		$s =~s/[tT]/U/g unless $dna;
		return $s;	
	} else {	
		my $s = ((Bio::Seq->new( -seq => $seq->subseq($start,$stop) , -verbose => -1))->revcom)->seq;
		$s =~s/[tT]/U/g unless $dna;
		return $s;						
	}	
}

sub chunk {
	my ($self) = @_;

	my $residues=0;
	for (@{$self->nheaders}){		
		$residues += ($self->db->fetch($_))->length;
	}
	$residues /= $self->parameter->threads;	
	
	my @chunks;
	my $c=0;
	my $outfile = catfile($self->parameter->tmp,$self->parameter->pid.'.'.$c);
	push @chunks, $outfile;
	open O , '>'.$outfile or die $!;
	my $written=0;
	my $pos=0;
	for my $h (@{$self->nheaders}){
		$pos=0;
		print O '>'.$pos.'.'.$h."\n";
		my $seqo = $self->db->fetch($h);		
		my @seqparts = unpack("(A80)*", $seqo->seq);
		for my $i (0..$#seqparts){
			my $seq=$seqparts[$i];
			print O $seq."\n";
			$written+=length($seq);
			$pos+=length($seq);
			if ($written > $residues){
				close O;
				$outfile = catfile($self->parameter->tmp,$self->parameter->pid.'.'.(++$c));
				push @chunks, $outfile;
				open O , '>'.$outfile or die $!;
				
				my @overlap;
				my $owritten=0;
				for (my $j=$i; $j>=0; $j--){
					unshift @overlap , $seqparts[$j];
					$pos-=length($seqparts[$j]);
					$owritten+=length($seqparts[$j]);
					last if $owritten > 1000;
				}

				print O '>'.$pos.'.'.$h."\n";				
				print O $_."\n" for @overlap;				
				$written=0;
			}
		}
	}
	close O;	
	
	return	\@chunks;
}

sub chunk_backmap {
	my ($self, $id, $start, $stop) = @_;
	
	my ($pos,$abbr,@orig) = split /\./ , $id;
	return (join(".",($abbr,@orig)),$start+$pos,$stop+$pos);
}

1;
