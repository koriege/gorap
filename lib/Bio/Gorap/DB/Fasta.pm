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

has 'headerToPos' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
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

has 'oheaders' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
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
			$residues += ($self->db->fetch($set_db_abbr.'.'.$_))->length;	
			$self->oheaders->{$_} = 1;
		}
		
		push @{$self->oheaderToDBsize} , [$oheader[0] , $residues];		
		#transform into cpucount many chunks 		
		do {push @{$self->chunks} , @{&chunk($self, $genome, $self->parameter->threads, 'g'.$i, \&_parse_id)} } if $self->do_chunks && $#{$self->parameter->queries} > -1;
		@oheader=();
	}
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
	my @id = split /\./ , $id;
	my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);

	my $seq=$self->db->fetch($abbr.'.'.$orig);	
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
	my ($self, $genome, $parts, $abbr, $idparser) = @_;

	my @files;
	my $genomesize = -s $genome;
	my $chunksize = $genomesize / $parts;
	my $outfile = catfile($self->parameter->tmp,$self->parameter->pid.'.'.$abbr.'c');
	my $oc = 0;
	my $ostop = 0;
	my $oheader;
	my $stop = 0;
	my $writtenbytes=0;
	
	open I , '<'.$genome or die $!;
	open O , '>'.$outfile.$oc.'.fa' or die $!;
	push @files , $outfile.$oc.'.fa';
	print O '>'.$abbr.'c'.$oc."\n";
	while(<I>){		
		if ( $writtenbytes >= $chunksize){
			#save position and start at 0
			close O;
			$writtenbytes=0;
			${$self->headerToPos->{$abbr.'c'.$oc}}[-1]->{'ostop'}=$ostop;
			${$self->headerToPos->{$abbr.'c'.$oc}}[-1]->{'stop'}=$stop;
			$stop = 0;
			$oc++;	
			push @{$self->headerToPos->{$abbr.'c'.$oc}} , { oid => $oheader, stop => $stop, ostop => $ostop};		
			open O , '>'.$outfile.$oc.'.fa' or die $!;			
			push @files , $outfile.$oc.'.fa';
			print O '>'.$abbr.'c'.$oc."\n";						
		}
		if ($_ =~ /^>/){			
			if ($#{$self->headerToPos->{$abbr.'c'.$oc}} > -1){
				${$self->headerToPos->{$abbr.'c'.$oc}}[-1]->{'ostop'}=$ostop;
				${$self->headerToPos->{$abbr.'c'.$oc}}[-1]->{'stop'}=$stop if $stop > 0;
			}
			$_=~ /^>\s*(\S+)/;
			$oheader = $idparser ? ${\&{$idparser}($_)} : $1;	
			$ostop = 0;
			push @{$self->headerToPos->{$abbr.'c'.$oc}} , { oid => $oheader, stop => $stop, ostop => $ostop};
		} else {
			chomp $_;
			$ostop += length $_;
			$stop += length $_;
			print O $_."\n";			
			$writtenbytes += length(Encode::encode_utf8($_."\n"));
		}
	}	
	close I;
	close O;
	if ($#{$self->headerToPos->{$abbr.'c'.$oc}} > -1){
		${$self->headerToPos->{$abbr.'c'.$oc}}[-1]->{'ostop'}=$ostop;
		${$self->headerToPos->{$abbr.'c'.$oc}}[-1]->{'stop'}=$stop if $stop > 0;
	}

	#overlaps	
	for $oc (0..$#files - 1){
		#print $abbr.'c'.$oc."\n";
		if(${$self->headerToPos->{$abbr.'c'.$oc}}[-1]->{'oid'} eq ${$self->headerToPos->{$abbr.'c'.($oc+1)}}[0]->{'oid'}){
			open O , '>>'.$outfile.$oc.'.fa' or die $!;
			open I , '<'.$outfile.($oc+1).'.fa' or die $!;
			my $stop = 0;
			while(<I>){
				next if $_ =~ /^>/;
				chomp $_;				
				if ($stop < 1000 && $stop < ${$self->headerToPos->{$abbr.'c'.($oc+1)}}[0]->{'stop'}){
					$stop += length $_;
					print O $_."\n";
				} else {
					last;
				}
			}
			close I;
			close O;
			${$self->headerToPos->{$abbr.'c'.$oc}}[-1]->{'stop'}+=$stop;
			${$self->headerToPos->{$abbr.'c'.$oc}}[-1]->{'ostop'}+=$stop;
		}
		#print ${$self->headerToPos->{$abbr.'c'.$oc}}[$_]->{'oid'}.' '.${$self->headerToPos->{$abbr.'c'.$oc}}[$_]->{'ostop'}.' '.${$self->headerToPos->{$abbr.'c'.$oc}}[$_]->{'stop'}."\n" for 0..$#{$self->headerToPos->{$abbr.'c'.$oc}};
	}
	#print $abbr.'c'.$oc."\n";
	#print ${$self->headerToPos->{$abbr.'c'.$oc}}[$_]->{'oid'}.' '.${$self->headerToPos->{$abbr.'c'.$oc}}[$_]->{'ostop'}.' '.${$self->headerToPos->{$abbr.'c'.$oc}}[$_]->{'stop'}."\n" for 0..$#{$self->headerToPos->{$abbr.'c'.$oc}};

	return \@files;
}

sub chunk_backmap {
	my ($self, $id, $start, $stop) = @_;

	my @e;	
	for (0..$#{$self->headerToPos->{$id}}){
		my $e = ${$self->headerToPos->{$id}}[$_];
		next if $start > $e->{'stop'};
		
		while ($stop > $e->{'stop'}){						
			push @e , [$e, $start, $e->{'stop'}];
			$start = $e->{'stop'} + 1;
			$e = ${$self->headerToPos->{$id}}[++$_];
		}		
		push @e , [$e,$start,$stop];

		@e = reverse sort {$$a[2]-$$a[1] <=> $$b[2]-$$b[1]} @e;
		
		my @ret;
		for (@e){
			# last if $$_[2]-$$_[1] < ($e[0][2] - $e[0][1]) * 0.9;
			my $ostop = $$_[0]->{'ostop'}-($$_[0]->{'stop'} - $$_[2]);
			my $ostart = $ostop - ($$_[2] - $$_[1]);
			push @ret , [$$_[0]->{'oid'}, $ostart, $ostop];
			# print $$_[0]->{'oid'}.' '.$ostart.' '.$ostop."\n";
		}

		return @ret;
	}
}

sub chunk_backmap_old {
	my ($self, $id, $start, $stop) = @_;

	for (0..$#{$self->headerToPos->{$id}}){
		my $e = ${$self->headerToPos->{$id}}[$_];
		next if $stop > $e->{'stop'};	

		#return bigger part if something overlaps
		if ($_ > 0 && $start <= ${$self->headerToPos->{$id}}[$_-1]->{'stop'}){			
			my $l1 = ${$self->headerToPos->{$id}}[$_-1]->{'stop'} - $start;
			my $l2 = $stop - ${$self->headerToPos->{$id}}[$_-1]->{'stop'};
			if ($l1 < $l2){
				$start = ${$self->headerToPos->{$id}}[$_-1]->{'stop'} + 1;
			} else {
				$e = ${$self->headerToPos->{$id}}[$_-1];
				$stop = $e->{'stop'};
			}
		}
		my $ostop = $e->{'ostop'}-($e->{'stop'}-$stop);
		my $ostart = $ostop - ($stop - $start);
		return ($e->{'oid'} , $ostart , $ostop );
	}
}
1;
