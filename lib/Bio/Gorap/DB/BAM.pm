package Bio::Gorap::DB::BAM;

use Moose;
use Bio::DB::Sam;
use List::Util qw(max);

has 'parameter' => (
	is => 'ro',
	isa => 'Bio::Gorap::Parameter',
	required => 1 ,
	trigger => \&_set_db 
);

has 'db' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub { [] }
);

has 'ids' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub { [] }
);

has 'sizes' => (
	is => 'rw',
    isa => 'ArrayRef',
    default => sub { [] }
);


sub _set_db {
	my ($self) = @_;	

	return unless $self->parameter->has_bams;
	print "Reading BAM files\n" if $self->parameter->verbose;
	for (@{$self->parameter->bams}){
		
		# push @{$self->ids} , {map { $_ => 1 } ${$self->db}[-1]->seq_ids};
		my $c = 0;
		my $map = {};
		my $next;
		for (`samtools idxstats $_`){			
			next if $_=~/^\*/;		
			$next = 1 if $_=~/fail/;
			last if $next;
			my ($header , $mapped, $unmapped) = split /\s+/ , $_;
			$c += $mapped + $unmapped;
			$map->{$header} = $mapped + $unmapped;
		}
		next if $next;		
		push @{$self->db} , Bio::DB::Sam->new(-bam => $_, -autoindex => 1 , -verbose => -1);
		push @{$self->sizes} , $c;
		push @{$self->ids} , $map;
	}	
}

sub rpkm(){
	my ($self,$id,$start,$stop) = @_;

	my $count = 0;
	# my $max = 0;
	my $libsize=0;
	my $ex=0;
	for (0..$#{$self->db}){
		next unless exists ${$self->ids}[$_]->{$id};
		my $bam = ${$self->db}[$_];		
		$ex=1;
		$libsize += ${$self->sizes}[$_];
		$count += scalar $bam->get_features_by_location(-seq_id => $id, -start => $start, -end => $stop);
		# my ($coverage) = $bam->features(-type=>'coverage', -seq_id => $id, -start => $start, -stop => $stop);
		# $max = max($max,@{$coverage->coverage()});		
	}
	unless($ex){
		return ('.','.');
	}	
	my $rpkm = $libsize ? ($count / ($libsize/10^6)) / (($stop-$start)/10^3) : 0;
	my $tpm = $libsize ? ($count / (($stop-$start)/10^3)) / ($libsize/10^6) : 0;

	return ($tpm == 0 ? 0 : sprintf("%.10f",$tpm), $rpkm == 0 ? 0 : sprintf("%.10f",$rpkm),$count);
}

1;