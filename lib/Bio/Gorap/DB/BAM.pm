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
	isa => 'HashRef',
	default => sub { {} }
);

has 'sizes' => (
	is => 'rw',
    isa => 'HashRef',
    default => sub { {} }
);

sub _set_db {
	my ($self) = @_;
	return unless $self->parameter->has_bams;
	print "Reading BAM files\n" if $self->parameter->verbose;

	my $sizes;
	for (0..$#{$self->parameter->bams}){
		my $abbr = ${$self->parameter->abbreviations}[$_];
		$sizes->{$abbr}=0;
		for my $f (@{${$self->parameter->bams}[$_]}){
			my $c = 0;
			my $map = {};
			my $next;
			for (`samtools idxstats $f`){
				next if $_=~/^\*/;
				$next = 1 if $_=~/fail/;
				last if $next;
				my ($header , $mapped, $unmapped) = split /\s+/ , $_;
				$c += $mapped + $unmapped;
				$map->{$header} = $mapped + $unmapped;
			}
			next if $next;
			push @{$self->db->{$abbr}} , Bio::DB::Sam->new(-bam => $f, -autindex => 1, -verbose => -1);
			# my ($foo) = ${$self->db->{$abbr}}[-1]->features();

			# exit;
			$sizes->{$abbr} += $c;
		}
	}
	$self->sizes($sizes);	
}

sub rpkm {
	#TODO strand specific
	my ($self,$abbr,$id,$start,$stop,$strand) = @_;

	my $count = 0;
	my $ex;
	for (0..$#{$self->db->{$abbr}}){		
		my $bam = ${$self->db->{$abbr}}[$_];		
		$ex=1;
		$count += scalar $bam->get_features_by_location(-type => 'read_pair', -seq_id => $id, -start => $start, -end => $stop);
	}
	
	return ('.','.') unless $ex;
	
	my $libsize=$self->sizes->{$abbr};
	my $rpkm = $libsize ? ($count / ($libsize/10^6)) / (($stop-$start)/10^3) : 0;
	my $tpm = $libsize ? ($count / (($stop-$start)/10^3)) / ($libsize/10^6) : 0;

	return ($tpm == 0 ? '.' : sprintf("%.10f",$tpm), $rpkm == 0 ? '.' : sprintf("%.10f",$rpkm), $count == 0 ? '.' : $count);
}

1;