package Bio::Gorap::DB::BAM;

use Moose;
use Bio::DB::Sam;
use List::Util qw(max);
use Try::Tiny;

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
			for (`samtools idxstats $f`){ #todo counts all fragments
				next if $_=~/^\*/;
				$next = 1 if $_=~/fail/;
				last if $next;
				my ($header , $size, $mapped, $unmapped) = split /\s+/ , $_;
				$c += $mapped;
				$map->{$header} = $mapped + $unmapped;
			}
			next if $next;
			push @{$self->db->{$abbr}} , Bio::DB::Sam->new(-bam => $f, -autindex => 1, -verbose => -1);
			$sizes->{$abbr} += $c;
		}
	}
	$self->sizes($sizes);
}

sub rpkm {	
	my ($self,$abbr,$id,$start,$stop,$strand) = @_;
	$strand = $strand=~/(\.|0)/ ? 0 : $strand=~/(\+|1)/ ? 1 : -1;

	my $count = 0;
	for (0..$#{$self->db->{$abbr}}){		
		my $bam = ${$self->db->{$abbr}}[$_];		

		my $readstrand = '+';
		try {
			if ($self->parameter->strandspec && $strand != 0){
				for ($bam->get_features_by_location(-type => 'read_pair', -seq_id => $id, -start => $start, -end => $stop)){					
					my @segments = $_->segments;
					my $rstrand = 1;
					my $unmapped = 0;
					for (@segments){
						my @flags = split /\|/ , $_->get_tag_values('FLAGS');
						for (@flags){
							$unmapped++ if $_ eq 'UNMAPPED'; #counts read also if only one mate maps
							$rstrand = -1 if $_ eq 'REVERSED'; #to exclude M_REVERSED
							# $mates++ if $_=~/FIRST/ || $_=~/SECOND/;
							# $matemapped = 0 if $_=~/M_UNMAPPED/;
							# $paired = 1 if $_=~/PAIRED/;
						}
					}
					$count++ if $unmapped < ($#segments +1) && $strand == $rstrand;
				}
			} else {
				for ($bam->get_features_by_location(-type => 'read_pair', -seq_id => $id, -start => $start, -end => $stop)){
					my @segments = $_->segments;
					my $unmapped = 0;
					for (@segments){
						for (split /\|/ , $_->get_tag_values('FLAGS')){
							$unmapped++ if $_ eq 'UNMAPPED';
						}
					}
					$count++ if $unmapped < ($#segments +1);
				}
			}
		} catch {
			
		};
	}	

	return ('.','.') unless $count;

	$count /= ($#{$self->db->{$abbr}} + 1);			
	my $libsize = $self->sizes->{$abbr} / ($#{$self->db->{$abbr}} + 1);

	my $rpkm = $libsize ? ($count / ($libsize/10^6)) / (($stop-$start)/10^3) : 0;
	my $tpm = $libsize ? ($count / (($stop-$start)/10^3)) / ($libsize/10^6) : 0;

	$tpm = $count/($stop-$start) * 

	return ($tpm == 0 ? '.' : sprintf("%.10f",$tpm), $rpkm == 0 ? '.' : sprintf("%.10f",$rpkm), $count == 0 ? '.' : $count);
}

1;