package Bio::Gorap::DB::BAM;

use Moose;
use Bio::DB::Sam;
use List::Util qw(max);
use Try::Tiny;
use IO::Select;
use IO::Pipe;

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

sub _set_db {
	my ($self) = @_;
	return unless $self->parameter->has_bams;
	print "Reading BAM files\n" if $self->parameter->verbose;

	my $sizes;
	for (0..$#{$self->parameter->bams}){
		my $abbr = ${$self->parameter->abbreviations}[$_];
		$sizes->{$abbr}=0;
		for my $f (@{${$self->parameter->bams}[$_]}){
			push @{$self->db->{$abbr}} , Bio::DB::Sam->new(-bam => $f, -autindex => 1, -verbose => -1, -expand_flags => 1);
			unless (-e $f.'.bai'){
				print "$f\n" if $self->parameter->verbose;
				system("samtools index $f");	
			}
		}
	}
}

sub calculate_tpm {
	my ($self,$features) = @_;

	my $select = IO::Select->new();
	my $thrs={};
	my @out;
	my $id2feature;
	for my $f (@$features){

		if (scalar(keys %{$thrs}) >= $self->parameter->threads){
			my $pid = wait();
			delete $thrs->{$pid};
			while( my @responses = $select->can_read(0) ){
				for my $pipe (@responses){
					push @out , $_ while <$pipe>;
					$select->remove( $pipe->fileno() );
				}
			}
		}
		my $fullid = $f->seq_id;
		my @id = split /\./,$f->seq_id;
		my $abbr = shift @id;
		my $count = pop @id;
		my $chr = join '.',@id;
		my $strand = $f->strand;
		my $start = $f->start;
		my $stop = $f->stop;
		my $l = $f->stop - $f->start;
		my $type = $f->primary_tag;
		my $source = $f->source;

		$id2feature->{$f->seq_id.".$type.$source"} = $f;

		my $pipe = IO::Pipe->new();
		if (my $pid = fork()) {
			$pipe->reader();
			$select->add( $pipe );
			$thrs->{$pid}++;
		} else {
			$pipe->writer();
			$pipe->autoflush(1);

			my $c = 0;
			for my $bam (@{$self->db->{$abbr}}){
				$bam->clone;
				for ($bam->get_features_by_location(-type => 'read_pair', -seq_id => $chr, -start => $start, -end => $stop)){
					my ($mate1,$mate2) = $_->get_SeqFeatures;
					if ($mate2){ #paired data
						#1: SE or FR, -1: RF
						if ($self->parameter->strandspec == 0 || $f->strand == 0){
							$c++ if ! $_->get_tag_values('UNMAPPED') && ! $mate1->get_tag_values('M_UNMAPPED');
						} elsif ($self->parameter->strandspec == 1) {
							$c++ if ! $_->get_tag_values('UNMAPPED') && $mate1->strand == $strand && ! $mate1->get_tag_values('M_UNMAPPED');
						} else {
							$c++ if ! $_->get_tag_values('UNMAPPED') && $mate2->strand == $strand && ! $mate1->get_tag_values('M_UNMAPPED');
						}
					} else {
						if ($self->parameter->strandspec == 0 || $f->strand == 0){
							$c++ if ! $_->get_tag_values('UNMAPPED');
						} else {
							$c++ if ! $_->get_tag_values('UNMAPPED') && $mate1->strand == $self->parameter->strandspec;
						}
					}
				}
			}
			print $pipe "$fullid.$type.$source"."\t".$c."\t",$c/($l/1000),"\n";
			exit;
		}
	}

	for (keys %{$thrs} ) {
		my $pid = wait;
		delete $thrs->{$pid};
		while( my @responses = $select->can_read(0) ){
			for my $pipe (@responses){
				push @out , $_ while <$pipe>;
				$select->remove( $pipe->fileno() );
			}
		}
	}

	my $libsize = 0;
	for (@out){
		my @l = split /\s+/,$_;
		$libsize += $l[-2];
	}
	$libsize /= 1000000;

	for (@out){
		my @l = split /\s+/,$_;
		my $tpm = pop @l;
		$tpm /= $libsize;
		my $reads = pop @l;
		my $f = $id2feature->{$l[0]};
		$f->remove_tag('tpm');
		$f->remove_tag('reads');
		$f->add_tag_value('tpm',$tpm);
		$f->add_tag_value('reads',$reads);
		$f->update;
	}
	
}

sub calculate_tpm_se {
	my ($self, $features) = @_;

	my $libsize = 0;
	my @tpms;
	my @counts;
	for my $f (@{$features}){
		my @id = split /\./,$f->seq_id;
		my $abbr = shift @id;
		my $count = pop @id;
		my $chr = join '.',@id;
		my $c = 0;
		my $l = $f->stop - $f->start;
		for my $bam (@{$self->db->{$abbr}}){
			for ($bam->get_features_by_location(-type => 'read_pair', -seq_id => $chr, -start => $f->start, -end => $f->stop)){
				my ($mate1,$mate2) = $_->get_SeqFeatures;
				if ($mate2){ #paired data
					#1: SE or FR, -1: RF
					if ($self->parameter->strandspec == 0 || $f->strand == 0){
						$c++ if ! $_->get_tag_values('UNMAPPED') && ! $mate1->get_tag_values('M_UNMAPPED');
					} elsif ($self->parameter->strandspec == 1) {
						$c++ if ! $_->get_tag_values('UNMAPPED') && $mate1->strand == $f->strand && ! $mate1->get_tag_values('M_UNMAPPED');
					} else {
						$c++ if ! $_->get_tag_values('UNMAPPED') && $mate2->strand == $f->strand && ! $mate1->get_tag_values('M_UNMAPPED');
					}
				} else {
					if ($self->parameter->strandspec == 0 || $f->strand == 0){
						$c++ if ! $_->get_tag_values('UNMAPPED');
					} else {
						$c++ if ! $_->get_tag_values('UNMAPPED') && $mate1->strand == $self->parameter->strandspec;
					}
				}
			}
		}

		push @tpms, $c/($l/1000);
		$libsize += $c;
		push @counts, $c;
	}

	$libsize /= 1000000;
	$_ /= $libsize for @tpms;

	return (\@tpms,\@counts);
}

1;