package Bio::Gorap::Tool::Piles;

use Moose; with 'Bio::Gorap::ToolI';

use List::MoreUtils qw(pairwise);
use IO::Select;
use IO::Pipe;
use List::Util qw(min max);

sub calc_features {
	my ($self) = @_;

	my $select = IO::Select->new();
	my $thrs={};
	my @out;

	for my $abbr (keys %{$self->fastadb->oheaders}){
		next if $#{$self->bamdb->db->{$abbr}} == -1;
		for my $chr (@{$self->fastadb->oheaders->{$abbr}}){
			my $start = 1;
			my $stop = 0;
			my $length = $self->fastadb->get_seq($abbr.'.'.$chr)->length;
			while($stop < $length){
				$stop = min($length,$start + 10000);

				if (scalar(keys %{$thrs}) >= $self->threads){
					my $pid = wait();
					delete $thrs->{$pid};
					while( my @responses = $select->can_read(0) ){
						for my $pipe (@responses){
							push @out , $_ while <$pipe>;
							$select->remove( $pipe->fileno() );
						}
					}
				}

				my $pipe = IO::Pipe->new();
				if (my $pid = fork()) {
					$pipe->reader();
					$select->add( $pipe );
					$thrs->{$pid}++;
				} else {
					$pipe->writer();
					$pipe->autoflush(1);
					my @coverage;
					for my $bam (@{$self->bamdb->db->{$abbr}}){
						$bam->clone;
						my ($c) = $bam->features(-type => 'coverage', -seq_id => $chr, -start => $start, -stop => $stop);
						my @c = $c->coverage;
						@coverage = $#coverage == -1 ? @c : pairwise { $a + $b } @coverage,@c;
					}
					for (my $j=0; $j<=$#coverage; $j++){
						my $sta=$j;
						$j++ while $coverage[$j] && $coverage[$j] >= $self->parameter->denovoheight;
						print $pipe $abbr.'.'.$chr,"\t",$start+$sta,"\t",$start+$j,"\n" unless $sta == $j;
					}
					exit;
				}
				$start = $stop + 1;				
			}
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

	@out = sort {my @a=split/\s+/,$a; my @b=split/\s+/,$b; $a[0] cmp $b[0] || $a[1] <=> $b[1]} @out;
	my $uid;
	for my $i (0..$#out-1){
		my @li = split /\s+/,$out[$i];
		my @lj = split /\s+/,$out[$i+1];
		my $chr = $li[0];
		my $start = $li[1];
		my $stop = $li[2];
		while($i < $#out && $li[0] eq $lj[0] && $li[2] == $lj[1]){
			$stop = $lj[2];
			$i++;
			@li = split /\s+/,$out[$i];
			@lj = split /\s+/,$out[$i+1];
		}
		if ($stop - $start >= $self->parameter->denovolength){
			$uid->{$chr}++;
			my @gff3entry = ($chr.'.denovo'.$uid->{$chr} , 'GORAPdenovo' , 'NEW_RNA' , $start , $stop, '.' , '.' , '.');
			$self->gffdb->add_gff3_entry(\@gff3entry,'');
		}
	}

}

1;
