package Bio::Gorap::Tool::Piles;

use Moose; with 'Bio::Gorap::ToolI';

use List::MoreUtils qw(pairwise);
use IO::Select;
use IO::Pipe;

sub calc_features {
	my ($self) = @_;

	my $select = IO::Select->new();
	my $thrs={};
	my @out;
	for (@{$self->fastadb->nheaders}){
		my $length = $self->fastadb->db->fetch($_)->length;
		my ($abbr,@orig) = split /\./ , $_;
		next unless exists $self->bamdb->db->{$abbr};
		my $orig = join '.' , @orig;

		my $chunk = sprintf("%.0f", $length / $self->threads);
		for (my $i=0; $i+$chunk<=$length; $i+=$chunk){
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
			my $start = $i+1;
			my $stop = $i+2*$chunk > $length ? $length : $i+$chunk;			
			
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
					#TODO iterate over feature alignments, take strand specific ones, store in new bam object and calc coverage
					my ($c) = $bam->features(-type => 'coverage', -strand => -1, -seq_id => $orig, -start => $start, -stop => $stop);
					my @c = $c->coverage;
					@coverage = $#coverage == -1 ? @c : pairwise { $a + $b } @coverage, @c;
				}
				for (my $j=0; $j<=$#coverage; $j++){
					my $sta=$j;
					$j++ while $coverage[$j] > $self->parameter->denovoheight;
					print $pipe "$abbr.$orig\t".($start+$sta)."\t".($start+$j)."\t.\n" if $j-$sta >= $self->parameter->denovolength || ($sta==0 && $start > 0) || ($j>=$#coverage && $stop < $length);
				}
				exit;
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

	my ($start,$stop,$uid);
	for my $e (sort {my @a=split/\s+/,$a; my @b=split/\s+/,$b; $a[0] cmp $b[0] || $a[1] <=> $b[1]} @out){
		my @l = split /\s+/ , $e;
		my $update;
		if ($stop && $l[1]-1 == $stop){
			$update=1;
			$e = $l[0]."\t".$start."\t".$l[2]."\t".$l[3]."\n" if $l[2] - $start >= $self->parameter->denovolength;
		}
		$start = $l[1];
		$stop = $l[2];		

		unless ($l[3] eq '.'){
			next if $self->parameter->check_overlaps && scalar $self->gffdb->get_user_overlaps($l[0],$update ? $start : $l[1],$l[2],$l[3]) > 0;
		}
		$uid->{$l[0]}++;
		my @gff3entry = ($l[0].'.'.$uid->{$l[0]} , 'GORAPde_novo' , 'NEW_RNA' , $update ? $start : $l[1] , $l[2], '.' , $l[3] , '.');
				
		my $seq = $l[3] eq '.' ? '' : $self->fastadb->get_gff3seq(\@gff3entry);
		$self->gffdb->add_gff3_entry(\@gff3entry,$seq);
	}
}

1;

