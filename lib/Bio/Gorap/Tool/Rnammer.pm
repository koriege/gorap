package Bio::Gorap::Tool::Rnammer;

use Moose; with 'Bio::Gorap::ToolI';
use IO::Select;
use IO::Pipe;
use IPC::Cmd qw(run);
use File::Spec::Functions;

sub calc_features {
	my ($self) = @_;

	#calculations and software calls
	#results are fetched and stored in DB structure
	for (0..$#{$self->parameter->genomes}){
		my $abbr = ${$self->parameter->abbreviations}[$_];
		#skip redundand calculations
		my @f = $self->gffdb->db->{$abbr}->features(-attributes => {source => 'GORAP'.$self->tool});
		return if $#f > -1;

		my @f = $self->gffdb->db->{$abbr}->features(-attributes => {source => $self->tool});
		for (@f){
			my $rfrna = $_->primary_tag;
			$self->stkdb->db->{$rfrna}->remove_seq($_) for $self->stkdb->db->{$rfrna}->get_seq_by_id($f->seq_id);
			$self->gffdb->db->{$abbr}->delete($f);
		}
	}

	my @kingdoms;

	push @kingdoms , 'arc' if exists $self->parameter->kingdoms->{'arc'};
	push @kingdoms , 'bac' if exists $self->parameter->kingdoms->{'bac'};
	push @kingdoms , 'euk' if exists $self->parameter->kingdoms->{'fungi'} || exists $self->parameter->kingdoms->{'euk'};

	my $select = IO::Select->new();
	my $thrs={};
	my @out;
	for my $kingdom (@kingdoms){

		for my $genome (@{$self->fastadb->chunks}){
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

				my $tmpfile = catfile($self->parameter->tmp,$$.'.rnammer');

				my $cmd = $self->cmd;
				$cmd =~ s/\$genome/$genome/;
				$cmd =~ s/\$cpus/$threads/;
				$cmd =~ s/\$kingdom/$kingdom/;
				$cmd =~ s/\$output/$tmpfile/;
				$cmd .= ' -T'.$self->parameter->tmp;
				my ($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => $cmd, verbose => 0 );

				open F ,'<'.$tmpfile or exit;
				while( <F> ) {
					chomp $_;
					$_ =~ s/^\s+|\s+$//g;
					next if $_=~/^#/;
					next if $_=~/^\s*$/;
					my @l = split /\s+/ , $_;
					next if $#l < 8;
					print $pipe $_."\t".$kingdom."\n";
				}
				close F;
				unlink $tmpfile;
				exit;
			}
		}
	}
	for (keys %{$thrs} ) {
		my $pid = wait();
		delete $thrs->{$pid};
		while( my @responses = $select->can_read(0) ){
			for my $pipe (@responses){
				push @out , $_ while <$pipe>;
				$select->remove( $pipe->fileno() );
			}
		}
	}

	my $uid;
	for (@out){
		my @l = split /\s+/, $_;
		my $kingdom = pop @l;

		my @gff3entry = &{$self->tool_parser}($kingdom,\@l);
		($gff3entry[0], $gff3entry[3], $gff3entry[4]) = $self->fastadb->chunk_backmap($gff3entry[0], $gff3entry[3], $gff3entry[4]);
		$gff3entry[0] .= '.'.(++$uid->{$gff3entry[0].'.'.$gff3entry[2]});

		my $seq = $self->fastadb->get_gff3seq(\@gff3entry);
		$self->gffdb->add_gff3_entry(\@gff3entry,$seq);
	}
}

1;