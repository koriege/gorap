package Bio::Gorap::Tool::Bcheck;

use Moose; with 'Bio::Gorap::ToolI';
use IPC::Cmd qw(run);
use File::Spec::Functions;
use IO::Select;
use IO::Pipe;

sub calc_features {
	my ($self) = @_;		

	#calculations and software calls
	#results are fetched and stored in DB structure
	for (0..$#{$self->parameter->genomes}){		
		my $abbr = ${$self->parameter->abbreviations}[$_];
		#skip redundand calculations
		my @f = $self->gffdb->db->{$abbr}->features(-attributes => {source => 'GORAP'.$self->tool});
		return if $#f > -1;

		for my $rfrna ( qw(RF00009_RNaseP_nuc RF00010_RNaseP_bact_a RF00011_RNaseP_bact_b RF00373_RNaseP_arch) ){
			for my $f ($self->gffdb->db->{$abbr}->features(-primary_tag => $rfrna , -attributes => {source => $self->tool})){
				$self->gffdb->db->{$abbr}->delete($f);
				if (exists $self->stkdb->db->{$rfrna}){								
					$self->stkdb->db->{$rfrna}->remove_seq($_) for $self->stkdb->db->{$rfrna}->get_seq_by_id($f->seq_id);								
				}
			}				
		}		
	}

	my @kingdoms;

	push @kingdoms , 'A' if exists $self->parameter->kingdoms->{'arc'};
	push @kingdoms , 'B' if exists $self->parameter->kingdoms->{'bac'};
	push @kingdoms , 'E' if exists $self->parameter->kingdoms->{'euk'};
	push @kingdoms , 'f' if exists $self->parameter->kingdoms->{'fungi'};

	my $select = IO::Select->new();
	my $thrs={};
	my @out;
	for my $kingdom (@kingdoms){
		$kingdom = uc substr $kingdom,0,1;
		$kingdom = 'f' if $kingdom eq 'F';

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
				my $tmpfile=catfile($self->parameter->tmp,$$);

				my $foo="/home/koriege/workspace-perl/gorap-dev/foo";
				for (@{$self->parameter->cfg->cmd}){
					$_ =~ s/\$genome/$genome/;
					$_ =~ s/\$kingdom/$kingdom/;
					$_ =~ s/\$output/$tmpfile/;
				}

				my @paths = ( glob(catdir($ENV{GORAP},'infernal-1.0','bin')) , glob(catdir($ENV{GORAP},'rnabob-*','bin')) );
				local $ENV{PATH} = $ENV{PATH} ? join(':',@paths).":$ENV{PATH}" : join(':',@paths);						

				my ($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => join(' ' , @{$self->parameter->cfg->cmd}), verbose => 0 );

				open F , '<'.$tmpfile.'_rnpB.ss' or exit;
				my @l;
				while( <F> ) {		
					chomp $_;
					$_ =~ s/^\s+|\s+$//g;
					next if $_=~/^#/;			
					next if $_=~/^\s*$/;
					if ($_=~/^>/){
						@l = ();
						push @l , $_;
						next;
					}
					next unless $_=~/^Score/;
					print $pipe join (' ' , (@l,$_,"\n"));
				}
				close F;
				unlink $tmpfile.'_rnpB.ss';
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
		my @gff3entry = &{$self->tool_parser}(\@l);
		($gff3entry[0], $gff3entry[3], $gff3entry[4]) = $self->fastadb->chunk_backmap($gff3entry[0], $gff3entry[3], $gff3entry[4]);
		
		my ($abbr,@orig) = split /\./ , $gff3entry[0];
		$uid->{$abbr.'.'.$gff3entry[2]}++;
		$gff3entry[0] = join('.',($abbr,@orig,$uid->{$abbr.'.'.$gff3entry[2]}));
		
		# due to overlapping chunks check for already annotated genes
		my $existingFeatures = $self->gffdb->get_overlapping_features(\@gff3entry);
		if ($#{$existingFeatures} > -1){
			$uid->{$abbr.'.'.$gff3entry[2]}--;	
			next;
		}		
		my $seq = $self->fastadb->get_gff3seq(\@gff3entry);			
		$self->gffdb->add_gff3_entry(\@gff3entry,$seq);
	}
}

1;