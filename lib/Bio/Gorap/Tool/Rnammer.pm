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

		for my $rfrna ( qw(RF00177_SSU_rRNA_bacteria RF01959_SSU_rRNA_archaea RF02540_LSU_rRNA_archaea RF02543_LSU_rRNA_eukarya RF01960_SSU_rRNA_eukarya RF02541_LSU_rRNA_bacteria RF00001_5S_rRNA RF02542_SSU_rRNA_microsporidia) ){
			for my $f ($self->gffdb->db->{$abbr}->features(-primary_tag => $rfrna , -attributes => {source => $self->tool})){
				$self->gffdb->db->{$abbr}->delete($f);
				if (exists $self->stkdb->db->{$rfrna}){								
					$self->stkdb->db->{$rfrna}->remove_seq($_) for $self->stkdb->db->{$rfrna}->get_seq_by_id($f->seq_id);								
				}
			}				
		}
	}
		
	my @kingdoms;
	
	push @kingdoms , 'arc' if exists $self->parameter->kingdoms->{'arc'};
	push @kingdoms , 'bac' if exists $self->parameter->kingdoms->{'bac'};
	push @kingdoms , 'euk' if exists $self->parameter->kingdoms->{'fungi'} || exists $self->parameter->kingdoms->{'euk'};

	splice @{$self->parameter->cfg->cmd}, 1, 0, '-T '.$self->parameter->tmp;

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
				for (@{$self->parameter->cfg->cmd}){
					$_ =~ s/\$genome/$genome/;
					$_ =~ s/\$kingdom/$kingdom/;
					$_ =~ s/\$output/$tmpfile/;
				}
				my ($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => join(' ' , @{$self->parameter->cfg->cmd}), verbose => 0 );
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