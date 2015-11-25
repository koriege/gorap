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

	my $kingdom;

	$kingdom = 'A' if exists $self->parameter->kingdoms->{'arc'};
	$kingdom = 'B' if exists $self->parameter->kingdoms->{'bac'};
	$kingdom = 'E' if exists $self->parameter->kingdoms->{'euk'};
	$kingdom = 'f' if exists $self->parameter->kingdoms->{'fungi'};

	my $select = IO::Select->new();
	my $thrs={};
	my @out;
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

	my $uid=0;
	for (@out){
		my @l = split /\s+/, $_;				

		@l = &{$self->tool_parser}('','',$self->parameter->cfg->rf_rna,\@l);
		$l[0] = (split /\./ , $l[0])[1];
		for ($self->fastadb->chunk_backmap($l[0], $l[3], $l[4])){
			($l[0],$l[3],$l[4]) = @{$_};							

			my ($abbr, @header) = split /\./,$l[0];
			#$uid->{$l[0]}++;
			$uid=0;
			$l[0] .='.'.$uid;
			
			my @gff3entry = @l;
			#due to overlapping chunks check for already annotated genes
			my $existingFeatures = $self->gffdb->get_overlapping_features(\@gff3entry,$abbr);
			next if $#{$existingFeatures} > -1;

			my $seq = $self->fastadb->get_gff3seq(\@gff3entry);
			$self->gffdb->add_gff3_entry(\@gff3entry,$seq,$abbr);		
		}
	}	
}

1;