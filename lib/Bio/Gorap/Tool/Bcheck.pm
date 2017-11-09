package Bio::Gorap::Tool::Bcheck;

use Moose; with 'Bio::Gorap::ToolI';
use IPC::Cmd qw(run);
use File::Spec::Functions;
use IO::Select;
use IO::Pipe;

sub calc_features {
	my ($self) = @_;

	my @kingdoms;

	push @kingdoms , 'A' if exists $self->parameter->kingdoms->{'arc'};
	push @kingdoms , 'B' if exists $self->parameter->kingdoms->{'bac'};
	push @kingdoms , 'E' if exists $self->parameter->kingdoms->{'euk'};
	push @kingdoms , 'f' if exists $self->parameter->kingdoms->{'fungi'};

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
				my $tmpfile=catfile($self->parameter->tmp,$$);

				my $cmd = $self->cmd;
				$cmd =~ s/\$genome/$genome/;
				$cmd =~ s/\$kingdom/$kingdom/;
				$cmd =~ s/\$output/$tmpfile/;

				local $ENV{PATH} = catdir($ENV{GORAP},'bin','infernal1').':'.$ENV{PATH};
				my ($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => $cmd, verbose => 0 );

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
		
		my $seq = $self->fastadb->get_gff3seq(\@gff3entry);	
		$self->gffdb->add_gff3_entry(\@gff3entry,$seq);
	}
}

1;