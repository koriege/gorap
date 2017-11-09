package Bio::Gorap::Tool::Crt;

use Moose; with 'Bio::Gorap::ToolI';

use POSIX qw(:sys_wait_h);				
use IPC::Cmd qw(run);
use File::Spec::Functions;
use Symbol qw(gensym);
use IPC::Open3;
use IO::Select;
use IO::Pipe;
use File::Basename;

sub calc_features {
	my ($self) = @_;
	
	#calculations and software calls
	#results are fetched and stored in DB structure
	
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

			my $cmd = $self->cmd;
			$cmd =~ s/\$genome/$genome/;
			my $pid = open3(gensym, \*READER, File::Spec->devnull, $cmd);

			my $id;	
			while( <READER> ) {	
				chomp $_;
				$_ =~ s/^\s+|\s+$//g;
				next if $_=~/^#/;
				next if $_=~/^\s*$/;
				if ($_=~/^ORGANISM:\s+(\S+)/){
					$id = $1;
				}
				next unless $id;
				if ($_=~/^(\d+)/){
					my @l = split /\s+/ , $_;
					my $start = $l[0];
					my $stop = $start + length($l[1]) - 1;
					print $pipe join(' ' , ($id, 'GORAPcrt' , 'CRISPR' , $start , $stop , '.' , '.' , '.', "\n"));
					$start = $stop;
				}
			}
			waitpid($pid, 0);
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

	my $uid;
	my $scorefile = catfile($self->parameter->tmp,$self->parameter->pid.'.score');
	for (@out){
		my @gff3entry = split /\s+/, $_;

		($gff3entry[0], $gff3entry[3], $gff3entry[4]) = $self->fastadb->chunk_backmap($gff3entry[0], $gff3entry[3], $gff3entry[4]);		
		$gff3entry[0].='.0';

		my @seqs;
		$gff3entry[6] = '+';
		push @seqs , $self->fastadb->get_gff3seq(\@gff3entry);
		$gff3entry[6] = '-';
		push @seqs , $self->fastadb->get_gff3seq(\@gff3entry);

		my $maxscore = -999999;
		my $maxfamily;
		my $strand;
		for my $cm (glob catfile($ENV{GORAP},'gorap','data','rfam','*','*CRISPR*.cm')){
			for my $i ( 0..1 ){
				my $seq = $seqs[$i];
				
				my ($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => "printf \"\>foo\\n$seq\" | cmalign --mxsize ".$self->parameter->mem." --noprob --sfile $scorefile --cpu ".$self->threads." $cm -", verbose => 0 );
		
				open S , '<'.$scorefile or die $!;
				while(<S>){
					chomp $_;
					$_ =~ s/^\s+|\s+$//g;
					next if $_=~/^#/;
					next if $_=~/^\s*$/;
					my $score = (split /\s+/ , $_)[6];
					if ($score > $maxscore){
						$maxscore = $score;
						$maxfamily = basename(dirname($cm));
						$strand = $i == 0 ? '+' : '-'; 
					}
				}
				close S;
				unlink $scorefile;
			}
		}

		if ($maxscore >= 10){					
			$gff3entry[2] = $maxfamily;
			$gff3entry[5] = $maxscore;
			$gff3entry[6] = $strand;
		} else {
			$gff3entry[6] = '.';
		}

		my ($abbr,@orig) = split /\./ , $gff3entry[0];
		pop @orig;
		$uid->{$abbr.'.'.$gff3entry[2]}++;
		$gff3entry[0] = join('.',($abbr,@orig,$uid->{$abbr.'.'.$gff3entry[2]}));

		my $seq = $gff3entry[6] eq '.' ? '' : $seqs[$gff3entry[6] eq '+' ? 0 : 1];
		$self->gffdb->add_gff3_entry(\@gff3entry,$seq);
	}		
}

1;
