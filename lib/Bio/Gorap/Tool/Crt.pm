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
use Bio::Seq;

sub calc_features {
	my ($self) = @_;

	return if $self->already_predicted;

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
					print $pipe $_."\t".$id."\n";
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

	return if $#out == -1;
	
	my @sequences;
	my $chr2gff;
	for (0..$#out){
		my @l = split /\s+/ , $out[$_];
		my $id = pop @l;
		my @gff3entry = &{$self->tool_parser}($self->tool,$id,\@l);
		($gff3entry[0], $gff3entry[3], $gff3entry[4]) = $self->fastadb->chunk_backmap($gff3entry[0], $gff3entry[3], $gff3entry[4]);

		$gff3entry[0] .= '.'.$_;
		$chr2gff->{$gff3entry[0]} = \@gff3entry;
		$gff3entry[6] = '+';
		push @sequences , Bio::Seq->new( -display_id => $gff3entry[0].'.'.$gff3entry[6] , -seq => $self->fastadb->get_gff3seq(\@gff3entry));
		$gff3entry[6] = '-';
		push @sequences , Bio::Seq->new( -display_id => $gff3entry[0].'.'.$gff3entry[6] , -seq => $self->fastadb->get_gff3seq(\@gff3entry));
	}

	my $chr2score;
	my $cfg = $self->parameter->cfg->cfg;
	for my $cm (glob catfile($ENV{GORAP},'gorap','data','rfam','*','*CRISPR*.cm')){
		my $rf_rna = basename(dirname($cm));

		$self->parameter->set_cfg(catfile($ENV{GORAP},'gorap','config',$rf_rna.'.cfg'));
		my ($scorefile,$stk) = $self->stkdb->align($rf_rna,\@sequences,$self->threads,$cm,catfile($self->parameter->tmp,$rf_rna.'.score'),1);

		open F,'<'.$scorefile or die $!;
		while(<F>){
			chomp $_ ;
			$_ =~ s/^\s+|\s+$//g;
			next if $_=~/^#/;
			my @l = split /\s+/,$_;
			my @chr = split /\./,$l[1];
			my $strand = pop @chr;
			$l[1] = join('.',@chr);
			if (exists $chr2score->{$l[1]}){
				if ($l[6] > ${$chr2score->{$l[1]}}[0]){
					$chr2score->{$l[1]} = [$l[6],$strand,$rf_rna];
				}
			} else {
				$chr2score->{$l[1]} = [$l[6],$strand,$rf_rna];
			}
			
		}
		close F;
	}
	$self->parameter->set_cfg($cfg);

	my $uid;
	my $types;
	while (my ($chr, $arr) = each %$chr2score) {
		my ($score,$strand,$rf_rna) = @{$arr};
		$types->{$rf_rna}=1;
		my @gff3entry = @{$chr2gff->{$chr}};
		my @chr = split /\./ , $gff3entry[0];
		if ($score >= 10){
			$chr[-1] = $rf_rna;
			$uid->{join('.',@chr)}++;
			$chr[-1] = 'crt'.$uid->{join('.',@chr)};
			$gff3entry[0] = join('.',@chr);
			$gff3entry[2] = $rf_rna;
			$gff3entry[5] = $score;
			$gff3entry[6] = $strand;
			my $seq = $self->fastadb->get_gff3seq(\@gff3entry);
			$self->gffdb->add_gff3_entry(\@gff3entry,$seq);
		} else {
			$chr[-1] = 'CRISPR';
			$uid->{join('.',@chr)}++;
			$chr[-1] = 'crt'.$uid->{join('.',@chr)};
			$gff3entry[0] = join('.',@chr);
			$gff3entry[5] = 0;
			$self->gffdb->add_gff3_entry(\@gff3entry,'');
		}
	}

	for(keys %$types){
		$self->gffdb->merge($_,$self->tool); #merge multi kingdoms and overlapping annotations du to genome chunks
	}
}

1;
