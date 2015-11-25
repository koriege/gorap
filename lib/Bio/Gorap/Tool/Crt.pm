package Bio::Gorap::Tool::Crt;

use Moose; with 'Bio::Gorap::ToolI';

use POSIX qw(:sys_wait_h);				
use IPC::Cmd qw(run);
use File::Spec::Functions;
use Symbol qw(gensym);
use IPC::Open3;
use IO::Select;
use IO::Pipe;

sub calc_features {
	my ($self) = @_;
	
	#calculations and software calls
	#results are fetched and stored in DB structure
	my @related_rf_rna = qw( RF01315_CRISPR-DR2
							RF01316_CRISPR-DR3
							RF01317_CRISPR-DR4
							RF01318_CRISPR-DR5
							RF01319_CRISPR-DR6
							RF01320_CRISPR-DR7
							RF01321_CRISPR-DR8
							RF01322_CRISPR-DR9
							RF01323_CRISPR-DR10
							RF01324_CRISPR-DR11
							RF01325_CRISPR-DR12
							RF01326_CRISPR-DR13
							RF01327_CRISPR-DR14
							RF01328_CRISPR-DR17
							RF01329_CRISPR-DR15
							RF01330_CRISPR-DR16
							RF01331_CRISPR-DR18
							RF01332_CRISPR-DR19
							RF01333_CRISPR-DR20
							RF01334_CRISPR-DR21
							RF01335_CRISPR-DR22
							RF01336_CRISPR-DR23
							RF01337_CRISPR-DR24
							RF01338_CRISPR-DR25
							RF01339_CRISPR-DR27
							RF01340_CRISPR-DR29
							RF01341_CRISPR-DR30
							RF01342_CRISPR-DR32
							RF01343_CRISPR-DR33
							RF01344_CRISPR-DR34
							RF01345_CRISPR-DR35
							RF01346_CRISPR-DR36
							RF01347_CRISPR-DR37
							RF01348_CRISPR-DR38
							RF01349_CRISPR-DR40
							RF01350_CRISPR-DR41
							RF01351_CRISPR-DR42
							RF01352_CRISPR-DR43
							RF01353_CRISPR-DR44
							RF01354_CRISPR-DR45
							RF01355_CRISPR-DR26
							RF01356_CRISPR-DR28
							RF01357_CRISPR-DR31
							RF01358_CRISPR-DR39
							RF01359_CRISPR-DR46
							RF01360_CRISPR-DR47
							RF01361_CRISPR-DR48
							RF01362_CRISPR-DR49
							RF01363_CRISPR-DR50
							RF01364_CRISPR-DR51
							RF01365_CRISPR-DR52
							RF01366_CRISPR-DR53
							RF01367_CRISPR-DR54
							RF01368_CRISPR-DR55
							RF01369_CRISPR-DR56
							RF01370_CRISPR-DR57
							RF01371_CRISPR-DR58
							RF01373_CRISPR-DR60
							RF01374_CRISPR-DR61
							RF01375_CRISPR-DR62
							RF01376_CRISPR-DR63
							RF01377_CRISPR-DR64
							RF01378_CRISPR-DR65
							RF01379_CRISPR-DR66 
	);

	for (0..$#{$self->parameter->genomes}){		
		my $abbr = ${$self->parameter->abbreviations}[$_];
		#skip redundand calculations
		my @f = $self->gffdb->db->{$abbr}->features(-attributes => {source => 'GORAP'.$self->tool});
		return if $#f > -1;
	
		for my $rf_rna ( (@related_rf_rna,'CRISPR') ){
			for my $f ($self->gffdb->db->{$abbr}->features(-primary_tag => $rf_rna , -attributes => {source => $self->tool})){
				$self->gffdb->db->{$abbr}->delete($f);
				if (exists $self->stkdb->db->{$rf_rna}){								
					$self->stkdb->db->{$rf_rna}->remove_seq($_) for $self->stkdb->db->{$rf_rna}->get_seq_by_id($f->seq_id);								
				}
			}				
		}
	}
		
	my $toolPath = `which CRT1.2-CLI.jar`;
	chomp $toolPath;
	
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

			for (@{$self->parameter->cfg->cmd}){
				$_ =~ s/crt/java -cp $toolPath crt/;
				$_ =~ s/\$genome/$genome/;				
			}		

			my $pid = open3(gensym, \*READER, File::Spec->devnull , join(' ' , @{$self->parameter->cfg->cmd}));

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
					print $pipe join(' ' , ($id, 'GORAPcrt' , 'CRISPR' , $start , $stop , '.' , '+' , '.', "\n"));
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

	my $uid=0;
	my $scorefile = catfile($self->parameter->tmp,$self->parameter->pid.'.score');
	for (@out){
		my @l = split /\s+/, $_;				
		for ($self->fastadb->chunk_backmap($l[0], $l[3], $l[4])){
			($l[0],$l[3],$l[4]) = @{$_};
			
			my ($abbr, @header) = split /\./,$l[0];
			# $uid->{$l[0]}++;
			$uid++;
			$l[0] = $l[0].'.'.$uid;
			my @gff3entry = @l;
			$l[6] = '-';
			my @seqs;
			
			#due to overlapping chunks check for already annotated genes
			my $existingFeatures = [@{$self->gffdb->get_overlapping_features(\@gff3entry,$abbr)},
									@{$self->gffdb->get_overlapping_features(\@l,$abbr)}];
			next if $#{$existingFeatures} > -1;
							
			push @seqs , $self->fastadb->get_gff3seq(\@gff3entry);					
			push @seqs , $self->fastadb->get_gff3seq(\@l);

			my $maxscore = -999999;
			my $maxfamily;
			my $strand;
			for my $rf_rna ( @related_rf_rna){
				my $cm = catfile($ENV{GORAP},'data','rfam',$rf_rna,$rf_rna.'.cm');						
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
							$maxfamily = $rf_rna;
							$strand = $i == 0 ? '+' : '-'; 
						}
					}
					close S;
				}
			}					

			if ($maxscore >= 10){					
				$gff3entry[2] = $maxfamily;
				$gff3entry[5] = $maxscore;
				$gff3entry[6] = $strand;
				$self->gffdb->add_gff3_entry(\@gff3entry,$seqs[$strand eq '+' ? 0 : 1],$abbr);					
			} else {
				$self->gffdb->add_gff3_entry(\@gff3entry,$seqs[0],$abbr);
			}
		}
	}		
}

1;
