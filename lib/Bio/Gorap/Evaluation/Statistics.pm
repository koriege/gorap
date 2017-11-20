package Bio::Gorap::Evaluation::Statistics;

use Moose;
use POSIX;
use List::Util qw(min max);
use Bio::SeqIO;
use File::Spec::Functions;

has 'parameter' => (
	is => 'ro',
	isa => 'Bio::Gorap::Parameter',
	required => 1
);

has 'taxdb' => (
	is => 'ro',
	isa => 'Bio::Gorap::DB::Taxonomy',
	required => 1
);

sub create_genomes {
	my ($self) = @_;

	my $sequences;
	my $lengths;

	for my $cfg (@{$self->parameter->queries}){
		$self->parameter->set_cfg($cfg);
		next if $self->parameter->cfg->rna=~/CRISPR/;
		next if $self->parameter->cfg->rf_rna =~/_rRNA/;
		my $testfasta = catfile($self->parameter->cfg->query_dir,$self->parameter->cfg->rf_rna.'.testing.fa');
		next unless -e $testfasta;

		my $seqio = Bio::SeqIO->new( '-format' => 'Fasta' , -file => $testfasta, -verbose => -1);
		while(my $seqobj = $seqio->next_seq()){
			my $name = (split /\./,$seqobj->id)[0];
			my $taxid = $self->taxdb->getIDfromName($name);
			my $kingdom;
			if ($taxid){
				my $nodes = $self->taxdb->getLineageNodes($taxid);
				if ($#{$nodes} > 0){
					$kingdom = ${$nodes}[1]->id;
					$kingdom = $kingdom == 2157 ? 'arc' : $kingdom == 2 ? 'bac' : 'euk';
				} else {
					$kingdom = 'euk';
				}
			} else {
				$kingdom = 'euk';
			}

			push @{$sequences->{$kingdom}} , { id => $self->parameter->cfg->rna.'.'.$seqobj->id, seq => $seqobj->seq};
			$lengths->{$kingdom} += length $seqobj->seq;
		}
	}

	my $genome;
	for my $kingdom (keys %$sequences){
		my $genomesize = 1000000 - $lengths->{$kingdom};
		for (@{$sequences->{$kingdom}}){
			my $chr = floor(rand 20);
			my $pos = floor(rand $genomesize);
			while (exists $genome->{$kingdom}->{$chr}->{$pos}){
				$pos = floor(rand $genomesize);
			}
			$genome->{$kingdom}->{$chr}->{$pos} = $_;
		}
	}

	my @genomes;
	my @abbreviations;
	for my $kingdom (keys %{$genome}){
		push @abbreviations , $kingdom;
		my $genomesize = 1000000 - $lengths->{$kingdom};

		push @genomes , catfile($self->parameter->output,'eval_'.$kingdom.'.fa');

		open GFF , '>'.catfile($self->parameter->output,'eval_'.$kingdom.'.gff') or die $!;
		open FA , '>'.$genomes[-1] or die $!;
		for my $chr (0..19){
			my $sl=0;
			my @arr = sort {$a <=> $b} keys %{$genome->{$kingdom}->{$chr}};
			my $in = $#arr > -1 ? shift @arr : $genomesize;
			print FA '>chr_'.$chr."\n";
			my $nl = 0;
			for my $pos (0..$genomesize-1){
				if ($nl > 0 && $nl % 80 == 0){
					$nl = 0;
					print FA "\n";
				}
				if( $pos < $in){
					my $num = floor(rand 4);
					print FA $num == 0 ? 'A' : $num == 1 ? 'C' : $num == 2 ? 'G' : 'T';
					$nl++;
					$sl++;
				} else {
					for(split // , $genome->{$kingdom}->{$chr}->{$pos}->{'seq'}){
						print FA $_;
						$nl++;
						$sl++;
						if ($nl % 80 == 0){
							$nl = 0;
							print FA "\n";
						}
					}
					print GFF 'chr_'.$chr."\tncRNA\t".$genome->{$kingdom}->{$chr}->{$pos}->{'id'}."\t",$sl-length($genome->{$kingdom}->{$chr}->{$pos}->{'seq'})+1,"\t",$sl,"\t.\t+\n";
					$in = $#arr > -1 ? shift @arr : $genomesize;
				}
			}
			print FA "\n";
		}
		close GFF;
		close FA;
	}

	return (\@genomes, \@abbreviations);
}

1;