package Bio::Gorap::Tool::Blast;

use Moose; with 'Bio::Gorap::ToolI';
use POSIX qw(:sys_wait_h);
use File::Spec::Functions;
use IPC::Cmd qw(run);
use IPC::Open3;
use List::Util qw(max any);
use Symbol qw(gensym);

sub calc_features {
	my ($self) = @_;

	#calculations and software calls
	#results are fetched and stored in DB structure
	for (0..$#{$self->parameter->genomes}){
		my $genome = ${$self->parameter->genomes}[$_];
		my $abbr = ${$self->parameter->abbreviations}[$_];
		my $threads = $self->threads;
		my $cmd = $self->cmd;
		$cmd=~s/\$genome/$genome/;
		$cmd=~s/\$cpus/$threads/;

		unless (-e $genome.'.nhr'){
			unlink $genome.$_ for qw( .nhr .nin .nnd .nni .nog .nsd .nsi .nsq);
			my ($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => "makeblastdb -in $genome -dbtype nucl -parse_seqids" , verbose => 0 );
		}

		my $pid = open3(gensym, \*READER, File::Spec->devnull , $cmd);
		my @tab;
		while( <READER> ) {
			chomp $_;
			$_ =~ s/^\s+|\s+$//g;
			next if $_=~/^#/;
			next if $_=~/^\s*$/;
			my @l = split /\s+/ , $_;
			if ($l[1]=~/^ref\|(.+)\|$/ && ! any {$_ eq $l[1]} @{$self->fastadb->oheaders->{$abbr}} ){
				$l[1]=$1;
			}
			if ($l[8]>$l[9]){
				push @tab, $l[1]."\tBlast\t".$self->parameter->cfg->rf_rna."\t$l[9]\t$l[8]\t$l[11]\t-\t$l[10]\t$l[0]\t$l[6]\t$l[7]";
			} else {
				push @tab, $l[1]."\tBlast\t".$self->parameter->cfg->rf_rna."\t$l[8]\t$l[9]\t$l[11]\t+\t$l[10]\t$l[0]\t$l[6]\t$l[7]";
			}
		}
		waitpid($pid, 0);

		my $uid;
		for (@{&merge_gff($self,\@tab)}){
			#tool_parser is set to blast_parser via Gorap.pl, static defined in Bio::Gorap::Functions::ToolParser
			chomp $_;
			my @l = split /\s+/ , $_;
			my $c = ++$uid->{join('.',($abbr,$l[0],$self->parameter->cfg->rf_rna))};
			my @gff3entry = &{$self->tool_parser}($c,$abbr,$self->parameter->cfg->rf_rna,\@l);
			$self->gffdb->add_gff3_entry(\@gff3entry,$self->fastadb->get_gff3seq(\@gff3entry));
		}
	}

	$self->gffdb->merge($self->parameter->cfg->rf_rna,$self->tool); #merge with infernal results
}

sub merge_gff {
	my ($self,$s) = @_;

	my @tab = @{$s};
	@tab = sort{my @a = split /\s+/,$a; my @b = split /\s+/,$b; $a[6] cmp $b[6] || $a[0] cmp $b[0] || $a[3] <=> $b[3] || $a[4] <=> $b[4]} @tab;

	my ($id, $strand, $sta, $sto, $q, $e);
	my @tmp;
	my @merge;
	for my $item (@tab){
		@tmp=split(/\s+/,$item);
		if (defined $id && $tmp[0] eq $id && $tmp[6] eq $strand && $tmp[3]-1<=$sto){
			$sto=max($sto,$tmp[4]);
			if ($tmp[5]>$e){
				$e=$tmp[5];
				$q=join("\t",@tmp[7..$#tmp]);
			}
		} else {
			push @merge, $id."\t".$tmp[1]."\t".$tmp[2]."\t".$sta."\t".$sto."\t".$e."\t".$strand."\t".$q if defined $id;
			$id=$tmp[0];
			$strand=$tmp[6];
			$sta=$tmp[3];
			$sto=$tmp[4];
			$e=$tmp[5];
			$q=join("\t",@tmp[7..$#tmp]);
		}
	}
	push @merge, $id."\t".$tmp[1]."\t".$tmp[2]."\t".$sta."\t".$sto."\t".$e."\t".$strand."\t".$q if defined $id;

	return \@merge;
}

1;