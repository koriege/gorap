package Bio::Gorap::Functions::CM;

use Moose;
use Symbol qw(gensym);
use IPC::Open3;

sub get_min_score {
	#open given covariance model and extract minimum of gathering-, trusted- and noise cutoff
	my ($self,$cm) = @_;
	
	open CM , '<'.$cm or die $!;
	my $score = 999999;
	while(<CM>){
		last if $_=~/^\/\//;
		if ($_=~/^(GA|TC|NC)/){
			my @line = split /\s+/ , $_;
			$score = $line[1] < $score ? $line[1] : $score;
		}
	}	
	close CM;
	return $score == 999999 ? -1 : $score;
}

sub get_rf_rna {
	my ($self,$cm) = @_;
	
	open CM , '<'.$cm or die $!;	
	my $rna='';
	my $rf='';
	while(<CM>){		
		$rna = $1 if $_=~/^NAME\s+(.+)/;
		$rf = $1 if $_=~/^ACC\s+(.+)/; 		
		do {close CM; return $rf.'_'.$rna} if $rna && $rf;
	}	
	close CM;
	return $rf.'_'.$rna;
}

sub compute {
	my ($self, $cpus, $stk, $cm);

	local *READER;
	my $pid = open3(gensym, \*READER, File::Spec->devnull , "cmbuild -F $cm $stk");
	waitpid($pid, 0);
	$pid = open3(gensym, \*READER, File::Spec->devnull , "cmcalibrate --cpu $cpus $cm");
	waitpid($pid, 0);	
}

1;