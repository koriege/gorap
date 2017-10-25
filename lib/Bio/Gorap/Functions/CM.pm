package Bio::Gorap::Functions::CM;

use Moose;
use Symbol qw(gensym);
use IPC::Open3;

sub compute {
	my ($self, $cpus, $stk, $cm);

	local *READER;
	my $pid = open3(gensym, \*READER, File::Spec->devnull , "cmbuild -F $cm $stk");
	waitpid($pid, 0);
	$pid = open3(gensym, \*READER, File::Spec->devnull , "cmcalibrate --cpu $cpus $cm");
	waitpid($pid, 0);	
}

1;