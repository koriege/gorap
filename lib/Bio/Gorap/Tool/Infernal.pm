package Bio::Gorap::Tool::Infernal;

use Moose; with 'Bio::Gorap::ToolI';

use POSIX qw(:sys_wait_h);
use IPC::Open3;
use File::Spec::Functions;
use Symbol qw(gensym);
use File::Temp;
use IPC::Cmd qw(run);

sub calc_features {
	my ($self) = @_;

	for (0..$#{$self->parameter->genomes}){
		my $genome = ${$self->parameter->genomes}[$_];
		my $abbr = ${$self->parameter->abbreviations}[$_];
		my $threads = $self->threads;
		my $cmd = $self->cmd;
		$cmd=~s/\$genome/$genome/;
		$cmd=~s/\$cpus/$threads/;

		my $tmpfile = File::Temp->new(DIR => $self->parameter->tmp)->filename;
		$cmd .= " > $tmpfile";
		my ($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => $cmd, verbose => 0 );
		open F,"<$tmpfile" or die $!;
		
		my $uid;
		#my $pid = open3(gensym, \*READER, File::Spec->devnull, $cmd);
		#while( <READER> ) {
		while(<F>){
			chomp $_;
			$_ =~ s/^\s+|\s+$//g;
			next if $_=~/^#/;
			next if $_=~/^\s*$/;
			my @l = split /\s+/ , $_;
			next if $#l < 11;

			my $c = ++$uid->{join('.',($abbr,$l[5],$self->parameter->cfg->rf_rna))};
			my @gff3entry = &{$self->tool_parser}($c,$abbr,$self->parameter->cfg->rf_rna,\@l);

			$self->gffdb->add_gff3_entry(\@gff3entry,$self->fastadb->get_gff3seq(\@gff3entry));
		}
		#waitpid($pid, 0);
		close F;
	}

	$self->gffdb->merge($self->parameter->cfg->rf_rna,$self->tool); #merge with blast results
}

1;