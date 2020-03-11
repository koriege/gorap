package Bio::Gorap::Tool::Barrnap;

use Moose; with 'Bio::Gorap::ToolI';
use IO::Select;
use IO::Pipe;
use IPC::Open3;
use IPC::Cmd qw(run);
use File::Spec::Functions;
use Symbol qw(gensym);
use File::Temp;

sub calc_features {
	my ($self) = @_;

	return if $self->already_predicted;

	my $types;
	for (0..$#{$self->parameter->genomes}){
		my $genome = ${$self->parameter->genomes}[$_];
		my $abbr = ${$self->parameter->abbreviations}[$_];
		my $threads = $self->threads;

		my @kingdoms;
		push @kingdoms , 'arc' if exists $self->parameter->kingdoms->{'arc'};
		push @kingdoms , 'bac' if exists $self->parameter->kingdoms->{'bac'};
		push @kingdoms , 'euk' if exists $self->parameter->kingdoms->{'fungi'} || exists $self->parameter->kingdoms->{'euk'};

		my $uid;
		for my $kingdom (@kingdoms){
			my $cmd = $self->cmd;
			$cmd =~ s/\$genome/$genome/;
			$cmd =~ s/\$cpus/$threads/;
			$cmd =~ s/\$kingdom/$kingdom/;

			my $tmpfile = File::Temp->new(DIR => $self->parameter->tmp)->filename;
			$cmd .= " > $tmpfile";
			my ($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => $cmd, verbose => 0 );
			open F,"<$tmpfile" or die $!;

			# my $pid = open3(gensym, \*READER, File::Spec->devnull , $cmd); replaced by run due to buffered output race conditions
			#while( <READER> ) {
			while(<F>) {
				chomp $_;
				$_ =~ s/^\s+|\s+$//g;
				next if $_=~/^#/;
				next if $_=~/^\s*$/;
				my @l = split /\s+/ , $_;
				my @gff3entry = &{$self->tool_parser}($self->tool,$kingdom,\@l);
				$types->{$gff3entry[2]} = 1;
				my @chr = ($abbr,$gff3entry[0],$gff3entry[2]);
				$chr[-1] = $self->tool.(++$uid->{join('.',@chr)});
				$gff3entry[0] = join('.',@chr);
				$self->gffdb->add_gff3_entry(\@gff3entry,$self->fastadb->get_gff3seq(\@gff3entry));
			}
			#waitpid($pid, 0);
			close F;
		}
	}

	for(keys %$types){
		$self->gffdb->merge($_,$self->tool); #merge multi kingdoms and overlapping annotations du to genome chunks
	}
}

1;