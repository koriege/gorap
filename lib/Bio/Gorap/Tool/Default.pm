package Bio::Gorap::Tool::Default;

use Moose; with 'Bio::Gorap::ToolI';

use IO::Select;
use IO::Pipe;
use POSIX qw(:sys_wait_h);
use IPC::Open3;
use IPC::Cmd qw(run);
use File::Spec::Functions;
use Symbol qw(gensym);
use List::MoreUtils qw(any);

sub calc_features {
	my ($self) = @_;

	#calculations and software calls
	#results are fetched and stored in DB structure
	if (grep {/\$cpus/} $self->cmd){
		for (0..$#{$self->parameter->genomes}){
			my $genome = ${$self->parameter->genomes}[$_];
			my $abbr = ${$self->parameter->abbreviations}[$_];
			my $uid = 0;
			my $tmpfile;
			my $threads = $self->threads;
			my $cmd = $self->cmd;
			$cmd =~ s/\$genome/$genome/;
			$cmd =~ s/\$cpus/$threads/;
			if ($cmd =~ /\$output/){
				$tmpfile = catfile($self->parameter->tmp,$self->parameter->pid.'.tmp');
				$cmd =~ s/\$output/$tmpfile/;
			}

			if ($tmpfile){
				my ($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => $cmd , verbose => 0 );
				open F , '<'.$tmpfile or next;
				while (<F>){
					chomp $_;
					$_ =~ s/^\s+|\s+$//g;
					next if $_=~/^#/;
					next if $_=~/^\s*$/;
					my @l = split /\s+/ , $_;
					#tool_parser is set to gff3_parser via Gorap.pl, static defined in Bio::Gorap::Functions::ToolParser
					my @gff3entry = &{$self->tool_parser}(++$uid,$abbr,$self->parameter->cfg->rf_rna,$self->tool,\@l);
					my $seq = $self->fastadb->get_gff3seq(\@gff3entry);
					$self->gffdb->add_gff3_entry(\@gff3entry,$seq);
				}
				close F;
			} else {
				my $pid = open3(gensym, \*READER, File::Spec->devnull , $cmd);
				while( <READER> ) {
					chomp $_;
					$_ =~ s/^\s+|\s+$//g;
					next if $_=~/^#/;
					next if $_=~/^\s*$/;
					my @l = split /\s+/ , $_;
					#tool_parser is set to gff3_parser via Gorap.pl, static defined in Bio::Gorap::Functions::ToolParser
					my @gff3entry = &{$self->tool_parser}(++$uid,$abbr,$self->parameter->cfg->rf_rna,$self->tool,\@l);
					my $seq = $self->fastadb->get_gff3seq(\@gff3entry);
					$self->gffdb->add_gff3_entry(\@gff3entry,$seq);
				}
				waitpid($pid, 0);
			}
		}
	} else {
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
				my $tmpfile;
				my $cmd = $self->cmd;
				$cmd =~ s/\$genome/$genome/;
				$cmd =~ s/\$cpus/$threads/;
				if ($cmd =~ /\$output/){
					$tmpfile = catfile($self->parameter->tmp,$self->parameter->pid.'.tmp');
					$cmd =~ s/\$output/$tmpfile/;
				}

				if ($tmpfile){
					my ($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => $cmd) , verbose => 0 );
					open F , '<'.$tmpfile or exit;
					while (<F>){
						chomp $_;
						$_ =~ s/^\s+|\s+$//g;
						next if $_=~/^#/;
						next if $_=~/^\s*$/;
						#tool_parser is set to gff3_parser via Gorap.pl, static defined in Bio::Gorap::Functions::ToolParser
						print $pipe $_."\n";
					}
					close F;
					unlink $tmpfile;
				} else {
					my $pid = open3(gensym, \*READER, File::Spec->devnull , $cmd);
					while( <READER> ) {
						chomp $_;
						$_ =~ s/^\s+|\s+$//g;
						next if $_=~/^#/;
						next if $_=~/^\s*$/;
						#tool_parser is set to gff3_parser via Gorap.pl, static defined in Bio::Gorap::Functions::ToolParser
						print $pipe $_."\n";
					}
					waitpid($pid, 0);
				}
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
		for(@out){
			my @gff3entry = split /\s+/, $_;

			($gff3entry[0], $gff3entry[3], $gff3entry[4]) = $self->fastadb->chunk_backmap($gff3entry[0], $gff3entry[3], $gff3entry[4]);
			my ($abbr,@orig) = split /\./ , $gff3entry[0];
			$gff3entry[0] = join '.' , @orig;
			@gff3entry = &{$self->tool_parser}(++$uid,$abbr,$self->parameter->cfg->rf_rna,$self->tool,\@gff3entry);

			my $seq = $self->fastadb->get_gff3seq(\@gff3entry);
			$self->gffdb->add_gff3_entry(\@gff3entry,$seq);
		}
	}
}

1;