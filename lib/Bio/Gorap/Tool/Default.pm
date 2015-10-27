package Bio::Gorap::Tool::Default;

use Moose; with 'Bio::Gorap::ToolI';

use IO::Select;
use IO::Pipe;
use POSIX qw(:sys_wait_h);				
use IPC::Open3;
use IPC::Cmd qw(run);
use File::Spec::Functions;
use Symbol qw(gensym);
use List::MoreUtils 'any';

sub calc_features {
	my ($self) = @_;
	
	#calculations and software calls
	#results are fetched and stored in DB structure
	if (any {/\$cpus/} @{$self->parameter->cfg->cmd}){
		for (0..$#{$self->parameter->genomes}){		
			my $genome = ${$self->parameter->genomes}[$_];
			my $abbr = ${$self->parameter->abbreviations}[$_];
			my $uid = 0;				
			my $tmpfile;
			my @cmd = @{$self->parameter->cfg->cmd};
			for (@cmd){
				$_ =~ s/\$genome/$genome/;
				$_ =~ s/\$cpus/$self->threads/;	
				if ($_ =~ /\$output/){
					$tmpfile = catfile($self->parameter->tmp,$self->parameter->pid.'.tmp');
					$_ =~ s/\$output/$tmpfile/;
				} 
			}

			if ($tmpfile){
				my ($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => join(' ' , @cmd) , verbose => 0 );
				open F , '<'.$tmpfile or next;
				while (<F>){
					chomp $_;
					$_ =~ s/^\s+|\s+$//g;
					next if $_=~/^#/;
					next if $_=~/^\s*$/;
					my @l = split /\s+/ , $_;
					#tool_parser is set to gff3_parser via Gorap.pl, static defined in Bio::Gorap::Functions::ToolParser	
					my @gff3entry = &{$self->tool_parser}(++$uid,$abbr,$self->parameter->cfg->rf_rna,\@l);						
					my $seq = $self->fastadb->get_gff3seq(\@gff3entry);
					$self->gffdb->add_gff3_entry(\@gff3entry,$seq,$abbr);
				}
				close F;
			} else {
				my $pid = open3(gensym, \*READER, File::Spec->devnull , join ' ' , @cmd);				
				while( <READER> ) {		
					chomp $_;
					$_ =~ s/^\s+|\s+$//g;
					next if $_=~/^#/;
					next if $_=~/^\s*$/;
					my @l = split /\s+/ , $_;			
					#tool_parser is set to gff3_parser via Gorap.pl, static defined in Bio::Gorap::Functions::ToolParser			
					my @gff3entry = &{$self->tool_parser}(++$uid,$abbr,$self->parameter->cfg->rf_rna,\@l);						
					my $seq = $self->fastadb->get_gff3seq(\@gff3entry);
					$self->gffdb->add_gff3_entry(\@gff3entry,$seq,$abbr);
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
				for (@{$self->parameter->cfg->cmd}){
					$_ =~ s/\$genome/$genome/;
					if ($_ =~ /\$output/){
						$tmpfile = catfile($self->parameter->tmp,$$.'.tmp');
						$_ =~ s/\$output/$tmpfile/;
					} 
				}		
				
				if ($tmpfile){
					my ($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => join(' ' , @{$self->parameter->cfg->cmd}) , verbose => 0 );
					open F , '<'.$tmpfile or exit;
					while (<F>){
						chomp $_;
						$_ =~ s/^\s+|\s+$//g;
						next if $_=~/^#/;
						next if $_=~/^\s*$/;
						#tool_parser is set to gff3_parser via Gorap.pl, static defined in Bio::Gorap::Functions::ToolParser	
						print $pipe $_;
					}
					close F;
					unlink $tmpfile;
				} else {
					my $pid = open3(gensym, \*READER, File::Spec->devnull , join(' ' , @{$self->parameter->cfg->cmd}));
					while( <READER> ) {		
						chomp $_;
						$_ =~ s/^\s+|\s+$//g;
						next if $_=~/^#/;
						next if $_=~/^\s*$/;
						#tool_parser is set to gff3_parser via Gorap.pl, static defined in Bio::Gorap::Functions::ToolParser	
						print $pipe $_;
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
			my @l = split /\s+/, $_;				
			for ($self->fastadb->chunk_backmap($l[0], $l[3], $l[4])){
				($l[0],$l[3],$l[4]) = @{$_};
																
				my ($abbr, @header) = split /\./,$l[0];
				$l[0] = join '.' , @header;
				$uid->{$l[0]}++;
				
				my @gff3entry = &{$self->tool_parser}($uid->{$l[0]},$abbr,$self->parameter->cfg->rf_rna,\@l);
				#due to overlapping chunks check for already annotated genes
				my $existingFeatures = $self->gffdb->get_overlapping_features(\@gff3entry,$abbr);
				next if $#{$existingFeatures} > -1;

				my $seq = $self->fastadb->get_gff3seq(\@gff3entry);
				$self->gffdb->add_gff3_entry(\@gff3entry,$seq,$abbr);			
			}
		}
	}
}