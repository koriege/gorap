package Bio::Gorap::Tool::Infernal;

use Moose; with 'Bio::Gorap::ToolI';

use POSIX qw(:sys_wait_h);				
use IPC::Open3;
use File::Spec::Functions;
use Symbol qw(gensym);

sub calc_features {
	my ($self) = @_;
	
	#calculations and software calls
	#results are fetched and stored in DB structure

	for (0..$#{$self->parameter->genomes}){		
		my $genome = ${$self->parameter->genomes}[$_];
		my $abbr = ${$self->parameter->abbreviations}[$_];
		my $uid = 0;				

		my @cmd = ('cmsearch' , '--noali' ,  '--cpu' , $self->threads , $self->parameter->cfg->cm , $genome);		
		my $pid = open3(gensym, \*READER, File::Spec->devnull , join(' ' , @cmd));
		
		while( <READER> ) {					
			chomp $_;
			$_ =~ s/^\s+|\s+$//g;
			next if $_=~/^#/;
			next if $_=~/^\s*$/;
			my @l = split /\s+/ , $_;
			next if $#l < 11;			
			
			#tool_parser is set to infernal_parser via Gorap.pl, static defined in Bio::Gorap::Functions::ToolParser			
			my @gff3entry = &{$self->tool_parser}(++$uid,$abbr,$self->parameter->cfg->rf_rna,\@l);						
			my $seq = $self->fastadb->get_gff3seq(\@gff3entry);			
			$self->gffdb->add_gff3_entry(\@gff3entry,$seq,$abbr);
		}
		waitpid($pid, 0);			
	}	
}