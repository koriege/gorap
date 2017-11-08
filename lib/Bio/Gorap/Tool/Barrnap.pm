package Bio::Gorap::Tool::Barrnap;

use Moose; with 'Bio::Gorap::ToolI';
use IO::Select;
use IO::Pipe;	
use IPC::Open3;
use File::Spec::Functions;
use Symbol qw(gensym);

sub calc_features {
	my ($self) = @_;

	for (0..$#{$self->parameter->genomes}){
		my $genome = ${$self->parameter->genomes}[$_];
		my $abbr = ${$self->parameter->abbreviations}[$_];
		#skip redundand calculations
		my @f = $self->gffdb->db->{$abbr}->features(-attributes => {source => 'GORAP'.$self->tool});
		next if $#f > -1;

		#remove all other previously families annotated by this tool
		my @f = $self->gffdb->db->{$abbr}->features(-attributes => {source => $self->tool}); 
		for (@f){
			my $rfrna = $_->primary_tag;
			$self->stkdb->db->{$rfrna}->remove_seq($_) for $self->stkdb->db->{$rfrna}->get_seq_by_id($f->seq_id);
			$self->gffdb->db->{$abbr}->delete($f);
		}

		my $uid;
		my $threads = $self->threads;

		my @kingdoms;
		push @kingdoms , 'arc' if exists $self->parameter->kingdoms->{'arc'};
		push @kingdoms , 'bac' if exists $self->parameter->kingdoms->{'bac'};
		push @kingdoms , 'euk' if exists $self->parameter->kingdoms->{'fungi'} || exists $self->parameter->kingdoms->{'euk'};

		for my $kingdom (@kingdoms){
			my @cmd = @{$self->parameter->cfg->cmd}; 
			for (@cmd){
				$_ =~ s/\$genome/$genome/;
				$_ =~ s/\$cpus/$threads/;
				$_ =~ s/\$kingdom/$kingdom/;
			}
			my $pid = open3(gensym, \*READER, File::Spec->devnull , join(' ' , @cmd));
			
			while( <READER> ) {
				chomp $_;
				$_ =~ s/^\s+|\s+$//g;
				next if $_=~/^#/;
				next if $_=~/^\s*$/;
				my @l = split /\s+/ , $_;
				#tool_parser is set to infernal_parser via Gorap.pl, static defined in Bio::Gorap::Functions::ToolParser			
				
				my @gff3entry = &{$self->tool_parser}($kingdom,\@l);
				$uid->{$abbr.'.'.$gff3entry[2]}++;
				$gff3entry[0] = join('.',($abbr,$gff3entry[0],$uid->{$abbr.'.'.$gff3entry[2]}));
				$self->gffdb->add_gff3_entry(\@gff3entry,$self->fastadb->get_gff3seq(\@gff3entry));
			}
			waitpid($pid, 0);
		}
	}
}

1;