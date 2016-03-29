package Bio::Gorap::Tool::Infernal;

use Moose; with 'Bio::Gorap::ToolI';

use POSIX qw(:sys_wait_h);				
use IPC::Open3;
use File::Spec::Functions;
use Symbol qw(gensym);
use List::Util qw(max);

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
			
			unless ($self->parameter->nofilter){
				next if $self->parameter->cfg->cs && $gff3entry[4]-$gff3entry[3] < length($self->parameter->cfg->cs)/2.5;

				if ($self->parameter->cfg->rf_rna=~/_mir/i || $self->parameter->cfg->rf_rna=~/_Afu/ || $self->parameter->cfg->rf_rna=~/_SNOR/ || $self->parameter->cfg->rf_rna=~/(-|_)sn?o?s?n?o?[A-WYZ]+[a-z]?-?\d/){
					my $existingFeatures = $self->gffdb->get_all_overlapping_features(\@gff3entry);
					my $snover=0;
					my $exscore = -999999;
					my @rmfeatures;
					for my $f (@{$existingFeatures}){
						if ($f->type=~/_mir/i || $f->type=~/_Afu/ || $f->type=~/_SNOR/ || $f->type=~/(-|_)sn?o?s?n?o?[A-WYZ]+[a-z]?-?\d/){
							$exscore = max($exscore,($f->get_tag_values('origscore'))[0]);
							push @rmfeatures , $f;
						}
					}
					if ($exscore < $gff3entry[5]){
						$self->gffdb->update_filter($_->seq_id,$_->primary_tag,"O") for @rmfeatures;
					} else {
						$uid--;
						next;
					}
				}
			}
			my $seq = $self->fastadb->get_gff3seq(\@gff3entry);
			$self->gffdb->add_gff3_entry(\@gff3entry,$seq,$abbr);
		}
		waitpid($pid, 0);
	}	
}