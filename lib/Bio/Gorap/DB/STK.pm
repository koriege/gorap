package Bio::Gorap::DB::STK;

use Moose;
use Bio::AlignIO;
use File::Basename;
use File::Spec::Functions;
use Bio::Gorap::Functions::STK;
use IPC::Cmd qw(run);
use List::Util qw(any);
use POSIX;
use List::Util qw(min max);

#uses the gorap parameter object to initialize the
#database of Bio::Align::Stockholm  objects
has 'parameter' => (
	is => 'ro',
	isa => 'Bio::Gorap::Parameter',
	required => 1 ,
	trigger => \&_set_db
);

has 'taxonomy' => (
	is => 'ro',
	isa => 'Bio::Gorap::DB::Taxonomy',
	predicate => 'has_taxonomy'
);

#genome file based hashmap of Bio::Align::AlignI databases
has 'db' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
);

#gorap specific mapping
has 'idToPath' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
);

#set up database taking existing alignments into account
sub _set_db {
	my ($self) = @_;

	print "Reading STK files from ".$self->parameter->output."\n" if $self->parameter->verbose;

	for (glob catfile($self->parameter->output,'alignments','*.stk')){
		#use gorap specific access ids: the rfam id and rna name
		$self->add_stk(substr(basename($_),0,-4),$_);
	}
}

#for manually adding files into this hashed database:
sub add_stk {
	my ($self,$id,$file) = @_;
	$self->idToPath->{$id} = $file;
	my $io = Bio::AlignIO->new(-format => 'stockholm', -file => $file, -verbose => -1);

	$self->db->{$id} = $io->next_aln;
	$self->db->{$id}->set_displayname_flat;
}

sub write_stk {
	my ($self,$stk,$file) = @_;

	my $out = '';
	open my $READER, '>', \$out;
	(Bio::AlignIO->new(-format  => 'stockholm', -fh => $READER, -verbose => -1 ))->write_aln($stk);
	close $READER;

	open STK , '>'.$file or die $!;
	for (split /\n/ , $out){
		next if $_=~/^\s*#=GS/;
		print STK $_."\n";
	}
	close STK;
}

sub store {
	my ($self,$id) = @_;

	if ($id){
		$self->write_stk($self->db->{$id},$self->idToPath->{$id});
	} else {
		for my $k (keys %{$self->db}){
			my $stk = $self->db->{$k};
			$stk = $self->taxonomy->sort_stk($stk) if $self->has_taxonomy && $self->parameter->sort;
			$self->remove_gap_columns_and_write($stk,$self->idToPath->{$k});
		}
	}
}

sub store_stk {
	my ($self,$stk, $file, $taxdb) = @_;

	$stk = $taxdb->sort_stk($stk) if $taxdb && $self->parameter->sort;
	$self->write_stk($stk,$file);
}

sub remove_gap_columns_and_write {
	my ($self, $stk, $file) = @_;
	my $tmpstk = $stk;

	my ($ss , $cs) = Bio::Gorap::Functions::STK->get_ss_cs_from_object($stk);
	my @ss = split // , $ss;
	my @cs = split // , $cs;

	my $del = 0 ;
	$tmpstk->map_chars('\.','-');
	my @gapcolm = @{$stk->gap_col_matrix};
	for (my $i=$#gapcolm; $i>=0; $i-- ){
		my $absent = 0;
		for my $k (keys %{$gapcolm[$i]}){
			$absent += $gapcolm[$i]->{$k};
		}

		if($absent == $tmpstk->num_sequences){
			$del = 1;
			$tmpstk = $tmpstk->remove_columns([$i,$i]);
			splice(@cs, $i, 1);
			splice(@ss, $i, 1);
		}
	}

	$tmpstk->set_displayname_flat;

	if ($del){
		my $out = '';
		open my $READER, '>', \$out;
		(Bio::AlignIO->new(-format  => 'stockholm', -fh => $READER, -verbose => -1 ))->write_aln($tmpstk);
		close $READER;
		my @lines = split /\n/ , $out;
		$lines[-1] = '#=GC RF '.join('',(' ') x ($tmpstk->maxdisplayname_length()-7)).' '.join('' , @cs);
		push @lines , '#=GC SS_cons '.join('',(' ') x ($tmpstk->maxdisplayname_length()-12)).' '.join('' , @ss);
		push @lines, '//';

		open STK , '>'.$file or die $!;
		for (@lines){
			next if $_=~/^\s*#=GS/;
			print STK $_."\n";
		}
		close STK;

		my $retaln = (Bio::AlignIO->new(-file => $file, -format => 'stockholm'))->next_aln;
		$retaln->set_displayname_flat;
		return $retaln;
	} else {
		$self->write_stk($stk,$file);
		return $stk;
	}

}

#align annotated sequences given as referenced array of Bio::Seq objects against the
#covariance model of gorap or a given one. resulting alignments are merged with those
#read into this hashreferenced database of Bio::Align::Stockholm objects
sub align {
	my ($self,$id,$sequences,$threads,$cm,$scorefile,$skipmerge) = @_;

	$scorefile = catfile($self->parameter->output,'meta',$id.'.scores') unless $scorefile;
	my $tmpfile1 = catfile($self->parameter->tmp,$self->parameter->pid.'.1.stk');
	my $tmpfile2 = catfile($self->parameter->tmp,$self->parameter->pid.'.2.stk');
	my $stkfile = catfile($self->parameter->output,'meta',$id.'.stk');
	my $fastafile = catfile($self->parameter->output,'meta',$id.'.fa');

	#print new annotated sequences into fasta file for aligning it via cmalign
	open FA , '>'.$fastafile or die $!;
	print FA '>'.$_->display_id."\n".$_->seq."\n" for @{$sequences};
	close FA;

	my ($cmd1, $cmd2, $success, $error_code, $full_buf, $stdout_buf, $stderr_buf);
	#if database was initialized with existing alignment, the new sequences are aligned in single, to merge 2 alignment files afterwards

	if ($skipmerge){
		$cmd1 = "cmalign --mxsize ".$self->parameter->mem." --noprob --sfile $scorefile --cpu $threads -o $stkfile ".$self->parameter->cfg->cm." $fastafile";
		$cmd1 = "cmalign --mxsize ".$self->parameter->mem." --noprob --sfile $scorefile --cpu $threads -o $stkfile $cm $fastafile" if $cm;
		($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => $cmd1 , verbose => 0 );
		return $scorefile;
	}

	$cm = $self->parameter->cfg->cm.".eval.orig" if -e $self->parameter->cfg->cm.".eval.orig";
	
	$cmd1 = "cmalign --mxsize ".$self->parameter->mem." --noprob --sfile $scorefile --cpu $threads -o $tmpfile1 ".$self->parameter->cfg->cm." $fastafile";
	$cmd1 = "cmalign --mxsize ".$self->parameter->mem." --noprob --sfile $scorefile --cpu $threads -o $tmpfile1 $cm $fastafile" if $cm;
	($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => $cmd1 , verbose => 0 );

	if (exists $self->db->{$id}){
		$self->store_stk($self->db->{$id},$tmpfile2);

		$cmd2="esl-alimerge --rna -o $stkfile $tmpfile2 $tmpfile1";
		($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => $cmd2 , verbose => 0 );

		if($success){
			my ($ss , $cs) = Bio::Gorap::Functions::STK->get_ss_cs_from_file($stkfile);
			$success = 0 if ! $ss || ! $cs;
		}

		unless($success){
			open FA , '>'.$fastafile or die $!;
			for ($self->db->{$id}->each_seq){
				my $s = $_->seq;
				$s=~s/\W//g;
				$s=uc($s);
				print FA '>'.$_->id."\n".$s."\n";
			}
			print FA '>'.$_->display_id."\n".$_->seq."\n" for @{$sequences};
			close FA;

			$cmd1 = "cmalign --mxsize ".$self->parameter->mem." --noprob --sfile $scorefile --cpu $threads -o $stkfile ".$self->parameter->cfg->cm." $fastafile";
			$cmd1 = "cmalign --mxsize ".$self->parameter->mem." --noprob --sfile $scorefile --cpu $threads -o $stkfile $cm $fastafile" if $cm;
			($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => $cmd1 , verbose => 0 );
		}
	} else {
		$cmd2 = "esl-alimerge --rna -o $stkfile ".$self->parameter->cfg->stk." $tmpfile1";
		($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => $cmd2 , verbose => 0 );

		if($success){
			my ($ss , $cs) = Bio::Gorap::Functions::STK->get_ss_cs_from_file($stkfile);
			$success = 0 if ! $ss || ! $cs; #TODO write own RF based merger
		}
		
		if (! $success && $self->parameter->cfg->fasta){
			open F , '<'.$self->parameter->cfg->fasta or die $!;
			open FA , '>'.$fastafile or die $!;
			while (<F>){
				if($_=~/^>/){
					print FA '>'.(split/\s+/,$_)[1]."\n"
				} else {
					print FA $_ ;
				}
			}
			close F;
			print FA '>'.$_->display_id."\n".$_->seq."\n" for @{$sequences};
			close FA;

			$cmd1 = "cmalign --mxsize ".$self->parameter->mem." --noprob --sfile $scorefile --cpu $threads -o $stkfile ".$self->parameter->cfg->cm." $fastafile";
			$cmd1 = "cmalign --mxsize ".$self->parameter->mem." --noprob --sfile $scorefile --cpu $threads -o $stkfile $cm $fastafile" if $cm;
			($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => $cmd1 , verbose => 0 );
		}
	}

	$self->add_stk($id,$stkfile);

	return ($scorefile, $self->db->{$id});
}

sub calculate_threshold {
	my ($self,$cpus) = @_;

	my $threshold=999999;
	if ($self->parameter->cfg->bitscore == $self->parameter->cfg->bitscore_cm){ #user didnt change anything

		my $scorefile;
		if ($self->parameter->taxonomy && (scalar keys %{$self->taxonomy->relatedSpeciesIDsToLineage} > 0 || scalar keys %{$self->taxonomy->relatedRankIDsToLineage} > 0)){

			my $fasta = Bio::SeqIO->new(-file => $self->parameter->cfg->fasta , -format => 'Fasta', -verbose => -1);

			my ($inRank,$inSpecies);
			while ( my $s = $fasta->next_seq() ) {
				if ($s->id=~/^(\d+)/){
					push @{$inSpecies} , $s if exists $self->taxonomy->relatedSpeciesIDsToLineage->{$1};
					push @{$inRank} , $s if exists $self->taxonomy->relatedRankIDsToLineage->{$1};
				}
			}

			if ($#{$inSpecies} > -1 ){
				$scorefile = $self->align(
					$self->parameter->cfg->rf_rna.'.tax',
					$inSpecies,
					$cpus,
					$self->parameter->cfg->cm,
					catfile($self->parameter->output,'meta',$self->parameter->cfg->rf_rna.'.tax.score'),
					1
				);
			}elsif($#{$inRank} > -1){
				$scorefile = $self->align(
					$self->parameter->cfg->rf_rna.'.tax',
					$inRank,
					$cpus,
					$self->parameter->cfg->cm,
					catfile($self->parameter->output,'meta',$self->parameter->cfg->rf_rna.'.tax.score'),
					1
				);
			}
		} else {
			$threshold = $self->parameter->cfg->bitscore * $self->parameter->thfactor;
			return ($threshold,0);
		}

		if ($scorefile){
			open S , '<'.$scorefile or die $!;
			while(<S>){
				chomp $_ ;
				$_ =~ s/^\s+|\s+$//g;
				next if $_=~/^#/;
				my @l = split /\s+/ , $_;
				$threshold = $l[6] if $l[6] < $threshold;
			}
			close S;

			if ($threshold > $self->parameter->cfg->bitscore * $self->parameter->thfactor){
				$threshold = floor( ($threshold - ($threshold - $self->parameter->cfg->bitscore * $self->parameter->thfactor)/2) * $self->parameter->thfactor);
			} else {
				$threshold = floor($threshold);
			}
			#returns threshold and nonTaxThreshold
			return ($threshold,$self->parameter->cfg->bitscore * $self->parameter->thfactor);
		} else { #no related species in seed alignment 
			if ($self->parameter->cmtaxbiascutoff > 0){
				my $fasta = Bio::SeqIO->new(-file => $self->parameter->cfg->fasta , -format => 'Fasta', -verbose => -1);
				my $ancestors;
				my $seqc=0;
				while ( my $s = $fasta->next_seq() ) {
					$seqc++;
					if ($s->id=~/^(\d+)/){
						my $ancestor;
						my @nodes;
						push @nodes, $_ for @{$self->taxonomy->getLineageNodes($1)};
						if ($#nodes != -1){
							unless ($nodes[0] == 28384 || $nodes[0] == 12908){ #label: other sequence(s)
								$ancestors->{$nodes[-1]}++; #e.g Genus epithet strand -> ancestor = Genus epithet 
								$ancestors->{$1}++; #e.g. Genus epithet -> ancestor = Genus => add also $1
								$ancestors->{$nodes[-2]}++; #e.g Genus epithet strand -> ancestor = Genus epithet => add also -2
							}
						}
					}
				}
				my ($overrepresented) = reverse sort { $ancestors->{$a} <=> $ancestors->{$b} } keys %$ancestors;
				return (999999,0) if $overrepresented && ($ancestors->{$overrepresented} == $seqc || $ancestors->{$overrepresented} > max(1,$seqc*$self->parameter->cmtaxbiascutoff));
			}
			$threshold = $self->parameter->cfg->bitscore * $self->parameter->thfactor;
			return ($threshold,0);
		}
	} else {
		#use user bitscore
		return ($self->parameter->cfg->bitscore,0);
	}
}

#gorap post stk filters
#gets stk as Bio::SimpleAlign, HashRef of Bio::SeqFeatures and threshold as float
sub filter_stk {
	my ($self, $id, $stk, $features, $threshold, $nonTaxThreshold, $taxdb, $gffdb) = @_;

	$stk = $taxdb->sort_stk($stk) if $taxdb && $self->parameter->sort;

	my $c=0;
	#necessary for deleten to map $c->feature , because 2 features can have same identifier if called print $object
	#due to seqfeature doesnt return memory references
	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';

	my @update;
	my $up;
	my $write;

	if ( ! $self->parameter->nofilter && ! $self->parameter->nobutkingsnofilter){
		($stk, $features, $up, $write) = Bio::Gorap::Functions::STK->length_filter($stk, $features, $self->parameter->cfg->cs);
		push @update , @{$up} if $up;
		return @update if scalar keys %{$features} == 0;
		$stk = $self->remove_gap_columns_and_write($stk,catfile($self->parameter->output,'meta',$id.'.L.stk'));# if $write;
	}

	if ( ! $self->parameter->nofilter || $self->parameter->nobutkingsnofilter){
		if ($self->parameter->cfg->types=~/CD-box/ && $self->parameter->cfg->userfilter){
			$threshold = 20 * $self->parameter->thfactor;
		}
	}
	($stk, $features, $up, $write) = Bio::Gorap::Functions::STK->score_filter($stk, $features, $threshold, $nonTaxThreshold);
	push @update , @{$up} if $up;
	return @update if scalar keys %{$features} == 0;
	$stk = $self->remove_gap_columns_and_write($stk,catfile($self->parameter->output,'meta',$id.'.B.stk'));# if $write;

	if ( ! $self->parameter->nofilter && ! $self->parameter->nobutkingsnofilter){
		($stk, $features, $up, $write) = Bio::Gorap::Functions::STK->structure_filter($stk, $features, $self->parameter->cfg->types);
		push @update , @{$up} if $up;
		return @update if scalar keys %{$features} == 0;
		$stk = $self->remove_gap_columns_and_write($stk,catfile($self->parameter->output,'meta',$id.'.S.stk'));# if $write;
	}

	if ($self->parameter->cfg->userfilter){
		unless ($self->parameter->nofilter){
			($stk, $features, $up, $write) = Bio::Gorap::Functions::STK->user_filter($stk, $features, $self->parameter->cfg->constrains, $self->parameter->cfg->cs, $self->parameter->cfg->stk, $self->parameter->cfg->types);
			push @update , @{$up} if $up;
			return @update if scalar keys %{$features} == 0;
			$stk = $self->remove_gap_columns_and_write($stk,catfile($self->parameter->output,'meta',$id.'.P.stk'));# if $write;
		}
	} else {
		if ( ! $self->parameter->nofilter && ! $self->parameter->nobutkingsnofilter){
			($stk, $features, $up, $write) = Bio::Gorap::Functions::STK->sequence_filter($stk, $features);
			push @update , @{$up} if $up;
			return @update if scalar keys %{$features} == 0;
			$stk = $self->remove_gap_columns_and_write($stk,catfile($self->parameter->output,'meta',$id.'.P.stk'));# if $write;
		}
	}

	if ( ! $self->parameter->nofilter && ! $self->parameter->nobutkingsnofilter){
		($stk, $features, $up, $write) = Bio::Gorap::Functions::STK->copy_filter($stk, $features, $self->parameter->cfg->pseudogenes);
		push @update , @{$up} if $up;
		return @update if scalar keys %{$features} == 0;
		$stk = $self->remove_gap_columns_and_write($stk,catfile($self->parameter->output,'meta',$id.'.C.stk'));# if $write;

		if ($self->parameter->check_overlaps){
			($stk, $features, $up, $write) = Bio::Gorap::Functions::STK->overlap_filter($stk, $features, $gffdb);
			push @update , @{$up} if $up;
			return @update if scalar keys %{$features} == 0;
			$stk = $self->remove_gap_columns_and_write($stk,catfile($self->parameter->output,'meta',$id.'.O.stk'));# if $write;
		}
	}

	$stk = $self->remove_gap_columns_and_write($stk,catfile($self->parameter->output,'alignments',$id.'.stk'),$taxdb);

	return @update;
}

sub rm_seq_from_stk {
	my ($self, $features, $type, $flag, $stk) = @_;

	if (-e catfile($self->parameter->output,'alignments',$type.'.stk') || $stk){
		$features = [$features] unless ref($features) eq 'ARRAY';

		$stk = (Bio::AlignIO->new(-format  => 'stockholm', -file => catfile($self->parameter->output,'alignments',$type.'.stk'), -verbose => -1 ))->next_aln unless $stk;

		for (@$features){
			my $seq = $stk->get_seq_by_id($_->seq_id);
			$stk->remove_seq($seq) if $seq;
		}
		$stk = $self->remove_gap_columns_and_write($stk,catfile($self->parameter->output,'meta',$type.'.'.$flag.'.stk'));
		$self->store_stk($stk,catfile($self->parameter->output,'alignments',$type.'.stk'));
	}

}

1;
