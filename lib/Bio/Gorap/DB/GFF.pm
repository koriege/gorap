package Bio::Gorap::DB::GFF;

use Moose;
use Bio::DB::SeqFeature::Store;
use Bio::SeqFeature::Generic;
use File::Spec::Functions;
use File::Basename;
use List::MoreUtils qw(any first_index);
use List::Util qw(min max);

#uses the gorap parameter object to initialize the
#database/hashmap of Bio::DB::SeqFeature objects
has 'parameter' => (
	is => 'ro',
	isa => 'Bio::Gorap::Parameter',
	required => 1,
	trigger => \&_set_db
);

#genome file based hashmap of Bio::DB::SeqFeature databases
has 'db' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
);

has 'rnas' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
);

has 'userdb' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
);

has 'bamdb' => (
	is => 'ro',
	isa => 'Bio::Gorap::DB::BAM',
	required => 1
);

has 'fastadb' => (
	is => 'ro',
	isa => 'Bio::Gorap::DB::Fasta',
	required => 1
);
#add a gff3 entry into this Bio::DB::SeqFeature database
sub add_gff3_entry {
	my ($self,$s,$seq,$flag) = @_;

	#id consists of abbreviation.original.copy
	if ($#{$s} > 6){
		my ($id,$source,$type,$start,$stop,$score,$strand,$phase,$attributes) = @{$s};
		my ($abbr,@orig) = split /\./ , $id;
		pop @orig;

		my @overlaps = ('.');
		my ($tpm,$reads,$filter,$notes) = ('.','.','!','.');

		$attributes.=';';
		$tpm = $1 ? $1 : '.' if $attributes=~/TPM=(.+?);/;
		$reads = $1 ? $1 : '.' if $attributes=~/Reads=(.+?);/;
		$filter = $1 ? $1 : '!' if $attributes=~/Filter=(.+?);/;
		$notes = $1 ? $1 : '!' if $attributes=~/Notes?=(.+?);/;
		@overlaps = split /,/ , $1 if $attributes=~/Overlaps=(.+?);/;

		$filter = $flag if $flag;

		$strand = $strand eq '+' ? 1 : $strand eq '-' ? -1 : 0;

		$self->rnas->{$type} = 1;
		$self->db->{$abbr}->new_feature(
			-start => $start,
			-stop => $stop,
			-strand => $strand,
			-seq_id => $id,
			-primary_tag => $type,
			-phase => $phase,
			#misused for updateable filter tag
			-display_name => $filter,
			-source => $source,
			#updateable score from search -> align
			-score => $score,
			-index => 1,
			-attributes => {
				chr => join('.',@orig),
				tpm => $tpm,
				reads =>  $reads,
				overlaps => \@overlaps,
				seq => $seq,
				source => $source,
				#search score
				origscore => $score,
				notes => $notes
			}
		);
	}
}

#gorap specific update filter tag into this Bio::DB::SeqFeature database
sub update_filter {
	my ($self,$id,$type,$filter) = @_;

	my $abbr = (split /\./ , $id)[0];
	return unless exists $self->db->{$abbr};
	my ($feature) = $self->db->{$abbr}->features(-seq_id => $id , -primary_tag => $type );

	$feature->display_name($filter);
	$feature->update;

	return $type;
}

#gorap specific update filter tag into this Bio::DB::SeqFeature database
sub update_score {
	my ($self,$s) = @_;

	my ($id,$type,$score) = split /\s+/ , $s;
	my $abbr = (split /\./ , $id)[0];

	return unless exists $self->db->{$abbr};
	my ($feature) = $self->db->{$abbr}->features(-seq_id => $id , -primary_tag => $type );

	return unless $feature;
	$feature->score($score);

	$feature->update;
}

#gorap specific update score into this Bio::DB::SeqFeature database
sub update_score_by_file {
	my ($self,$rf_rna,$file) = @_;

	open S , '<'.$file or die $!;
	while(<S>){
		chomp $_ ;
		$_ =~ s/^\s+|\s+$//g;
		next if $_=~/^#/;
		my @l = split /\s+/ , $_;
		$self->update_score(join(' ',($l[1],$rf_rna,$l[6])));
	}
	close S;
}

sub get_all_sequences {
	my ($self,$type) = @_;

	return $self->get_sequences($type,[keys %{$self->db}]);
}

sub get_sequences {
	my ($self,$type,$abbreviations) = @_;

	my @sequences;
	for my $abbr (@{$abbreviations}){
		push @sequences , Bio::Seq->new( -display_id => $_->seq_id , -seq => ($_->get_tag_values('seq'))[0] , -verbose => -1) for $self->db->{$abbr}->features(-primary_tag => $type);
	}
	return \@sequences;
}

sub get_seq {
	my ($self,$id,$type,$source,$abbr) = @_;

	return (($self->db->{$abbr}->features(-seq_id => $id , -primary_tag => $type , -attributes => {source => $source}))[0]->get_tag_values('seq'))[0];
}

sub get_features {
	my ($self,$type,$abbreviations,$filter) = @_;

	my @features;
	if ($filter) {
		for (@{$abbreviations}){
			next unless exists $self->db->{$_};
			if ($type){
				for ($self->db->{$_}->features(-primary_tag => $type)){
					push @features , $_ if $_->display_name eq $filter;
				}
			} else {
				for ($self->db->{$_}->features){
					push @features , $_ if $_->display_name eq $filter;
				}
			}
		}
	} else {
		if ($type){
			push @features , $self->db->{$_}->features(-primary_tag => $type) for @{$abbreviations};
		} else {
			push @features , $self->db->{$_}->features for @{$abbreviations};
		}
	}
	return \@features;
}

sub get_all_features {
	my ($self,$type,$filter) = @_;

	return $self->get_features($type,[keys %{$self->db}],$filter);
}

sub get_filtered_features {
	my ($self,$filter) = @_;

	my $types;
	for (keys %{$self->db}){
		for ($self->db->{$_}->features()){
			push @{$types->{$_->primary_tag}} , $_ if $_->display_name eq $filter;
		}
	}

	return $types;
}

sub get_feature {
	my ($self,$id,$type,$source,$abbr) = @_;

	my ($f) = $self->db->{$abbr}->features(-seq_id => $id , -primary_tag => $type);

	return $f;
}

sub get_all_overlapping_features { #bad name but returns all features (type independent) overlapping at least 50% in both direction
	my ($self,$s) = @_;

	my ($f,$id,$source,$type,$start,$stop,$score,$strand,$phase,$attributes);

	if (ref($s) eq 'ARRAY'){
	#id consists of abbreviation.ori.g.inal.copy
		($id,$source,$type,$start,$stop,$score,$strand,$phase,$attributes) = @{$s};
		$strand = $strand eq '+' ? 1 : $strand eq '-' ? -1 : 0;
		$f = Bio::SeqFeature::Generic->new( -start => $start, -end => $stop, -strand => $strand);
	} else {
		($id,$type,$source,$start,$stop,$strand) = ($s->seq_id, $s->primary_tag, $s->source, $s->start, $s->stop, $s->strand);
		$f = $s;
	}
	my ($abbr,@orig) = split /\./ , $id;
	pop @orig;
	$id = join '.' , ($abbr,@orig);

	my @features;

	for ($self->db->{$abbr}->features()){
		my @tmp = split /\./, $_->seq_id;
		pop @tmp;
		next unless join('.' , @tmp) eq $id;
		next unless $_->strand == $strand;
		next unless $_->display_name eq '!';
		next if $_->primary_tag eq $type;

		my ($istart, $istop, $istrand) = $f->intersection($_);
		push @features , $_ if $istart && $istop - $istart >= ($stop - $start) * 0.5 && $istop - $istart >= ($_->stop - $_->start) * 0.5;
	}

	return \@features;
}

sub get_user_overlaps {
	my ($self,$id,$start,$stop,$strand) = @_;

	my ($abbr,@orig) = split /\./ , $id;
	my $orig = join '.' , @orig;

	my @features = $self->userdb->{$abbr}->features(-seq_id => $orig, -start => $start , -strand => $strand, -stop => $stop , -range_type => 'overlaps');
	return \@features;
}

#gorap specific set up databases for all genome files from parameter object and
#fill in existing data in defined output directory
sub _set_db {
	my ($self) = @_;

	for my $abbr (@{$self->parameter->abbreviations}){
		$self->db->{$abbr} = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -verbose => -1 );
		$self->db->{$abbr}->init_database([1]);
		$self->_read($abbr);
	}

	return unless $self->parameter->has_gffs;
	print "Reading GFF files\n" if $self->parameter->verbose;

	my $idc=0;
	for (0..$#{$self->parameter->gffs}){
		my $abbr = ${$self->parameter->abbreviations}[$_];
		for my $f (@{${$self->parameter->gffs}[$_]}){
			$self->userdb->{$abbr} = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -verbose => -1 ) unless $self->userdb->{$abbr};
			print "$f\n" if -e $f && $self->parameter->verbose;
			open GFF, '<'.$f or die $!;
			$self->add_user_entry($_,$abbr,(++$idc)) while <GFF>;
			close GFF;
		}
	}
}

sub add_user_entry {
	my ($self,$line,$abbr,$idc) = @_;

	chomp $line;
	return unless $line || $line=~/^#/ || $line=~/^\s*$/;
	my @line = split /\t+/ , $line;

	return if $#line < 7;
	return if $line[2]=~/region/i;

	push @line , '.' if $#line == 7;

	$line[8].=';';
	$line[8]=~/ID=(.+?);/i;

	my $id;
	my $attr;
	if ($1){
		$id = $1;
		$attr = substr($line[8],0,$-[0]).substr($line[8],$+[0],-1);
	} else{
		$id = 'feature'.$idc;
		if ($line[8] eq '.;'){
			$attr = '';
		} else {
			$attr = length($line[8]) < 3 ? 'Note='.substr($line[8],0,-1) : substr($line[8],0,-1);
		}
	}

	return if $attr=~/Note=(\?)/ || $attr=~/Filter=(\?)/;

	my ($tpm,$reads) = ('.','.');

	$self->userdb->{$abbr}->new_feature(
		-start => $line[3],
		-stop => $line[4],
		-seq_id => $line[0],
		-display_name => '.',
		-strand => $line[6],
		-primary_tag => $line[2],
		-source => $line[1],
		-score => $line[5],
		-index => 1,
		-attributes => {
			ID => $id,
			attr => $attr,
			tpm => $tpm,
			reads => $reads
		}
	);
}

sub add_db {
	my ($self,$abbr) = @_;

	unless (exists $self->db->{$abbr}){
		$self->db->{$abbr} = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -verbose => -1 );
		$self->db->{$abbr}->init_database([1]);
	}
}

sub update_attribute {
	my ($self,$f,$attribute,$value) = @_;

	if (ref($f) eq 'ARRAY'){
		my ($abbr) = split /\./, $$f[0];
		($f) = $self->db->{$abbr}->features(-seq_id => $$f[0] , -primary_tag => $$f[2]);
	}
	
	$f->remove_tag($attribute);
	$f->add_tag_value($attribute,$value);
	$f->update;
}


sub get_features_by_source {
	my ($self,$abbr,$source,$type) = @_;

	my @features;
	if ($type){
		@features = $self->db->{$abbr}->features(-primary_tag => $type , -attributes => {source => $source});
	}  else {
		@features = $self->db->{$abbr}->features(-attributes => {source => $source});
	}

	return \@features;
}

#gorap way to fill gff database with existing data in defined output directory
sub _read {
	my ($self,$abbr) = @_;

	print "Reading previous annotations for $abbr\n" if $self->parameter->verbose;

	my @fastas = (catfile($self->parameter->output,'annotations',$abbr.'.fa'), catfile($self->parameter->output,'annotations',$abbr.'.passed.fa'));
	my @gffs = (catfile($self->parameter->output,'annotations',$abbr.'.gff'),catfile($self->parameter->output,'annotations',$abbr.'.passed.gff'));

	my $headerMapSeq={};
	for my $file (@fastas){
		next unless -e $file;
		open FA , '<'.$file or die $!;
		my $seqid;
		my $seq;
		while(<FA>){
			chomp $_;
			if ($_=~/^>\s*(\S+)/){
				#temp store of sequences
				$headerMapSeq->{$seqid}=$seq if $seqid;
				$seqid = $1;
				$seq = '';
			} else {
				$seq .= $_;
			}
		}
		close FA;
		$headerMapSeq->{$seqid}=$seq if $seqid;
	}

	for my $file (@gffs){
		next unless -e $file;
		open GFF , '<'.$file or die $!;
		while(<GFF>){
			chomp $_;
			my @l = split /\s+/ , $_;
			next unless $l[1];
			#build header like for fasta: abbreviation.original.type.tool.copy out of abbreviation.original.copy
			my @id = split /\./ , $l[0];
			my $abbr = shift @id;
			my $copy = pop @id;
			my $orig = join '.',@id;
			my $id = join '.' , ($abbr,$orig,$l[2],$l[1],$copy);
			$l[1] =~ s/\W//g;
			$l[1] = lc $l[1];
			if (exists $headerMapSeq->{$id}){
				$self->add_gff3_entry(\@l,$headerMapSeq->{$id});
			} else {
				$self->add_gff3_entry(\@l,$self->fastadb->get_gff3seq(\@l));
			}
		}
		close GFF;
	}
}

sub store {
	my ($self,$type) = @_;

	for my $abbr (keys %{$self->db}){
		my @features;
		if ($type){
			open GFF , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.gff') or die $!;
			open GFFO , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.orig.gff') or die $!;
			open FA , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.fa') or die $!;
			open FAO , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.orig.fa') or die $!;
			open GFFP , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.passed.gff') or die $!;
			open GFFPO , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.passed.orig.gff') or die $!;
			open FAP , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.passed.fa') or die $!;
			open FAPO , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.passed.orig.fa') or die $!;
			@features = $self->db->{$abbr}->features(-primary_tag => $type);
		} else {
			open GFF , '>'.catfile($self->parameter->output,'annotations',$abbr.'.gff') or die $!;
			open GFFO , '>'.catfile($self->parameter->output,'annotations',$abbr.'.orig.gff') or die $!;
			open FA , '>'.catfile($self->parameter->output,'annotations',$abbr.'.fa') or die $!;
			open FAO , '>'.catfile($self->parameter->output,'annotations',$abbr.'.orig.fa') or die $!;
			open GFFP , '>'.catfile($self->parameter->output,'annotations',$abbr.'.passed.gff') or die $!;
			open GFFPO , '>'.catfile($self->parameter->output,'annotations',$abbr.'.passed.orig.gff') or die $!;
			open FAP , '>'.catfile($self->parameter->output,'annotations',$abbr.'.passed.fa') or die $!;
			open FAPO , '>'.catfile($self->parameter->output,'annotations',$abbr.'.passed.orig.fa') or die $!;
			@features = $self->db->{$abbr}->features();
		}
		for (@features){
			my @id = split /\./,$_->seq_id;
			my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
			if ( $_->display_name eq '!' ){
				print GFFP $self->to_gff3_string($_);
				print FAP $self->to_fasta_string($_);
				print GFFPO $self->to_gff3_string($_,{},$orig);
				print FAPO $self->to_fasta_string($_,join('.',($orig,$copy)));
			} else {
				print GFF $self->to_gff3_string($_);
				print FA $self->to_fasta_string($_);
				print GFFO $self->to_gff3_string($_,{},$orig);
				print FAO $self->to_fasta_string($_,join('.',($orig,$copy)));
			}
		}
		close GFF;
		close FA;
		close GFFP;
		close FAP;
		close GFFO;
		close FAO;
		close GFFPO;
		close FAPO;
	}
}

sub get_user_entries {
	my ($self,$abbr) = @_;

	my @r;
	for my $f ($self->userdb->{$abbr}->features()){
		my @out = split /\s+/ , $f->gff_string;
		splice @out , 8 , 2;
		my @attributes;
		for (sort {$a cmp $b} $f->get_all_tags()){
			my ($v) = $f->get_tag_values($_);
			push @attributes , $_ eq 'attr' ? $v: $_."=".$v if $v;
		}
		push @out , (join ';' , @attributes);

		push @r, join "\t" , @out;
	}

	return @r;
}

sub merge {
	my ($self, $type, $tool) = @_;
	$tool = 'GORAP'.lc($tool);
	for my $abbr (@{$self->parameter->abbreviations}){
		my @features;
		if ($tool=~/blast/ || $tool=~/infernal/){
			@features = sort {($a->get_tag_values('chr'))[0] cmp ($b->get_tag_values('chr'))[0] || $a->strand <=> $b->strand || $a->start <=> $b->start || $a->stop <=> $b->stop} ($self->db->{$abbr}->features(-primary_tag => $type , -attributes => {source => 'GORAPinfernal'}),$self->db->{$abbr}->features(-primary_tag => $type , -attributes => {source => 'GORAPblast'}));
		} else {
			@features = sort {($a->get_tag_values('chr'))[0] cmp ($b->get_tag_values('chr'))[0] || $a->strand <=> $b->strand || $a->start <=> $b->start || $a->stop <=> $b->stop} $self->db->{$abbr}->features(-primary_tag => $type , -attributes => {source => $tool});
		}
		my $toremove;
		for (my $i=0; $i<=$#features; $i++) {
			next if exists $toremove->{$i};
			my $score = $features[$i]->source=~/blast/ ? 0 : $features[$i]->score;
			my $j = $i+1;
			while(defined $features[$j] && ($features[$i]->get_tag_values('chr'))[0] eq ($features[$j]->get_tag_values('chr'))[0] && $features[$i]->strand == $features[$j]->strand && $features[$i]->stop >= $features[$j]->start){
				#merge overlapping annotations of same type, due to multi tool screening or hits in overlapping regions (see fasta split);
				my $thisscore = $features[$j]->source=~/blast/ ? 0 : $features[$j]->score;
				if ($thisscore > $score){
					$toremove->{$i} = 1;
				} else {
					$toremove->{$j} = 1;
				}
				$j++;
			}
		}
		my @rm;
		push @rm, $features[$_] for keys %$toremove;
		$self->db->{$abbr}->delete(@rm);
	}
}

sub store_overlaps { #final store function which includes overlapping infos and tpms
	my ($self) = @_;

	for my $abbr (keys %{$self->db}){
		open GFF , '>'.catfile($self->parameter->output,'annotations',$abbr.'.gff') or die $!;
		open GFFO , '>'.catfile($self->parameter->output,'annotations',$abbr.'.orig.gff') or die $!;
		open FA , '>'.catfile($self->parameter->output,'annotations',$abbr.'.fa') or die $!;
		open FAO , '>'.catfile($self->parameter->output,'annotations',$abbr.'.orig.fa') or die $!;
		open GFFP , '>'.catfile($self->parameter->output,'annotations',$abbr.'.passed.gff') or die $!;
		open GFFPO , '>'.catfile($self->parameter->output,'annotations',$abbr.'.passed.orig.gff') or die $!;
		open FAP , '>'.catfile($self->parameter->output,'annotations',$abbr.'.passed.fa') or die $!;
		open FAPO , '>'.catfile($self->parameter->output,'annotations',$abbr.'.passed.orig.fa') or die $!;

		my @features = sort {($a->get_tag_values('chr'))[0] cmp ($b->get_tag_values('chr'))[0] || $a->strand <=> $b->strand || $a->start <=> $b->start || $a->stop <=> $b->stop} $self->db->{$abbr}->features();

		$self->bamdb->calculate_tpm(\@features) if ! $self->parameter->notpm && $self->parameter->has_bams;

		my $overlaps;
		my $overlapspassed;
		for (my $i=0; $i<=$#features; $i++) {
			my $stop = $features[$i]->stop;
			my $j = $i+1;
			while(defined $features[$j] && ($features[$i]->get_tag_values('chr'))[0] eq ($features[$j]->get_tag_values('chr'))[0] && $features[$i]->strand == $features[$j]->strand && $features[$i]->stop >= $features[$j]->start){
				if ($features[$i]->display_name eq '!' && $features[$j]->display_name eq '!'){
					$overlapspassed->{$features[$i]->seq_id.'.'.$features[$i]->primary_tag}->{$features[$j]->primary_tag} = 1;
					$overlapspassed->{$features[$j]->seq_id.'.'.$features[$j]->primary_tag}->{$features[$i]->primary_tag} = 1;
				} elsif ($features[$i]->display_name ne '!' && $features[$j]->display_name ne '!'){
					$overlaps->{$features[$i]->seq_id.'.'.$features[$i]->primary_tag}->{$features[$j]->primary_tag} = 1;
					$overlaps->{$features[$j]->seq_id.'.'.$features[$j]->primary_tag}->{$features[$i]->primary_tag} = 1;
				}
				$j++;
			}
			my @orig = split /\./,$features[$i]->seq_id;
			shift @orig;
			my $copy = pop @orig;
			if(exists $self->userdb->{$abbr}){
				for($self->userdb->{$abbr}->features(-seq_id => join('.',@orig), -start => $features[$i]->start , -strand => $features[$i]->strand, -stop => $features[$i]->stop , -range_type => 'overlaps')){
					my ($userid) = $_->get_tag_values('ID');
					next unless $userid;
					$overlapspassed->{$features[$i]->seq_id}->{$userid} = 1;
					$overlaps->{$features[$i]->seq_id}->{$userid} = 1;
				}
			}				
			if ($features[$i]->display_name eq '!'){
				print GFFP $self->to_gff3_string($features[$i],$overlapspassed);
				print FAP $self->to_fasta_string($features[$i]);
				print GFFPO $self->to_gff3_string($features[$i],$overlapspassed,join('.',@orig));
				print FAPO $self->to_fasta_string($features[$i],join('.',(@orig,$copy)));
			} else {
				print GFF $self->to_gff3_string($features[$i],$overlaps);
				print FA $self->to_fasta_string($features[$i]);
				print GFFO $self->to_gff3_string($features[$i],$overlaps,join('.',@orig));
				print FAO $self->to_fasta_string($features[$i],join('.',(@orig,$copy)));
			}
		}
	}
}

sub to_gff3_string {
	my ($self, $feature, $overlaps, $chr) = @_;
	$chr = $feature->seq_id unless $chr;

	my $o = join(',',keys %{$overlaps->{$feature->seq_id.'.'.$feature->primary_tag}});
	$o='.' unless $o;

	my $source = $feature->source;
	$source =~ s/GORAP//;

	my ($reads) = $feature->get_tag_values('reads');
	my ($tpm) = $feature->get_tag_values('tpm');
	$tpm = sprintf("%.2f", $tpm) unless $tpm eq '.';

	my $score = $feature->source=~/blast/ || $feature->source=~/denovo/ ? $feature->score : max($feature->score,($feature->get_tag_values('origscore'))[0]);

	my $attributes = 'TPM='.$tpm.';Reads='.$reads.';Filter='.$feature->display_name.';Note='.($feature->get_tag_values('notes'))[0].';Overlaps='.$o;

	return $chr."\t".$source."\t".$feature->primary_tag."\t".$feature->start."\t".$feature->stop."\t".$score."\t",$feature->strand ? $feature->strand > 0 ? '+' : '-' : '.',"\t".$feature->phase."\t".$attributes."\n";
}

sub to_fasta_string {
	my ($self, $feature, $chr) = @_;
	$chr = $feature->seq_id unless $chr;

	my $source = $feature->source;
	$source =~ s/GORAP//;

	my $seq = $feature->get_tag_values('seq');
	my @chr = split /\./,$chr;
	my $copy = pop @chr;
	
	return $seq ? '>'.join('.',(@chr,$feature->primary_tag,$source,$copy))."\n".$seq."\n" : '';
}

1;