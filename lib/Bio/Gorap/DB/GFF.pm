package Bio::Gorap::DB::GFF;

use Moose;
use Bio::DB::SeqFeature::Store;
use Bio::SeqFeature::Generic;
use File::Spec::Functions;
use File::Basename;
use List::MoreUtils qw(any first_index);

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

has 'userdb' => (
	is => 'rw',
    isa => 'HashRef',
    default => sub { {} }
);

has 'bamdb' => (
	is => 'ro',
	isa => 'Bio::Gorap::DB::BAM'
);

#add a gff3 entry into this Bio::DB::SeqFeature database
sub add_gff3_entry {	
	my ($self,$s,$seq) = @_;	

	#id consists of abbreviation.original.copy	
	if ($#{$s} > 6){
		my ($id,$source,$type,$start,$stop,$score,$strand,$phase,$attributes) = @{$s};		
		my ($abbr,@orig) = split /\./ , $id;
		pop @orig;
		
		my @overlaps = ('.');
		my ($tpm,$rpkm,$reads,$filter,$notes) = ('.','.','.','!','.');
		
		$attributes.=';';
		$tpm = $1 ? $1 : '.' if $attributes=~/TPM=(.+?);/;
		$rpkm = $1 ? $1 : '.' if $attributes=~/\wPKM=(.+?);/;
		$reads = $1 ? $1 : '.' if $attributes=~/Reads=(.+?);/;
		$filter = $1 ? $1 : '!' if $attributes=~/Filter=(.+?);/;
		$notes = $1 ? $1 : '!' if $attributes=~/Notes?=(.+?);/;
		@overlaps = split /,/ , $1 if $attributes=~/Overlaps=(.+?);/;

		if (! $self->parameter->notpm && $self->parameter->has_bams){
			($tpm,$rpkm,$reads) = $self->bamdb->rpkm($abbr,join('.',@orig),$start,$stop,$strand);
		}

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
	        	rpkm => $rpkm, 
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
	return unless $feature;
	$feature->display_name($filter);
	$feature->update;

	return $type;
	# $self->db->{$abbr}->store($feature);
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
		&update_score($self,join(' ',($l[1],$rf_rna,$l[6])));
	}
	close S;		
}

sub get_all_sequences {
	my ($self,$type) = @_;	

	return &get_sequences($self,$type,[keys %{$self->db}]);
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
			for ($self->db->{$_}->features(-primary_tag => $type)){
				push @features , $_ if $_->display_name eq $filter;
			}			
		}		
	} else {
		push @features , $self->db->{$_}->features(-primary_tag => $type) for @{$abbreviations};
	}	
	return \@features;
}

sub get_all_features {
	my ($self,$type,$filter) = @_;

	return &get_features($self,$type,[keys %{$self->db}],$filter);
}

sub get_old_features {
	my ($self,$type,$abbreviations) = @_;

	my @features;	
	for ([keys %{$self->db}]){		
		next if any { /$_/ } @{$abbreviations};
		push @features , $self->db->{$_}->features(-primary_tag => $type);
	}

	return \@features;
}

sub get_feature {
	my ($self,$id,$type,$source,$abbr) = @_;
	
	my @features = $self->db->{$abbr}->features(-seq_id => $id , -primary_tag => $type , -attributes => {source => $source});

	return $#features > -1 ? $features[0] : undef;
}

sub get_overlapping_features { #bad name but returns all overlapping features of same type, strand independent
	my ($self,$s) = @_;
	
	my ($f,$id,$source,$type,$start,$stop,$score,$strand,$phase,$attributes);

	if (ref($s) eq 'ARRAY'){
	#id consists of abbreviation.ori.g.inal.copy
		($id,$source,$type,$start,$stop,$score,$strand,$phase,$attributes) = @{$s};		
		$f = Bio::SeqFeature::Generic->new( -start => $start, -end => $stop, -strand => $strand);
	} else {
		($id,$type,$source,$start,$stop,$strand) = ($s->seq_id, $s->type, $s->source, $s->start, $s->stop, $s->strand);		
		$f = $s;
	}
	my ($abbr,@orig) = split /\./ , $id;
	pop @orig;
	$id = join '.' , ($abbr,@orig);

	my @features;
	#without strand info due to a) crt and b) same tool wont annotate same type of gene at both strands

	for ($self->db->{$abbr}->features(-primary_tag => $type , -attributes => {source => $source})){		
		my @tmp = split /\./, $_->seq_id;
		pop @tmp;
		next unless join('.' , @tmp) eq $id;
		my ($start, $stop, $strand) = $f->intersection($_);
		push @features , $_ if $start && ($stop - $start) > 0;
	}
	return \@features;
}

sub get_all_overlapping_features { #bad name but returns all features (type independent) overlapping at least 70% in both direction
	my ($self,$s) = @_;
	
	my ($f,$id,$source,$type,$start,$stop,$score,$strand,$phase,$attributes);

	if (ref($s) eq 'ARRAY'){
	#id consists of abbreviation.ori.g.inal.copy
		($id,$source,$type,$start,$stop,$score,$strand,$phase,$attributes) = @{$s};		
		$f = Bio::SeqFeature::Generic->new( -start => $start, -end => $stop, -strand => $strand);
	} else {
		($id,$type,$source,$start,$stop,$strand) = ($s->seq_id, $s->type, $s->source, $s->start, $s->stop, $s->strand);		
		$f = $s;
	}
	my ($abbr,@orig) = split /\./ , $id;
	pop @orig;
	$id = join '.' , ($abbr,@orig);

	my @features;
	
	for ($self->db->{$abbr}->features(-attributes => {source => $source})){		
		my @tmp = split /\./, $_->seq_id;
		pop @tmp;
		next unless join('.' , @tmp) eq $id;
		my ($istart, $istop, $istrand) = $f->intersection($_);
		push @features , $_ if $istart && $istop - $istart >= ($stop - $start) * 0.7 && $istop - $istart >= ($_->stop - $_->start) * 0.7;
	}
	
	return \@features;
}

sub get_user_overlaps {
	my ($self,$id,$start,$stop,$strand) = @_;
	return () unless $strand || $strand eq '.';

	my ($abbr,@orig) = split /\./ , $id;
	my $orig = join '.' , @orig;

	return $self->userdb->{$abbr}->features(-seq_id => $orig, -start => $start , -strand => $strand eq '+' ? 1 : -1, -stop => $stop , -range_type => 'overlaps');
}

#gorap specific set up databases for all genome files from parameter object and 
#fill in existing data in defined output directory
sub _set_db {
	my ($self) = @_;

	my @abbr;
	push @abbr , basename($_) for glob catfile($self->parameter->output,'annotations','*.final.gff');	
	$_=~s/\.final\.gff// for @abbr;

	for ( @{$self->parameter->abbreviations},@abbr ){		
		$self->db->{$_} = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -verbose => -1 );
		$self->db->{$_}->init_database([1]);
	}		

	&_read;

	return unless $self->parameter->has_gffs;
	if ($self->parameter->verbose){
		print ! $self->parameter->has_bams && $self->parameter->notpm ? "Reading GFF files\n" : "Reading GFF files and calculating mapping statistics\n";
	}
	
	my $idc=0;
	for (0..$#{$self->parameter->gffs}){
		my $abbr = ${$self->parameter->abbreviations}[$_];		
		for my $f (@{${$self->parameter->gffs}[$_]}){
			$self->userdb->{$abbr} = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -verbose => -1 ) unless $self->userdb->{$abbr};
			print $f."\n";
			open GFF, '<'.$f or die $!;
			&add_user_entry($self,$self->userdb->{$abbr},$_,$abbr,(++$idc)) while <GFF>;
			close GFF;
		}		
	}
}

sub add_user_entry(){
	my ($self,$db,$line,$abbr,$idc) = @_;

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

	my ($tpm,$rpkm,$reads) = ('.','.','.');
	if (! $self->parameter->notpm && $self->parameter->has_bams && ($line[2]=~/exon/i || $attr!~/parent=/i)) {
		($tpm,$rpkm,$reads) = $self->bamdb->rpkm($abbr,$line[0],$line[3],$line[4],$line[6]);		
	}	
	
	$db->new_feature(
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
			RPKM => $rpkm,
			TPM => $tpm,
			Reads => $reads
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

sub add_attribute {
	my ($self,$gffentry,$attribute,$value) = @_;

	my ($id,$source,$type,$start,$stop,$score,$strand,$phase,$attributes) = @$gffentry;
	my ($abbr,@orig) = split /\./, $id;
	pop @orig;

	my ($feature) = $self->db->{$abbr}->features(-seq_id => $id , -primary_tag => $type , -attributes => {source => $source});
	$feature->remove_tag($attribute);
	$feature->add_tag_value($attribute,$value);
	$feature->update;
}


sub add_seq {
	my ($self,$seq,$id,$type,$source,$abbr) = @_;

	my ($feature) = $self->db->{$abbr}->features(-seq_id => $id , -primary_tag => $type , -attributes => {source => $source});
	$feature->remove_tag('seq');
	$feature->add_tag_value('seq',$seq);
	$feature->update;
}

#gorap way to fill database with existing data in defined output directory
sub _read {
	my ($self) = @_;

	my @fastas = glob catfile($self->parameter->output,'annotations','*.fa');
	if ($self->parameter->verbose && $#fastas > -1 ){
		print "Reading previous annotations from ".$self->parameter->output." : ";
		print "!!!do not kill GORAP here!!!\n";
	}	

	for my $file (@fastas){		
		my $headerMapSeq={};		
		open FA , '<'.$file or die $!;
		my $seqid;
		my $seq;
		while(<FA>){
			chomp $_;
			if ($_=~/^>\s*(\S+)/){		
				#temp store of sequences
				#id consists of abbreviation.original.type.tool.copy									
				$headerMapSeq->{$seqid}=$seq if $seqid;				
				$seqid = $1;
				$seq = '';									
			} else {
				$seq .= $_;
			}
		}
		close FA;		
		unlink $file;
			
		$headerMapSeq->{$seqid}=$seq if $seqid;

		$file=~s/\.fa$/\.gff/;

		open GFF , '<'.$file or die $!;
		while(<GFF>){	
			chomp $_;
			my @l = split /\s+/ , $_;		
			next unless $l[1];
			#build header like for fasta: abbreviation.original.type.tool.copy out of abbreviation.original.copy
			my @id = split /\./ , $l[0];										
			my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
			my $id = join '.' , ($abbr,$orig,$l[2],$l[1],$copy);				
			$l[1] =~ s/\W//g;
			$l[1] = lc $l[1];
			&add_gff3_entry($self,\@l, exists $headerMapSeq->{$id} ? $headerMapSeq->{$id} : '',$abbr);			
		}
		close GFF;			
		unlink $file;
	}	
}

sub store {
	my ($self,$type) = @_;
	
	if ($type){
		for my $abbr (keys %{$self->db}){			
			my @features = sort {$a->seq_id cmp $b->seq_id || $a->start <=> $b->start || $a->stop <=> $b->stop} $self->db->{$abbr}->features(-primary_tag => $type);
			open GFF , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.gff') or die $!;
			open FA , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.fa') or die $!;
			open GFFF , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.final.gff') or die $!;
			open FAF , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.final.fa') or die $!;
			for my $i (0..$#features){
				my $f1 = $features[$i];

				my $source = $f1->source;	
				$source =~ s/GORAP//;

				my @id = split /\./ , $f1->seq_id;
				my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
				my $faid = join '.' , ($abbr,$orig,$f1->primary_tag,$source,$copy);			

				my $reads = ($f1->get_tag_values('reads'))[0];
				$reads = 0 unless $reads;
				my $rpkm = ($f1->get_tag_values('rpkm'))[0];
				$rpkm = 0 unless $rpkm;
				my $tpm = ($f1->get_tag_values('tpm'))[0];
				$tpm = 0 unless $tpm;

				my $attributes = 'TPM='.$tpm.';FPKM='.$rpkm.';Reads='.$reads.';Filter='.$f1->display_name.';Note='.($f1->get_tag_values('notes'))[0];

				if ( $f1->display_name eq '!' ){
					print GFFF $f1->seq_id."\t".$source."\t".$f1->primary_tag."\t".$f1->start."\t".$f1->stop."\t".$f1->score."\t",$f1->strand ? $f1->strand > 0 ? '+' : '-' : '.',"\t".$f1->phase."\t".$attributes."\n";
					my ($seq) = $f1->get_tag_values('seq');
					print FAF '>'.$faid."\n".$seq."\n" if $seq;
				} else {
					print GFF $f1->seq_id."\t".$source."\t".$f1->primary_tag."\t".$f1->start."\t".$f1->stop."\t".$f1->score."\t",$f1->strand ? $f1->strand > 0 ? '+' : '-' : '.',"\t".$f1->phase."\t".$attributes."\n";
					my ($seq) = $f1->get_tag_values('seq');
					print FA '>'.$faid."\n".$seq."\n" if $seq;
				}
			}
			close GFF;
			close FA;
			close GFFF;
			close FAF;
		}
	} else {

		for my $abbr (keys %{$self->db}){
			my @features = sort {$a->seq_id cmp $b->seq_id || $a->start <=> $b->start || $a->stop <=> $b->stop} $self->db->{$abbr}->features();
			
			open GFF , '>'.catfile($self->parameter->output,'annotations',$abbr.'.gff') or die $!;
			open FA , '>'.catfile($self->parameter->output,'annotations',$abbr.'.fa') or die $!;
			open GFFF , '>'.catfile($self->parameter->output,'annotations',$abbr.'.final.gff') or die $!;
			open FAF , '>'.catfile($self->parameter->output,'annotations',$abbr.'.final.fa') or die $!;
			for my $i (0..$#features){
				my $f1 = $features[$i];

				my $source = $f1->source;	
				$source =~ s/GORAP//;

				my @id = split /\./ , $f1->seq_id;
				my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
				my $faid = join '.' , ($abbr,$orig,$f1->primary_tag,$source,$copy);			

				my $reads = ($f1->get_tag_values('reads'))[0];
				$reads = 0 unless $reads;
				my $rpkm = ($f1->get_tag_values('rpkm'))[0];
				$rpkm = 0 unless $rpkm;
				my $tpm = ($f1->get_tag_values('tpm'))[0];
				$tpm = 0 unless $tpm;

				my $attributes = 'TPM='.$tpm.';FPKM='.$rpkm.';Reads='.$reads.';Filter='.$f1->display_name.';Note='.($f1->get_tag_values('notes'))[0];

				if ( $f1->display_name eq '!' ){
					print GFFF $f1->seq_id."\t".$source."\t".$f1->primary_tag."\t".$f1->start."\t".$f1->stop."\t".$f1->score."\t",$f1->strand ? $f1->strand > 0 ? '+' : '-' : '.',"\t".$f1->phase."\t".$attributes."\n";						
					my $seq = $f1->get_tag_values('seq');
					print FAF '>'.$faid."\n".$seq."\n" if $seq;
				} else {
					print GFF $f1->seq_id."\t".$source."\t".$f1->primary_tag."\t".$f1->start."\t".$f1->stop."\t".$f1->score."\t",$f1->strand ? $f1->strand > 0 ? '+' : '-' : '.',"\t".$f1->phase."\t".$attributes."\n";
					my $seq = $f1->get_tag_values('seq');
					print FA '>'.$faid."\n".$seq."\n" if $seq;
				}
			}
			close GFF;
			close FA;
			close GFFF;
			close FAF;
		}
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

sub store_overlaps {
	my ($self) = @_;
	
	for my $abbr (keys %{$self->db}){		
		my @features = sort {my @id = split /\./ , $a->seq_id; my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
								my @id2 = split /\./ , $b->seq_id; my ($abbr2,$orig2,$copy2) = ($id2[0] , join('.',@id2[1..($#id2-1)]) , $id2[-1]);
								$orig cmp $orig2 || $a->start <=> $b->start || $a->stop <=> $b->stop} $self->db->{$abbr}->features();
		my $overlaps;
		open GFF , '>'.catfile($self->parameter->output,'annotations',$abbr.'.gff') or die $!;
		if (exists $self->userdb->{$abbr}){
			open GFFP , '>'.catfile($self->parameter->output,'annotations',$abbr.'.external.gff') or die $! ;			
			print GFFP $_."\n" for &get_user_entries($self,$abbr);			
			close GFFP if exists $self->userdb->{$abbr};
		}		

		open FA , '>'.catfile($self->parameter->output,'annotations',$abbr.'.fa') or die $!;
		open GFFF , '>'.catfile($self->parameter->output,'annotations',$abbr.'.final.gff') or die $!;
		open FAF , '>'.catfile($self->parameter->output,'annotations',$abbr.'.final.fa') or die $!;
		for my $i (0..$#features){
			my $f1 = $features[$i];

			my $source = $f1->source;	
			$source =~ s/GORAP//;
			
			my @id = split /\./ , $f1->seq_id;
			my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
			my $faid = join '.' , ($abbr,$orig,$f1->primary_tag,$source,$copy);			
			my $f1id = $abbr.'.'.$orig;

			if ($f1->display_name eq '!'){
				for ($i+1..$#features){					
					my $f2 = $features[$_];					
					if ($f1->strand && $f2->strand){
						next if $f1->strand != $f2->strand;
					}
					last if $f2->start > $f1->stop;
					next unless $f2->display_name eq '!';
					my @id2 = split /\./ , $f2->seq_id;
					my ($abbr2,$orig2,$copy2) = ($id2[0] , join('.',@id2[1..($#id2-1)]) , $id2[-1]);
					next unless $orig eq $orig2;
					$overlaps->{$f1->seq_id}->{$f2->primary_tag}=1;
					$overlaps->{$f2->seq_id}->{$f1->primary_tag}=1;
				}
			}

			#with given gff
			if(exists $self->userdb->{$abbr}){
				for($self->userdb->{$abbr}->features(-seq_id => $orig, -start => $f1->start , -strand => $f1->strand, -stop => $f1->stop , -range_type => 'overlaps')){
					my ($userid) = $_->get_tag_values('ID');
					next unless $userid;
					$overlaps->{$f1->seq_id}->{$userid}=1;
				}
			}
			
			my $o = join(',',keys %{$overlaps->{$f1->seq_id}});
			$o='.' unless $o;

			my $reads = ($f1->get_tag_values('reads'))[0];
			$reads = 0 unless $reads;
			my $rpkm = ($f1->get_tag_values('rpkm'))[0];
			$rpkm = 0 unless $rpkm;
			my $tpm = ($f1->get_tag_values('tpm'))[0];
			$tpm = 0 unless $tpm;

			my $attributes = 'TPM='.$tpm.';FPKM='.$rpkm.';Reads='.$reads.';Filter='.$f1->display_name.';Note='.($f1->get_tag_values('notes'))[0].';Overlaps='.$o;

			if ( $f1->display_name eq '!' ){
				print GFFF $f1->seq_id."\t".$source."\t".$f1->primary_tag."\t".$f1->start."\t".$f1->stop."\t".$f1->score."\t",$f1->strand ? $f1->strand > 0 ? '+' : '-' : '.',"\t".$f1->phase."\t".$attributes."\n";
				my $seq = $f1->get_tag_values('seq');
				print FAF '>'.$faid."\n".$seq."\n" if $seq;
			} else {
				print GFF $f1->seq_id."\t".$source."\t".$f1->primary_tag."\t".$f1->start."\t".$f1->stop."\t".$f1->score."\t",$f1->strand ? $f1->strand > 0 ? '+' : '-' : '.',"\t".$f1->phase."\t".$attributes."\n";
				my $seq = $f1->get_tag_values('seq');
				print FA '>'.$faid."\n".$seq."\n" if $seq;
			}
		}
		close GFF;
		close FA;
		close GFFF;
		close FAF;
	}
}

1;