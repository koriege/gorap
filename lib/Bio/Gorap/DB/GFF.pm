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

has 'bamdb' => (
	is => 'ro',
	isa => 'Bio::Gorap::DB::BAM'
);

#add a gff3 entry into this Bio::DB::SeqFeature database
sub add_gff3_entry {	
	my ($self,$s,$seq,$abbr) = @_;
	
	#id consists of abbreviation.original.copy	
	if ($#{$s} > 6 && $seq){
		my ($id,$source,$type,$start,$stop,$score,$strand,$phase,$attributes) = @{$s};
		
		my @overlaps = ('.');
		my ($rpkm,$reads,$filter) = ('.','.','!');

		if ($attributes && $#{$self->parameter->queries} > -1){
			$attributes.=';';
			$rpkm = $1 ? $1 : '.' if $attributes=~/RPKM=(.+?);/;
			$reads = $1 ? $1 : '.' if $attributes=~/Reads=(.+?);/;
			$filter = $1 ? $1 : '!' if $attributes=~/Filter=(.+?);/;		
			@overlaps = split /,/ , $1 if $attributes=~/Overlaps=(.+?);/;
		} else {
			$attributes.=';' if $attributes;
			$rpkm = $1 ? $1 : '.' if $attributes && $attributes=~/RPKM=(.+?);/;
			$reads = $1 ? $1 : '.' if $attributes && $attributes=~/Reads=(.+?);/;
			$filter = $1 ? $1 : '!' if $attributes && $attributes=~/Filter=(.+?);/;		
			@overlaps = split /,/ , $1 if $attributes && $attributes=~/Overlaps=(.+?);/;
			if ($self->parameter->has_bams){
				my @id = split /\./ , $id;										
				my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
				($rpkm,$reads) = $self->bamdb->rpkm($orig,$start,$stop);
			}
		}

		$self->db->{$abbr}->new_feature(
		 	-start => $start, 
	        -stop => $stop,
	        -strand => $strand eq '+' ? +1 : -1,         
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
	        	rpkm => $rpkm, 
	        	reads =>  $reads,        	
	        	overlaps => \@overlaps,
	        	seq => $seq,
	        	source => $source,
	        	#search score
	        	origscore => $score
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

sub get_overlapping_features {
	my ($self,$s,$abbr) = @_;
	
	#id consists of abbreviation.ori.g.inal.copy
	my ($id,$source,$type,$start,$stop,$score,$strand,$phase,$attributes) = @{$s};
	my @tmp = split /\./, $id;
	pop @tmp;
	$id = join '.' , @tmp;

	my @features;
	my $f = Bio::SeqFeature::Generic->new( -start => $start, -end => $stop, -strand => $strand eq '+' ? 1 : -1);
	for ($self->db->{$abbr}->features(-primary_tag => $type , -attributes => {source => $source})){
		next unless $_->strand == $f->strand;
		@tmp = split /\./, $_->seq_id;
		pop @tmp;
		next unless join '.' , @tmp eq $id;
		my ($start, $stop, $strand) = $f->intersection($_);
		push @features , $_ if $start && ($stop - $start) > 0;
	}
	return \@features;
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
}

sub add_db {
	my ($self,$abbr) = @_;

	unless (exists $self->db->{$abbr}){
		$self->db->{$abbr} = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -verbose => -1 );
		$self->db->{$abbr}->init_database([1]);	
	}
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

	print "Reading GFF files from ".$self->parameter->output."\n" if $self->parameter->verbose;
	for my $file (glob catfile($self->parameter->output,'annotations','*.fa')){		
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
			&add_gff3_entry($self,\@l, $headerMapSeq->{$id},$abbr);			
		}
		close GFF;			
	}	
}

sub store {
	my ($self,$type) = @_;
	
	# if ($type){
	# 	for my $abbr (keys %{$self->db}){			
	# 		my @features = sort {$a->seq_id cmp $b->seq_id || $a->strand <=> $b->strand || $a->start <=> $b->start || $a->stop <=> $b->stop} $self->db->{$abbr}->features(-primary_tag => $type);
	# 		open GFF , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.gff') or die $!;
	# 		open FA , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.fa') or die $!;
	# 		open GFFF , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.final.gff') or die $!;
	# 		open FAF , '>>'.catfile($self->parameter->output,'annotations',$abbr.'.final.fa') or die $!;
	# 		for my $i (0..$#features){
	# 			my $f1 = $features[$i];

	# 			my $source = $f1->source;	
	# 			$source =~ s/GORAP//;

	# 			my @id = split /\./ , $f1->seq_id;
	# 			my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
	# 			my $faid = join '.' , ($abbr,$orig,$f1->primary_tag,$source,$copy);			

	# 			my $attributes = 'RPKM='.($f1->get_tag_values('rpkm'))[0].';Reads='.($f1->get_tag_values('reads'))[0].';Filter='.$f1->display_name;

	# 			if ( $f1->display_name eq '!' ){
	# 				print GFFF $f1->seq_id."\t".$source."\t".$f1->primary_tag."\t".$f1->start."\t".$f1->stop."\t".$f1->score."\t",$f1->strand > 0 ? '+' : '-',"\t".$f1->phase."\t".$attributes."\n";						
	# 				print FAF '>'.$faid."\n".($f1->get_tag_values('seq'))[0]."\n";
	# 			} else {
	# 				print GFF $f1->seq_id."\t".$source."\t".$f1->primary_tag."\t".$f1->start."\t".$f1->stop."\t".$f1->score."\t",$f1->strand > 0 ? '+' : '-',"\t".$f1->phase."\t".$attributes."\n";
	# 				print FA '>'.$faid."\n".($f1->get_tag_values('seq'))[0]."\n";
	# 			}
	# 		}
	# 		close GFF;
	# 		close FA;
	# 		close GFFF;
	# 		close FAF;
	# 	}
	# } else {

		for my $abbr (keys %{$self->db}){
			my @features = sort {$a->seq_id cmp $b->seq_id || $a->strand <=> $b->strand || $a->start <=> $b->start || $a->stop <=> $b->stop} $self->db->{$abbr}->features();
			
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

				my $attributes = 'RPKM='.($f1->get_tag_values('rpkm'))[0].';Reads='.($f1->get_tag_values('reads'))[0].';Filter='.$f1->display_name;

				if ( $f1->display_name eq '!' ){
					print GFFF $f1->seq_id."\t".$source."\t".$f1->primary_tag."\t".$f1->start."\t".$f1->stop."\t".$f1->score."\t",$f1->strand > 0 ? '+' : '-',"\t".$f1->phase."\t".$attributes."\n";						
					print FAF '>'.$faid."\n".($f1->get_tag_values('seq'))[0]."\n";
				} else {
					print GFF $f1->seq_id."\t".$source."\t".$f1->primary_tag."\t".$f1->start."\t".$f1->stop."\t".$f1->score."\t",$f1->strand > 0 ? '+' : '-',"\t".$f1->phase."\t".$attributes."\n";
					print FA '>'.$faid."\n".($f1->get_tag_values('seq'))[0]."\n";
				}
			}
			close GFF;
			close FA;
			close GFFF;
			close FAF;
		}
	# }
}

sub store_overlaps {
	my ($self) = @_;
	
	for my $abbr (keys %{$self->db}){		
		my @features = sort {my @id = split /\./ , $a->seq_id; my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
								my @id2 = split /\./ , $b->seq_id; my ($abbr2,$orig2,$copy2) = ($id2[0] , join('.',@id2[1..($#id2-1)]) , $id2[-1]);
								$orig cmp $orig2 || $a->strand <=> $b->strand || $a->start <=> $b->start || $a->stop <=> $b->stop
							} $self->db->{$abbr}->features();
		my $overlaps;
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
			my $f1id = $abbr.'.'.$orig;

			if ($f1->display_name eq '!'){
				for ($i+1..$#features){					
					my $f2 = $features[$_];
					last if $f1->strand != $f2->strand;
					last if $f2->start > $f1->stop;
					next unless $f2->display_name eq '!';					
					my @id2 = split /\./ , $f2->seq_id;
					my ($abbr2,$orig2,$copy2) = ($id2[0] , join('.',@id2[1..($#id2-1)]) , $id2[-1]);
					next unless $orig eq $orig2;
					$overlaps->{$f1->seq_id}->{$f2->primary_tag}=1;
					$overlaps->{$f2->seq_id}->{$f1->primary_tag}=1;
				}
			}
			
			my $o = join(',',keys %{$overlaps->{$f1->seq_id}});
			$o='.' unless $o;			
			my $attributes = 'RPKM='.($f1->get_tag_values('rpkm'))[0].';Reads='.($f1->get_tag_values('reads'))[0].';Filter='.$f1->display_name.';Overlaps='.$o;

			if ( $f1->display_name eq '!' ){
				print GFFF $f1->seq_id."\t".$source."\t".$f1->primary_tag."\t".$f1->start."\t".$f1->stop."\t".$f1->score."\t",$f1->strand > 0 ? '+' : '-',"\t".$f1->phase."\t".$attributes."\n";				
				print FAF '>'.$faid."\n".($f1->get_tag_values('seq'))[0]."\n";
			} else {
				print GFF $f1->seq_id."\t".$source."\t".$f1->primary_tag."\t".$f1->start."\t".$f1->stop."\t".$f1->score."\t",$f1->strand > 0 ? '+' : '-',"\t".$f1->phase."\t".$attributes."\n";
				print FA '>'.$faid."\n".($f1->get_tag_values('seq'))[0]."\n";
			}
		}
		close GFF;
		close FA;
		close GFFF;
		close FAF;
	}
}

1;