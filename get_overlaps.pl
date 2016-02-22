#! /usr/bin/perl
use strict;
use warnings;
use Bio::DB::SeqFeature::Store;
use List::Util qw(min max);
use Getopt::Long;
use File::Basename;
use Bio::Index::Fasta;
use Bio::Seq;
use IPC::Cmd qw(run);
use Sys::MemInfo qw(totalmem);
use File::Spec::Functions;
use File::Path qw(make_path);

(Getopt::Long::Parser->new)->getoptions (
	'gff1|gff1=s' => \my $gff1,
	'g1c|g1constrain=s' => \my $g1constrain,
	'g1cn|g1constrneg=s' => \my $g1constrneg,
	'gff2|gff2=s' => \my $gff2,
	'g2c|g2constrain=s' => \my $g2constrain,
	'g2cn|g2constrneg=s' => \my $g2constrneg,
	'o|overlap=f' => \my $overlap,
	't|threads=i' => \my $threads,
	'mem|memory=i' => \my $memory,
	'g|genome=s' => \my $genome,
	'cm|cm=s' => \my $cms,
	'c|combine' => \my $combine,
	'v|verbose' => \my $verbose,
	'h|help' => \my $help
) or &usage();

&usage() if $help || ! $gff1 || ! $gff2;

sub usage(){
	print <<EOF;
	-gff1 | --gff1        => (reqiered) <FILE> reference GFF
	-g1c  | --g1constrain => <STRING> e.g. 'gene' to compare type=gene only
	-g1cn | --g1constrneg => <STRING> e.g. 'gene' to compare type!=gene only
	-gff2 | --gff2        => (reqiered) <FILE> GFF to compare against reference GFF
	-g2c  | --g2constrain => <STRING> e.g. 'gene' to compare type=gene only
	-g2cn | --g2constrneg => <STRING> e.g. 'gene' to compare type!=gene only
	-o    | --overlap     => <FLOAT> percentage of necessary overlap (0.75)
	-t    | --threads     => <INT> number of cpu cores to use (1)
	-mem  | --memory      => <INT> amount of RAM to use (all)
	-g    | --genome      => <FILE> GFF related FASTA file to assign Rfam ids to genes of -gff2
	-cm   | --cm          => (reqiered for -g) <STRING> wildcarded path to Rfam covariance models, e.g. '/dir/*.cm'
	-c    | --combine     => combines duplicated annotations of -gff2 to count them only once
	-v    | --verbose     => print to STDOUT instead of files: <gff2>.known and <gff2>.unknown
	-h    | --help        => this message
EOF
	exit;
}

my @cm;
my $fadb;
my $temp;
$overlap = 0.75 unless $overlap;
if ($genome){
	&usage() unless $cms;
	@cm = glob $cms;
	$threads = 1 unless $threads;
	unless ($memory){
		$memory = sprintf("%.0f",Sys::MemInfo::totalmem()/1024/1024 * 0.9);	
		$memory = 4000 unless $memory=~/^\d+$/;		
	}	
	if ($ENV{TMPDIR}){
		$temp = $ENV{TMPDIR};
	} elsif (-e catdir(rootdir,'tmp')){
		$temp = catdir(rootdir,'tmp');
	} else {
		make_path(catdir($ENV{PWD},'tmp'));
		$temp = catdir($ENV{PWD},'tmp');
	}
	$fadb = Bio::Index::Fasta->new(-filename => catfile($temp,"$$.idx"), -write_flag => 1 , -verbose => -1);
	$fadb->make_index($genome);
}

my $gffdb1 = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -verbose => -1 );
my $gffdb2 = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -verbose => -1 );
$gffdb1->init_database([1]);
$gffdb2->init_database([1]);

open G1, '<'.$gff1 or die $!;
open G2, '<'.$gff2 or die $!;
my $idc=0;
&add_entry($gffdb1,$_,$g1constrain,$g1constrneg) while <G1>;
$idc=0;
&add_entry($gffdb2,$_,$g2constrain,$g2constrneg) while <G2>;
close G1;
close G2;

sub add_entry(){
	my ($db,$line,$constrain,$constraineg) = @_;

	chomp $line;
	return unless $line || $line=~/^#/ || $line=~/^\s*$/;
	my @line = split /\t+/ , $line;
	
	return if $#line < 7;
	
	push @line , '.' if $#line == 7;
	if ($constrain){
		return unless $line[2]=~/$constrain/i;
	}
	if ($constraineg){
		return if $line[2]=~/$constraineg/i;
	}

	$line[8].=';';
	$line[8]=~/ID=(.+?);/i;

	my $id;
	my $attr;
	if ($1){		
		$id = $1;
		$attr = substr($line[8],0,$-[0]).substr($line[8],$+[0],-1);
	} else{		
		$id = 'feature'.($idc++);
		if ($line[8] eq '.;'){
			$attr = '';
		} else {
			$attr = length($line[8]) < 3 ? 'Note='.substr($line[8],0,-1) : substr($line[8],0,-1);	
		}		
	}

	return if $attr=~/Note=(\?)/ || $attr=~/Filter=(\?)/;	
	
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
			attr => $attr
		}
	);	
}

sub entry_tostring(){
	my ($f) = @_;	

	my @out = split /\s+/ , $f->gff_string;
	splice @out , 8 , 2;
	my @attributes;
	for (sort {$a cmp $b} $f->get_all_tags()){
		my ($v) = $f->get_tag_values($_);
		push @attributes , $_ eq 'attr' ? $v: $_."=".$v if $v;
	}	
	push @out , (join ';' , @attributes);

	return join "\t" , @out;
}

my $detected = 0;
my $total = scalar $gffdb2->features();
my $known = 0;

open KNOW, '>'.$gff2.'known' or die $! unless $verbose;
open UNDE, '>'.$gff2.'undetected' or die $! unless $verbose;
my @features = $gffdb1->features();
my $c=0;
for my $f1 (@features){
	print ''.($c++).' of '.$#features."\r" unless $verbose;
	my @overlaps = $gffdb2->features(-seq_id => $f1->seq_id, -start => $f1->start , -stop => $f1->stop , -range_type => 'overlaps');
		
	my @tp;
	for my $f2 (@overlaps){
		if ($f2->strand && $f2->primary_tag!~/PK-G12/i){
			next unless $f1->strand == $f2->strand;
		
			my ($start, $stop, $strand) = $f1->intersection($f2);		
			if ( $stop - $start >= ($f1->stop - $f1->start) * $overlap && $stop - $start >= ($f2->stop - $f2->start) * $overlap ){
				push @tp , $f2;
				$known++;			
			}	
		} else {
			my ($start, $stop, $strand) = $f1->intersection($f2);		
			if ( ($stop - $start >= ($f1->stop - $f1->start) * $overlap && $stop - $start >= ($f2->stop - $f2->start) * $overlap) ||
				 ($stop - $start ==  $f1->stop - $f1->start) ||
				 ($stop - $start ==  $f2->stop - $f2->start) ){
				push @tp , $f2;
				$known++;
			}
		}
	}
	if($#tp > -1){
		$detected++;		
		unless ($verbose){
			print KNOW &entry_tostring($f1)."\n";
			print KNOW "\t".&entry_tostring($_)."\n" for @tp;
		} else {
			print &entry_tostring($f1)."\n";
			print "\t".&entry_tostring($_)."\n" for @tp;
		}
		$gffdb2->delete(@tp);
	} else {
		unless ($verbose){
			print UNDE &entry_tostring($f1)."\n";			
		} else {
			print &entry_tostring($f1)."\n";			
		}
	}
}

unless ($verbose) {	
	print KNOW "total $total\n";
	print KNOW "known $known\n";
	print KNOW "unknown ".($total - $known)."\n";
	print KNOW "undetected ".($#features + 1 - $detected)."\n";
	close KNOW;
	close UNDE;
} else {
	print "total $total\n";
	print "known $known\n";
	print "unknown ".($total - $known)."\n";
	print "undetected ".($#features + 1 - $detected)."\n";
}

my $unknown = 0;
open NOKNOW, '>'.$gff2.'unknown' or die $! unless $verbose;
while( my ($f1) = $gffdb2->features()){

	print ''.(scalar $gffdb2->features())." to go \r" unless $verbose;

	my @overlaps = $gffdb2->features(-seq_id => $f1->seq_id, -start => $f1->start , -stop => $f1->stop , -range_type => 'overlaps');
	my @tp;	

	for my $f2 (@overlaps){
		if ($f2->strand){
			next unless $f1->strand == $f2->strand;
		}
		my ($start, $stop, $strand) = $f1->intersection($f2);
		if ( $stop - $start > ($f1->stop - $f1->start) * $overlap && $stop - $start > ($f2->stop - $f2->start) * $overlap ){
			push @tp , $f2;
		} 
	}

	$unknown++;	
		
	unless($verbose){
		if ($combine){			
			if ($genome){
				my @f = reverse sort {$a->score <=> $b->score} grep {$_->score ne '.'} @tp;
				my $f;
				if ($#f == -1){
					$f = $f1;
				} else {
					$f = $f[0]->score > 1 ? $f[0] : $f[-1];
				}
				my ($newtag, $score, $seq) = &align($f);				
				$f->primary_tag($newtag);
				$f->score($score);
				$f->update;				
				print NOKNOW &entry_tostring($f)."\t$seq\n";
			} else {
				my @f = reverse sort {$a->score <=> $b->score} grep {$_->score ne '.'} @tp;
				push @f , $f1 if $#f == -1;
				my $f = $f[0]->score > 1 ? $f[0] : $f[-1];
				print NOKNOW &entry_tostring($f)."\n";
			}			
		} else {
			if ($genome){
				for(@tp){
					my ($newtag, $score, $seq) = &align($_);					
					$_->primary_tag($newtag);
					$_->score($score);
					$_->update;				
					print NOKNOW &entry_tostring($_)."\t$seq\n";
				}
			} else {				
				print NOKNOW &entry_tostring($_)."\n" for @tp;
			}			
		}
	} else {
		if ($combine){			
			if ($genome){
				my @f = reverse sort {$a->score <=> $b->score} grep {$_->score ne '.'} @tp;
				my $f;
				if ($#f == -1){
					$f = $f1;
				} else {
					$f = $f[0]->score > 1 ? $f[0] : $f[-1];
				}
				my ($newtag, $score, $seq) = &align($f);				
				$f->primary_tag($newtag);
				$f->score($score);
				$f->update;				
				print &entry_tostring($f)."\t$seq\n";
			} else {
				my @f = reverse sort {$a->score <=> $b->score} grep {$_->score ne '.'} @tp;
				push @f , $f1 if $#f == -1;
				my $f = $f[0]->score > 1 ? $f[0] : $f[-1];
				print &entry_tostring($f)."\n";
			}			
		} else {
			if ($genome){
				for(@tp){
					my ($newtag, $score, $seq) = &align($_);					
					$_->primary_tag($newtag);
					$_->score($score);
					$_->update;					
					print &entry_tostring($_)."\t$seq\n";
				}
			} else {				
				print &entry_tostring($_)."\n" for @tp;
			}			
		}
	}
	$gffdb2->delete(@tp);
}

unless($verbose){
	print NOKNOW "combined:\n";
	print NOKNOW "total ".($unknown + $detected)."\n";
	print NOKNOW "known ".($detected)."\n";
	print NOKNOW "unknown $unknown\n";
	print NOKNOW "undetected ".($#features + 1 - $detected)."\n";
	close NOKNOW;
} else {
	print "combined:\n";
	print "total ".($unknown + $detected)."\n";
	print "known ".($detected)."\n";
	print "unknown $unknown\n";
	print "undetected ".($#features + 1 - $detected)."\n";
}

unlink $_ for glob catfile($temp,"$$*");

sub align(){
	my ($f) = @_;

	my $rna=$f->primary_tag;		
	
	return ($f->primary_tag, $f->score, '.') if ! $f->strand || $f->primary_tag eq '.';
	my $seq = &get_subseq($f->seq_id, $f->start, $f->stop, $f->strand);
	return ($f->primary_tag, $f->score, $seq) if ($f->stop - $f->start) > 3000;

	my $max=-999999999999;
	my $maxm;

	my @m = grep { $_=~/_$rna\.cm/i } @cm;
	@m = grep { $_=~/_$rna/i } @cm if $#m == -1;

	if ($#m == -1){
		($rna) = split /_/ , $rna;
		@m = grep { $_=~/_$rna\.cm/i } @cm;
		@m = grep { $_=~/_$rna/i } @cm if $#m == -1;
	}
	if ($#m == -1){
		($rna) = split /-/ , $rna;
		@m = grep { $_=~/_$rna\.cm/i } @cm;
		@m = grep { $_=~/_$rna/i } @cm if $#m == -1;
	}
	if ($#m == -1){		
		my ($v) = $f->get_tag_values('attr');
		
		return ($f->primary_tag, $f->score, $seq) unless $v;
		$v=~/rfam-id=(.+?)(;|$)/;
		$v=~/product=(.+?)(;|$)/ unless $1;
		$v=~/gene_biotype=(.+?)(;|$)/ unless $1;		
		
		return ($f->primary_tag, $f->score, $seq) unless $1;
		$v = $1;		

		$rna = $v;
		@m = grep { $_=~/_$rna\.cm/i } @cm;
		@m = grep { $_=~/_$rna/i } @cm if $#m == -1;

		($rna) = split /_/ , $v;

		if ($#m == -1){
			($rna) = split /-/ , $v;
			@m = grep { $_=~/_$rna\.cm/i } @cm;
			@m = grep { $_=~/_$rna/i } @cm if $#m == -1;
		}
		if ($#m == -1){
			($rna) = split /\s+/ , $v;
			@m = grep { $_=~/_$rna\.cm/i } @cm;
			@m = grep { $_=~/_$rna/i } @cm if $#m == -1;
		}
	}

	if ($#m > -1){
		if ($#m > 0){			
			for my $m (@m){				
				my ($success, $error_code, $full_buf, $stdout_buf, $stderr_buf) = run( command => "printf \"\>foo\\n$seq\" | cmalign --mxsize $memory --sfile ".catfile($temp,"$$.tmp")." --noprob --cpu $threads $m -", verbose => 0 );
				my $score = &get_score();
				if($score > $max){
					$max = $score;
					$maxm = $m;
				}			
			}
			return ($f->primary_tag.':'.basename($maxm), $max, $seq);
		} else {
			return ($f->primary_tag.':'.basename($m[0]), $f->score, $seq);
		}
	} else {
		return ($f->primary_tag, $f->score, $seq);
	}	
}

sub get_subseq(){
	my ($chr,$sta,$sto,$strand) = @_;

	my $seq = $fadb->fetch($chr);
	return $strand eq "+" || $strand > 0 ? $seq->subseq($sta,$sto) : ((Bio::Seq->new( -seq => $seq->subseq($sta,$sto) , -verbose => -1))->revcom)->seq;	
}

sub get_score(){	

	open S , '<'.catfile($temp,"$$.tmp") or die $!;
	my $l;
	while(<S>){
		chomp $_;
		$_ =~ s/^\s+|\s+$//g;
		next if $_=~/^#/;
		next if $_=~/^\s*$/;
		$l=$_;
	}
	close S;

	return (split /\s+/ , $l)[6];
}