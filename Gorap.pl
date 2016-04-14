#! /usr/bin/perl

use strict;
use warnings;
use sigtrap qw(handler SIGABORT normal-signals);

use Bio::Gorap::ThrListener;
use Bio::Gorap::Parameter;
use Bio::Gorap::DB::GFF;
use Bio::Gorap::DB::STK;
use Bio::Gorap::DB::Fasta;
use Bio::Gorap::DB::BAM;
use Bio::Gorap::DB::Taxonomy;
use Bio::Gorap::Evaluation::HTML;
use Bio::Gorap::Tool::Piles;
use Bio::Gorap::CFG;
use Bio::Gorap::Functions::ToolParser;
use File::Spec::Functions;
use File::Basename;
use Try::Tiny;
use Cwd 'abs_path';
use Bio::Tree::Draw::Cladogram;
use Bio::TreeIO;
use List::Util qw(any);
use Hash::Merge qw(merge);

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$year = $year + 1900;
$mon += 1;
print "$mday.$mon.$year-$hour:$min:$sec\n";
my $stamp = "$mday.$mon.$year-$hour:$min:$sec";

print "\nFor help run Gorap.pl -h\n\n";

#push gorap tools to $PATH
my $PATHtools;
for(reverse glob(catdir($ENV{GORAP},'*','bin'))){	
	$PATHtools .= $PATHtools ? ":$_" : $_;
}
local $ENV{PATH} = $ENV{PATH} ? "$PATHtools:$ENV{PATH}" : $PATHtools;
my ($trnalib) = reverse glob catdir($ENV{GORAP},'tRNAscan-SE*');
local $ENV{PERL5LIB} = $ENV{PERL5LIB} ? "$trnalib:$ENV{PERL5LIB}" : $trnalib if $trnalib;

#set defaults and read CLI parameter or parameter file
my $parameter = Bio::Gorap::Parameter->new(
	#pwd => dirname(abs_path(__FILE__)),
	pwd => $ENV{PWD},
	pid => $$,
	commandline => 1
);

my $taxdb;
my $stkdb;
if ($parameter->taxonomy){
	$taxdb = Bio::Gorap::DB::Taxonomy->new(
		parameter => $parameter
	);
	$taxdb->findRelatedSpecies;
	$stkdb = Bio::Gorap::DB::STK->new(
		parameter => $parameter,
		taxonomy => $taxdb
	);
} else {
	$stkdb = Bio::Gorap::DB::STK->new(
		parameter => $parameter
	);
}

my $fastadb = Bio::Gorap::DB::Fasta->new(
	parameter => $parameter
);

my $bamdb = Bio::Gorap::DB::BAM->new(
	parameter => $parameter
);

#gorap gff3 storage database initialization without loss of existing data in output directory 
my $gffdb = Bio::Gorap::DB::GFF->new(
	parameter => $parameter,
	bamdb => $bamdb
); 

#starts a necessary stdout listener for forked jobs using io::select and io::pipe
#with related storage object and a codeRef for parsing/storing a string of information
#fetched during background calculations
my $thrListener = Bio::Gorap::ThrListener->new(
	threads => $parameter->threads,		
	storage => $gffdb, 
	#gorap background jobs filter and mark specific entries
	storage_saver => \&Bio::Gorap::DB::GFF::update_filter
);

&run();
#stops the thread listener and waits for remaining background jobs to be finished
$thrListener->stop;
#remove overlap marked sequences
my ($type_map_features) = $gffdb->get_filtered_features('O');
for (keys %$type_map_features){
	$stkdb->rm_seq_from_stk($type_map_features->{$_},$_,'O');
}
#store final annotation results
$gffdb->store_overlaps;
Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb->idToPath,$stamp);

if ($parameter->has_outgroups){	
	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
	$year = $year + 1900;
	$mon += 1;
	print "$mday.$mon.$year-$hour:$min:$sec\n";

	print "Preparing phylogeny reconstruction\n" if $parameter->verbose;
	my @newQ;	
	my @oldqueries = @{$parameter->queries};
	my @genomes = @{$parameter->genomes};
	my @abbres = @{$parameter->abbreviations};

	$parameter->set_queries();
	for my $cfg (@{$parameter->queries}){		
		$parameter->set_cfg($cfg);
		push @newQ , $parameter->cfg->rf if $#{$gffdb->get_all_features($parameter->cfg->rf_rna , '!')} > -1;
	}

	if ($#newQ > -1){	
		my $outdir = catdir($parameter->output,'phylogeny');

		unlink $_ for glob catfile($outdir,'RAxML_*');			

		$parameter->set_genomes($parameter->outgroups,$parameter->ogabbreviations);

		push @genomes, @{$parameter->genomes};
		push @abbres, @{$parameter->abbreviations};	

		$parameter->set_queries(\@newQ);
		
		unless ($parameter->skip_comp){

			print "\nAnnotation of outgroups for phylogeny reconstruction\n" if $parameter->verbose;
			
			$gffdb->add_db($_) for @{$parameter->abbreviations};			

			$fastadb = Bio::Gorap::DB::Fasta->new(
				parameter => $parameter
			);
			$thrListener = Bio::Gorap::ThrListener->new(
				threads => $parameter->threads,		
				storage => $gffdb, 
				#gorap background jobs filter and mark specific entries
				storage_saver => \&Bio::Gorap::DB::GFF::update_filter
			);			
			&run();	
			#stops the thread listener and waits for remaining background jobs to be finished
			$thrListener->stop;
			my ($type_map_features) = $gffdb->get_filtered_features('O');
			for (keys %$type_map_features){
				$stkdb->rm_seq_from_stk($type_map_features->{$_},$_,'O');
			}
			#store final annotation results
			$gffdb->store_overlaps;
			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb->idToPath,$stamp.'-outgroup');
		}
				
		my ($speciesSSU,$coreFeatures,$stkFeatures,$stkCoreFeatures) = &get_phylo_features(\@abbres);		

		if ( (any { exists $speciesSSU->{$_} } @{$parameter->abbreviations}) && scalar keys %$speciesSSU > 3){
			open FA , '>'.catfile($outdir,'SSU.fasta') or die $!;	
			for my $k (keys %$speciesSSU ){
				print FA '>'.$k."\n";
				$speciesSSU->{$k}=~s/(\W|_)/-/g;				
				print FA $_."\n" for unpack("(a80)*",$speciesSSU->{$k});				
			}			
			close FA;			

			my $ex = system('mafft --localpair --maxiterate 1000 --thread '.$parameter->threads.' '.catfile($outdir,'SSU.fasta').' > '.catfile($outdir,'SSU.mafft'));
			&ABORT("mafft not found") unless $ex == 0;
			$ex = system('raxml -T '.$parameter->threads.' -f a -# 1000 -x 1234 -p 1234 -s '.catfile($outdir,'SSU.mafft').' -w '.$outdir.' -n SSU.mafft.tree -m GTRGAMMA -o '.join(',',grep { exists $speciesSSU->{$_} } @{$parameter->abbreviations}));
			&ABORT("raxml not found") unless $ex == 0;

			my $obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($outdir,'RAxML_bipartitions.SSU.mafft.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($outdir,'SSU.mafft.eps'));	

			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb->idToPath,$stamp.'-outgroup');			
		} 		

		if ( (any { exists $coreFeatures->{$_} } @{$parameter->abbreviations}) && scalar keys %$coreFeatures > 3){
			open FA , '>'.catfile($outdir,'coreRNome.fasta') or die $!;	
			for my $k (keys %$coreFeatures ){
				print FA '>'.$k."\n";
				$coreFeatures->{$k}=~s/(\W|_)/-/g;				
				print FA $_."\n" for unpack("(a80)*",$coreFeatures->{$k});				
			}
			close FA;
			open FA , '>'.catfile($outdir,'coreRNome.stkfa') or die $!;
			for my $k (keys %$stkCoreFeatures ){
				print FA '>'.$k."\n";
				$stkCoreFeatures->{$k}=~s/(\W|_)/-/g;				
				print FA $_."\n" for unpack("(a80)*",$stkCoreFeatures->{$k});				
			}			
			close FA;

			my $ex = system('mafft --localpair --maxiterate 1000 --thread '.$parameter->threads.' '.catfile($outdir,'coreRNome.fasta').' > '.catfile($outdir,'coreRNome.mafft'));
			&ABORT("mafft not found") unless $ex == 0;
			$ex = system('raxml -T '.$parameter->threads.' -f a -# 1000 -x 1234 -p 1234 -s '.catfile($outdir,'coreRNome.mafft').' -w '.$outdir.' -n coreRNome.mafft.tree -m GTRGAMMA -o '.join(',',grep { exists $coreFeatures->{$_} } @{$parameter->abbreviations}));
			&ABORT("raxml not found") unless $ex == 0;
			$ex = system('raxml -T '.$parameter->threads.' -f a -# 1000 -x 1234 -p 1234 -s '.catfile($outdir,'coreRNome.stkfa').' -w '.$outdir.' -n coreRNome.stk.tree -m GTRGAMMA -o '.join(',',grep { exists $coreFeatures->{$_} } @{$parameter->abbreviations}));
			&ABORT("raxml not found") unless $ex == 0;
			
			my $obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($outdir,'RAxML_bipartitions.coreRNome.mafft.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($outdir,'coreRNome.mafft.eps'));	
			$obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($outdir,'RAxML_bipartitions.coreRNome.stk.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($outdir,'coreRNome.stk.eps'));	

			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb->idToPath,$stamp.'-outgroup');			
		} 

		if ( (any { exists $stkFeatures->{$_} } @{$parameter->abbreviations}) && scalar keys %$stkFeatures > 3){
			open FA , '>'.catfile($outdir,'RNome.stkfa') or die $!;				
			for my $k (keys %$stkFeatures ){
				print FA '>'.$k."\n";
				$stkFeatures->{$k}=~s/(\W|_)/-/g;				
				print FA $_."\n" for unpack("(a80)*",$stkFeatures->{$k});				
			}
			close FA;
			
			my $ex = system('raxml -T '.$parameter->threads.' -f a -# 1000 -x 1234 -p 1234 -s '.catfile($outdir,'RNome.stkfa').' -w '.$outdir.' -n RNome.stk.tree -m GTRGAMMA -o '.join(',',grep { exists $stkFeatures->{$_} } @{$parameter->abbreviations}));
			&ABORT("raxml not found") unless $ex == 0;
			
			my $obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($outdir,'RAxML_bipartitions.RNome.stk.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($outdir,'RNome.stk.eps'));	
			
			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb->idToPath,$stamp.'-outgroup');			
		}				
	}
} 

#remove temp files
unlink $_ for glob catfile($parameter->tmp,'*');
system("rm -rf ".$parameter->tmp);

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$year = $year + 1900;
$mon += 1;
print "$mday.$mon.$year-$hour:$min:$sec\n";

sub run {
	print "Writing to ".$parameter->output."\n" if $parameter->verbose;

	my $c=0;
	for my $cfg (@{$parameter->queries}){		
		#parse the query related cfg file and store it into parameter object
		$parameter->set_cfg($cfg);		
		$c++;
		print $c.' of ',$#{$parameter->queries}+1,' - '.$parameter->cfg->rf."\n" if $parameter->verbose;

		#start screening if rfam query belongs to a kingdom of interest
		my $next=1;
		for (keys %{$parameter->kingdoms}){
			$next=0 if exists $parameter->cfg->kingdoms->{$_};	
		}
		if($next && ! $parameter->nofilter){			
			#print "- skipped due to kingdom restriction\n";
			next;
		}		
		
		#reverse sort to process infernal before blast
		my $thcalc=0;
		my ($threshold,$nonTaxThreshold);
		for my $tool (reverse sort @{$parameter->cfg->tools}){
			$tool = 'infernal' if $parameter->nofilter;
			$tool=~s/[\W\d_]//g;
			$tool = lc $tool;
			next if $tool eq 'blast' && $parameter->noblast;			
			$thcalc++ if $tool eq 'infernal' || $tool eq 'blast';
			$tool = ucfirst $tool;
			#print '- by '.$tool."\n";

			#dynamically initialize toolname dependent modules
			my $toolparser = lc $tool.'_parser';
			my $f = catfile('Bio','Gorap','Tool',$tool.'.pm');	
			my $class = 'Bio::Gorap::Tool::'.$tool;
			try {
		    	require $f;
			} catch {
				if (/^Can't locate .*?\.pm in \@INC/){
					print 'Unknown tool '.$tool.': using Default.pm'."\n";
					$toolparser = 'gff3_parser';
					$f = catfile('Bio','Gorap','Tool','Default.pm');
					$class = 'Bio::Gorap::Tool::Default';				
					require $f;
				} else {					
					&SIGABORT($_);
				}
			};
			$class->import;		
			
			#create an instance with related parser for change tool output into gff3 
			#with gorap attributes to store in gff database		
			my $obj = $class->new(
				threads => $parameter->threads - $thrListener->get_workload +1,
				parameter => $parameter,
				tool_parser => Bio::Gorap::Functions::ToolParser->can($toolparser) ? \&{'Bio::Gorap::Functions::ToolParser::'.$toolparser} : \&Bio::Gorap::Functions::ToolParser::gff3_parser,
				gffdb => $gffdb,
				fastadb => $fastadb,
				stkdb => $stkdb,
				bamdb => $bamdb,
				tool => $tool
			);
			
			#waits for resources
			$thrListener->push_obj($obj);

			#run software, use parser, store new gff3 entries 
			($threshold,$nonTaxThreshold) = $stkdb->calculate_threshold(($parameter->threads - $thrListener->get_workload)) if $thcalc == 1;			
			last if $threshold && $threshold == 999999;

			$obj->calc_features;
			last if $parameter->nofilter;
		}		
		next if $threshold && $threshold == 999999;
		next if $parameter->cfg->rf_rna=~/SU_rRNA/;
		my $sequences = $gffdb->get_sequences($parameter->cfg->rf_rna,$parameter->abbreviations);
		next if $#{$sequences} == -1;
		#print "align\n";
		if ($parameter->cfg->cm){
			my ($scorefile,$stk);
			try{
				($scorefile,$stk) = $stkdb->align(
					$parameter->cfg->rf_rna,
					$sequences,
					($parameter->threads - $thrListener->get_workload)
				);		
			} catch {
				&SIGABORT($_);
			};
			$gffdb->update_score_by_file($parameter->cfg->rf_rna,$scorefile);

			#start of time consuming single threaded background job with a subroutine reference,
			#which returns an array of String, parsable by pipe_parser
			#print "getfeatures\n";
			my $features = $gffdb->get_features($parameter->cfg->rf_rna,$parameter->abbreviations);	

			next if $#{$features} == -1;
			#print "bgjob\n";
			if ($thcalc){				
				$thrListener->calc_background(sub {$stkdb->filter_stk($parameter->cfg->rf_rna,$stk,$features,$threshold,$nonTaxThreshold,$taxdb,$gffdb)});
			} else {				
				$stkdb->store_stk($stk,catfile($parameter->output,'alignments',$parameter->cfg->rf_rna.'.stk'),$taxdb);
			}
		}
		
		#store annotations already in the database in case of errors		
		while($#{$thrListener->finished}>-1){
			$gffdb->store(shift @{$thrListener->finished});
			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb->idToPath,"$mday.$mon.$year-$hour:$min:$sec");
		}				
	}	
	
	$parameter->cfg( Bio::Gorap::CFG->new(rf_rna => 'NEW_RNA'));
	
	(Bio::Gorap::Tool::Piles->new(
		threads => $parameter->threads - $thrListener->get_workload,
		parameter => $parameter,
		tool_parser => \&Bio::Gorap::Functions::ToolParser::gff3_parser,
		gffdb => $gffdb,
		fastadb => $fastadb,
		stkdb => $stkdb,
		bamdb => $bamdb,
		tool => 'de_novo'
	))->calc_features;
}

sub get_phylo_features {
	my ($abbres) = @_;	

	my ($speciesSSU,$coreFeatures,$stkFeatures,$stkCoreFeatures);	
	my ($ssuToAbbr,$coreToAbbr,$rnomeToAbbr);

	for my $cfg (@{$parameter->queries}){
		$parameter->set_cfg($cfg);

		my $featureScore;	
		if ($cfg=~/_SSU_/){			
			for my $f ( @{$gffdb->get_features($parameter->cfg->rf_rna , $abbres, '!')}){				
				my @id = split /\./ , $f->seq_id;
				my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
				$ssuToAbbr->{$abbr} = 1;
				if (exists $featureScore->{$abbr}){
					my $s = ($f->get_tag_values('seq'))[0];
					if ($f->score > $featureScore->{$abbr}){
						$speciesSSU->{ $abbr } = $s;
						$featureScore->{$abbr} = $f->score;
					} elsif ($f->score == $featureScore->{$abbr} && length($s) > length($speciesSSU->{ $abbr })){
						$speciesSSU->{ $abbr } = $s;
						$featureScore->{$abbr} = $f->score;
					}					
				} else {
					$speciesSSU->{ $abbr } = ($f->get_tag_values('seq'))[0];
					$featureScore->{$abbr} = $f->score;					
				}
			}
		}		
		next if $cfg=~/_tRNA/ || $cfg=~/_rRNA/ || $cfg=~/CRISPR/;

		my $speciesFeature;	
		my $speciesSTKseq;	
		$featureScore={};

		for (@{$gffdb->get_features($parameter->cfg->rf_rna,$abbres,'!')}){			
			my @id = split /\./ , $_->seq_id;										
			my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
			if (exists $featureScore->{$abbr}){
				if ($_->score > $featureScore->{$abbr}){
					$speciesFeature->{ $abbr } = ($_->get_tag_values('seq'))[0];
					$speciesSTKseq->{$abbr} = ($stkdb->db->{$parameter->cfg->rf_rna}->get_seq_by_id($_->seq_id))->seq if exists $stkdb->db->{$parameter->cfg->rf_rna};
					$featureScore->{$abbr} = $_->score;
				}
			} else {
				$speciesFeature->{ $abbr } = ($_->get_tag_values('seq'))[0];
				$featureScore->{$abbr} = $_->score;
				$speciesSTKseq->{$abbr} = ($stkdb->db->{$parameter->cfg->rf_rna}->get_seq_by_id($_->seq_id))->seq if exists $stkdb->db->{$parameter->cfg->rf_rna};
			}	
		}

		next if scalar keys %$speciesFeature == 0;		
		
		my $l = length($speciesSTKseq->{(keys %$speciesSTKseq)[0]});
		for(@$abbres){			
			if(exists $speciesSTKseq->{$_}){
				$rnomeToAbbr->{$_}->{$parameter->cfg->rf_rna} = 1;
				$stkFeatures->{$_} .= $speciesSTKseq->{$_};
			} else {				
				$stkFeatures->{$_} .= join('',('?') x $l);
			}
		}

		next if scalar keys %$speciesFeature < ($#{$abbres} + 1);		

		for(@$abbres){
			$coreToAbbr->{$_}->{$parameter->cfg->rf_rna} = 1;
			$coreFeatures->{$_} .= $speciesFeature->{$_};			
			$stkCoreFeatures->{$_} .= $speciesSTKseq->{$_};
		}		
	}

	open TXT , '>'.catfile($parameter->output,'phylogeny','INFO') or die $!;	
	print TXT "SSU\n";
	print TXT $_."\n" for sort keys %$ssuToAbbr;
	print TXT "RNome\n";
	print TXT $_."\t".join("\t",sort keys %{$rnomeToAbbr->{$_}})."\n" for sort keys %$rnomeToAbbr;
	print TXT "coreRNome\n";
	print TXT $_."\t".join("\t",sort keys %{$coreToAbbr->{$_}})."\n" for sort keys %$coreToAbbr;
	close TXT;
	return ($speciesSSU,$coreFeatures,$stkFeatures,$stkCoreFeatures);
}

sub ABORT {
	print $_[0]."\n";
	unlink $_ for glob catfile($parameter->tmp,'*');
	system("rm -rf ".$parameter->tmp);
	exit 1;
}

sub SIGABORT {	
	$thrListener->stop if $thrListener;
	print "\n".'Safety store in progress..'."\n";
	$gffdb->store_overlaps if $gffdb;
	&ABORT;
}

__END__

=head1 NAME

GORAP - Genomewide ncRNA Annotation Pipeline

=head1 SYNOPSIS

Gorap.pl [OPTION]...
  
example: Gorap.pl -a x,y -i 1.fa,2.fa,*.fa -g x1.gff:x2.gff,,3.gff -c 1 -k bac -q 1:20,169,1852: -r 543 -s 'species name'

=head1 DESCRIPTION

For more parameters check also RNA family specific configuration files $GORAP/config/*.cfg
Please read the manual for improved annotations by the use of specific:
queries, 
covariance models, 
thresholds, 
structural properties
and how to add own software or scripts

B<-h>, B<--help>	

	this (help) message

B<-example>, B<--example>	

	Apply GORAP on $GORAP/example/ecoli.fa
	
B<-update>, B<--update>=I<all,rfam,ncbi,silva,cfg>	

	(optional, default: all) 
	updates internal used databases (Rfam, NCBI, Silva)
	!!! updating the configuration files will cause loss of all mismatch/indel settings
	
B<-file>, B<--file>=F<FILE> 
	
	(optional)
	(example file at: $GORAP/parameter/parameter.txt) 
	run GORAP with a parameter file
	following parameters will override settings from the parameter file
	
	
-------------------------
	
		
B<-i>, B<--fastas>=F<FILE>,...

	(optional, default $GORAP/example/ecoli.fa) 
	path(s) of comma separated species FASTA file(s) - wildcards allowed
		
B<-a>, B<--abbreviations>=I<abbreviation,> ...	

	(optional, default: build from FASTA file name(s)) 
	list of comma separated abbreviations as unique identifiers
	note: list lengths of B<-i> and B<-a> must be equal

B<-q>, B<--queries>=I<RF00001:RF00005,RF00008,> ...

	(default: all Rfam families)
	list of comma separated Rfam ids/numbers or ranges by ':'
	to enable additional features only and skip annotation set B<-q> 0

B<-k>, B<--kingdom>=I<bac,arc,euk,fungi,virus>

	(optional, default: all) 
	list of comma separated kingdoms to screen for kingdom specific ncRNAs only

B<-r>, B<--rank>=[I<INT/STRING>]	

	(optional) 
	NCBI taxonomy matching scientific name or taxonomy id of rank/genus/... for given species
	please escape white spaces i.e. Bacillus\\ subtilis\\ group or write 'Bacillus subtilis group'
	
B<-s>, B<--species>=[I<INT/STRING>]

	(optional) 
	NCBI taxonomy matching scientific name or taxonomy id of given species
	please escape white spaces i.e. Bacillus\\ subtilis or write 'Bacillus subtilis'
								
B<-og>, B<--outgroup>=F<FILE>

	(optional) 
	path to an additional species FASTA file to be used as outgroup for
	SSU rRNA based phylogeny reconstruction including all SSU rRNA annotations of B<-i>

B<-oga>, B<--ogabbreviations>=F<FILE>
	
	(optional)
	(default: build from FASTA file name of B<-og> 
	list of comma separated abbreviations as unique identifiers
	
B<-b>, B<--bams>=F<FILE>,...

	(optional) 
	paths(s) of comma separated list of colon separated indexed BAM file(s) - wildcards allowed

B<minl>, B<--minlength>=I<INT>
	
	(optional, default 50)
	minimum length of de novo predicted genes

B<minh>, B<--minheigth>=I<INT>

	(optional, default 1000)
	minimum number of reads at the same locus for de novo gene prediction

B<-g>, B<--gffs>=F<FILE>,...

	(optional) 
	paths(s) of comma separated list of colon separated GFF3 file(s) - wildcards allowed
	
B<-o>, B<--output>=F<PATH>	

	(optional, default: <working directory>/gorap_out) 
	output directory

B<-c>, B<--cpu>=I<INT>	

	(optional, default: 1) 
	number of threads (cpu cores to use)

B<-t>, B<--tmp>=F<PATH>	

	(optional, default: $TMPDIR or /tmp or $GORAP/tmp) 
	set the temporary directory - will be removed afterwards
	
B<-sort>, B<--sort>
	
	(optional)
	enable resulting alignments to be sorted taxonomical by given rank or species 

B<-notax>, B<--notaxonomy>
	
	(optional)
	disables taxonomic sorting and filter using rank and species information

B<-noo>, B<--nooverlap>
	
	(optional)
	disables deletion of de novo predictions if they overlap with a given GFF3 file

B<-notpm>, B<--notpm>
	
	(optional)
	disables FPKM and TPM calculation

B<-nofi>, B<--nofilter>
	
	(optional)
	disables GORAP specific sequence and structure filter

=head1 AUTHOR

Konstantin Riege, E<lt>konstantin.riege@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 by Konstantin Riege

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.14.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 SEE ALSO

The full description of all in-/outputs and parameters is maintained as PDF manual. 
You can access it at L<www.rna.uni-jena.de/software.php>.

=cut


