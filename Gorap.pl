#! /usr/bin/perl

use strict;
use warnings;

#TODO remove 
use lib 'lib';

use Bio::Gorap::ThrListener;
use Bio::Gorap::Parameter;
use Bio::Gorap::DB::GFF;
use Bio::Gorap::DB::STK;
use Bio::Gorap::DB::Fasta;
use Bio::Gorap::DB::BAM;
use Bio::Gorap::DB::Taxonomy;
use Bio::Gorap::Evaluation::HTML;
use Bio::Gorap::Functions::ToolParser;
use File::Spec::Functions;
use File::Basename;
use Try::Tiny;
use Cwd 'abs_path';
use Bio::Tree::Draw::Cladogram;
use Bio::TreeIO;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$year = $year + 1900;
$mon += 1;
print "$mday.$mon.$year-$hour:$min:$sec\n";
my $stamp = "$mday.$mon.$year-$hour:$min:$sec";

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

my $bamdb = Bio::Gorap::DB::BAM->new(
	parameter => $parameter
);

#gorap gff3 storage database initialization without loss of existing data in output directory 
my $gffdb = Bio::Gorap::DB::GFF->new(
	parameter => $parameter,
	bamdb => $bamdb
); 

my $fastadb = Bio::Gorap::DB::Fasta->new(
	parameter => $parameter
);

my $stkdb = Bio::Gorap::DB::STK->new(
	parameter => $parameter
);

my $taxdb = Bio::Gorap::DB::Taxonomy->new(
	parameter => $parameter
);
$taxdb->findRelatedSpecies;

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
#store final annotation results
$gffdb->store;
Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$fastadb->oheaderToDBsize,$stkdb->idToPath,$stamp);

if ($parameter->has_outgroups){	
		
	my @newQ;	
	my @oldqueries = @{$parameter->queries};
	my @genomes = @{$parameter->genomes};
	my @abbres = @{$parameter->abbreviations};

	$parameter->set_queries();
	for my $cfg (@{$parameter->queries}){
		$parameter->set_cfg($cfg);
		if ($cfg=~/_SSU_/ || ! ($cfg=~/_tRNA/ || $cfg=~/_rRNA/ || $cfg=~/CRISPR/)){
			push @newQ , $parameter->cfg->rf if $#{$gffdb->get_all_features($parameter->cfg->rf_rna , '!')} > -1;
		}
	}
	
	if ($#newQ > 1){	
		
		unlink $_ for glob catfile($parameter->output,'RAxML_*');			

		$parameter->set_genomes($parameter->outgroups,$parameter->ogabbreviations);

		push @genomes, @{$parameter->genomes};
		push @abbres, @{$parameter->abbreviations};	

		$parameter->set_queries(\@newQ);				
		
		if ($#oldqueries > -1){
			print "\nAnnotation of outgroups for phylogeny reconstruction\n";
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
			#store final annotation results
			$gffdb->store;
			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$fastadb->oheaderToDBsize,$stkdb->idToPath,$stamp.'-outgroup');
		}
				
		my ($speciesSSU,$coreFeatures,$stkSSU,$stkFeatures,$stkCoreFeatures) = &get_phylo_features(\@abbres);

		if (any { exists $speciesSSU->{$_} } @{$parameter->abbreviations} && scalar keys %$speciesSSU > 3){
			open FA , '>'.catfile($parameter->output,'SSU.fasta') or die $!;	
			print FA '>'.$_."\n".$speciesSSU->{$_}."\n" for keys %$speciesSSU;
			close FA;
			open FA , '>'.catfile($parameter->output,'SSU.stkfa') or die $!;	
			print FA '>'.$_."\n".$stkSSU->{$_}."\n" for keys %$stkSSU;
			close FA;

			system('mafft --localpair --maxiterate 1000 --thread '.$parameter->threads.' '.catfile($parameter->output,'SSU.fasta').' > '.catfile($parameter->output,'SSU.mafft'));
			system('raxml -T '.$parameter->threads.' -f a -# 100 -x 1234 -p 1234 -s '.catfile($parameter->output,'SSU.mafft').' -w '.$parameter->output.' -n SSU.mafft.tree -m GTRGAMMA -o '.join(',',grep { exists $speciesSSU->{$_} } @{$parameter->abbreviations}));
			system('raxml -T '.$parameter->threads.' -f a -# 100 -x 1234 -p 1234 -s '.catfile($parameter->output,'SSU.stkfa').' -w '.$parameter->output.' -n SSU.stk.tree -m GTRGAMMA -o '.join(',',grep { exists $speciesSSU->{$_} } @{$parameter->abbreviations}));
			
			my $obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($parameter->output,'RAxML_bipartitions.SSU.mafft.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($parameter->output,'SSU.mafft.eps'));	
			$obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($parameter->output,'RAxML_bipartitions.SSU.stk.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($parameter->output,'SSU.stk.eps'));	
			my $test = `which convert`;
			system('convert -density 300 '.catfile($parameter->output,'SSU.mafft.eps').' '.catfile($parameter->output,'SSU.mafft.png')) if $test;
			system('convert -density 300 '.catfile($parameter->output,'SSU.stk.eps').' '.catfile($parameter->output,'SSU.stk.png')) if $test;

			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$fastadb->oheaderToDBsize,$stkdb->idToPath,$stamp);			
		} 

		if (any { exists $coreFeatures->{$_} } @{$parameter->abbreviations} && scalar keys %$coreFeatures > 3){
			open FA , '>'.catfile($parameter->output,'coreRNome.fasta') or die $!;	
			print FA '>'.$_."\n".$coreFeatures->{$_}."\n" for keys %$coreFeatures;
			close FA;
			open FA , '>'.catfile($parameter->output,'coreRNome.stkfa') or die $!;				
			print FA '>'.$_."\n".$stkCoreFeatures->{$_}."\n" for keys %$stkCoreFeatures;			
			close FA;

			system('mafft --localpair --maxiterate 1000 --thread '.$parameter->threads.' '.catfile($parameter->output,'coreRNome.fasta').' > '.catfile($parameter->output,'coreRNome.mafft'));
			system('raxml -T '.$parameter->threads.' -f a -# 100 -x 1234 -p 1234 -s '.catfile($parameter->output,'coreRNome.mafft').' -w '.$parameter->output.' -n coreRNome.mafft.tree -m GTRGAMMA -o '.join(',',grep { exists $coreFeatures->{$_} } @{$parameter->abbreviations}));
			system('raxml -T '.$parameter->threads.' -f a -# 100 -x 1234 -p 1234 -s '.catfile($parameter->output,'coreRNome.stkfa').' -w '.$parameter->output.' -n coreRNome.stk.tree -m GTRGAMMA -o '.join(',',grep { exists $coreFeatures->{$_} } @{$parameter->abbreviations}));
			
			my $obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($parameter->output,'RAxML_bipartitions.coreRNome.mafft.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($parameter->output,'coreRNome.mafft.eps'));	
			$obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($parameter->output,'RAxML_bipartitions.coreRNome.stk.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($parameter->output,'coreRNome.stk.eps'));	
			my $test = `which convert`;
			system('convert -density 300 '.catfile($parameter->output,'coreRNome.mafft.eps').' '.catfile($parameter->output,'coreRNome.mafft.png')) if $test;
			system('convert -density 300 '.catfile($parameter->output,'coreRNome.stk.eps').' '.catfile($parameter->output,'coreRNome.stk.png')) if $test;

			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$fastadb->oheaderToDBsize,$stkdb->idToPath,$stamp);			
		} 

		if (any { exists $stkFeatures->{$_} } @{$parameter->abbreviations} && scalar keys %$stkFeatures > 3){
			open FA , '>'.catfile($parameter->output,'RNome.stkfa') or die $!;	
			print FA '>'.$_."\n".$stkFeatures->{$_}."\n" for keys %$stkFeatures;
			close FA;
			
			system('raxml -T '.$parameter->threads.' -f a -# 100 -x 1234 -p 1234 -s '.catfile($parameter->output,'RNome.stkfa').' -w '.$parameter->output.' -n RNome.stk.tree -m GTRGAMMA -o '.join(',',grep { exists $stkFeatures->{$_} } @{$parameter->abbreviations}));
			
			my $obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($parameter->output,'RAxML_bipartitions.RNome.stk.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($parameter->output,'RNome.stk.eps'));	
			my $test = `which convert`;			
			system('convert -density 300 '.catfile($parameter->output,'RNome.stk.eps').' '.catfile($parameter->output,'coreRNome.stk.png')) if $test;

			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$fastadb->oheaderToDBsize,$stkdb->idToPath,$stamp);			
		}				
	}
} 

#remove temp files
unlink $_ for glob catfile($parameter->tmp,$parameter->pid.'.*');

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$year = $year + 1900;
$mon += 1;
print "$mday.$mon.$year-$hour:$min:$sec\n";

sub run {
	my $c=0;
	for my $cfg (@{$parameter->queries}){
		#parse the query related cfg file and store it into parameter object
		$parameter->set_cfg($cfg);	
		$c++;
		print $c.' of ',$#{$parameter->queries}+1,' - '.$parameter->cfg->rf."\n";

		#start screening if rfam query belongs to a kingdom of interest
		my $next=1;
		for (keys %{$parameter->kingdoms}){
			$next=0 if exists $parameter->cfg->kingdoms->{$_};	
		}
		if($next){
			#print "- skipped due to kingdom restriction\n";
			next;
		}
		
		#reverse sort to process infernal before blast
		my $thcalc;
		for my $tool (reverse sort @{$parameter->cfg->tools}){
			$thcalc = 1 if $tool eq 'infernal' || $tool eq 'blast';
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
					$toolparser = 'gff3_parser';
					$f = catfile('Bio','Gorap','Tool','Default.pm');
					$class = 'Bio::Gorap::Tool::Default';				
					require $f;
				} else {
					print $_;
				}
			};
			$class->import;		
			
			#create an instance with related parser for change tool output into gff3 
			#with gorap attributes to store in gff database		
			my $obj = $class->new(
				threads => $parameter->threads - $thrListener->get_workload,
				parameter => $parameter,
				tool_parser => Bio::Gorap::Functions::ToolParser->can($toolparser) ? \&{'Bio::Gorap::Functions::ToolParser::'.$toolparser} : \&Bio::Gorap::Functions::ToolParser::gff3_parser,
				gffdb => $gffdb,
				fastadb => $fastadb,
				stkdb => $stkdb,				
				tool => $tool
			);		
			
			#waits for resources
			$thrListener->push_obj($obj);

			#run software, use parser, store new gff3 entries 
			# print "calc\n";
			$obj->calc_features;		
		}
		#print "getseqs\n";
		next if $parameter->cfg->rf_rna=~/SU_rRNA/;
		my $sequences = $gffdb->get_sequences($parameter->cfg->rf_rna,$parameter->abbreviations);
		next if $#{$sequences} == -1;
		#print "align\n";
		my ($scorefile,$stk) = $stkdb->align(
			$parameter->cfg->rf_rna,
			$sequences,
			($parameter->threads - $thrListener->get_workload)
		);
		#print "updatescores\n";
		$gffdb->update_score_by_file($parameter->cfg->rf_rna,$scorefile);

		#start of time consuming single threaded background job with a subroutine reference,
		#which returns an array of String, parsable by pipe_parser
		#print "getfeatures\n";
		my $features = $gffdb->get_features($parameter->cfg->rf_rna,$parameter->abbreviations);	
		next if $#{$features} == -1;
		#print "bgjob\n";
		if ($thcalc){
			my ($threshold,$nonTaxThreshold) = $stkdb->calculate_threshold(($parameter->threads - $thrListener->get_workload),$taxdb->relatedRankIDsToLineage,$taxdb->relatedSpeciesIDsToLineage);	
			$thrListener->calc_background(sub {$stkdb->filter_stk($parameter->cfg->rf_rna,$stk,$features,$threshold,$nonTaxThreshold,$taxdb)});	
		} else {
			$thrListener->calc_background(sub {$stkdb->scorefilter_stk($parameter->cfg->rf_rna,$stk,$features,0)});
			# $stkdb->store_stk($stk,catfile($parameter->output,'alignments',$parameter->cfg->rf_rna.'.stk'),$taxdb);
		}
		#print "store\n";
		#store annotations already in the database in case of errors
		$gffdb->store;
		#print "eval\n";
		Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$fastadb->oheaderToDBsize,$stkdb->idToPath,"$mday.$mon.$year-$hour:$min:$sec");
	}
}

sub get_phylo_features {
	my ($abbres) = @_;

	my ($speciesSSU,$coreFeatures,$stkSSU,$stkFeatures,$stkCoreFeatures);	

	for my $cfg (@{$parameter->queries}){
		$parameter->set_cfg($cfg);
		my $featureScore;	
		if ($cfg=~/_SSU_/){
			
			for (@{$gffdb->get_features($parameter->cfg->rf_rna , $abbres, '!')}){
				my @id = split /\./ , $_->seq_id;										
				my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
				if (exists $featureScore->{$abbr}){
					if ($_->score > $featureScore->{$abbr}){
						$speciesSSU->{ $abbr } = ($_->get_tag_values('seq'))[0];
						$stkSSU->{$abbr} = ($stkdb->{$parameter->cfg->rf_rna}->get_seq_by_id($_->seq_id))->seq;	
					}					
				} else {
					$speciesSSU->{ $abbr } = ($_->get_tag_values('seq'))[0];
					$featureScore->{$abbr} = $_->score;
					$stkSSU->{$abbr} = ($stkdb->{$parameter->cfg->rf_rna}->get_seq_by_id($_->seq_id))->seq;	
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
					$speciesSTKseq->{$abbr} = ($stkdb->{$parameter->cfg->rf_rna}->get_seq_by_id($_->seq_id))->seq;	
				}
			} else {
				$speciesFeature->{ $abbr } = ($_->get_tag_values('seq'))[0];
				$featureScore->{$abbr} = $_->score;
				$speciesSTKseq->{$abbr} = ($stkdb->{$parameter->cfg->rf_rna}->get_seq_by_id($_->seq_id))->seq;
			}	
		}

		next if scalar keys %$speciesFeature == 0;		
		
		my $l = length($speciesSTKseq->{(keys %$$speciesSTKseq)[0]});
		for(@$abbres){
			if(exists $speciesSTKseq->{$_}){					
				$stkFeatures->{$_} .= $speciesSTKseq->{$_};
			} else {			
				$stkFeatures->{$_} .= join('',('?') x $l);
			}
		}

		next if scalar keys %$speciesFeature < ($#{$abbres} + 1);		

		for(@$abbres){
			$coreFeatures->{$_} .= $speciesFeature->{$_};			
			$stkCoreFeatures->{$_} .= $speciesSTKseq->{$_};
		}		
	}

	return ($speciesSSU,$coreFeatures,$stkSSU,$stkFeatures,$stkCoreFeatures);
}


__END__

=head1 NAME

GORAP - Genomewide ncRNA Annotation Pipeline

=head1 SYNOPSIS

Gorap.pl [OPTION]...
  
example: Gorap.pl -a ecoli -i $GORAP/example/ecoli.fa -o results -c 1 -k bac -q RF00001:RF00020,RF00169,RF001854: -r 543 -s Escherichia\\ coli

=head1 DESCRIPTION

For more parameters check also RNA family specific configuration files $GORAP/parameter/families/*.cfg
Please read the manual for improved annotations by the use of specific:
queries, 
covariance models, 
thresholds, 
structural properties
and how to add own software or scripts

B<-h>, B<--help>	

	this (help) message
	
B<-update>, B<--update>=I<all,rfam,ncbi,silva,cfg>	

	(optional, default: all) 
	updates internal used databases (Rfam, NCBI, Silva)
	
B<-file>, B<--file>=F<FILE> 
	
	(optional)
	(default: E. coli example with parameter file at: $GORAP/parameter/parameter.txt) 
	run GORAP with a parameter file
	
	
-------------------------
	
		
B<-i>, B<--fastas>=F<FILE>,...

	(requierd) 
	(regex) path(s) of comma separated species FASTA file(s)
		
B<-a>, B<--abbreviations>=I<abbreviation,> ...	

	(default: build from FASTA file name(s)) 
	list of comma separated abbreviations as unique identifiers
	
B<-o>, B<--output>=F<PATH>	

	(default: $PWD/gorap_out) 
	output directory

B<-t>, B<--tmp>=F<PATH>	

	(default: $TMPDIR or /tmp or $GORAP/tmp) 
	set the temporary directory
	
B<-c>, B<--cpu>=I<INT>	

	(default: 1) 
	count of cpu cores to use
	
B<-k>, B<--kingdom>=I<bac,arc,euk,fungi,virus>

	(default: all) 
	list of comma separated kingdoms to screen for kingdom specific ncRNAs only
	
B<-q>, B<--queries>=I<RF00001:RF00005,RF00008,> ...

	(default: all Rfam families)
	list of comma separated Rfam ids/numbers or ranges by ':'
	to enable additional features only, skip annotation with 0
	 
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
	SSU rRNA based phylogeny reconstruction including all given FASTA files of B<-i>

B<-oga>, B<--ogabbreviation>=F<FILE>
	
	(optional)
	(default: build from FASTA file name of B<-og> 
	abbreviation as unique identifiers
	
B<-b>, B<--bam>=F<FILE>,...

	(optional) 
	(regex) paths(s) of comma separated mapping results as SAM/BAM file(s)	

-sort, --sort
	
	(optional)
	enable resulting alignments to be sorted taxonomical by given rank or species 

-force, --force
	
	(optional)
	avoid parameter restrictions for debugging and own implementations

=head1 AUTHOR

Konstantin Riege, E<lt>konstantin.riege@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 by Konstantin Riege

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.14.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 SEE ALSO

The full description of all in-/outputs and parameters is maintained as PDF manual. 
You can access it on L<www.rna.uni-jena.de/software.php>.

=cut


