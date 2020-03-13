#!/usr/bin/env perl
use lib './lib'; # devel

use v5.10;
use strict;
use warnings;
use sigtrap qw(handler safety_store normal-signals);
use Cwd qw(abs_path);
use File::Spec::Functions;
use File::Basename;
use File::Copy qw(copy);
use File::Path qw(make_path remove_tree);
use List::Util qw(any);
use Try::Tiny;
use Hash::Merge qw(merge);

use Bio::Gorap::Parameter;
use Bio::Gorap::ThrListener;
use Bio::Gorap::DB::GFF;
use Bio::Gorap::DB::STK;
use Bio::Gorap::DB::Fasta;
use Bio::Gorap::DB::BAM;
use Bio::Gorap::DB::Taxonomy;
use Bio::Gorap::Evaluation::HTML;
use Bio::Gorap::Tool::Piles;
use Bio::Gorap::CFG;
use Bio::Gorap::Functions::ToolParser;
use Bio::Tree::Draw::Cladogram;
use Bio::TreeIO;

BEGIN {
	if ( ! $ENV{GORAP} ){
		say "Export GORAP environment variable pointing towards installation directory and try again!";
		exit 1;
	}
	my @path;
	push @path,abs_path($_) for glob("$ENV{GORAP}/*/*/bin");

	if($#path == -1){
		say "Please ensure Gorap is correctly installed and GORAP environment variable points towards installation directory!";
		exit 1;
	}
	
	$ENV{PATH} = $ENV{PATH} ? join(":",($ENV{PATH},@path)) : join(":",@path);
}

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$year += 1900;
$mon += 1;
my $stamp = "$mday.$mon.$year-$hour:$min:$sec";

my $parameter = Bio::Gorap::Parameter->new(
	pwd => $ENV{PWD},
	pid => $$,
	commandline => 1,
	label => $stamp,
);

say "Gorap started - $stamp";

my $taxdb;
my $stkdb;

if ($parameter->refresh){
	print "Updating GFF and FASTA files from Stockholm alignments\n";

	my $fastadb = Bio::Gorap::DB::Fasta->new(
		parameter => $parameter,
		do_chunks => 0
	);

	$stkdb = Bio::Gorap::DB::STK->new(
		parameter => $parameter
	);

	my $bamdb = Bio::Gorap::DB::BAM->new(
		parameter => $parameter
	);

	my $gffdb = Bio::Gorap::DB::GFF->new(
		parameter => $parameter,
		bamdb => $bamdb,
		fastadb => $fastadb
	);

	for my $cfg (@{$parameter->queries}){
		$parameter->set_cfg($cfg);
		my $type = $parameter->cfg->rf_rna;
		next if $type=~/(L|S)SU_rRNA/;
		my $hold = 0;
		for my $f (@{$gffdb->get_all_features($type)}){
			my $seq;
			$seq = $stkdb->db->{$type}->get_seq_by_id($f->seq_id) if exists $stkdb->db->{$type};
			if ($seq){
				$hold=1;
				$gffdb->update_filter($f->seq_id,$type,'!');
			} else {
				$gffdb->update_filter($f->seq_id,$type,'X') if $f->display_name eq '!';
			}
		}
		unlink $stkdb->idToPath->{$type} if ! $hold && exists $stkdb->db->{$type};
	}
	print "Storing changes\n";
	$gffdb->store_overlaps;

	Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);
	print "\nResults stored with label: ".$parameter->label."\n";

	remove_tree($parameter->tmp,{verbose => 0, error  => \my $err});

	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
	$year = $year + 1900;
	$mon += 1;
	print "$mday.$mon.$year-$hour:$min:$sec\n";

	exit 0;
}

my $fastadb = Bio::Gorap::DB::Fasta->new(
	parameter => $parameter
);

if ($parameter->taxonomy){
	$taxdb = Bio::Gorap::DB::Taxonomy->new(
		parameter => $parameter
	);
	$taxdb->findRelatedSpecies unless $parameter->skip_comp;
	$stkdb = Bio::Gorap::DB::STK->new(
		parameter => $parameter,
		taxonomy => $taxdb
	);
} else {
	$stkdb = Bio::Gorap::DB::STK->new(
		parameter => $parameter
	);
}

my $bamdb = Bio::Gorap::DB::BAM->new(
	parameter => $parameter,
);

#gorap gff3 storage database initialization without loss of existing data in output directory
my $gffdb = Bio::Gorap::DB::GFF->new(
	parameter => $parameter,
	bamdb => $bamdb,
	fastadb => $fastadb
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

unless ($parameter->skip_comp){
	&run();
} else {
	$stkdb->store if $parameter->sort;
}

if ($parameter->has_outgroups){

	print "Preparing phylogeny reconstruction\n" if $parameter->verbose;
	my @newQ;
	my @oldqueries = @{$parameter->queries};
	my @oldgenomes = @{$parameter->genomes};
	my @oldabbres = @{$parameter->abbreviations};
	my @newgenomes = @{$parameter->outgroups};
	my @newabbres = @{$parameter->ogabbreviations};

	$parameter->set_queries();
	for my $cfg (@{$parameter->queries}){
		$parameter->set_cfg($cfg);
		push @newQ , $parameter->cfg->rf if $#{$gffdb->get_all_features($parameter->cfg->rf_rna , '!')} > -1;
	}

	if ($#newQ > -1){
		my $outdir = catdir($parameter->output,'phylogeny-'.$parameter->label);
		make_path($outdir);
		unlink $_ for glob catfile($outdir,'RAxML_*');

		$parameter->set_queries(\@newQ);

		for my $a (@newabbres){ # do not use $_ - somewhere in the code it will be overridden and $_ entries becomes unitialized
			$gffdb->add_db($a);
		}

		unless ($parameter->skip_comp){
			print "\nAnnotating outgroups for phylogeny reconstruction\n" if $parameter->verbose;

			$parameter->set_genomes(\@newgenomes,\@newabbres); #new referenfce, otherwise og will be unset during set_genomes

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

			$parameter->set_genomes(\@oldgenomes,\@oldabbres);
			$parameter->set_queries(\@oldqueries); # reset for html output
		}

		my @allabbres = (@oldabbres,@newabbres);
		my ($speciesSSU,$coreFeatures,$stkFeatures,$stkCoreFeatures,$stk50Features) = &get_phylo_features(\@allabbres,$outdir);

		Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);

		if ( (any { exists $speciesSSU->{$_} } @newabbres) && scalar keys %$speciesSSU > 3){
			say "raxml on SSU";
			open FA , '>'.catfile($outdir,'SSU.fasta') or die $!;
			for my $k (keys %$speciesSSU ){
				print FA '>'.$k."\n";
				$speciesSSU->{$k}=~s/(\W|_)/-/g;
				print FA $_."\n" for unpack("(a80)*",$speciesSSU->{$k});
			}
			close FA;

			system('mafft --localpair --maxiterate 1000 --thread '.$parameter->threads.' '.catfile($outdir,'SSU.fasta').' > '.catfile($outdir,'SSU.mafft'));
			system('raxmlHPC-PTHREADS-AVX2 -T '.$parameter->threads.' -f a -# 1000 -x 1234 -p 1234 -s '.catfile($outdir,'SSU.mafft').' -w '.$outdir.' -n SSU.mafft.tree -m GTRGAMMA -o '.join(',',grep { exists $speciesSSU->{$_} } @newabbres));

			my $obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($outdir,'RAxML_bipartitionsBranchLabels.SSU.mafft.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($outdir,'SSU.mafft.eps'));
			copy catfile($outdir,'RAxML_bipartitionsBranchLabels.SSU.mafft.tree'), catfile($outdir,'SSU.mafft.tree');
			system('newicktopdf -pc 1 -boot -notitle '.catfile($outdir,'RAxML_bipartitionsBranchLabels.SSU.mafft.tree'));
			copy catfile($outdir,'RAxML_bipartitionsBranchLabels.pdf'), catfile($outdir,'SSU.mafft.pdf');

			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);
		}

		if ( (any { exists $coreFeatures->{$_} } @newabbres) && scalar keys %$coreFeatures > 3){
			say "raxml on coreRNome";
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

			system('mafft --localpair --maxiterate 1000 --thread '.$parameter->threads.' '.catfile($outdir,'coreRNome.fasta').' > '.catfile($outdir,'coreRNome.mafft'));
			system('raxmlHPC-PTHREADS-AVX2 -T '.$parameter->threads.' -f a -# 1000 -x 1234 -p 1234 -s '.catfile($outdir,'coreRNome.mafft').' -w '.$outdir.' -n coreRNome.mafft.tree -m GTRGAMMA -o '.join(',',grep { exists $coreFeatures->{$_} } @newabbres));
			system('raxmlHPC-PTHREADS-AVX2 -T '.$parameter->threads.' -f a -# 1000 -x 1234 -p 1234 -s '.catfile($outdir,'coreRNome.stkfa').' -w '.$outdir.' -n coreRNome.stk.tree -m GTRGAMMA -o '.join(',',grep { exists $coreFeatures->{$_} } @newabbres));

			my $obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($outdir,'RAxML_bipartitionsBranchLabels.coreRNome.mafft.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($outdir,'coreRNome.mafft.eps'));
			copy catfile($outdir,'RAxML_bipartitionsBranchLabels.coreRNome.mafft.tree'), catfile($outdir,'coreRNome.mafft.tree');
			system('newicktopdf -pc 1 -boot -notitle '.catfile($outdir,'RAxML_bipartitionsBranchLabels.coreRNome.mafft.tree'));
			copy catfile($outdir,'RAxML_bipartitionsBranchLabels.pdf'), catfile($outdir,'coreRNome.mafft.pdf');

			$obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($outdir,'RAxML_bipartitionsBranchLabels.coreRNome.stk.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($outdir,'coreRNome.stk.eps'));
			copy catfile($outdir,'RAxML_bipartitionsBranchLabels.coreRNome.stk.tree'), catfile($outdir,'coreRNome.stk.tree');
			system('newicktopdf -pc 1 -boot -notitle '.catfile($outdir,'RAxML_bipartitionsBranchLabels.coreRNome.stk.tree'));
			copy catfile($outdir,'RAxML_bipartitionsBranchLabels.pdf'), catfile($outdir,'coreRNome.stk.pdf');

			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);
		}

		if ( (any { exists $stkFeatures->{$_} } @newabbres) && scalar keys %$stkFeatures > 3){
			say "raxml on RNome";
			open FA , '>'.catfile($outdir,'RNome.stkfa') or die $!;
			for my $k (keys %$stkFeatures ){
				print FA '>'.$k."\n";
				$stkFeatures->{$k}=~s/(\W|_)/-/g;
				print FA $_."\n" for unpack("(a80)*",$stkFeatures->{$k});
			}
			close FA;

			system('raxmlHPC-PTHREADS-AVX2 -T '.$parameter->threads.' -f a -# 1000 -x 1234 -p 1234 -s '.catfile($outdir,'RNome.stkfa').' -w '.$outdir.' -n RNome.stk.tree -m GTRGAMMA -o '.join(',',grep { exists $stkFeatures->{$_} } @newabbres));

			my $obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($outdir,'RAxML_bipartitionsBranchLabels.RNome.stk.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($outdir,'RNome.stk.eps'));
			copy catfile($outdir,'RAxML_bipartitionsBranchLabels.RNome.stk.tree'), catfile($outdir,'RNome.stk.tree');
			system('newicktopdf -pc 1 -boot -notitle '.catfile($outdir,'RAxML_bipartitionsBranchLabels.RNome.stk.tree'));
			copy catfile($outdir,'RAxML_bipartitionsBranchLabels.pdf'), catfile($outdir,'RNome.stk.pdf');

			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);
		}

		if ( (any { exists $stk50Features->{$_} } @newabbres) && scalar keys %$stk50Features > 3){
			open FA , '>'.catfile($outdir,'core50RNome.stkfa') or die $!;
			for my $k (keys %$stk50Features ){
				print FA '>'.$k."\n";
				$stk50Features->{$k}=~s/(\W|_)/-/g;
				print FA $_."\n" for unpack("(a80)*",$stk50Features->{$k});
			}
			close FA;

			system('raxmlHPC-PTHREADS-AVX2 -T '.$parameter->threads.' -f a -# 1000 -x 1234 -p 1234 -s '.catfile($outdir,'core50RNome.stkfa').' -w '.$outdir.' -n core50RNome.stk.tree -m GTRGAMMA -o '.join(',',grep { exists $stk50Features->{$_} } @newabbres));

			my $obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($outdir,'RAxML_bipartitionsBranchLabels.core50RNome.stk.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($outdir,'core50RNome.stk.eps'));
			copy catfile($outdir,'RAxML_bipartitionsBranchLabels.core50RNome.stk.tree'), catfile($outdir,'core50RNome.stk.tree');
			system('newicktopdf -pc 1 -boot -notitle '.catfile($outdir,'RAxML_bipartitionsBranchLabels.core50RNome.stk.tree'));
			copy catfile($outdir,'RAxML_bipartitionsBranchLabels.pdf'), catfile($outdir,'core50RNome.stk.pdf');

			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);
		}
	}
}

sub run {
	print "Writing to ".$parameter->output."\n" if $parameter->verbose;

	my $c=0;
	for my $cfg (@{$parameter->queries}){
		#parse the query related cfg file and store it into parameter object
		$parameter->set_cfg($cfg);
		$c++;
		print $c.' of ',$#{$parameter->queries}+1,' - '.$parameter->cfg->rf_rna if $parameter->verbose;
		#start screening if rfam query belongs to a kingdom of interest
		my $next=1;
		for (keys %{$parameter->kingdoms}){
			$next=0 if exists $parameter->cfg->kingdoms->{$_};
		}
		if($next && ! $parameter->nofilter){
			say " - skipped due to kingdom restriction" if $parameter->verbose;
			next;
		}

		my ($threshold,$nonTaxThreshold);
		my @tools = @{$parameter->cfg->tools};
		@tools = ('infernal') if $parameter->rfamscan || $parameter->pureinfernal;
		for my $tool (@tools){
			my $cmd = $parameter->cfg->cmd->{$tool};
			$tool=~s/[\W\d_]//g;
			$tool = lc $tool;
			$tool = ucfirst $tool;
			next if $tool eq 'Blast' && $parameter->noblast;

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
					&safety_store($_);
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
				tool => $tool,
				cmd => $cmd
			);

			if ($tool eq 'Infernal' || $tool eq 'Blast'){
				($threshold,$nonTaxThreshold) = $stkdb->calculate_threshold($parameter->threads - $thrListener->get_workload);
				last if $threshold == 999999; #cm is biases by species not related to input
			}
			$obj->calc_features;
		}
		if ($threshold && $threshold == 999999){  #cm is biases by species not related to input
			say " - skipped due to heavy bias of unrelated species" if $parameter->verbose;
			next;
		}
		say "" if $parameter->verbose;
		next if $parameter->cfg->rf_rna=~/SU_rRNA/; #do not align long RNAs
		my $sequences = $gffdb->get_sequences($parameter->cfg->rf_rna,$parameter->abbreviations);
		next if $#{$sequences} == -1;

		# print "align\n";
		if ($parameter->cfg->cm){
			my ($scorefile,$stk);
			try{
				($scorefile,$stk) = $stkdb->align(
					$parameter->cfg->rf_rna,
					$sequences,
					($parameter->threads - $thrListener->get_workload)
				);
			} catch {
				&safety_store($_);
			};
			$gffdb->update_score_by_file($parameter->cfg->rf_rna,$scorefile);

			#start of time consuming single threaded background job with a subroutine reference,
			#which returns an array of String, parsable by thredlisteners storage_saver - the gff objects function to update the filter tag
			if ($threshold){
				my $features = $gffdb->get_features($parameter->cfg->rf_rna,$parameter->abbreviations);
				next if $#{$features} == -1;
				$thrListener->start(sub {$stkdb->filter_stk($parameter->cfg->rf_rna,$stk,$features,$threshold,$nonTaxThreshold,$taxdb,$gffdb)});
			} else {
				$stkdb->remove_gap_columns_and_write($stk,catfile($parameter->output,'alignments',$parameter->cfg->rf_rna.'.stk'),$taxdb);
				$gffdb->store($parameter->cfg->rf_rna);
			}
		}
		#store annotations already in the database in case of errors
		while($#{$thrListener->finished}>-1){
			$gffdb->store(shift @{$thrListener->finished});
			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);
		}
	}

	if ($parameter->has_bams){
		$parameter->cfg( Bio::Gorap::CFG->new(rf_rna => 'NEW_RNA'));
		(Bio::Gorap::Tool::Piles->new(
			threads => $parameter->threads - $thrListener->get_workload +1,
			parameter => $parameter,
			tool_parser => \&Bio::Gorap::Functions::ToolParser::gff3_parser,
			gffdb => $gffdb,
			fastadb => $fastadb,
			stkdb => $stkdb,
			bamdb => $bamdb,
			tool => 'denovo',
			cmd => ''
		))->calc_features;
	}

	#stops the thread listener and waits for remaining background jobs to be finished
	$thrListener->stop;
	#store final annotation results
	$gffdb->store_overlaps;
	Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);
	remove_tree($parameter->tmp,{verbose => 0, error  => \my $err});
	print "\nResults stored with label ".$parameter->label."\n";
	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
	$year = $year + 1900;
	$mon += 1;
	print "$mday.$mon.$year-$hour:$min:$sec\n";
}

sub get_phylo_features {
	my ($abbres,$outdir) = @_;

	my ($speciesSSU,$coreFeatures,$stkFeatures,$stk50Features,$stkCoreFeatures,);
	my ($ssuToAbbr,$coreToAbbr,$core50ToAbbr,$rnomeToAbbr);

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
					$speciesFeature->{$abbr} = ($_->get_tag_values('seq'))[0];
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

		next if scalar keys %$speciesFeature < ($#{$abbres} + 1)/2;

		for(@$abbres){
			if(exists $speciesSTKseq->{$_}){
				$core50ToAbbr->{$_}->{$parameter->cfg->rf_rna} = 1;
				$stk50Features->{$_} .= $speciesSTKseq->{$_};
			} else {
				$stk50Features->{$_} .= join('',('?') x $l);
			}
		}

		next if scalar keys %$speciesFeature < ($#{$abbres} + 1);

		for(@$abbres){
			$coreToAbbr->{$_}->{$parameter->cfg->rf_rna} = 1;
			$coreFeatures->{$_} .= $speciesFeature->{$_};
			$stkCoreFeatures->{$_} .= $speciesSTKseq->{$_};
		}
	}

	open TXT , '>'.catfile($outdir,'INFO_SSU.txt') or die $!;
	print TXT "SSU\n";
	print TXT $_."\n" for sort keys %$ssuToAbbr;
	close TXT;
	open TXT , '>'.catfile($outdir,'INFO_RNome.txt') or die $!;
	print TXT "RNome\n";
	print TXT $_."\t".join("\t",sort keys %{$rnomeToAbbr->{$_}})."\n" for sort keys %$rnomeToAbbr;
	close TXT;
	open TXT , '>'.catfile($outdir,'INFO_coreRNome.txt') or die $!;
	print TXT "coreRNome\n";
	print TXT $_."\t".join("\t",sort keys %{$coreToAbbr->{$_}})."\n" for sort keys %$coreToAbbr;
	close TXT;
	open TXT , '>'.catfile($outdir,'INFO_core50RNome.txt') or die $!;
	print TXT "core50RNome\n";
	print TXT $_."\t".join("\t",sort keys %{$core50ToAbbr->{$_}})."\n" for sort keys %$core50ToAbbr;
	close TXT;

	return ($speciesSSU,$coreFeatures,$stkFeatures,$stkCoreFeatures,$stk50Features);
}

sub safety_store {
	my ($e) = @_;
	say $e;
	say ":ERROR: Safety store in progress";
	$thrListener->stop if $thrListener;
	$gffdb->store if $gffdb;
	remove_tree($parameter->tmp,{verbose => 0, error  => \my $err});
	exit 1;
}

__END__

=encoding utf8

=head1 NAME

Gorap - Genomewide ncRNA Annotation Pipeline

=head1 SYNOPSIS

Gorap.pl [OPTIONS]

help: Gorap.pl -h

example: Gorap.pl -i s1.fa,s2.fa -a sa,sb  -g s1.gff,s2.gff -c 4 -k bac -q 1:20,169,1852: -r <taxid> -s 'species name' -sort


=head1 DESCRIPTION

B<Gorap> screens genomic sequences for non-coding RNAs present in the Rfam database using either a generalized strategy 
applying multiple filters on Infernal and Blast results or utilizes specialized software like tRNA-scan, Bcheck, Barrnap, 
CRT, RNAmmer. Gorap provides ncRNA based reconstruction of phylogenetic trees and is able to perform de novo predictions 
including TPM calculations from RNA-Seq experiments. RNA family specific screening options, threshold and constrains can 
be easily amended and completed by custom queries.

=head1 OPTIONS

=over 4

=item B<-h>, B<--help>

this help message

=item B<-v>, B<--version>

this Gorap version

=item B<-update>, B<--update>=[all,rfam,ncbi,silva]

updates internal databases (Rfam, NCBI, Silva)

=item B<-file>, B<--file>=FILE

run Gorap with a parameter file - see provided example for detailed information.
note: command line parameters priorize parameter file settings

=item B<-l>, B<--label>=STRING

(! recommended) - a label for this run

=item B<-i>, B<--fastas>=FILE,...

(! required) - comma separated paths to input FASTA files

=item B<-a>, B<--abbreviations>=STRING,...

(default: filenames, in equal order and number to option -i) - unique list of comma separated abbreviations

=item B<-q>, B<--queries>=1:5,8,...

(default: 1:) - comma separated list of single Rfam ids e.g. RF00001 (or short: 1) and/or colon defined ranges

=item B<-k>, B<--kingdom>=[bac,arc,euk,fungi,virus]

(default: bac,arc,euk,fungi,virus) - list of comma separated kingdoms to screen for kingdom specific ncRNAs

=item B<-r>, B<--rank>=[INT/STRING]

NCBI taxonomy matching id or scientific name of a rank like class/order/genus/... matching all input sequences. please quote inputs with white spaces.

=item B<-s>, B<--species>=[INT/STRING]

NCBI taxonomy matching id or scientific name of a species matching all input sequences. please quote inputs with white spaces.

=item B<-o>, B<--output>=PATH

(default: gorap_out) - output directory

=item B<-c>, B<--cpu>=INT

(default: 1) - number of threads to use

=item B<-t>, B<--tmp>=PATH

(default: output directory) - set the temporary directory

=item B<-sort>, B<--sort>

taxonomically sort resulting alignments

=item B<-nobl>, B<--noblast>

disable Blast search

=item B<-nofi>, B<--nofilter>

disables Gorap specific filters for length, identitiy, secondary structure, overlapping ncRNA families

=item B<-nodel>, B<--nooverlapdeletion>

keep all overlapping ncRNA annotations

=item B<-rfamscan>, B<--rfamscan>

behave like multithreaded rfam_scan on latest Rfam, i.e. use Infernal only, use Rfam thresholds and disable all Gorap filters

=item B<-g>, B<--gffs>=FILE,...

(in equal order and number to option -i) - comma separated paths to known annotations in GFF3 format to mark overlapping ncRNA predictions

=back

=item B<-skip>, B<--skipanno>

skip ncRNA prediction and fastforward to downstream analysis

=item B<-refresh>, B<--refresh>

(! requires option -l) - update label related HTML page and annotation files (GFF and FASTA) according to curated Stockholm alignments

=back

=head1 OPTIONS FOR PHYLOGENY RECONSTRUCTION

=over 4

=item B<-og>, B<--outgroups>=FILE

comma separated paths to outgroup FASTA files, to trigger RNome and SSU rRNA based phylogeny reconstruction

=item B<-oga>, B<--ogabbreviations>=FILE

(default: filenames, in equal order and number to option -og) - unique list of comma separated abbreviations

=head1 EXPERIMENTAL OPTIONS FOR DE NOVO NCRNA PREDICTION

=over 4

=item B<-b>, B<--bams>=FILE,...

(in equal order and number to -i) - comma separated paths to sorted mapping results from RNA-Seq experiments in BAM format to trigger TPM calculation and de novo ncRNA prediction

=item B<-strand>, B<--strandspecific>=[-1,1]

specify if mapping data resulted from strand specific library preparation: (1) forward stranded or (-1) reverse stranded

=item B<-notpm>, B<--notpm>

skip TPM calculation and fastforward to de novo ncRNA prediction

=item B<-minl>, B<--minlength>=INT

(default: 50) - minimum length to report a de novo predicted ncRNA

=item B<-minh>, B<--minheight>=INT

(default: 1000) - minimum observed nucleotide coverage to report a de novo predicted ncRNA

=back

=head1 AUTHOR

Konstantin Riege, E<lt>konstantin.riege@leibniz-fli.deE<gt>

=head1 COPYRIGHT AND LICENSE

MIT License

Copyright (c) 2017 Konstantin Riege

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

=head1 SEE ALSO

L<www.rna.uni-jena.de>

=cut

