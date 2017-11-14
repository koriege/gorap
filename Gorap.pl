#! /usr/bin/env perl

use v5.10;
use strict;
use warnings;
use sigtrap qw(handler safety_store normal-signals);

use File::Spec::Functions;
use Cwd qw(abs_path);

BEGIN {
	if ( ! $ENV{GORAP} ){
		say "Install and export Gorap environment variable - see README";
		exit 1;
	}

	my @path;
	push @path,abs_path($_) for glob("$ENV{GORAP}/bin/*");
	$ENV{PATH} = $ENV{PATH} ? join(":",(@path,$ENV{PATH})) : join(":",@path);

	if ( $ENV{PERL5LIB} ){
		$ENV{PERL5LIB} = join(":",
			glob(catdir($ENV{GORAP},"gorap","*","perl5","x86_64*")),
			glob(catdir($ENV{GORAP},"gorap","*","perl5")),
			$ENV{PERL5LIB}
		);
	} else {
		$ENV{PERL5LIB} = join(":",
			glob(catdir($ENV{GORAP},"gorap","*","perl5","x86_64*")),
			glob(catdir($ENV{GORAP},"gorap","*","perl5"))
		);
	}

	unshift(@INC,
		glob(catdir($ENV{GORAP},"gorap","*","perl5","x86_64*")),
		glob(catdir($ENV{GORAP},"gorap","*","perl5"))
	);

	if (`tRNAscan-SE -h &> /dev/null; echo \$?` != 0){
		say ":ERROR: check failed - try: setup -i trnascan";
		exit 1;
	}
	if (`java -h &> /dev/null; echo \$?` != 0){
		say ":ERROR: check failed - try: setup -i java";
		exit 1;
	}
	if (`raxml -h &> /dev/null; echo \$?` != 0){
		say ":ERROR: check failed - try: setup -i raxml";
		exit 1;
	}
	if (`newicktopdf -h &> /dev/null; echo \$?` != 0){
		say ":ERROR: check failed - try: setup -i newicktopdf";
		exit 1;
	}
	if (`hmmsearch -h -h &> /dev/null; echo \$?` != 0){
		say ":ERROR: check failed - try: setup -i hmmer";
		exit 1;
	}
	if (`rnabob -h &> /dev/null; echo \$?` != 0){
		say ":ERROR: check failed - try: setup -i rnabob";
		exit 1;
	}
	if (`PATH=\$GORAP/bin/infernal1:\$PATH && Bcheck -h &> /dev/null; echo \$?` != 0){
		say ":ERROR: check failed - try: setup -i bcheck";
		exit 1;
	}
	if (`cmsearch -h &> /dev/null; echo \$?` != 0){
		say ":ERROR: check failed - try: setup -i infernal";
		exit 1;
	}
	if (`blastn -h &> /dev/null; echo \$?` != 0){
		say ":ERROR: check failed - try: setup -i bast";
		exit 1;
	}
	if (`barrnap 2>&1 | grep -Fc Torsten` == 0){
		say ":ERROR: check failed - try: setup -i barrnap";
		exit 1;
	}
	if (`samtools 2>&1 | grep -Fc Version` == 0){
		say ":ERROR: check failed - try: setup -i samtools";
		exit 1;
	}
	if (`crt 2>&1 | grep -Fc OPTIONS` == 0){
		say ":ERROR: check failed - try: setup -i crt";
		exit 1;
	}
	if (`mafft -h 2>&1 | grep -Fc MAFFT` == 0){
		say ":ERROR: check failed - try: setup -i mafft";
		exit 1;
	}
}

use File::Basename;
use File::Copy qw(copy);
use File::Path qw(make_path rmtree);
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

	if ($parameter->taxonomy){
		$taxdb = Bio::Gorap::DB::Taxonomy->new(
			parameter => $parameter
		);
		$stkdb = Bio::Gorap::DB::STK->new(
			parameter => $parameter,
			taxonomy => $taxdb
		);
	} else {
		$stkdb = Bio::Gorap::DB::STK->new(
			parameter => $parameter
		);
	}

	my $gffdb = Bio::Gorap::DB::GFF->new(
		parameter => $parameter
	);

	for my $cfg (@{$parameter->queries}){
		$parameter->set_cfg($cfg);
		my $type = $parameter->cfg->rf_rna;
		my $hold = 0;
		for my $f (@{$gffdb->get_all_features($type)}){
			my $seq;
			$seq = $stkdb->db->{$type}->get_seq_by_id($f->seq_id) if exists $stkdb->db->{$type};
			if ($seq){
				$hold=1;
				$gffdb->update_filter($f->seq_id,$type,'!');
			} else {
				$gffdb->update_filter($f->seq_id,$type,'X') if $f->display_name eq '!' && $f->primary_tag!~/SU_rRNA/;
			}
		}
		unlink $stkdb->idToPath->{$type} if ! $hold && exists $stkdb->db->{$type};
	}
	print "Storing changes\n";
	$stkdb->store;
	$gffdb->store_overlaps;

	Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);
	print "\nResults stored with label: ".$parameter->label."\n";

	rmtree($parameter->tmp);

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

unless ($parameter->skip_comp){
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

	Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);
}

if ($parameter->has_outgroups){

	print "Preparing phylogeny reconstruction\n" if $parameter->verbose;
	my @newQ;
	my @oldqueries = @{$parameter->queries};
	my @oldgenomes = @{$parameter->genomes};
	my @oldabbres = @{$parameter->abbreviations};

	$parameter->set_queries();
	for my $cfg (@{$parameter->queries}){
		$parameter->set_cfg($cfg);
		push @newQ , $parameter->cfg->rf if $#{$gffdb->get_all_features($parameter->cfg->rf_rna , '!')} > -1;
	}

	if ($#newQ > -1){
		my $outdir = catdir($parameter->output,'phylogeny-'.$parameter->label);
		make_path($outdir);
		unlink $_ for glob catfile($outdir,'RAxML_*');

		my @newgenomes = @{$parameter->outgroups};
		my @newabbres = @{$parameter->ogabbreviations};
		$parameter->set_genomes(\@newgenomes,\@newabbres); #new referenfce, otherwise og will be unset during set_genomes

		$parameter->set_queries(\@newQ);

		$gffdb->add_db($_) for @{$parameter->ogabbreviations};

		unless ($parameter->skip_comp){

			print "\nAnnotation of outgroups for phylogeny reconstruction\n" if $parameter->verbose;

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
		}

		my @allabbres = (@oldabbres,@{$parameter->ogabbreviations});
		my ($speciesSSU,$coreFeatures,$stkFeatures,$stkCoreFeatures,$stk50Features) = &get_phylo_features(\@allabbres,$outdir);

		$parameter->set_genomes(\@oldgenomes,\@oldabbres);
		$parameter->set_queries(\@oldqueries);

		Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);

		if ( (any { exists $speciesSSU->{$_} } @{$parameter->ogabbreviations}) && scalar keys %$speciesSSU > 3){
			open FA , '>'.catfile($outdir,'SSU.fasta') or die $!;
			for my $k (keys %$speciesSSU ){
				print FA '>'.$k."\n";
				$speciesSSU->{$k}=~s/(\W|_)/-/g;
				print FA $_."\n" for unpack("(a80)*",$speciesSSU->{$k});
			}
			close FA;

			system('mafft --localpair --maxiterate 1000 --thread '.$parameter->threads.' '.catfile($outdir,'SSU.fasta').' > '.catfile($outdir,'SSU.mafft'));
			system('raxml -T '.$parameter->threads.' -f a -# 1000 -x 1234 -p 1234 -s '.catfile($outdir,'SSU.mafft').' -w '.$outdir.' -n SSU.mafft.tree -m GTRGAMMA -o '.join(',',grep { exists $speciesSSU->{$_} } @{$parameter->ogabbreviations}));

			my $obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($outdir,'RAxML_bipartitions.SSU.mafft.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($outdir,'SSU.mafft.eps'));
			copy catfile($outdir,'RAxML_bipartitions.SSU.mafft.tree'), catfile($outdir,'SSU.mafft.tree');
			system('newicktopdf -pc 1 -boot -notitle '.catfile($outdir,'RAxML_bipartitions.SSU.mafft.tree'));
			copy catfile($outdir,'RAxML_bipartitions.pdf'), catfile($outdir,'SSU.mafft.pdf');

			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);
		}

		if ( (any { exists $coreFeatures->{$_} } @{$parameter->ogabbreviations}) && scalar keys %$coreFeatures > 3){
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
			system('raxml -T '.$parameter->threads.' -f a -# 1000 -x 1234 -p 1234 -s '.catfile($outdir,'coreRNome.mafft').' -w '.$outdir.' -n coreRNome.mafft.tree -m GTRGAMMA -o '.join(',',grep { exists $coreFeatures->{$_} } @{$parameter->ogabbreviations}));
			system('raxml -T '.$parameter->threads.' -f a -# 1000 -x 1234 -p 1234 -s '.catfile($outdir,'coreRNome.stkfa').' -w '.$outdir.' -n coreRNome.stk.tree -m GTRGAMMA -o '.join(',',grep { exists $coreFeatures->{$_} } @{$parameter->ogabbreviations}));

			my $obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($outdir,'RAxML_bipartitions.coreRNome.mafft.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($outdir,'coreRNome.mafft.eps'));
			copy catfile($outdir,'RAxML_bipartitions.coreRNome.mafft.tree'), catfile($outdir,'coreRNome.mafft.tree');
			system('newicktopdf -pc 1 -boot -notitle '.catfile($outdir,'RAxML_bipartitions.coreRNome.mafft.tree'));
			copy catfile($outdir,'RAxML_bipartitions.pdf'), catfile($outdir,'coreRNome.mafft.pdf');

			$obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($outdir,'RAxML_bipartitions.coreRNome.stk.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($outdir,'coreRNome.stk.eps'));
			copy catfile($outdir,'RAxML_bipartitions.coreRNome.stk.tree'), catfile($outdir,'coreRNome.stk.tree');
			system('newicktopdf -pc 1 -boot -notitle '.catfile($outdir,'RAxML_bipartitions.coreRNome.stk.tree'));
			copy catfile($outdir,'RAxML_bipartitions.pdf'), catfile($outdir,'coreRNome.stk.pdf');

			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);
		}

		if ( (any { exists $stkFeatures->{$_} } @{$parameter->ogabbreviations}) && scalar keys %$stkFeatures > 3){
			open FA , '>'.catfile($outdir,'RNome.stkfa') or die $!;
			for my $k (keys %$stkFeatures ){
				print FA '>'.$k."\n";
				$stkFeatures->{$k}=~s/(\W|_)/-/g;
				print FA $_."\n" for unpack("(a80)*",$stkFeatures->{$k});
			}
			close FA;

			system('raxml -T '.$parameter->threads.' -f a -# 1000 -x 1234 -p 1234 -s '.catfile($outdir,'RNome.stkfa').' -w '.$outdir.' -n RNome.stk.tree -m GTRGAMMA -o '.join(',',grep { exists $stkFeatures->{$_} } @{$parameter->ogabbreviations}));

			my $obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($outdir,'RAxML_bipartitions.RNome.stk.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($outdir,'RNome.stk.eps'));
			copy catfile($outdir,'RAxML_bipartitions.RNome.stk.tree'), catfile($outdir,'RNome.stk.tree');
			system('newicktopdf -pc 1 -boot -notitle '.catfile($outdir,'RAxML_bipartitions.RNome.stk.tree'));
			copy catfile($outdir,'RAxML_bipartitions.pdf'), catfile($outdir,'RNome.stk.pdf');

			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);
		}

		if ( (any { exists $stk50Features->{$_} } @{$parameter->ogabbreviations}) && scalar keys %$stk50Features > 3){
			open FA , '>'.catfile($outdir,'core50RNome.stkfa') or die $!;
			for my $k (keys %$stk50Features ){
				print FA '>'.$k."\n";
				$stk50Features->{$k}=~s/(\W|_)/-/g;
				print FA $_."\n" for unpack("(a80)*",$stk50Features->{$k});
			}
			close FA;

			system('raxml -T '.$parameter->threads.' -f a -# 1000 -x 1234 -p 1234 -s '.catfile($outdir,'core50RNome.stkfa').' -w '.$outdir.' -n core50RNome.stk.tree -m GTRGAMMA -o '.join(',',grep { exists $stk50Features->{$_} } @{$parameter->ogabbreviations}));

			my $obj = Bio::Tree::Draw::Cladogram->new(-tree => (Bio::TreeIO->new(-format => 'newick', '-file' => catfile($outdir,'RAxML_bipartitions.core50RNome.stk.tree')))->next_tree , -bootstrap => 1 , -size => 4, -tip => 4 );
			$obj->print(-file => catfile($outdir,'core50RNome.stk.eps'));
			copy catfile($outdir,'RAxML_bipartitions.core50RNome.stk.tree'), catfile($outdir,'core50RNome.stk.tree');
			system('newicktopdf -pc 1 -boot -notitle '.catfile($outdir,'RAxML_bipartitions.core50RNome.stk.tree'));
			copy catfile($outdir,'RAxML_bipartitions.pdf'), catfile($outdir,'core50RNome.stk.pdf');

			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);
		}
	}
}

#remove temp files
rmtree($parameter->tmp);

$stkdb->store if $parameter->skip_comp && $parameter->sort;
print "\nResults stored with label ".$parameter->label."\n";

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
		print $c.' of ',$#{$parameter->queries}+1,' - '.$parameter->cfg->rf if $parameter->verbose;
		#start screening if rfam query belongs to a kingdom of interest
		my $next=1;
		for (keys %{$parameter->kingdoms}){
			$next=0 if exists $parameter->cfg->kingdoms->{$_};
		}
		if($next && ! $parameter->nofilter){
			print " - skipped due to kingdom restriction" if $parameter->verbose;
			next;
		}
		say "" if $parameter->verbose;

		my ($threshold,$nonTaxThreshold);
		for my $tool (@{$parameter->cfg->tools}){
			my $cmd = $parameter->cfg->cmd->{$tool};
			$tool=~s/[\W\d_]//g;
			$tool = lc $tool;
			$tool = ucfirst $tool;
			$tool = 'Infernal' if $parameter->nofilter;
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
			$thrListener->push_obj($obj);

			if ($tool eq 'Infernal' || $tool eq 'Blast'){
				($threshold,$nonTaxThreshold) = $stkdb->calculate_threshold($parameter->threads - $thrListener->get_workload);
				last if $threshold && $threshold == 999999;
			}

			$obj->calc_features;
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
				&safety_store($_);
			};
			$gffdb->update_score_by_file($parameter->cfg->rf_rna,$scorefile);

			#start of time consuming single threaded background job with a subroutine reference,
			#which returns an array of String, parsable by pipe_parser
			if ($threshold){
				my $features = $gffdb->get_features($parameter->cfg->rf_rna,$parameter->abbreviations);
				next if $#{$features} == -1;
				$thrListener->calc_background(sub {$stkdb->filter_stk($parameter->cfg->rf_rna,$stk,$features,$threshold,$nonTaxThreshold,$taxdb,$gffdb)});
			} else {
				$stkdb->remove_gap_columns_and_write($stk,catfile($parameter->output,'alignments',$parameter->cfg->rf_rna.'.stk'),$taxdb);
			}
		}
		#store annotations already in the database in case of errors
		while($#{$thrListener->finished}>-1){
			$gffdb->store(shift @{$thrListener->finished});
			Bio::Gorap::Evaluation::HTML->create($parameter,$gffdb,$stkdb,$gffdb->rnas,$parameter->label);
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
		tool => 'de_novo',
		cmd => ''
	))->calc_features;
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
	# $gffdb->store_overlaps if $gffdb;
	# rmtree($parameter->tmp);
	exit 1;
}

__END__

=encoding utf8

=head1 NAME

GORAP - Genomewide ncRNA Annotation Pipeline

=head1 SYNOPSIS

Gorap.pl [OPTION]...

example: Gorap.pl -a sa,sb -i s1.fa,s2.fa -g s1.gff,s2.gff -c 4 -k bac -q 1:20,169,1852: -r 123 -s 'species name' -sort

=head1 DESCRIPTION

B<GORAP> will read given input sequences and screen them for all non-coding RNAs present in the Rfam database with
generalized or specialized software by use of developed filtering strategies. Furthermore the pipeline is able to
reconstruct phylogenetic trees and perform de novo predictions as well as TPM/FPKM calculations from RNA-Seq experiments.
Screening options are defined in RNA family configuration files at $GORAP/config/*.cfg, setting up software and
parameters, queries, covariance models, thresholds and sequence constrains. These files can be easily amended and completed by own queries.

=head1 OPTIONS

=over 4

=item B<-h>, B<--help>

this (help) message

=item B<-example>, B<--example>

apply GORAP on $GORAP/example/ecoli.fa

=item B<-update>, B<--update>=[all,rfam,ncbi,silva]

updates internal used databases (Rfam, NCBI, Silva)	!!! check your edited configuration files afterwards

=item B<-file>, B<--file>=FILE

run GORAP with a parameter file - see provided example for detailed information.
note: command line parameters priorize parameter file settings

=item B<-l>, B<--label>=STRING

a label for this RUN - useful to find it in the HTML output and for an output refesh
note: see -refresh

=item B<-i>, B<--fastas>=FILE,...

(required)
paths of comma separated FASTA files

=item B<-a>, B<--abbreviations>=STRING,...

(default: build from FASTA file names)
list of comma separated abbreviations as unique identifiers.
note: in equal order and list size to -i (otherwise use parameter file)

=item B<-q>, B<--queries>=1:5,8,...

list of comma separated, single Rfam ids (e.g. RF00001) or numbers (e.g. 1) and ranges defined by ':' (default: all queries)

=item B<-k>, B<--kingdom>=[bac,arc,euk,fungi,virus]

list of comma separated kingdoms to screen for kingdom specific ncRNAs (default: all)

=item B<-r>, B<--rank>=[INT/STRING]

NCBI taxonomy matching id or scientific name of rank/genus/... for given sequences.	please use quotas if using names e.g. 'Bacillus subtilis group'

=item B<-s>, B<--species>=[INT/STRING]

NCBI taxonomy matching scientific name or taxonomy id of given species. please use quotas if using names e.g. 'Bacillus subtilis'

=item B<-o>, B<--output>=PATH

(default: <working directory>/gorap_out)
output directory

=item B<-c>, B<--cpu>=INT

(default: 1)
number of threads to use

=item B<-t>, B<--tmp>=PATH

(default: $TMPDIR or /tmp or $GORAP/tmp)
set the temporary directory - will be removed after successful GORAP run

=item B<-sort>, B<--sort>

enable resulting alignments to be sorted in taxonomic order

=item B<-notax>, B<--notaxonomy>

disables taxonomic sorting and filters based on given rank/species information useful to skip time consuming initializations for testing purposes

=item B<-nofi>, B<--nofilter>

disables GORAP specific sequence and structure filter

=back

=head1 ADDITIONAL OPTIONS

=over 4

=item B<-skip>, B<--skipanno>

disables annotation process - useful for e.g. additional phylogeny reconstruction

=item B<-refresh>, B<--refresh>

disables any calculations - updates HTML page and annotation files (GFF and FASTA) according to manual changes in Stockholm alignment files
note: see -l

=item B<-og>, B<--outgroups>=FILE

path to additional FASTA files as trigger for starting phylogeny reconstructions based on RNome and SSU rRNA annotations in given sequences

=item B<-oga>, B<--ogabbreviations>=FILE

(default: build from FASTA file name of -og)
list of comma separated abbreviations as unique identifiers.
note: in equal order and list size to -og (otherwise use parameter file)

=item B<-g>, B<--gffs>=FILE,...

comma separated paths of known annotations in GFF3 format with necessary ID tag for overlap filter and to calculate TMP/FPKM values
note1: in equal order and list size to -i (otherwise use parameter file).
note2: separate multiple BAMs related to one input FASTA by colons
e.g. s1g1.gff:s1g2.gff,s2.gff

=item B<-noc>, B<--nooverlapcheck>

disables deletion of predictions even if they overlap with a given GFF3 file

=item B<-b>, B<--bams>=FILE,...

comma separated paths of mapping results in BAM format, triggers TPM/FPKM calculation and de novo prediction.
note1: in equal order and list size to -i (otherwise use parameter file).
note2: separate multiple BAMs related to one input FASTA by colons
e.g. s1g1.bam:s1g2.bam,s2.bam

=item B<-strand>, B<--strandspecific>

mapping data (BAM files) resulted from strand specific library preparation

=item B<-notpm>, B<--notpm>

disables TPM/FPKM calculation

=item B<-minl>, B<--minlength>=INT

(default: 50)
minimum length for de novo gene prediction

=item B<-minh>, B<--minheigth>=INT

(default: 1000)
minimum nucleotide coverage for de novo gene prediction

=back

=head1 AUTHOR

Konstantin Riege, E<lt>konstantin.riege@uni-jena.deE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 by Konstantin Riege

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10 or,
at your option, any later version of Perl 5 you may have available.

=head1 SEE ALSO

L<www.rna.uni-jena.de>

=cut

