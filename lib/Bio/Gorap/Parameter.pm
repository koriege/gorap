package Bio::Gorap::Parameter;

use Moose; 
use Getopt::Long;
use Pod::Usage;
use Bio::Gorap::CFG;
use Bio::Gorap::Update;
use Bio::Gorap::Gorap;
use File::Spec::Functions;
use Switch;
use File::Basename;
use File::Path qw(make_path);
use Try::Tiny;
use Sys::MemInfo qw(totalmem);
use List::Util qw(min max);
use List::MoreUtils qw(uniq);
use Config::IniFiles;

has 'label' => (
	is => 'rw',
	isa => 'Str',
);

has 'mem' => (
	is => 'ro',
	isa => 'Int',
	default => sub { min(40000,sprintf("%.0f",Sys::MemInfo::totalmem()/1024/1024 * 0.9)) }
);

has 'pwd' => (
	is => 'ro',
	isa => 'Str',
	required => 1
);

has 'pid' => (
	is => 'ro',
	isa => 'Int',
	required => 1
);

has 'cfg' => (
	is => 'rw',
	isa => 'Bio::Gorap::CFG'
);

has 'genomes' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub { [catfile($ENV{GORAP},'gorap','example','ecoli.fa')] }
);

has 'threads' => (
	is => 'rw',
	isa => 'Int',		
	default => 1
);

has 'abbreviations' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub { [ 'ecoli' ] }
);

has 'tmp' => (
	is => 'rw',
	isa => 'Str',
	default => catdir($ENV{GORAP},'tmp')
);

has commandline => (
	is => 'rw',
	isa => 'Bool',
	required => 1
);

has ['taxonomy', 'verbose', 'check_overlaps'] => (
	is => 'rw',
	isa => 'Bool',
	default => 1
);

has ['sort', 'skip_comp', 'notpm', 'noblast', 'nofilter', 'nobutkingsnofilter', 'strandspec', 'notax', 'refresh'] => (
	is => 'rw',
	isa => 'Bool',		
	default => 0
);

has 'kingdoms' => (
	is => 'rw',
	isa => 'HashRef',		
	default => sub { {bac => 1 , arc => 1 , euk => 1 , fungi => 1 , virus => 1} }
);

has 'queries' => (
	is => 'rw',
	isa => 'ArrayRef',		
	builder => '_set_queries'
);

has 'output' => (
	is => 'rw',
	isa => 'Str',		
	lazy => 1,
	default => sub { my $self = shift; 
		make_path(catdir($self->pwd,'gorap_out','alignments'));
		make_path(catdir($self->pwd,'gorap_out','annotations'));
		make_path(catdir($self->pwd,'gorap_out','meta'));
		make_path(catdir($self->pwd,'gorap_out','html'));
		return catdir($self->pwd,'gorap_out'); 
	},
	trigger => \&_make_paths
);

has 'rank' => (
	is => 'rw',
	isa => 'Str',
	predicate => 'has_rank'
);

has 'species' => (
	is => 'rw',
	isa => 'Str',
	predicate => 'has_species'		    
);

has 'bams' => (
	is => 'rw',
	isa => 'ArrayRef',
	predicate => 'has_bams',
);

has 'gffs' => (
	is => 'rw',
	isa => 'ArrayRef',
	predicate => 'has_gffs',
);

has 'outgroups' => (
	is => 'rw',
	isa => 'ArrayRef',
	predicate => 'has_outgroups'
);

has 'ogabbreviations' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub { [] }
);

has 'querystring' => (
	is => 'rw',
	isa => 'Str',
	default => ''
);

has 'denovolength' => (
	is => 'rw',
	isa => 'Int',
	default => 50
);

has 'denovoheigth' => (
	is => 'rw',
	isa => 'Int',
	default => 1000
);

has 'thfactor' => (
	is => 'rw',
	isa => 'Num',
	default => 0.9
);

has 'cmtaxbiascutoff' => (
	is => 'rw',
	isa => 'Num',
	default => 0.333
);

has 'file' => (
	is => 'rw',
	isa => 'Str',
	predicate => 'has_file',
);

sub BUILD {
	my ($self) = @_;
	my $file='x';

	(Getopt::Long::Parser->new)->getoptions (
		'l|label=s' => \my $label, 
		'i|fastas=s' => \my $genomes, 
		'o|output=s' => \my $output, 
		'c|cpu=i' => \my $threads,
		'k|kingdom=s' => \my $kingdoms, 
		'q|queries=s' => \my $queries,
		'a|abbreviations=s' => \my $abbreviations,
		'r|rank=s' => \my $rank,
		's|species=s' => \my $species,
		'v|version' => \my $version,
		'og|outgroups=s' => \my $outgroups,
		'oga|outgroupabbreviations=s' => \my $ogabbreviations,
		'b|bams=s' => \my $bams,	
		'g|gffs=s' => \my $gffs,
		'update|update=s' => \my $update,
 		'file|file:s' => \$file,
		'h|help' => \my $help,
		't|tmp:s' => \my $tmp,
		'notax|notaxonomy' => \my $notax,
		'notpm|notpm' => \my $notpm,
		'nobl|noblast' => \my $noblast,
		'sort|sort' => \my $sort,
		'refresh|refresh' => \my $refresh,
		'skip|skipanno' => \my $skipanno,
		'example|example' => \my $example,
		'noc|nooverlapcheck' => \my $nooverlaps,
		'minl|minlength=i' => \my $denovolength,
		'minh|minheigth=i' => \my $denovoheigth,
		'nobutkingsnofi|nobutkingsnofi' => \my $nobutkingsnofilter, #hidden dev option
		'nofi|nofilter' => \my $nofilter, 
		'strand|strandspecific' => \my $strandspec,
		'thfactor|thresholdfactor=f' => \my $thfactor, #hidden dev option
		'biasco|biascutoff=f' => \my $taxbiascutoff #hidden dev option
	) or pod2usage(-exitval => 1, -verbose => 0) if $self->commandline;

	&read_parameter($self,$file) if $file && $file ne 'x';

	if ($self->commandline){
		if($update){
			#downloading and parsing latest databases
			switch (lc $update) {
				case 'ncbi' {
					print "Updating NCBI Taxonomy\n";
					my $taxdb = Bio::Gorap::DB::Taxonomy->new(
						parameter => $self
					);
					Bio::Gorap::Update->dl_ncbi($self,$taxdb);
					exit;
				}
				case 'silva' {
					print "Updating Silva tree\n";
					my $taxdb = Bio::Gorap::DB::Taxonomy->new(
						parameter => $self
					);
					Bio::Gorap::Update->dl_silva($self,$taxdb);
					exit;
				}
				case 'rfam' {
					my $taxdb = Bio::Gorap::DB::Taxonomy->new(
						parameter => $self
					);
					print "Updating Rfam database\n";
					Bio::Gorap::Update->dl_rfam($self,$taxdb);
					print "Updating configuration files\n";
					Bio::Gorap::Update->create_cfgs($self);
					exit;
				}
				case 'cfg' {
					print "Updating configuration files\n";
					my $taxdb = Bio::Gorap::DB::Taxonomy->new(
						parameter => $self
					);
					Bio::Gorap::Update->create_cfgs($self,$taxdb);
					exit;
				}
				else {
					print "Updating all databases\n";																	
					my $taxdb = Bio::Gorap::DB::Taxonomy->new(
						parameter => $self
					);
					Bio::Gorap::Update->dl_ncbi($self,$taxdb);
					Bio::Gorap::Update->dl_silva($self,$taxdb);
					Bio::Gorap::Update->dl_rfam($self,$taxdb);	
					Bio::Gorap::Update->create_cfgs($self,$taxdb);											
					exit;
				}
			}
		}

		pod2usage(-exitval => 0, -verbose => 0, -message => "Version:\n    ".Bio::Gorap::Gorap->VERSION) if $version;
		pod2usage(-exitval => 0, -verbose => -1, -message => "Version:\n    ".Bio::Gorap::Gorap->VERSION) if $file eq 'x' && $help;
		pod2usage(-exitval => 1, -verbose => 0, -message => "Option i requiered. Use option h to get help\n") if $file eq 'x' && ! $genomes && ! $example;
	}

	if ($label){
		$label=~s/\s+/_/g;
		$self->label($label);
	}
	$self->skip_comp(1) if $skipanno;
	$self->notax(1) if $notax;
	$self->refresh(1) if $refresh;

	#store arguments into data structure
	$self->threads($threads) if $threads;
	if ($genomes){
		my @g;
		for (split(/\s*,\s*/,$genomes)){
			for (glob $_){
				pod2usage(-exitval => 1, -verbose => 0, -message => ":ERROR: Option i. File does not exists") unless -e $_;
				push @g , $_;
			}
		}		
		&set_genomes($self, \@g ,[split(/\s*,\s*/,$abbreviations ? $abbreviations : "")]);
	}

	my @ogg;		
	do { push @ogg , glob $_ for split(/\s*,\s*/,$outgroups); $self->outgroups(\@ogg) } if $outgroups;	
	if ($ogabbreviations){			
		$self->ogabbreviations([split(/\s*,\s*/,$ogabbreviations)]);
	} elsif ($outgroups) {
		my @ogabbre;
		for(@ogg){
			my $abbr = basename($_);
			my @abbr = split /\./ , $abbr;
			pop @abbr if $#abbr > 0;
			$abbr = join '' , @abbr;				
			$abbr=~s/\W//g;				
			push @ogabbre , $abbr;
		}
		$self->ogabbreviations(\@ogabbre);
	}		
	if($kingdoms){
		$self->kingdoms({});
		for (split(/\s*,\s*/,$kingdoms)) {
			my $s = lc $_;
			pod2usage(-exitval => 1, -verbose => 0, -message => ":ERROR: Option k. Unknown kingdom definition") unless $s =~ /^(bac|arc|euk|fungi|virus)$/;
			$self->kingdoms->{$s}=1
		}
	}		
			
	if (defined $queries){
		$self->querystring($queries);
		&set_queries($self,[split(/\s*,\s*/,$queries)]);
	}
	
	if ($bams){
		my @bams;
		my $c=-1;
		for (split(/\s*,\s*/,$bams)){
			$c++;
			next unless $_;
			for (split(/\s*:\s*/,$_)) {
				pod2usage(-exitval => 1, -verbose => 0, -message => ":ERROR: Option b. File does not exists") unless -e $_;
				push @{$bams[$c]} , $_;
			}
		}		
		$self->bams(\@bams);
	}		
	if ($gffs){
		my @gffs;
		my $c=-1;
		for (split(/\s*,\s*/,$gffs)){
			$c++;
			next unless $_;
			for (split(/\s*:\s*/,$_)) {
				pod2usage(-exitval => 1, -verbose => 0, -message => ":ERROR: Option g. File does not exists") unless -e $_;
				push @{$gffs[$c]} , $_;
			}
		}		
		$self->gffs(\@gffs);
	}
	if ($output){
	    try {
			make_path(catdir(rootdir, $output)); 
			$self->output($output);
	    } catch {
			$self->output(catdir($self->pwd, $output));
	    };
		
	}
	$self->rank($rank) if $rank;
	$self->species($species) if $species;

	$self->tmp($tmp) if $tmp;
	$self->sort(1) if ! $notax && $sort;
	$self->taxonomy(0) if (! $sort && $refresh) || $notax || (!($self->has_rank || $self->has_species) && ! $sort) || $#{$self->queries} == -1 || $skipanno;
	$self->notpm(1) if $notpm;
	$self->check_overlaps(0) if $nooverlaps;
	$self->denovoheigth($denovoheigth) if $denovoheigth;
	$self->denovolength($denovolength) if $denovolength;
	$self->noblast(1) if $noblast;
	$self->nofilter(1) if $nofilter;
	$self->nobutkingsnofilter(1) if $nobutkingsnofilter;
	$self->thfactor($thfactor) if $thfactor;
	$self->cmtaxbiascutoff($taxbiascutoff) if $taxbiascutoff;
	$self->strandspec(1) if $strandspec;
		
	make_path(catdir($self->tmp,$self->pid));
	$self->tmp(catdir($self->tmp,$self->pid));
}

sub _make_paths {
	my ($self) = @_;
	
	make_path(catdir($self->output,'alignments'));
	make_path(catdir($self->output,'annotations'));
	make_path(catdir($self->output,'meta'));
	make_path(catdir($self->output,'html'));
}

sub set_genomes {
	my ($self, $genomes, $abbreviations) = @_;

	undef @{$self->genomes};
	undef @{$self->abbreviations};
	$self->genomes($genomes);

	if ($#{$abbreviations}>-1){		
		$self->abbreviations($abbreviations);
	} else {
		for(@{$self->genomes}){
			my $abbr = basename($_);
			my @abbr = split /\./ , $abbr;
			pop @abbr if $#abbr > 0;
			$abbr = join '' , @abbr;				
			$abbr=~s/\W//g;				
			push @{$self->abbreviations} , $abbr;
		}
	}
}

#store parsed rfam query related gorap configuration file of interest
sub set_cfg {
	my ($self,$cfg) = @_;	
	
	$self->cfg(Bio::Gorap::CFG->new(
		cfg => $cfg
	));
}

sub set_queries {
	my ($self,$queries) = @_;

	$self->queries(&_set_queries()) , return unless defined $queries;
	#get the rfam queries related gorap configuration files
	my @queries = ();
	for (@{$queries}){
		if ($_ eq '0'){
			@queries = ();
			$self->skip_comp(1);
			# $self->taxonomy(0); #dont do this e.g. if someone want to sort available alignments wo calculation
			last;
		}
		if ($_=~/R?F?0*(\d+)\s*:\s*R?F?0*(\d+)/){
			for ($1..$2){
				my ($q) = glob(catfile($ENV{GORAP},'gorap','config','RF'.((0) x (5-length($_))).$_.'*.cfg'));
				push @queries , $q if $q;
			}
			
		}elsif($_=~/R?F?0*(\d+)\s*:\s*/) {
			my $nr1 = $1;
			my @q = glob(catfile($ENV{GORAP},'gorap','config','*.cfg'));
			basename($q[$#q])=~/R?F?0*(\d+)/;
			my $nr2=$1;			
			for ($nr1..$nr2){
				my ($q) = glob(catfile($ENV{GORAP},'gorap','config','RF'.((0) x (5-length($_))).$_.'*.cfg'));
				push @queries , $q if $q;
			}										
		} elsif ($_=~/R?F?0*(\d+)/){
			my ($q) = glob(catfile($ENV{GORAP},'gorap','config','RF'.((0) x (5-length($1))).$1.'*.cfg'));
			push @queries , $q if $q;
		} 
	}

	@queries = uniq @queries;
	$self->queries(\@queries);
}

sub _set_queries {
	my ($self) = @_;

	my @queries;
	push @queries , glob(catfile($ENV{GORAP},'gorap','config','*.cfg'));

	return \@queries;
}

sub read_parameter {
	my ($self,$file) = @_;
	$self->file($file);

	my $cfg = Config::IniFiles->new( -file => $file , -nomultiline => 1, -handle_trailing_comment => 1);
	my @genomes;
	my @abbreviations;
	my @ogenomes;
	my @ogabbreviations;	
	my $assignment;

	my $c=1;
	while( my @g = glob $cfg->val( 'input', 'genome'.$c )){		
		for (@g){
			pod2usage(-exitval => 1, -verbose => 0, -message => ":ERROR: Option i. File does not exists") unless -e $_;
			push @genomes , $_;
			push @{$assignment->{'genome'.$c}} , $#genomes;
			my $abbr = $cfg->GetParameterTrailingComment('input', 'genome'.$c);			
			unless ($abbr) {
				$abbr = basename($_);
				my @abbr = split /\./ , $abbr;
				pop @abbr if $#abbr > 0;
				$abbr = join '' , @abbr;
				$abbr=~s/\W//g;
			}			
			push @abbreviations , $abbr;
		}		
		$c++;
	}	
	$c=1;
	while( my @g = glob $cfg->val( 'addons', 'genome'.$c )){
		for (@g){
			pod2usage(-exitval => 1, -verbose => 0, -message => ":ERROR: Option og: File does not exists") unless -e $_;
			push @ogenomes , $_;
			my $abbr = $cfg->GetParameterTrailingComment('addons', 'genome'.$c);
			unless ($abbr) {
				$abbr = basename($_);
				my @abbr = split /\./ , $abbr;
				pop @abbr if $#abbr > 0;
				$abbr = join '' , @abbr;
				$abbr=~s/\W//g;
			}
			push @ogabbreviations , $abbr;
		}
		$c++;
	}
	my $v = $cfg->val('query','kingdom');
	if ($v){
		my %h = map { if ($_ =~ /^(bac|arc|euk|fungi|virus)$/) { $_ => 1 } else { pod2usage(-exitval => 1, -verbose => 0, -message => ":ERROR: Option k. Unknown kingdom definition")} } split /\n/ , $cfg->val('query','kingdom');
		$self->kingdoms(\%h) if scalar keys %h > 0;
	}
	$v = $cfg->val('query','rfam');	
	if (defined $v){
		&set_queries($self,[split /\n/ , $cfg->val('query','rfam')]);
		$self->querystring(join(",",split(/\n/ , $cfg->val('query','rfam'))));
	}

	$v = $cfg->val('system','threads');
	$self->threads($v) if $v;
	$v = $cfg->val('output','out');
	$self->output($v) if $v;
	$v = $cfg->val('system','temp');
	$self->tmp($v) if $v;

	$v = $cfg->val('query','family');
	$self->rank($v) if $v;
	$v = $cfg->val('query','species');
	$self->species($v) if $v;

	$self->genomes(\@genomes) if $#genomes > -1;
	$self->abbreviations(\@abbreviations) if $#abbreviations > -1;
	$self->outgroups(\@ogenomes) if $#ogenomes > -1;
	$self->ogabbreviations(\@ogabbreviations) if $#ogabbreviations > -1;
	
	my @bams;
	my @gffs;
	$#bams=$#genomes;
	$#gffs=$#genomes;
	$c=1;	
	while( my @g = glob $cfg->val( 'addons', 'bam'.$c )){
		for (@g){
			pod2usage(-exitval => 1, -verbose => 0, -message => ":ERROR: Option b. File does not exists") unless -e $_;
		}
		for (@{$assignment->{$cfg->GetParameterTrailingComment('addons','bam'.$c)}}){			
			push @{$bams[$_]} , @g;
		}
		$c++;
	}
	$self->bams(\@bams) if $c > 1;
	$self->strandspec(1) if $cfg->val('addons','strandspecific');
	$c=1;	
	while( my @g = glob $cfg->val( 'addons', 'gff'.$c )){
		for (@g){
			pod2usage(-exitval => 1, -verbose => 0, -message => ":ERROR: Option g. File does not exists") unless -e $_;
		}
		for (@{$assignment->{$cfg->GetParameterTrailingComment('addons','gff'.$c)}}){
			push @{$gffs[$_]} , @g;	
		}		
		$c++;
	}	
	$self->gffs(\@gffs) if $c > 1;
	$v = $cfg->val('addons','label');
	if ($v){
		$v=~s/\s+/_/g;
		$self->label($v);
	}
}

1;