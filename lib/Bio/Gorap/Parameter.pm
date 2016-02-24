package Bio::Gorap::Parameter;

use Moose; 
use Getopt::Long;
use Pod::Usage;
use Bio::Gorap 2.2;
use Bio::Gorap::CFG;
use Bio::Gorap::Update;
use File::Spec::Functions;
use Switch;
use File::Basename;
use File::Path qw(make_path);
use Try::Tiny;
use Sys::MemInfo qw(totalmem);
use List::Util qw(min max);
use Config::IniFiles;

has 'mem' => (
	is => 'ro',
	isa => 'Int',
	default => sub { min(40000,sprintf("%.0f",Sys::MemInfo::totalmem()/1024/1024 * 0.95)) }
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

has 'commandline' => (
	is => 'ro',
    isa => 'Bool',
    required => 1
);

has 'cfg' => (
	is => 'rw',
	isa => 'Bio::Gorap::CFG'
);

has 'genomes' => (
	is => 'rw',
    isa => 'ArrayRef',
    default => sub { [catfile($ENV{GORAP},'example','ecoli.fa')] }
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
    default => sub { 
		if ($ENV{TMPDIR}){
			return $ENV{TMPDIR};
		} elsif (-e catdir(rootdir,'tmp')){
			return catdir(rootdir,'tmp');
		} else {
			make_path(catdir($ENV{GORAP},'tmp'));
			return catdir($ENV{GORAP},'tmp');
		}
    }
);

has 'skip_comp' => (
	is => 'rw',
    isa => 'Bool',		
    default => 0
);

has 'taxonomy' => (
	is => 'rw',
    isa => 'Bool',		
    default => 1
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

has 'verbose' => (
	is => 'rw',
	isa => 'Bool',
	default => 1
);

has 'sort' => (
	is => 'rw',
	isa => 'Bool',
	default => 0
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

has 'force' => (
	is => 'rw',
	isa => 'Bool',
	default => 0
);

has 'notpm' => (
	is => 'rw',
	isa => 'Bool',
	default => 0
);

has 'noblast' => (
	is => 'rw',
	isa => 'Bool',
	default => 0
);

has 'querystring' => (
	is => 'rw',
	isa => 'Str',
	default => ''
);

has 'check_overlaps' => (
	is => 'rw',
	isa => 'Bool',
	default => 1
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

sub BUILD {
	my ($self) = @_;
	my $file='x';

	(Getopt::Long::Parser->new)->getoptions (
		'i|fastas=s' => \my $genomes, 
		'o|output=s' => \my $output, 
		'c|cpu=i' => \my $threads,
		'k|kingdom=s' => \my $kingdoms, 
		'q|queries=s' => \my $queries,
		'a|abbreviations=s' => \my $abbreviations,
		'r|rank=s' => \my $rank,
		's|species=s' => \my $species,		
		'og|outgroups=s' => \my $outgroups,
		'oga|outgroupabbreviations=s' => \my $ogabbreviations,
		'b|bams=s' => \my $bams,	
		'g|gffs=s' => \my $gffs,
		'update|update=s' => \my $update,
 		'file|file:s' => \$file,
		'h|help' => \my $help,
		'force|force' => \my $force,
		't|tmp:s' => \my $tmp,
		'notax|notaxonomy' => \my $notax,
		'notpm|notpm' => \my $notpm,
		'nobl|noblast' => \my $noblast,
		'sort|sort' => \my $sort,
		'example|example' => \my $example,
		'no|nooverlap' => \my $nooverlaps,
		'minl|minlength=i' => \my $denovolength,
		'minh|minheigth=i' => \my $denovoheigth
	) or pod2usage(-exitval => 1, -verbose => 3) if $self->commandline;

	
	&read_parameter($self,$file) if $file && $file ne 'x';

	unless ($force || ! $self->commandline){
		if($update){
			#downloading and parsing latest databases
			switch (lc $update) {
				case 'ncbi' {
					print "Updating NCBI Taxonomy\n";
					Bio::Gorap::Update->dl_ncbi();
					exit;
				}
				case 'silva' {
					print "Updating Silva tree\n";
					Bio::Gorap::Update->dl_silva($self);
					exit;
				}
				case 'rfam' {
					print "Updating Rfam database\n";
					Bio::Gorap::Update->dl_rfam($self);
					exit;
				}
				case 'cfg' {
					print "Updating configuration files\n";
					Bio::Gorap::Update->create_cfgs($self);
					exit;
				}
				else {
					print "Updating all databases\n";																	
					Bio::Gorap::Update->dl_ncbi();
					my $taxdb = Bio::Gorap::DB::Taxonomy->new(
						parameter => $self
					);
					Bio::Gorap::Update->dl_silva($self,$taxdb);
					Bio::Gorap::Update->dl_rfam($self,$taxdb);	
					Bio::Gorap::Update->create_cfgs($self,$taxdb);											
					exit;
				}
			}
		}					
		pod2usage(-exitval => 1, -verbose => 3) if ($file eq 'x' && $help) || ($file eq 'x' && ! $genomes && ! $example);
	}

	$self->force(1) if $force;

	#store arguments into data structure
	$self->threads($threads) if $threads;
	if ($genomes){
		my @g;
		push @g , glob $_ for split(/\s*,\s*/,$genomes);
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
			pod2usage(-exitval => 1, -verbose => 3) unless $s =~ /^(bac|arc|euk|fungi|virus)$/;
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
			push @{$bams[$c]} , split(/\s*:\s*/,$_);
		}		
		$self->bams(\@bams);
	}		
	if ($gffs){
		my @gffs;
		my $c=-1;
		for (split(/\s*,\s*/,$gffs)){
			$c++;
			next unless $_;
			push @{$gffs[$c]} , split(/\s*:\s*/,$_);
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
	$self->sort(1) if ! $notax && $sort && ($self->has_rank || $self->has_species);
	$self->taxonomy(0) if $notax || ! ($self->has_rank || $self->has_species);
	$self->notpm(1) if $notpm;
	$self->check_overlaps(0) if $nooverlaps;
	$self->denovoheigth($denovoheigth) if $denovoheigth;
	$self->denovolength($denovolength) if $denovolength;
	$self->noblast(1) if $noblast;
		
	make_path(catdir($self->tmp,$self->pid));
	$self->tmp(catdir($self->tmp,$self->pid));
}

sub _make_paths {
	my ($self) = @_;
	
	make_path(catdir($self->output,'alignments'));
	make_path(catdir($self->output,'annotations'));
	make_path(catdir($self->output,'meta'));
	make_path(catdir($self->output,'html'));
	make_path(catdir($self->output,'phylogeny'));
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
			last;
		}
		if ($_=~/R?F?0*(\d+)\s*:\s*R?F?0*(\d+)/){
			for ($1..$2){
				my ($q) = glob(catfile($ENV{GORAP},'config','RF'.((0) x (5-length($_))).$_.'*.cfg'));
				push @queries , $q if $q;
			}
			
		}elsif($_=~/R?F?0*(\d+)\s*:\s*/) {
			my $nr1 = $1;
			my @q = glob(catfile($ENV{GORAP},'config','*.cfg'));
			basename($q[$#q])=~/R?F?0*(\d+)/;
			my $nr2=$1;			
			for ($nr1..$nr2){
				my ($q) = glob(catfile($ENV{GORAP},'config','RF'.((0) x (5-length($_))).$_.'*.cfg'));
				push @queries , $q if $q;
			}										
		} elsif ($_=~/R?F?0*(\d+)/){
			my ($q) = glob(catfile($ENV{GORAP},'config','RF'.((0) x (5-length($1))).$1.'*.cfg'));
			push @queries , $q if $q;
		} 
	}

	$self->queries(\@queries);
}

sub _set_queries {
	my ($self) = @_;

	my @queries;
	push @queries , glob(catfile($ENV{GORAP},'config','*.cfg'));

	return \@queries;
}

sub read_parameter {
	my ($self,$file) = @_;

	my $cfg = Config::IniFiles->new( -file => $file , -nomultiline => 1, -handle_trailing_comment => 1);
	my @genomes;
	my @abbreviations;
	my @ogenomes;
	my @ogabbreviations;	
	my $assignment;

	my $c=1;
	while( my @g = glob $cfg->val( 'input', 'genome'.$c )){		
		if ($#g > 0){
			for (@g){				
				push @genomes , $_;
				push @{$assignment->{'genome'.$c}} , $#genomes;
				my $abbr = basename($_);
				my @abbr = split /\./ , $abbr;
				pop @abbr if $#abbr > 0;
				$abbr = join '' , @abbr;
				$abbr=~s/\W//g;
				push @abbreviations , $abbr;
			}
		} else {
			push @genomes , $g[0];
			push @{$assignment->{'genome'.$c}} , $#genomes;
			my $abbr = $cfg->GetParameterTrailingComment('input', 'genome'.$c);
			unless ($abbr){
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
		if ($#g > 0){
			for (@g){
				push @ogenomes , $_;
				my $abbr = basename($_);
				my @abbr = split /\./ , $abbr;
				pop @abbr if $#abbr > 0;
				$abbr = join '' , @abbr;
				$abbr=~s/\W//g;
				push @ogabbreviations , $abbr;
			}
		} else {
			my $abbr = $cfg->GetParameterTrailingComment('input', 'genome'.$c);
			push @ogenomes , $g[0];
			unless ($abbr){
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
	my %h = map { $_ => 1 } split /\n/ , $cfg->val('query','kingdom');
	$self->kingdoms(\%h) if scalar keys %h > 0;
	&set_queries($self,[split /\n/ , $cfg->val('query','rfam')]);
	$self->querystring(join(",",split(/\n/ , $cfg->val('query','rfam'))));

	my $v = $cfg->val('system','threads');
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
		for (@{$assignment->{$cfg->GetParameterTrailingComment('addons','bam'.$c)}}){			
			push @{$bams[$_]} , @g;
		}
		$c++;
	}
	
	$c=1;	
	while( my @g = glob $cfg->val( 'addons', 'gff'.$c )){		
		for (@{$assignment->{$cfg->GetParameterTrailingComment('addons','gff'.$c)}}){
			push @{$gffs[$_]} , @g;	
		}		
		$c++;
	}
	$self->bams(\@bams) if $#bams > -1;
	$self->gffs(\@gffs) if $#gffs > -1;	
}

1;