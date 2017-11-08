package Bio::Gorap::CFG;

use Moose;
use Switch;
use File::Spec::Functions;
use Bio::Gorap::Functions::CM;
use List::Util qw(min max);
use Config::IniFiles;

has 'cfg' => (
	is => 'ro',
    isa => 'Str',
    trigger => \&_set
);

has ['rf' , 'rna' , 'rf_rna' , 'query_dir' , 'fasta' , 'stk' , 'cm', 'types'] => (
	is => 'rw',
    isa => 'Str',
    default => ''
);

has 'cmd' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub { [] }
);

has 'tools' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub { ['infernal','blast'] }
);

has ['bitscore', 'bitscore_cm'] => (
	is => 'rw',
	isa => 'Num',
	default => 20
);

has 'evalue' => (
	is => 'rw',
	isa => 'Num',
	default => 1e-3
);

has 'pseudogenes' => (
	is => 'rw',
	isa => 'Int',
	default => 999999
);

has 'kingdoms' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {bac => 1 , arc => 1 , euk => 1 , fungi => 1 , virus => 1} }
);

has 'constrains' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub { [] }
);

has 'cs' => (
	is => 'rw',
	isa => 'Str'
);

has 'ss' => (
	is => 'rw',
	isa => 'Str'
);

has 'userfilter' => (
	is => 'rw',
    isa => 'Bool',
    default => 0
);

sub _set {
	my ($self) = @_;	
	
	my $cfg = Config::IniFiles->new( -file => $self->cfg , -nomultiline => 1, -handle_trailing_comment => 1);
	
	$self->rf($cfg->val('family','id')) or die 'Check your parameter file '.$self->cfg;
	$self->rna($cfg->val('family','name')) or die 'Check your parameter file '.$self->cfg;
	$self->rf_rna($self->rf.'_'.$self->rna);
	my ($querydir) = glob(catdir($ENV{GORAP},'gorap','data','rfam',$self->rf.'*'));
	my ($fastafile) = glob(catfile($querydir,$self->rf.'*.fa'));
	my ($cmfile) = glob(catfile($querydir,$self->rf.'*.cm'));
	my ($stkfile) = glob(catfile($querydir,$self->rf.'*.stk'));
	$self->query_dir($querydir);

	my $v = $cfg->val('cmd','tool') or die 'Check your parameter file '.$self->cfg;
	if ($v){		
		$self->tools([split /\n/ , $v]);
		$v = $cfg->val('cmd','parameter');
		if ($v){
			my @parameter = (${$self->tools}[-1] , split(/\n/ , $v));	
			$self->cmd(\@parameter);
		}
	}

	$v = $cfg->val('thresholds','evalue');
	$self->evalue($v) if $v;
	$v = $cfg->val('thresholds','bitscore');	
	$self->bitscore($v) if $v;

	$v = $cfg->val('query','fasta');
	if (-e $v){
		$self->fasta($v);
	} else {
		die 'Check fasta definition in parameter file '.$self->cfg unless -e $fastafile;
		$self->fasta($fastafile);
	}
	$v = $cfg->val('query','stk');
	if (-e $v){
		$self->stk($v);
	} else {
		die 'Check stk definition in parameter file '.$self->cfg unless -e $stkfile;
		$self->stk($stkfile);
	}
	$v = $cfg->val('query','cm');
	if (-e $v){
		$self->cm($v);
	} else {
		die 'Check cm definition in parameter file '.$self->cfg unless -e $cmfile;
		$self->cm($cmfile);
	}

	$self->bitscore_cm(Bio::Gorap::Functions::STK->get_min_score($self->stk));

	my $types = Bio::Gorap::Functions::STK->get_rna_types($self->stk);
	$types=~s/:*CD-box// if $self->rf_rna =~ /_U[0-9](1|2|atac|_|$)/;
	$self->types($types);

	$v = $cfg->val('query','pseudogenes');
	$self->pseudogenes($v) if $v && $v=~/^\d+$/;
	$v = $cfg->val('query','kingdom');
	if ($v){
		my %h = map { $_ => 1} split(/\n/ , $v);
		$self->kingdoms(\%h);
	}

	$v = $cfg->val('constrains','structure');
	$self->ss($v) if $v;
	$v = $cfg->val('constrains','conserved');
	$self->cs($v) if $v;
	$v = $cfg->val('constrains','constrain');
	if ($v){
		$self->userfilter(1);		
		for (split/\n/ , $v){
			while ($_=~/\|(\.*\d+\.*)\|/g){					
				my ($sta,$sto) = ($-[0]+1,$+[0]-1);
				$1=~/(\d+)/;
				push @{$self->constrains} , [$sta+1,$sto,$1,substr($self->cs,$sta,$sto-$sta),$_];
			}
		}
	}
}

1;
