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

has ['rf' , 'rna' , 'rf_rna' , 'query_dir' , 'fasta' , 'stk' , 'cm'] => (
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
	default => 0
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
	$self->query_dir(catdir($ENV{GORAP},'data','rfam',$self->rf_rna));

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
	if ($v){
		if(catfile($self->query_dir,$self->rf_rna.'.fa')=~/$v/){
			$self->fasta(catfile($self->query_dir,$self->rf_rna.'.fa'));
		} else {
			$self->fasta($v);
		}
	}
	$v = $cfg->val('query','stk');
	if ($v){
		if(catfile($self->query_dir,$self->rf_rna.'.stk')=~/$v/){
			$self->stk(catfile($self->query_dir,$self->rf_rna.'.stk'));
		} else {
			$self->stk($v);
		}
	}
	$v = $cfg->val('query','cm');
	if ($v){
		if(catfile($self->query_dir,$self->rf_rna.'.cm')=~/$v/){
			$self->cm(catfile($self->query_dir,$self->rf_rna.'.cm'));
		} else {
			$self->cm($v);
		}
		my $score = Bio::Gorap::Functions::CM->get_min_score($self->cm);		
		$self->bitscore_cm($score) if $score;
	}
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
