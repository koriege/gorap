package Bio::Gorap::CFG;

use Moose;
use Switch;
use File::Spec::Functions;
use Bio::Gorap::Functions::CM;
use List::Util qw(min max);

has 'cfg' => (
	is => 'ro',
    isa => 'Str',
    required => 1
);

has ['rf' , 'rna' , 'rf_rna' , 'query_dir' , 'fasta' , 'stk' , 'cm'] => (
	is => 'rw',
    isa => 'Str'
);

has 'cmd' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub { [] }
);

has 'tools' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub { [] }
);

has ['evalue' , 'bitscore', 'bitscore_cm'] => (
	is => 'rw',
	isa => 'Num'
);

has 'pseudogenes' => (
	is => 'rw',
	isa => 'Int'
);

has 'kingdoms' => (
	is => 'rw',
	isa => 'HashRef',
	default => sub { {} }
);

has 'cssep' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub { [] }
);

has 'csindels' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub { [] }
);

has 'userfilter' => (
	is => 'rw',
    isa => 'Bool',
    default => 0
);

sub BUILD {
	my ($self) = @_;
	#parse rfam query related gorap configuration file of interest into data structure
	open PARAM , '<'.$self->cfg or die $!;
	my $c=-1;	
	my $command;
	my @csparts;
	while(<PARAM>){
		chomp($_);
		$_ =~ s/^\s+|\s+$//g if $c < 9;
		next if $_=~/^\s*$/;
		if ($_=~/^#/){
			$c++;
			if ($c==0){
				substr($_,1)=~/(RF\d+)_(.+)/;
				$self->rf($1);
				$self->rna($2);
				$self->rf_rna($1.'_'.$2);
				$self->query_dir(catdir($ENV{GORAP},'data','rfam',$self->rf_rna));
			} elsif($c==9){	
				for(split(/\|/,substr($_,1))){
					push @csparts , length($_);
					my @matches = $_=~/(\w)/g;					
					push @{$self->cssep} , $#matches > -1 ? $#matches+1 : 0;
					push @{$self->csindels} , -1;									
				}
			}			
			next;
		}		
		switch($c){
			case 1 {				
				push @{$self->cmd} , $_ ;									
			}
			case 2 {				
				if(catfile($self->query_dir,$self->rf_rna.'.fa')=~/$_$/){
					$self->fasta(catfile($self->query_dir,$self->rf_rna.'.fa'));					
				} else {
					$self->fasta($_);					
				}
			}
			case 3 {				
				if(catfile($self->query_dir,$self->rf_rna.'.stk')=~/$_$/){
					$self->stk(catfile($self->query_dir,$self->rf_rna.'.stk'));					
				} else {
					$self->stk($_);					
				}
			}
			case 4 {				
				if(catfile($self->query_dir,$self->rf_rna.'.cm')=~/$_$/){
					$self->cm(catfile($self->query_dir,$self->rf_rna.'.cm'));					
				} else {
					$self->cm($_);					
				}
			}
			case 5 {
				$self->evalue($_)
			}
			case 6 {
				$self->bitscore($_)
			}
			case 7 {
				$self->pseudogenes(lc $_ eq 'all' ? 999999 : $_)
			}
			case 8 {								
				$self->kingdoms->{$_}=1 for split(/,/,$_);								
			}
			case 10 {					
				$self->userfilter(1);
				my $l = substr($_,1);				
				my $start=0;
				for (0..$#csparts){																											
					if (substr($l,$start,$csparts[$_]) =~ /(\d+)/){						
						${$self->csindels}[$_] = $1;
					}
					$start += $csparts[$_]+1;
					last if length($l) <= $start;
				}											
			}
			
			else {
				die 'Check your parameter file!'
			}
		}
	}
	close PARAM;

	my $tool = lc ${$self->cmd}[0];
	$tool =~ s/\W//g;
	push @{$self->tools} , $tool;
	push @{$self->tools} , lc ${$self->cmd}[1] if $#{$self->cmd} > 0 && ${$self->cmd}[1]=~/(B|b)last|(I|i)nfernal/ ;

	my $score = Bio::Gorap::Functions::CM->get_min_score(catfile($self->query_dir,$self->rf_rna.'.cm'));		
	$self->bitscore_cm($score == -1 ? $self->bitscore : $score);
}

1;