package Bio::Gorap::ToolI;

use Moose::Role; requires qw(calc_features);

has 'threads' => (
	is => 'ro',
	isa => 'Int',	
	default => 1
);

has 'parameter' => (
	is => 'ro',
	isa => 'Bio::Gorap::Parameter',
	required => 1
);

has 'tool_parser' => (
	is => 'ro',
	isa => 'CodeRef',
	required => 1
);

has 'gffdb' => (
	is => 'ro',
	isa => 'Bio::Gorap::DB::GFF',
	required => 1
);

has 'fastadb' => (
	is => 'ro',
	isa => 'Bio::Gorap::DB::Fasta',
	required => 1
);

has 'stkdb' => (
	is => 'ro',
	isa => 'Bio::Gorap::DB::STK',
	required => 1
);

has 'tool' => (
	is => 'ro',
	isa => 'Str',
	required => 1
);

#remove entries from database, which will be recomputed now
sub BUILD {
	my ($self) = @_;
	my $abbres;
	for (0..$#{$self->parameter->genomes}){		
		my $genome = ${$self->parameter->genomes}[$_];
		my $abbr = ${$self->parameter->abbreviations}[$_];
		$abbres->{$abbr}=1; #removes all old entries, i.e. String don't start with GORAP+toolname
		for ($self->gffdb->db->{$abbr}->features(-primary_tag => $self->parameter->cfg->rf_rna , -attributes => {source => $self->tool})){												
			$self->gffdb->db->{$abbr}->delete($_);			
		}		
	}
	if (exists $self->stkdb->db->{$self->parameter->cfg->rf_rna}){	
		for ($self->stkdb->db->{$self->parameter->cfg->rf_rna}->each_seq){
			my @id = split /\./ , $_->id;
			next if $#id < 2;
			my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
			$self->stkdb->db->{$self->parameter->cfg->rf_rna}->remove_seq($_) if exists $abbres->{$abbr} && exists $self->fastadb->oheaders->{$orig};
		}
	}	
}

1;