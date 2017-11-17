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

has 'bamdb' => (
	is => 'ro',
	isa => 'Bio::Gorap::DB::BAM',
	required => 1
);

has 'tool' => (
	is => 'ro',
	isa => 'Str',
	required => 1
);

has ['already_predicted'] => (
	is => 'rw',
	isa => 'Bool',
	default => 0
);

has 'cmd' => (
	is => 'rw',
	isa => 'Str',
	required => 1
);

#remove entries from database, which will be recomputed now
sub BUILD {
	my ($self) = @_;

	for (0..$#{$self->parameter->genomes}){
		my $abbr = ${$self->parameter->abbreviations}[$_];

		#remove all previously annotated data
		for my $f ($self->gffdb->db->{$abbr}->features(-primary_tag => $self->parameter->cfg->rf_rna , -attributes => {source => $self->tool})){
			$self->gffdb->db->{$abbr}->delete($f);
			next unless exists $self->stkdb->db->{$f->primary_tag};
			$self->stkdb->db->{$f->primary_tag}->remove_seq($_) for $self->stkdb->db->{$f->primary_tag}->get_seq_by_id($f->seq_id);
		}

		unless ($self->tool eq "Infernal" || $self->tool eq "Blast"){ # if tool covers multi families and was executed previously, dont execute again
			my @f = $self->gffdb->db->{$abbr}->features(-attributes => {source => 'GORAP'.$self->tool}); #GORAP needs to be added via ToolParser.pm
			$self->already_predicted(1) if $#f > -1;
		}
	}
}

1;