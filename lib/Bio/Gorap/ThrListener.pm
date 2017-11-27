package Bio::Gorap::ThrListener;

use Moose;
use IO::Select;
use IO::Pipe;
use POSIX qw(:sys_wait_h);

has 'thrList' => (
	is => 'ro',
    isa => 'HashRef',
    init_arg => undef,
	default => sub { {} }
);

has 'storage' => (
	is => 'ro',
	isa => 'Object',
	required => 1
);

has 'threads' => (
	is => 'rw',
    isa => 'Int',
	default => 1,
);

has 'select' => (
	is => 'ro',
    isa => 'IO::Select',
    init_arg => undef,
	default => sub { IO::Select->new }
);

has 'storage_saver' => (
	is => 'ro',
	isa => 'CodeRef',
	required => 1
);

has 'finished' => (
	is => 'rw',
	isa => 'ArrayRef',
	default => sub {[]}
);

#returns running background threads
sub get_workload {
	my ($self) = @_;

	return scalar keys %{$self->thrList};
}

#terminates thread listener and all forked long time computations
sub stop {
	my ($self) = @_;

	for (keys %{$self->thrList}){
		waitpid($_,0);
		delete $self->thrList->{$_};
		push @{$self->finished} , $self->_read();
	}
}

#check for resources
sub _enqueue {
	my ($self) = @_;

	for (keys %{$self->thrList}){
		waitpid($_, &WNOHANG);
		if (WIFEXITED($?)){
			delete $self->thrList->{$_};
			push @{$self->finished} , $self->_read();
		}
	}

	while($self->get_workload >= 0.6*$self->threads){
		sleep(0.1);
		my ($key) = keys %{$self->thrList};
		waitpid($key,0);
		delete $self->thrList->{$key};
		push @{$self->finished} , $self->_read();
	}
}

#prints final results from IPC into a datastructure by code referenced subroutine
sub _read {
	my ($self) = @_;
	my $type;
	while( my @responses = $self->select->can_read(0) ){
		for my $pipe (@responses){
			while(<$pipe>){
				$type = &{$self->storage_saver}($self->storage,split(/\s+/,$_));
			}
			$self->select->remove($pipe->fileno);
		}
	}

	return $type;
}

sub start {
	#subroutine reference, which returns array of strings to print into pipe
	my ($self,$sub) = @_;

	$self->_enqueue;

	my $pipe = IO::Pipe->new;
	#start single threaded filter calculations
	my $pid = fork();
	unless ($pid) {
		$pipe->writer();
		$pipe->autoflush(1);
		#background calculation starts here
		#returned array of strings of given subroutine is written to IO::Pipe
		for (&$sub){
			chomp $_;
			print $pipe $_."\n";
		}
		exit;
	} else {
		$pipe->reader();
		$self->select->add($pipe);
	}

	$self->thrList->{$pid}=1;
}

1;