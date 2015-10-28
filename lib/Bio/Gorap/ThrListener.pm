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
	default => sub { IO::Select->new },	
);

has 'storage_saver' => (
	is => 'ro',
	isa => 'CodeRef',
	required => 1
);

#returns running background threads
sub get_workload {
	my ($self) = @_;

	return scalar keys %{$self->thrList};
}

#terminates thread listener and all forked long time computations
sub stop {
	my ($self) = @_;

	print "Waiting for background jobs to be finished\n" if scalar keys %{$self->thrList} > 0;
	for (keys %{$self->thrList}){
		waitpid($_,0);	
		delete $self->thrList->{$_};					
		&_read($self);
	}	
}

#check for resources and extends an object by inter process communication
sub push_obj {
	my ($self,$obj) = @_;	

	for (keys %{$self->thrList}){
		waitpid($_, &WNOHANG);		
		if (WIFEXITED($?)){			
			delete $self->thrList->{$_};	
			&_read($self);							
		} 
	}

	while($self->get_workload >= $self->threads/2){
		sleep(0.1);
		my ($key) = keys %{$self->thrList};
		waitpid($key,0);	
		delete $self->thrList->{$key};			
		&_read($self);
	}	
}

#prints final results from IPC into a datastructure by code referenced subroutine
sub _read {	
	my ($self) = @_;
	
	while( my @responses = $self->select->can_read(0) ){
		for my $pipe (@responses){			
			while(<$pipe>){	
				&{$self->storage_saver}($self->storage,split(/\s+/,$_));				
			}			
			$self->select->remove($pipe->fileno);
		}
	}
}

sub calc_background {
	#subroutine reference, which returns array of strings to print into pipe
	my ($self,$sub) = @_;

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