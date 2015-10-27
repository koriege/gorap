package Bio::Gorap::Evaluation::Clustering;

use Moose;
use Bio::AlignIO;
use POSIX;
use Bio::SeqIO;
use Tree::Simple;
use List::Util qw(min max);
use File::Spec::Functions;

has 'parameter' => (
	is => 'ro',
	isa => 'Bio::Gorap::Parameter',
	required => 1
);

sub write_iid_60_sets {
	my ($self) = @_;

	my $idToLength;
	my $stk = (Bio::AlignIO->new(-file => $self->parameter->cfg->stk, -format => 'Stockholm', -verbose => -1))->next_aln;

	my @lengths = 0;
	for my $seqObj ($stk->each_seq){
		my $seq = $seqObj->seq;
		$seq =~ s/\W//g;
		$seq = uc $seq;		
		$seq=~y/U/T/;
		$idToLength->{$seqObj->id} = length $seq;
		push @lengths , length $seq;		
	}

	my $medianlength = $lengths[ceil(($#lengths+1)/2)]*0.7;   
	for (keys %{$idToLength}){
		if ($idToLength->{$_} < $medianlength){			
			$stk->remove_seq($stk->get_seq_by_id($_))
		}
	}	
	
	#compute distace matrix
	my @D;	
	my $max = {'max' => 0 , 'maxi' => 1 , 'maxj' => 1};	
	my @names;
	for my $i (1..$stk->num_sequences){
		push @names , ($stk->get_seq_by_pos($i))->id;
		for my $j (1..$stk->num_sequences){
			if ($j<=$i){
				$D[$i-1][$j-1] = 0;
				# print '0  ';
				next;
			}
			my $twoseqaln = $stk->select_noncont($i, $j);
			$D[$i-1][$j-1] = floor($twoseqaln->percentage_identity);
			# print $D[$i-1][$j-1] >= 10 ? $D[$i-1][$j-1].' ' : $D[$i-1][$j-1].'  ';
			if ($D[$i-1][$j-1]> $max->{'max'}){
				$max = {'max' => $D[$i-1][$j-1] , 'maxi' => $i , 'maxj' => $j};				
			}		
		}
		# print "\n";
	}

	my $trees={};
	
	my $node = Tree::Simple->new($max->{'max'});
	my $n1 = Tree::Simple->new($max->{'maxi'});
	my $n2 = Tree::Simple->new($max->{'maxj'});	
	$node->addChildren($n1,$n2);	
	#compute upgma tree with node lables of max sequence iid
	my $c=0;	
	$trees->{0} = $node;	
	while ($max->{'max'}>0){		
		my $maxi = $max->{'maxi'};
		my $maxj = $max->{'maxj'};
		
#		#upate
		for my $i (1..$maxj){			
			$D[$i-1][$maxj-1] = 0;			
		}
		for my $j ($maxj..$stk->num_sequences){						
			$D[$maxj-1][$j-1] = 0;
		}
		for my $i (1..$maxi){
			next if $D[$i-1][$maxi-1] == 0;
			my $twoseqaln1 = $stk->select_noncont($maxi, $i);
			my $twoseqaln2 = $stk->select_noncont($maxj, $i);			
			$D[$i-1][$maxi-1] = floor(max($twoseqaln1->percentage_identity,$twoseqaln2->percentage_identity));
		}
		for my $j ($maxi+1..$stk->num_sequences){
			next if $D[$maxi-1][$j-1] == 0 ;
			my $twoseqaln1 = $stk->select_noncont($maxi, $j);
			my $twoseqaln2 = $stk->select_noncont($maxj, $j);
			$D[$maxi-1][$j-1] = floor(max($twoseqaln1->percentage_identity,$twoseqaln2->percentage_identity));			
		}	
		
		$max = {'max' => 0 , 'maxi' => 1 , 'maxj' => 1};
		#find max in distance matrix
		for my $i (0..$stk->num_sequences-1){
			for my $j (0..$stk->num_sequences-1){						
				if ($D[$i][$j]> $max->{'max'}){
					$max = {'max' => $D[$i][$j] , 'maxi' => $i+1 , 'maxj' => $j+1};				
				}
			}	
		}

		my $in=0;
		for my $x (sort {$a <=> $b} keys %$trees){
			my ($i,$j) = (0,0);
			$trees->{$x}->traverse(sub {
		    	my ($_tree) = @_;
		    	if ($_tree->isLeaf){
		    		$in++ , $i=1 if $_tree->getNodeValue()==$max->{'maxi'};
		    		$in++ , $j=1 if $_tree->getNodeValue()==$max->{'maxj'};
		    	} 
		    });
			
			if( $in==1 ){
				for my $y (sort {$a <=> $b} keys %$trees){
					next if $y<=$x;
					$trees->{$x}->traverse(sub {
				    	my ($_tree) = @_;
				    	$in++ if $_tree->isLeaf && ($_tree->getNodeValue()==$max->{'maxi'} || $_tree->getNodeValue()==$max->{'maxj'} );
				    });
				    
					if( $in==2 ){
						$node = Tree::Simple->new($max->{'max'});
						$node->addChildren($trees->{$x},$trees->{$y});
						$trees->{$x}=$node;
						delete $trees->{$y};	
						last;
					}
				}
				if($in==1){
					$node = Tree::Simple->new($max->{'max'});
					#leaf node stores stk position instead of iid
					$n1 = Tree::Simple->new($i ? $max->{'maxj'} : $max->{'maxi'});				
					$node->addChildren($trees->{$x},$n1);
					$trees->{$x}=$node;
				}
				last;
			}  
		}
		if ( $in==0 ){
			$node = Tree::Simple->new($max->{'max'});
			my $n1 = Tree::Simple->new($max->{'maxi'});			
			my $n2 = Tree::Simple->new($max->{'maxj'});			
			$node->addChildren($n1,$n2);				
			$trees->{++$c}=$node;
		}					
	}	

	#print "\n";	
	# $trees->{(keys %$trees)[0]}->traverse(sub {
	#         my ($_tree) = @_;
	#     	print ((' ' x $_tree->getDepth()), '<',  $_tree->getNodeValue(),'>',"\n");
	#     },
	# 	sub {
	#     	my ($_tree) = @_;
	#     	print ((' ' x $_tree->getDepth()), '</', $_tree->getNodeValue(),'>',"\n");
	# 	}
	# );

	#split tree into sets 
	my $index={}; #stores visited node UID
	my $idx = (sort keys %$trees)[0];
	my $sets={};
	my $setC=0;
	#initialize all sequences as remaining ones - unchecked, if assigned to a set
	my @remainingTMP = (1) x $stk->num_sequences;	
	$trees->{$idx}->traverse(sub {
	        my ($_tree) = @_;
	        if (not exists $index->{$_tree->getUID}){ 	        	    	
		    	if($_tree->getNodeValue()>=60){ #mark node as visited and store leaf nodes = sequences by traversing in sets
		    		$index->{$_tree->getUID}=1; 		    				    		
		        	$_tree->traverse(sub {
		        			my ($_tree) = @_;
		        			if (not exists $index->{$_tree->getUID}){		        			
			        			$index->{$_tree->getUID}=1;			        			
			        			if ($_tree->isLeaf){
			        				push @{$sets->{$setC}} , $_tree->getNodeValue();
			        				#uncheck
			        				$remainingTMP[$_tree->getNodeValue()-1]=0;
			        			}			        			
		        			}
		        		}
		        	);		        		
		        	$setC++;	        	
		        }
	        }	        
	    }
	);
	
	#hold remaining single sequences with iid < 60% 
	my @remaining;
	for (0..$#remainingTMP){		
		push @remaining , $_+1 if $remainingTMP[$_];
	}		
	
	# remove overrepresented sequences in all sets
	for my $key (keys %$sets){
		my $del={};
		for my $i (0..$#{$sets->{$key}}){
			for my $j ($i+1..$#{$sets->{$key}}){
				my @seq1 = split(//,($stk->get_seq_by_pos(${$sets->{$key}}[$i]))->seq);
				my @seq2 = split(//,($stk->get_seq_by_pos(${$sets->{$key}}[$j]))->seq);
				my $ident=0;
				my $length=0;
				for(0..$#seq1){
					next if $seq1[$_]=~/(\.|-)/ && $seq2[$_]=~/(\.|-)/;
					$length++;
					$ident++ if $seq1[$_] eq $seq2[$_];
				}
				# mark overrepresented sequences 
				$del->{$j}=1 if $ident/$length > 0.9;	
			}													
		}
		# remove overrepresented sequences 
		for (reverse sort {$a <=> $b} keys %$del){
			splice @{$sets->{$key}}, $_, 1;
		}
	}	
	
	# remove overrepresented sequences in remaining set (just to be save - should not find anything per definition)
	my $del={};
	for my $i (0..$#remaining){
		for my $j ($i+1..$#remaining){			
			my @seq1 = split(//,($stk->get_seq_by_pos($remaining[$i]))->seq);
			my @seq2 = split(//,($stk->get_seq_by_pos($remaining[$j]))->seq);
			my $ident=0;
			my $length=0;
			for(0..$#seq1){
				next if $seq1[$_]=~/(\.|-)/ && $seq2[$_]=~/(\.|-)/;
				$length++;
				$ident++ if $seq1[$_] eq $seq2[$_];
			}
			
			$del->{$j}=1 if $ident/$length > 0.9;
		}
	}
	for (reverse sort {$a <=> $b} keys %$del){
		splice @remaining, $_, 1;
	}

	#print sets and remaining sequences to fasta files
	my $counter=0;
	my $o =  catfile($self->parameter->cfg->query_dir,$self->parameter->cfg->rf_rna);
	unlink $_ for glob $o.'.set.*.fa';
	unlink $o.'.remaining.fa';
	for my $key (reverse sort {$#{$sets->{$a}} <=> $#{$sets->{$b}}} keys %$sets){		
		$counter++;		
		print $o.'.set.'.$counter.'.fa'."\n";
		open SET , '>'.$o.'.set.'.$counter.'.fa' or die $!;
		for (@{$sets->{$key}}){
			my $seq = $stk->get_seq_by_pos($_);
			print SET '>'.$seq->id."\n";				
			my $s = $seq->seq;
			$s=~s/[._-]//g;
			$s=~y/acgtUu/ACGTTT/;						
			print SET $s."\n";
		}
		close SET;
	}
	if ($#remaining > -1){
		open SET , '>'.$o.'.remaining.fa' or die $!;
		for (@remaining){
			my $seq = $stk->get_seq_by_pos($_);
			print SET '>'.$seq->id."\n";
			my $s = $seq->seq;
			$s=~s/[._-]//g;
			$s=~y/acgtUu/ACGTTT/;			
			print SET $s."\n";
		}
		close SET;
	}	
}

sub build_cm_from_sets {
	my ($self) = @_;

	my @sets;
	#read and sort all sets of sequences with idd < 60% between sets and > 60% in each set into array of hashrefs
	my $o = catfile($self->parameter->cfg->query_dir,$self->parameter->cfg->rf_rna);	
	for (glob $o.'.set.*.fa'){
		my $set;
		my $seqio = Bio::SeqIO->new( '-format' => 'Fasta' , -file => $_, -verbose => -1);
		while(my $seqobj = $seqio->next_seq()){
			$set->{$seqobj->id}=$seqobj->seq;
		}
		push @sets , $set;
	}	
	@sets = sort {scalar keys %$a <=> scalar keys %$b} @sets;
	
	#read remaining single sequences with iid < 60 % to all others sequences into hashref
	my $remaining={};
	if (-e $o.'.remaining.fa'){ 
		my $seqio = Bio::SeqIO->new( -format => 'Fasta' , -file => $o.'.remaining.fa', -verbose => -1);
		while(my $seqobj = $seqio->next_seq()){
			$remaining->{$seqobj->id}=$seqobj->seq;
		}
	}
	
	my $training;
	my $testing;
	
	#constrain for training set is #sequences >=5 
	if($#sets==-1){
		#use remaining only 
		if(scalar keys %$remaining > 5){
			for (keys %$remaining){
				if (scalar keys %$training < 5){
					$training->{$_}=$remaining->{$_};							
				} else {
					$testing->{$_}=$remaining->{$_};					
				}
			}				
		}					
	} else {
		#fill training set with sequences, beginning with smallest set
		my $usedSets;
		my $setseqC = 0;
		for my $i (0..$#sets){ 
			$setseqC += scalar keys %{$sets[$i]};			
			if($setseqC + (scalar keys %$remaining) > 4){				
				for (my $j=$i; $j>=0; $j--){
					$usedSets->{$j} = 1;
					$training->{$_} = $sets[$j]->{$_} for keys %{$sets[$j]};
					last if (scalar keys %$training) + (scalar keys %$remaining) > 4
				}
				last;
			}
		}
		for my $i (0..$#sets){
			next if exists $usedSets->{$i};
			$testing->{$_} = $sets[$i]->{$_} for keys %{$sets[$i]};			
		}
		for (keys %$remaining){
			if (scalar keys %$training < 5){
				$training->{$_} = $remaining->{$_};			
			} else {
				$testing->{$_} = $remaining->{$_};
			}
		}
	}
		
	if(scalar keys %$testing > 0 && scalar keys %$training > 4){		
		open FA , '>'.$o.'.training.fa' or die $!;
		for (keys %$training){
			print FA '>'.$_."\n";
			print FA $training->{$_}."\n";
		}
		close FA;
		open FA , '>'.$o.'.testing.fa' or die $!;
		for (keys %$testing){
			print FA '>'.$_."\n";
			print FA $testing->{$_}."\n";
		}
		close FA;
		unlink $o.'.training.stk';
		unlink $o.'.training.cm';
				
		system("cmalign --mxsize 3000 --noprob --cpu ".$self->parameter->threads." -o ".$o.".training.stk ".$self->parameter->cfg->cm." ".$o.".training.fa" );		
		system("cmbuild -F ".$o.".training.cm ".$o.".training.stk");
		system("cmcalibrate --cpu ".$self->parameter->threads." ".$o.".training.cm");		
	} 
}

1;