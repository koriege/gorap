package Bio::Gorap::Functions::STK;

use Bio::AlignIO;
use Bio::SimpleAlign;
use POSIX;
use Switch;

sub score_filter {
	my ($self, $stk, $features, $threshold, $nonTaxThreshold) = @_;

	my $c=0;
	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';
	my $write;
	my @update;

	if ($nonTaxThreshold){		
		my $tmpfeatures;
		for (keys %{$features}){		
			my $f = $features->{$_};
			next if $f->score eq '.';

			my @id = split /\./ , $f->seq_id;
			my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
			$tmpfeatures->{$abbr}++;

			if (($f->get_tag_values('source'))[0] =~ /infernal/ || ($f->get_tag_values('source'))[0] =~ /blast/){						
				if (($f->get_tag_values('origscore'))[0] < $threshold){
					$tmpfeatures->{$abbr}--;
				}
			} else {
				if ($f->score < $threshold){
					$tmpfeatures->{$abbr}--;
				}
			}
		}

		for (keys %{$features}){
			my $f = $features->{$_};
			next if $f->score eq '.';

			my @id = split /\./ , $f->seq_id;
			my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
			if ($tmpfeatures->{$abbr}==0){
				if (($f->get_tag_values('source'))[0] =~ /infernal/ || ($f->get_tag_values('source'))[0] =~ /blast/){						
					#print $f->seq_id." ".$f->score." ".ceil(($f->get_tag_values('origscore'))[0]).' '.$threshold."\n";
					if (($f->get_tag_values('origscore'))[0] < $nonTaxThreshold){
						delete $features->{$_};
						$write = 1;
						$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
						push @update , $f->seq_id.' '.$f->primary_tag.' B';
					}
				} else {
					if ($f->score < $nonTaxThreshold){
						delete $features->{$_};
						$write = 1;
						$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
						push @update , $f->seq_id.' '.$f->primary_tag.' B';
					}
				}
			} else {
				if (($f->get_tag_values('source'))[0] =~ /infernal/ || ($f->get_tag_values('source'))[0] =~ /blast/){						
					#print $f->seq_id." ".$f->score." ".ceil(($f->get_tag_values('origscore'))[0]).' '.$threshold."\n";
					if (($f->get_tag_values('origscore'))[0] < $threshold){
						delete $features->{$_};
						$write = 1;
						$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
						push @update , $f->seq_id.' '.$f->primary_tag.' B';
					}
				} else {
					if ($f->score < $threshold){
						delete $features->{$_};
						$write = 1;
						$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
						push @update , $f->seq_id.' '.$f->primary_tag.' B';
					}
				}
			}
		}

	} else {
		for (keys %{$features}){		
			my $f = $features->{$_};
			next if $f->score eq '.';
			
			if (($f->get_tag_values('source'))[0] =~ /infernal/ || ($f->get_tag_values('source'))[0] =~ /blast/){						
				#print $f->seq_id." ".$f->score." ".ceil(($f->get_tag_values('origscore'))[0]).' '.$threshold."\n";
				if (($f->get_tag_values('origscore'))[0] < $threshold){
					delete $features->{$_};
					$write = 1;
					$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
					push @update , $f->seq_id.' '.$f->primary_tag.' B';
				}
			} else {
				if ($f->score < $threshold){
					delete $features->{$_};
					$write = 1;
					$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
					push @update , $f->seq_id.' '.$f->primary_tag.' B';
				}
			}
		}
	}

	return ($stk , $features, \@update , $write);
}

sub structure_filter(){
	my ($self, $stk, $features, $minstructures) = @_;

	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';
	my @update;
	my $write;

	my @ss;

	my $out;
	open my $READER, '>', \$out;
	(Bio::AlignIO->new(-format  => 'stockholm', -fh => $READER, -verbose => -1 ))->write_aln($stk);
	close $READER;
	my @lines = split /\n/ , $out;
	for(@lines){		
		if ($_=~/^#=GC\s+SS_cons\s+(.+)/){			
			my $line=$1;
			chomp $line;
			$line=~s/[{\[<]/(/g;
			$line=~s/[}\]>]/)/g;
			push @ss , split(// , $line);
		}
	}

	my @so;
	my @ssAreas;
	
	for (my $i=0 ; $i < $#ss ; $i++){
		if ($ss[$i] eq '('){
			push @so , $i;	
		} else {
			next if $ss[$i] ne ')'; 	
			
			my $pos = pop @so;
			$pos = pop @so while(not defined $pos);
				 
			my @areaStart = ($pos);
			my @areaStop= ($i); 		
			while($i < $#ss && $ss[++$i] ne '(' && defined $so[$#so]){			
				next if $ss[$i] ne ')';			 			
				push @areaStart , pop @so;
				push @areaStop , $i;
			} 
			$i--;		
			push @ssAreas , {$areaStart[$#areaStart] => $areaStart[0]};
			push @ssAreas , {$areaStop[0] => $areaStop[$#areaStop]};
			push @so , undef;				
		}
	}

	return ($stk , $features, \@update , $write) if $#ssAreas==-1;
	
	my $ssPresentCount;

	for (@ssAreas){
		my ($sta, $sto) = each(%$_);
		my $tmpStk = $stk->slice($sta+1, $sto+1);

		for (keys %{$features}){		
			my $f = $features->{$_};
			$ssPresentCount->{$f->seq_id}++ if $tmpStk->get_seq_by_id($f->seq_id);						
		}
	}
	
	$minstructures = ($#ssAreas +1)/2 unless $minstructures;	
	for (keys %{$features}){
		my $f = $features->{$_};
		if ($#ssAreas == 1){
			if ($ssPresentCount->{$f->seq_id} < 1){
				delete $features->{$_};
				$write = 1;
				$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
				push @update , $f->seq_id.' '.$f->primary_tag.' S';
			}
		} else {
			if ($ssPresentCount->{$f->seq_id} < $minstructures){
				delete $features->{$_};
				$write = 1;
				$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
				push @update , $f->seq_id.' '.$f->primary_tag.' S';
			}
		}
	}
	
	return ($stk , $features, \@update , $write);
}

sub sequence_filter { #$_!~/_snora/i && $_=~/_sn?o?(r|z|u|y)d?\d/i 
	my ($self, $stk, $features) = @_;

	my $c=0;
	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';
	my @update;
	my $write;

	my $tmpStk = $stk->select(1, $stk->num_sequences);
	for (keys %{$features}){
		my $f = $features->{$_};		
		$tmpStk->remove_seq($tmpStk->get_seq_by_id($f->seq_id));
	}
	my @consensus = split //, $tmpStk->consensus_string(95);

	for (keys %{$features}){
		my $f = $features->{$_};		
		my @seq = split // , ($stk->get_seq_by_id($f->seq_id))->seq;	
		# print $f->seq_id."\n";
		# print join '' , @consensus; print "\n";
		# print join '' , @seq; 
		my $consC=0;
		my $allCons=scalar(@seq); 
		for my $i (0..$#seq){
			if($consensus[$i] eq '?'){
				$allCons--;
			} else {				
				if($seq[$i]!~/[acgtuACGTU]/){
					$allCons--;
				} else {
					$consC++ if $seq[$i] eq $consensus[$i];	
				}				
			}
		}
		#print "\n",$consC," ",scalar(@seq)," ",$allCons," ",($consC/$allCons)," ",($allCons/scalar(@seq))," ------------\n";
		if($allCons>0 && $consC > 9){			
			if ($consC/$allCons < 0.7 ){
				delete $features->{$_};
				$write = 1;
				$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
				push @update , $f->seq_id.' '.$f->primary_tag.' P';
			}
		}		
	}
	
	return ($stk , $features, \@update , $write);
}

sub user_filter {
	my ($self, $stk, $features, $constrains, $cs, $seedstk) = @_;	

	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';
	my @update;
	my $write;

	my $seedseqo = (Bio::AlignIO->new(-format  => 'stockholm', -file => $seedstk, -verbose => -1 ))->next_aln->get_seq_by_pos(1);	
    my $seq = $seedseqo->seq;    
    $seq=~s/\./-/g;
 
	my @seedseq = split // , $seq;
    $seq = ($stk->get_seq_by_id($seedseqo->id))->seq;
    $seq=~s/\./-/g;
 
	my @newseedseq = split // , $seq;
	my @seedcs;
	open STK , '<'.$seedstk or die $!;
	while(<STK>){
		chomp $_;
		if ($_=~/^#=GC\s+RF\s+(.+)/){			
			push @seedcs , split( // , $1);
		}
	}
	close STK;

	my $out;
	open my $READER, '>', \$out;
	(Bio::AlignIO->new(-format  => 'stockholm', -fh => $READER, -verbose => -1 ))->write_aln($stk);
	close $READER;
	my @lines = split /\n/ , $out;
	my $newcs;
	for(@lines){
		chomp $_;
		if ($_=~/^#=GC\s+RF\s+(.+)/){
			$newcs.=$1;
		}
	}

	my @seedpos;
	my @cs = split(//,$cs);
	my $c=0;
	for my $i (0..$#cs){
		$c++ while lc($seedcs[$i+$c]) !~ /\w/;
		push @seedpos , $i+$c;
	}

	my @csNewpos;
	$c=0;
	for my $i (0..$#cs){
		$c++ while lc($seedseq[$seedpos[$i]]) ne lc($newseedseq[$seedpos[$i]+$c]);
		push @csNewpos , $seedpos[$i]+$c;
	}
	
	for my $k (keys %{$features}){				
    	my $f = $features->{$k};

		# my $foo = ($stk->get_seq_by_id($f->seq_id))->seq;
		my $presentu;
		for (0..$#{$constrains}){			

			my ($sta,$sto,$mm,$query) = @{$$constrains[$_]};
			my @query = split // , $query;
			$sta = $csNewpos[$sta-1]+1;
			$sto = $csNewpos[$sto-1]+1;			
			my @subseq = split // , ($stk->get_seq_by_id($f->seq_id))->subseq($sta,$sto);
			my $cssubseq = substr($newcs,$sta-1,$sto-$sta+1);			
			$subseq =~ s/\W//g;	
			if($f->type=~/_Afu/ || $f->type=~/_SNORD/ || $f->type=~/_sn?o?s?n?o?[A-WYZ]+[a-z]?\d/){
				$presentu++ if $_ == 0 && join('',$subseq[1..$#subseq])=~/[uU]/;
				if ($_ == 1){
					$presentu++ if join('',$subseq)=~/[uU]/;
					unless ($presentu){
						delete $features->{$k};
						$write = 1;
						$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
						push @update , $f->seq_id.' '.$f->primary_tag.' P';
						last;
					}				
				}
			}
			my $i=-1;
            my $j=-1;
			for(split // , $cssubseq){
                $i++;
                if($_=~/\w/){
                  $j++;
				} else {
					$mm-- if $subseq[$i]=~/\w/;
					next;
				}
				switch(lc($query[$j])){
					case /[acgtu]/ {
						$mm-- unless lc($subseq[$i]) eq lc($query[$j]);
					}
					case "r" {
						$mm-- unless lc($subseq[$i])=~/[ag]/;
					}
					case "y" {
						$mm-- unless lc($subseq[$i])=~/[ctu]/;	
					}
					case "s" {
						$mm-- unless lc($subseq[$i])=~/[gc]/;	
					}
					case "w" {
						$mm-- unless lc($subseq[$i])=~/[atu]/;	
					}
					case "k" {
						$mm-- unless lc($subseq[$i])=~/[gtu]/;	
					}
					case "m" {
						$mm-- unless lc($subseq[$i])=~/[ac]/;	
					}
					case "b" {
						$mm-- unless lc($subseq[$i])=~/[cgtu]/;	
					}
					case "d" {
						$mm-- unless lc($subseq[$i])=~/[agtu]/;	
					}
					case "h" {
						$mm-- unless lc($subseq[$i])=~/[actu]/;	
					}
					case "v" {
						$mm-- unless lc($subseq[$i])=~/[acg]/;	
					}
					case "n" {
						$mm-- unless lc($subseq[$i])=~/[acgtu]/;	
					}
					else {}
				}						
			}
			# print $mm." ".$cssubseq." ".$query." ".join("",@subseq)."\n";			
			if ($mm < 0){
				delete $features->{$k};
				$write = 1;
				$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
				push @update , $f->seq_id.' '.$f->primary_tag.' P';
				last;
			}
		}				
	}

	return ($stk , $features, \@update , $write);
}

sub get_ss_cs_from_file {
	my ($self,$stk) = @_;

	my $s='';
	my $c='';
	open STK, '<'.$stk or die $!;
	while(<STK>){
		$s.=$1 if $_=~/^#=GC\sSS_cons\s+(.+)/;
		$c.=$1 if $_=~/^#=GC\sRF\s+(.+)/;
	}
	close STK;

	return ($s,$c);
}

sub get_ss_cs_from_object {
	my ($self,$stk) = @_;

	open my $READER, '>', \$out;
	(Bio::AlignIO->new(-format  => 'stockholm', -fh => $READER, -verbose => -1 ))->write_aln($stk);
	close $READER;
	my @lines = split /\n/ , $out;
	my $s='';
	my $c='';
	for(@lines){		
		$s.=$1 if $_=~/^#=GC\sSS_cons\s+(.+)/;
		$c.=$1 if $_=~/^#=GC\sRF\s+(.+)/;
	} 

	return ($s,$c);
}

1;
