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
	my ($self, $stk, $features) = @_;

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
			if ($ssPresentCount->{$f->seq_id} < ($#ssAreas +1)/2){
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
	my ($self, $stk, $features, $cssep , $csindels) = @_;	

	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';

	my @update;
	my $write;

	my $out;
	open my $READER, '>', \$out;
	(Bio::AlignIO->new(-format  => 'stockholm', -fh => $READER, -verbose => -1 ))->write_aln($stk);
	close $READER;
	my @lines = split /\n/ , $out;
	my @cs;
	for(@lines){		
		if ($_=~/^#=GC\s+RF\s+(.+)/){			
			my $line=$1;
			chomp $line;
			push @cs , split(// , $line);
		}
	}
	
	my $csposstart;
	my $csposend = 0;
	my $rmfeature;
	for my $i (0..$#{$cssep}){		
		my $c = 0;
		for((defined $csposstart ? $csposstart : 0)..$#cs){
			if ($cs[$_]=~/\w/){
				$c++;
				$csposstart = $_ unless defined $csposstart;
			}
			if ($c == ${$cssep}[$i] || $_ == $#cs){
				$csposend = $_;
				last;
			}
		}	

		if (${$csindels}[$i]>-1){
			#for each feature seq get csstart csstop subseq to compare with cs[csstart..csend];
			for my $k (keys %{$features}){
				next if exists $rmfeature->{$k};
				my $f = $features->{$k};		
				my @seq = split // , ($stk->get_seq_by_id($f->seq_id))->subseq($csposstart+1,$csposend+1);	
				my $mm = 0;
				for (0..$#seq){	
					if ($cs[$csposstart+$_]=~/\w/){	
						switch(lc($cs[$csposstart+$_])){
							case /[acgtu]/ {
								$mm++ unless lc($seq[$_]) eq lc($cs[$csposstart+$_]);
							}
							case "r" {
								$mm++ unless lc($seq[$_])=~/[ag]/;
							}
							case "y" {
								$mm++ unless lc($seq[$_])=~/[ctu]/;	
							}
							case "s" {
								$mm++ unless lc($seq[$_])=~/[gc]/;	
							}
							case "w" {
								$mm++ unless lc($seq[$_])=~/[atu]/;	
							}
							case "k" {
								$mm++ unless lc($seq[$_])=~/[gtu]/;	
							}
							case "m" {
								$mm++ unless lc($seq[$_])=~/[ac]/;	
							}
							case "b" {
								$mm++ unless lc($seq[$_])=~/[cgtu]/;	
							}
							case "d" {
								$mm++ unless lc($seq[$_])=~/[agtu]/;	
							}
							case "h" {
								$mm++ unless lc($seq[$_])=~/[actu]/;	
							}
							case "v" {
								$mm++ unless lc($seq[$_])=~/[acg]/;	
							}
							case "n" {
								$mm++ unless lc($seq[$_])=~/[acgtu]/;	
							}
							else {}
						}						
					} else {
						$mm++ if $seq[$_]=~/\w/;
					}
				}	
				$rmfeature->{$k}=1 if $mm > ${$csindels}[$i];
			}			
		}
		$csposstart=$csposend+1;
	}

	for(keys %$rmfeature){
		my $f = $features->{$_};				
		$write = 1;
		$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
		push @update , $f->seq_id.' '.$f->primary_tag.' P';				
		delete $features->{$_};
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
