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

	if ($nonTaxThreshold){ #gorap was startet with taxonomy - $threshold is taxonomy based
		my $tmpfeatures;
		for (keys %{$features}){		
			my $f = $features->{$_};
			next if $f->score eq '.';

			my @id = split /\./ , $f->seq_id;
			my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
			$tmpfeatures->{$abbr}++;

			if (($f->get_tag_values('source'))[0] =~ /infernal/ || ($f->get_tag_values('source'))[0] =~ /blast/){						
				if (($f->get_tag_values('origscore'))[0] < $threshold){ #check for hits witch cmsearch score - usually higher than cmalign scores
					$tmpfeatures->{$abbr}--;
				}
			} else {
				if ($f->score < $threshold){
					$tmpfeatures->{$abbr}--;
				}
			}
		}

		$threshold = $threshold-$threshold*0.1;

		for (keys %{$features}){
			my $f = $features->{$_};
			next if $f->score eq '.';

			my @id = split /\./ , $f->seq_id;
			my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
			if ($tmpfeatures->{$abbr}==0){ #if hits, but not in this species, lower threshold again 
				if (($f->get_tag_values('source'))[0] =~ /infernal/ || ($f->get_tag_values('source'))[0] =~ /blast/){						
					#print $f->seq_id." ".$f->score." ".ceil(($f->get_tag_values('origscore'))[0]).' '.$threshold."\n";
					if (($f->get_tag_values('origscore'))[0] < $nonTaxThreshold && ($f->get_tag_values('origscore'))[0] < $threshold){
						delete $features->{$_};
						$write = 1;
						$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
						push @update , $f->seq_id.' '.$f->primary_tag.' B';
					}
				} else {
					if ($f->score < $nonTaxThreshold && $f->score < $threshold){
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
	$cs=~s/\W/-/g;

	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';
	my @update;
	my $write;

	my $seedseqo = (Bio::AlignIO->new(-format  => 'stockholm', -file => $seedstk, -verbose => -1 ))->next_aln->get_seq_by_pos(1);	
    my $seq = lc($seedseqo->seq);
    $seq=~s/\W/-/g;

	my @seedseq = split // , lc($seq);
    $seq = ($stk->get_seq_by_id($seedseqo->id))->seq;
    $seq=~s/\W/-/g;

	my @stkseedseq = split // , lc($seq);
	my @seedcs;
	open STK , '<'.$seedstk or die $!;
	while(<STK>){
		chomp $_;
		if ($_=~/^#=GC\s+RF\s+(.+)/){
			my $s = lc($1);
			$s=~s/\W/-/g;
			push @seedcs , split( // , $s);
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
			my $s = lc($1);
			$s=~s/\W/-/g;
			$newcs.=$s;
		}
	}

	#mapping of cfg cs to seed stk pos
	my @cs_pos_in_seed;
	my @cs = split(//,$cs);
	my $c=0;
	for my $i (0..$#cs){
		$c++ while $i+$c <= $#seedcs && $seedcs[$i+$c] !~ /\w/;
		push @cs_pos_in_seed , $i+$c;
	}

	#mapping of seed stk seq to new stk seq
	my @seedseq_pos_in_stk;	
	$c=0;	
	my $j=0;
	# print ''.join('',@seedseq)."\n";
	# print ''.join('',@stkseedseq)."\n";
	for (my $i=0; $i+$j<=$#seedseq; $i++){
		if ($seedseq[$i+$j]=~/\W/ && $stkseedseq[$i+$c]=~/\W/){
			# print $seedseq[$i+$j]." 1 ".$stkseedseq[$i+$c]."\n";
			push @seedseq_pos_in_stk , $i+$c;
		} elsif($seedseq[$i+$j]=~/\w/) {
			# print $seedseq[$i+$j]." 2 ".$stkseedseq[$i+$c]." before\n";
			$c++ while $i+$c <= $#stkseedseq && $seedseq[$i+$j] ne $stkseedseq[$i+$c];
			# print $seedseq[$i+$j]." 2 ".$stkseedseq[$i+$c]."\n";
			push @seedseq_pos_in_stk , $i+$c;
		} else {			
			while ($i+$j <= $#seedseq && $seedseq[$i+$j] ne $stkseedseq[$i+$c]){				
				# print $seedseq[$i+$j]." 3 ".$stkseedseq[$i+$c]."\n";
				push @seedseq_pos_in_stk , $i+$c;
				$j++;
			}
			# print $seedseq[$i+$j]." 3 ".$stkseedseq[$j+$c]."\n";
			push @seedseq_pos_in_stk , $j+$c;			
		}		
	}
	
	for my $k (keys %{$features}){				
    	my $f = $features->{$k};

		# my $foo = ($stk->get_seq_by_id($f->seq_id))->seq;

		my $presentu;
		my @uga_ug;
		for my $ci (0..$#{$constrains}){			

			my ($sta,$sto,$mm,$query) = @{$$constrains[$ci]};
			my @query = split // , lc($query);
			$sta = $seedseq_pos_in_stk[$cs_pos_in_seed[$sta-1]]+1;
			$sto = $seedseq_pos_in_stk[$cs_pos_in_seed[$sto-1]]+1;
			
			my @subseq = split // , lc(($stk->get_seq_by_id($f->seq_id))->subseq($sta,$sto));
			my $cssubseq = substr($newcs,$sta-1,$sto-$sta+1);			
			$subseq =~ s/\W//g;	
			if($f->type=~/_Afu/ || $f->type=~/_SNOR.?D/ || $f->type=~/_sn?o?s?n?o?[A-WYZ]+[a-z]?\d/){
				$presentu++ if $ci == 0 && join('',$subseq[1..$#subseq])=~/[uU]/;
				if ($ci == 1){
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
            my $last;
			for(split // , $cssubseq){
                $i++;
                if($_=~/\w/){
                  $j++;
				} else {
					$mm-- if $subseq[$i]=~/\w/;
					next;
				}
				#check covalent bp of uga(ug)a and (cu)ga
				if ($f->type=~/_Afu/ || $f->type=~/_SNOR.?D/ || $f->type=~/_sn?o?s?n?o?[A-WYZ]+[a-z]?\d/){
					if ($ci==0){						
						push @uga_ug , $subseq[$i] if $j==3 || $j==4;							
					} elsif ($ci==1){
						# print "uga(ug)a ".join('',@uga_ug)." <> " if $j==0;
						# print $subseq[$i] if $j==0;
						# print $subseq[$i]." (cu)ga\n" if $j==1;
						if ($j==0 && $subseq[$i]=~/\w/ && $uga_ug[1]=~/\w/){
							unless (($subseq[$i] eq 'a' && $uga_ug[1] eq 'u') || 
									($subseq[$i] eq 'c' && $uga_ug[1] eq 'g') || 
									($subseq[$i] eq 'g' && $uga_ug[1] eq 'u') || 
									($subseq[$i] eq 'g' && $uga_ug[1] eq 'c') ||
									($subseq[$i] eq 'u' && $uga_ug[1] eq 'a') ){
								delete $features->{$k};
								$write = 1;
								$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
								push @update , $f->seq_id.' '.$f->primary_tag.' P';
								$last = 1;
							}
						}elsif($j==1 && $subseq[$i]=~/\w/ && $uga_ug[0]=~/\w/){
							unless (($subseq[$i] eq 'a' && $uga_ug[0] eq 'u') || 
									($subseq[$i] eq 'c' && $uga_ug[0] eq 'g') || 
									($subseq[$i] eq 'g' && $uga_ug[0] eq 'u') || 
									($subseq[$i] eq 'g' && $uga_ug[0] eq 'c') ||
									($subseq[$i] eq 'u' && $uga_ug[0] eq 'a') ||
									($subseq[$i] eq 'u' && $uga_ug[0] eq 'u') ){
								delete $features->{$k};
								$write = 1;
								$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
								push @update , $f->seq_id.' '.$f->primary_tag.' P';
								$last = 1;
							}
						}
					}					
				}
				last if $last;

				switch($query[$j]){
					case /[acgtu]/ {
						$mm-- unless $subseq[$i] eq $query[$j];
					}
					case "r" {
						$mm-- unless $subseq[$i]=~/[ag]/;
					}
					case "y" {
						$mm-- unless $subseq[$i]=~/[ctu]/;	
					}
					case "s" {
						$mm-- unless $subseq[$i]=~/[gc]/;	
					}
					case "w" {
						$mm-- unless $subseq[$i]=~/[atu]/;	
					}
					case "k" {
						$mm-- unless $subseq[$i]=~/[gtu]/;	
					}
					case "m" {
						$mm-- unless $subseq[$i]=~/[ac]/;	
					}
					case "b" {
						$mm-- unless $subseq[$i]=~/[cgtu]/;	
					}
					case "d" {
						$mm-- unless $subseq[$i]=~/[agtu]/;	
					}
					case "h" {
						$mm-- unless $subseq[$i]=~/[actu]/;	
					}
					case "v" {
						$mm-- unless $subseq[$i]=~/[acg]/;	
					}
					case "n" {
						$mm-- unless $subseq[$i]=~/[acgtu]/;	
					}
					else {}
				}						
			}
			last if $last;
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
