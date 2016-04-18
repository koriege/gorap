package Bio::Gorap::Functions::STK;

use Bio::AlignIO;
use Bio::SimpleAlign;
use POSIX;
use Switch;
use List::Util qw(min max);
use File::Spec::Functions;

sub length_filter {
	my ($self, $stk, $features, $seedcs) = @_;

	my $c=0;
	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';
	my $write;
	my @update;	
	for (keys %{$features}){
		my $f = $features->{$_};
		if ( ($f->source =~ /blast/ ? $f->score : max($f->score,($f->get_tag_values('origscore'))[0])) < 40 && ($f->stop - $f->start < length($seedcs)/2.5)){
			delete $features->{$_};
			$write = 1;
			$stk->remove_seq($stk->get_seq_by_id($f->seq_id));			
			push @update , $f->seq_id.' '.$f->primary_tag.' L';
		}
	}

	return ($stk , $features, \@update , $write);
}

sub score_filter {
	my ($self, $nofilter, $userfilter, $stk, $features, $threshold, $nonTaxThreshold) = @_;

	my $c=0;
	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';
	my $write;
	my @update;	
	my $type = $features->{(keys %{$features})[0]}->primary_tag;
	
	#HACAS: _SNORA _snopsi
	if ( ! $nofilter && $userfilter && ($type=~/_Afu/ || $type=~/_SNOR[ND\d]/ || $type=~/_sn\d/ || $type=~/_sno[A-Z]/ || $type=~/(-|_)sn?o?s?n?o?[A-WYZ]+[a-z]?-?\d/)){
		for (keys %{$features}){		
			my $f = $features->{$_};
			next if $f->score eq '.';

			if ( ($f->source =~ /blast/ ? $f->score : max($f->score,($f->get_tag_values('origscore'))[0])) < 8){
				delete $features->{$_};
				$write = 1;
				$stk->remove_seq($stk->get_seq_by_id($f->seq_id));
				push @update , $f->seq_id.' '.$f->primary_tag.' B';
			}
		}

		return ($stk , $features, \@update , $write);
	}	

	if ($nonTaxThreshold){ #gorap was startet with taxonomy - $threshold is taxonomy based
		if ($type=~/_sn/i){
			for (keys %{$features}){
				my $f = $features->{$_};
				next if $f->score eq '.';

				my @id = split /\./ , $f->seq_id;
				my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
				
				if ($f->source =~ /infernal/ || $f->source =~ /blast/){
					my $score = $f->source =~ /blast/ ? $f->score : max($f->score,($f->get_tag_values('origscore'))[0]);
					if ($score < min($nonTaxThreshold,$threshold)){
						delete $features->{$_};
						$write = 1;
						$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
						push @update , $f->seq_id.' '.$f->primary_tag.' B';
					}
				} else {
					if ($f->score < 8 || $f->score < min($nonTaxThreshold,$threshold) ){
						delete $features->{$_};
						$write = 1;
						$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
						push @update , $f->seq_id.' '.$f->primary_tag.' B';
					}
				}			
			}
		} else {	
			for (keys %{$features}){
				my $f = $features->{$_};
				next if $f->score eq '.';

				my @id = split /\./ , $f->seq_id;
				my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
				
				if ($f->source =~ /infernal/ || $f->source =~ /blast/){
					my $score = $f->source =~ /blast/ ? $f->score : max($f->score,($f->get_tag_values('origscore'))[0]);
					if ($score < $threshold){
						delete $features->{$_};
						$write = 1;
						$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
						push @update , $f->seq_id.' '.$f->primary_tag.' B';
					}
				} else {
					if ($f->score < 8 || $f->score < min($nonTaxThreshold,$threshold) ){
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
			
			if ($f->source =~ /infernal/ || $f->source =~ /blast/){
				my $score = $f->source =~ /blast/ ? $f->score : max($f->score,($f->get_tag_values('origscore'))[0]);
				if ($score < $threshold){
					delete $features->{$_};
					$write = 1;
					$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
					push @update , $f->seq_id.' '.$f->primary_tag.' B';
				}
			} else {
				if ($f->score < 8 || $f->score < $threshold){
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

	my ($ss , $cs) = &get_ss_cs_from_object($self,$stk);
	$ss=~s/[{\[<]/(/g;
	$ss=~s/[}\]>]/)/g;
	my @ss = split // , $ss;

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
			my $seqo = $tmpStk->get_seq_by_id($f->seq_id);
			if ($seqo){
				my $ncc = $seqo->seq;				
				$ncc=~s/\W//g;
				my $minl = $sto - $sta + 1 > 4 ? 1 : 0;				
				$ssPresentCount->{$f->seq_id}++ if length $ncc > $minl;
			}			
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

sub sequence_filter {
	my ($self, $stk, $features) = @_;

	my $c=0;
	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';
	my $type = $features->{(keys %{$features})[0]}->primary_tag;
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

		if($allCons>0 && $consC > 9){			
			# if (($f->get_tag_values('source'))[0] =~ /infernal/ || ($f->get_tag_values('source'))[0] =~ /blast/){
			# 	if (max($f->score,($f->get_tag_values('origscore'))[0]) < 30 && ($consC/$allCons < 0.7 || ($type=~/_mir/i && $consC/$allCons < 0.9))){
			# 		delete $features->{$_};
			# 		$write = 1;
			# 		$stk->remove_seq($stk->get_seq_by_id($f->seq_id));
			# 		push @update , $f->seq_id.' '.$f->primary_tag.' P';	
			# 	}
			# } elsif ($consC/$allCons < 0.7 || ($type=~/_mir/i && $consC/$allCons < 0.9)){
			if ($consC/$allCons < 0.7 || ($type=~/_mir/i && $consC/$allCons < 0.9)){
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
	my ($seedss , $seedcs) = &get_ss_cs_from_file($self,$seedstk);
	my @seedcs = split // , $seedcs;
	my ($newss , $newcs) = &get_ss_cs_from_object($self,$stk);

	#mapping of cfg cs to seed stk pos
	my @cs_pos_in_seed;
	my @cs = split(//,$cs);
	my $c=0;
	for my $i (0..$#cs){
		$c++ while $i+$c <= $#seedcs && $seedcs[$i+$c] !~ /\w/;
		push @cs_pos_in_seed , $i+$c;
	}

	# mapping of seed stk seq to new stk seq
	my @seedseq_pos_in_stk;		
	my $j=0;
	for (my $i=0; $i<=$#seedseq; $i++){
		$c = $i;
		$i++ while $seedseq[$i]=~/\W/;		
		#print $seedseq[$i].' '.$stkseedseq[$j]." $i  $j\n";
		
		if ($seedseq[$i] eq $stkseedseq[$j]){
			$j++;
		} else {			
			$j++ while $seedseq[$i] ne $stkseedseq[$j];	
			$i--;
		}

		push @seedseq_pos_in_stk , ($j)x($i-$c+1);
		# print ($j)x($i-$c+1);
		# exit if $i >9;
	}
	# print $#seedseq." ".$#seedseq_pos_in_stk."\n";

	# print ''.join('',@cs)."\n";
	# print ''.join('',@seedseq)."\n";
	# for (0..$#{$constrains}){
	# 	my ($sta,$sto,$mm,$query) = @{$$constrains[$_]};
	# 	print $sta." ".$sto." ".join('',@cs[$sta-1..$sto-1])."\n";
	# 	print $cs_pos_in_seed[$sta-1]." ".$cs_pos_in_seed[$sto-1].' '.join('',@seedseq[$cs_pos_in_seed[$sta-1]..$cs_pos_in_seed[$sto-1]])."\n";
	# }

	# print ''.join('',@stkseedseq)."\n";
	# # print ''.join('',@seedseq_pos_in_stk)."\n";
	# for (0..$#{$constrains}){
	# 	my ($sta,$sto,$mm,$query) = @{$$constrains[$_]};
	# 	print $sta." ".$sto." ".join('',@cs[$sta-1..$sto-1])."\n";
	# 	print $seedseq_pos_in_stk[$cs_pos_in_seed[$sta-1]-1]." ".$seedseq_pos_in_stk[$cs_pos_in_seed[$sto-1]-1].' '.join('',@stkseedseq[$seedseq_pos_in_stk[$cs_pos_in_seed[$sta-1]-1]..$seedseq_pos_in_stk[$cs_pos_in_seed[$sto-1]-1]])."\n";
	# }	
	
	for my $k (keys %{$features}){				
		my $f = $features->{$k};
		# print $f->seq_id."\n";
		# print ''.join('',@stkseedseq)."\n";
		# print ''.($stk->get_seq_by_id($f->seq_id))[0]->seq."\n";
		# for (0..$#{$constrains}){
		# 	my ($sta,$sto,$mm,$query) = @{$$constrains[$_]};
		# 	print $sta." ".$sto." ".join('',@cs[$sta-1..$sto-1])."\n";
		# 	print $seedseq_pos_in_stk[$cs_pos_in_seed[$sta-1]-1]." ".$seedseq_pos_in_stk[$cs_pos_in_seed[$sto-1]-1].' '.lc(($stk->get_seq_by_id($f->seq_id))[0]->subseq($seedseq_pos_in_stk[$cs_pos_in_seed[$sta-1]-1]+1,$seedseq_pos_in_stk[$cs_pos_in_seed[$sto-1]-1]+1))."\n";
		# }		

		next if $f->score ne '.' && $f->score > 30 && ($f->primary_tag=~/_sno[A-Z]/ || $f->primary_tag=~/_Afu/ || $f->primary_tag=~/_SNOR[ND\d]/ || $f->primary_tag=~/_sn\d/ || $f->primary_tag=~/(-|_)sn?o?s?n?o?[A-WYZ]+[a-z]?-?\d/);

		my $hold=1;
		my @uga_ug;
		my @cu_ga;
		for my $c (0..$#{$constrains}){
			my ($sta,$sto,$mm,$query) = @{$$constrains[$c]};
			my @query = split // , lc($query);			
			my @stkseq = split // , lc(($stk->get_seq_by_id($f->seq_id))[0]->subseq($seedseq_pos_in_stk[$cs_pos_in_seed[$sta-1]-1]+1,$seedseq_pos_in_stk[$cs_pos_in_seed[$sto-1]-1]+1));
			# print @query;
			# print ' vs ';
			# print @stkseq;
			# print "\n";
			my $j=0;
			
			for (my $i=0; $i<=$#query; $i++){
				$j++ while $i+$j <= $#stkseq && $stkseq[$i+$j]!~/\w/;				

				push @uga_ug , $stkseq[$i+$j] if $c==0 && ($i==3 || $i==4);
				push @cu_ga , $stkseq[$i+$j] if $c==1 && ($i==0 || $i==1);					
				
				switch($query[$i]){
					case /[acgtu]/ {
						$mm-- unless $stkseq[$i+$j] eq $query[$i];
					}
					case "r" {
						$mm-- unless $stkseq[$i+$j]=~/[ag]/;
					}
					case "y" {
						$mm-- unless $stkseq[$i+$j]=~/[ctu]/;	
					}
					case "s" {
						$mm-- unless $stkseq[$i+$j]=~/[gc]/;	
					}
					case "w" {
						$mm-- unless $stkseq[$i+$j]=~/[atu]/;	
					}
					case "k" {
						$mm-- unless $stkseq[$i+$j]=~/[gtu]/;	
					}
					case "m" {
						$mm-- unless $stkseq[$i+$j]=~/[ac]/;	
					}
					case "b" {
						$mm-- unless $stkseq[$i+$j]=~/[cgtu]/;	
					}
					case "d" {
						$mm-- unless $stkseq[$i+$j]=~/[agtu]/;	
					}
					case "h" {
						$mm-- unless $stkseq[$i+$j]=~/[actu]/;	
					}
					case "v" {
						$mm-- unless $stkseq[$i+$j]=~/[acg]/;	
					}
					case "n" {
						$mm-- unless $stkseq[$i+$j]=~/[acgtu]/;	
					}
					else {}
				}					
			}
			# print $mm." ".join('',@query)." ".join('',@stkseq)."\n";
			if ($mm < 0){
				$hold=0;
				last;
			}
		}

		if ($hold && $f->score ne '.' && $f->score < 25 && ($f->primary_tag=~/_sno[A-Z]/ || $f->primary_tag=~/_Afu/ || $f->primary_tag=~/_SNOR[ND\d]/ || $f->primary_tag=~/_sn\d/ || $f->primary_tag=~/(-|_)sn?o?s?n?o?[A-WYZ]+[a-z]?-?\d/)){
			# print ''.join('',@uga_ug)." ".join('',@cu_ga)."\n";
			my $bpmm=0;
			switch ($uga_ug[0]){
				case "a" {
					$bpmm++ unless $cu_ga[1]=~/[tu]/;
				}
				case "c" {
					$bpmm++ unless $cu_ga[1] eq 'g';
				}
				case "g" {
					$bpmm++ unless $cu_ga[1]=~/[cu]/;
				}
				case /[tu]/ {
					$bpmm++ unless $cu_ga[1]=~/[au]/;
				}
			}
			switch ($uga_ug[1]){
				case "a" {
					$bpmm++ unless $cu_ga[0]=~/[tu]/;
				}
				case "c" {
					$bpmm++ unless $cu_ga[0] eq 'g';
				}
				case "g" {
					$bpmm++ unless $cu_ga[0]=~/[cu]/;
				}
				case /[tu]/ {
					$bpmm++ unless $cu_ga[0] eq 'a';
				}
			}			
			$hold=0 if $bpmm > 1;
			$hold=0 if $uga_ug[0]!~/[ut]/ && $cu_ga[1]!~/[ut]/;
		}		

		unless ($hold){
			# print $features->{$k}->seq_id."\n";
			delete $features->{$k};
			$write = 1;
			$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
			push @update , $f->seq_id.' '.$f->primary_tag.' P';			
		}		
	}

	return ($stk , $features, \@update , $write);
}

sub copy_filter {
	my ($self, $stk, $features, $copynumber) = @_;

	my $c=0;
	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';	
	my @update;
	my $write;

	if (scalar keys %{$features} > $copynumber){
		my @featureKeys = reverse sort { max($features->{$a}->score,($features->{$a}->get_tag_values('origscore'))[0]) <=> max($features->{$b}->score,($features->{$b}->get_tag_values('origscore'))[0]) } keys %{$features};
		for ($copynumber + 1 .. scalar keys %{$features}){
			my $f = $features->{$featureKeys[$_]};
			delete $features->{$featureKeys[$_]};
			$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
			push @update , $f->seq_id.' '.$f->primary_tag.' C';
		}			
	}

	return ($stk , $features, \@update , $write);
}

sub overlap_filter {
	my ($self, $stk, $features, $gffdb) = @_;

	my $c=0;
	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';
	my $type = $features->{(keys %{$features})[0]}->primary_tag;	
	my @update;
	my $write;	

	for my $i (keys %{$features}){
		my $f = $features->{$i};		
		my $existingFeatures = $gffdb->get_all_overlapping_features($f);
		my $exscore = -999999;
		my @rmfeatures;
		for my $f2 (@{$existingFeatures}){
			next if $f2->primary_tag=~/rRNA/ && $type=~/rRNA/;
			$exscore = max($exscore,  $f2->source=~/blast/ ?  $f2->score : max($f2->score,($f2->get_tag_values('origscore'))[0]) );
			push @rmfeatures , $f2;			 
		}	
		if ($exscore < ($f->source=~/blast/ ? $f->score : max($f->score,($f->get_tag_values('origscore'))[0])) ){
			push @update , $_->seq_id.' '.$_->primary_tag.' O' for @rmfeatures; 
		} else {			
			delete $features->{$i};
			$stk->remove_seq($stk->get_seq_by_id($f->seq_id));
			push @update , $f->seq_id.' '.$f->primary_tag.' O';
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
