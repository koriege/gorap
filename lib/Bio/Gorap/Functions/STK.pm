package Bio::Gorap::Functions::STK;

use Bio::AlignIO;
use Bio::SimpleAlign;
use POSIX;
use Switch;

sub score_filter {
	my ($self, $nofilter, $stk, $features, $threshold, $nonTaxThreshold) = @_;

	my $c=0;
	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';
	my $write;
	my @update;	
	my $type = $features->{(keys %{$features})[0]}->type;
	
	if ( ! $nofilter && ($type=~/_Afu/ || $type=~/_SNOR/ || $type=~/(-|_)sn?o?s?n?o?[A-WYZ]+[a-z]?-?\d/)){
		for (keys %{$features}){		
			my $f = $features->{$_};
			next if $f->score eq '.';

			if ($f->score < 10){
				delete $features->{$_};
				$write = 1;
				$stk->remove_seq($stk->get_seq_by_id($f->seq_id));
				push @update , $f->seq_id.' '.$f->primary_tag.' B';
			}
		}

		return ($stk , $features, \@update , $write);
	}	

	if ($nonTaxThreshold){ #gorap was startet with taxonomy - $threshold is taxonomy based

		for (keys %{$features}){
			my $f = $features->{$_};
			next if $f->score eq '.';

			my @id = split /\./ , $f->seq_id;
			my ($abbr,$orig,$copy) = ($id[0] , join('.',@id[1..($#id-1)]) , $id[-1]);
			
			if (($f->get_tag_values('source'))[0] =~ /infernal/ || ($f->get_tag_values('source'))[0] =~ /blast/){
				#print $f->seq_id." ".$f->score." ".ceil(($f->get_tag_values('origscore'))[0]).' '.$threshold."\n";
				if (($f->get_tag_values('origscore'))[0] < $threshold || $f->score < 10){
					delete $features->{$_};
					$write = 1;
					$stk->remove_seq($stk->get_seq_by_id($f->seq_id)); 
					push @update , $f->seq_id.' '.$f->primary_tag.' B';
				}
			} else {
				if ($f->score < 10 || $f->score < $nonTaxThreshold && $f->score < $threshold){
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

sub sequence_filter {
	my ($self, $stk, $features) = @_;

	my $c=0;
	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';
	my $type = $features->{(keys %{$features})[0]}->type;
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
	my $j=0;
	for (my $i=0; $i<=$#seedseq; $i++){
		$c = $i;
		$i++ while $seedseq[$i]=~/\W/;
		$j++ while $seedseq[$i] ne $stkseedseq[$j];

		push @seedseq_pos_in_stk , ($j)x($i-$c+1);		
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
	# for (0..$#{$constrains}){
	# 	my ($sta,$sto,$mm,$query) = @{$$constrains[$_]};
	# 	print $sta." ".$sto." ".join('',@cs[$sta-1..$sto-1])."\n";
	# 	print $seedseq_pos_in_stk[$cs_pos_in_seed[$sta-1]]." ".$seedseq_pos_in_stk[$cs_pos_in_seed[$sto-1]].' '.join('',@stkseedseq[$seedseq_pos_in_stk[$cs_pos_in_seed[$sta-1]]..$seedseq_pos_in_stk[$cs_pos_in_seed[$sto-1]]])."\n";
	# }
	
	for my $k (keys %{$features}){				
		my $f = $features->{$k};
		# print ''.join('',@stkseedseq)."\n";
		# print ''.($stk->get_seq_by_id($f->seq_id))[0]->seq."\n";
		# for (0..$#{$constrains}){
		# 	my ($sta,$sto,$mm,$query) = @{$$constrains[$_]};
		# 	print $sta." ".$sto." ".join('',@cs[$sta-1..$sto-1])."\n";
		# 	print $seedseq_pos_in_stk[$cs_pos_in_seed[$sta-1]]." ".$seedseq_pos_in_stk[$cs_pos_in_seed[$sto-1]].' '.lc(($stk->get_seq_by_id($f->seq_id))[0]->subseq($seedseq_pos_in_stk[$cs_pos_in_seed[$sta-1]]+1,$seedseq_pos_in_stk[$cs_pos_in_seed[$sto-1]]+1))."\n";
		# }

		my $hold=1;
		my @uga_ug;
		my @cu_ga;
		for my $c (0..$#{$constrains}){
			my ($sta,$sto,$mm,$query) = @{$$constrains[$c]};
			my @query = split // , lc($query);			
			my @stkseq = split // , lc(($stk->get_seq_by_id($f->seq_id))[0]->subseq($seedseq_pos_in_stk[$cs_pos_in_seed[$sta-1]]+1,$seedseq_pos_in_stk[$cs_pos_in_seed[$sto-1]]+1));
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
		if ($hold && $f->score ne '.' && $f->score < 35 && ($f->type=~/_Afu/ || $f->type=~/_SNOR/ || $f->type=~/(-|_)sn?o?s?n?o?[A-WYZ]+[a-z]?-?\d/)){
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
					$bpmm++ unless $cu_ga[1]=~/[gu]/;
				}
				case /[tu]/ {
					$bpmm++ unless $cu_ga[0]=~/[au]/;
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
					$bpmm++ unless $cu_ga[0]=~/[gu]/;
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
