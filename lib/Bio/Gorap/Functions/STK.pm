package Bio::Gorap::Functions::STK;

use Bio::AlignIO;
use Bio::SimpleAlign;
use POSIX;
use Switch;
use List::Util qw(min max);
use File::Spec::Functions;

sub get_min_score {
	#open given covariance model and extract minimum of gathering-, trusted- and noise cutoff
	my ($self,$stk) = @_;
	open STK , '<'.$stk or die $!;
	my $score = 999999;
	while(<STK>){
		chomp;
		if ($_=~/^#=GF\s+(GA|TC|NC)\s+(.+)/){
			$score = $2 < $score ? $2 : $score;
		}
	}
	close STK;
	return $score == 999999 ? 0 : $score;
}

sub get_rf_rna {
	my ($self,$stk) = @_;

	my $rna='';
	my $rf='';
	open STK , '<'.$stk or die $!;
	while(<STK>){
		chomp;
		if ($_=~/^#=GF\s+AC\s+(.+)/){
			$rf=$1;
		}
		if ($_=~/^#=GF\s+ID\s+(.+)/){
			$rna=$1;
		}
	}
	close STK;
	$rna=~s/(^[\s\._]+|[\s\._]+$)//g;
	$rna=~s/\s+/_/g;
	$rna=~s/\.+/\./g;
	$rna=~s/__+/_/g;
	$rna=~s/[^a-zA-Z0-9_\.]*//g;

	return $rf.'_'.$rna;
}

sub get_rna_types {
	my ($self,$stk) = @_;

	my $types = '';
	open STK , '<'.$stk or die $!;
	while(<STK>){
		chomp;
		if ($_=~/^#=GF\s+TP\s+(.+);/){
			$types = $1;
			$types =~ s/;\s+/:/g;
		}
	}
	close STK;

	return $types;
}

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
	my ($self, $stk, $features, $threshold, $nonTaxThreshold) = @_;

	my $c=0;
	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';
	my $write;
	my @update;

	$nonTaxThreshold = 999999 unless $nontaxthreshold; #else: gorap was startet with taxonomy - $threshold is taxonomy based, $nontaxthreshold comes from cfg
		
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
			if ($f->score < min($nonTaxThreshold,$threshold) ){
				delete $features->{$_};
				$write = 1;
				$stk->remove_seq($stk->get_seq_by_id($f->seq_id));
				push @update , $f->seq_id.' '.$f->primary_tag.' B';
			}
		}
	}
	
	return ($stk , $features, \@update , $write);
}

sub structure_filter {
	my ($self, $stk, $features, $rna_types, $minstructures) = @_;

	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';

	my @update;
	my $write;

	my ($ss , $cs) = $self->get_ss_cs_from_object($stk);
	$ss=~y/\{\<\[\}\>\]/\(\(\(\)\)\)/;
	my @ss = split // , $ss;
	my @cs = split // , $cs;

	my @open;
	my @close;
	my $x = 0;
	for (my $i=0 ; $i <= $#ss ; $i++){
		if ($ss[$i] eq '('){
			push @open , $i;
			$x = 0;
		}elsif($ss[$i] eq ')'){
			$#close = $#open;
			if ($close[-1*(1+$x)]){
				$x++ while $close[-1*(1+$x)];
				$close[-1*(1+$x)] = $i;
			} else {
				$close[-1*(1+$x)] = $i;
				$x++;
			}
		}
	}
	my %c2o = map {$close[$_] => $open[$_]} (0..$#open);
	# for (0..$#open){
	# 	say $open[$_].' '.$close[$_];
	# }
	my @stack;
	$#stack = 0;
	my $bracket='(';
	for (0..$#ss){
		next unless $ss[$_]=~/[\(\)]/;
		if ($ss[$_] eq $bracket){
			push @{$stack[-1]} , {$ss[$_] => $_};
		} else {
			$bracket = $bracket eq ')' ? '(' : ')';
			$#stack++;
			push @{$stack[-1]} , {$ss[$_] => $_};
		}
	}
	# for (@stack){
	# 	for (@$_){
	# 		($bracket) = keys %$_;
	# 		print $bracket.$_->{$bracket};
	# 	}
	# 	say '';
	# }
	my ($o,$c,@areas);
	for (my $i = 1; $i<=$#stack; $i+=2){
		my $j = $i - 1;
		$j-=2 while $#{$stack[$j]} == -1;
		#area open
		($bracket) = keys %{$stack[$i][0]};
		$o = $stack[$i][0]->{$bracket};
		while(my $h = shift @{$stack[$i]}){
			if ($#{$stack[$j]} == -1){
				#area close
				push @areas , ([$o,$c],[$c2o{$c},$c2o{$o}]);
				# say "$o-$c $c2o{$c}-$c2o{$o}";
				$j-=2 while $#{$stack[$j]} == -1;
				#area open
				($bracket) = keys %$h;
				$o = $h->{$bracket};
				pop @{$stack[$j]};
			} else {
				pop @{$stack[$j]};
				($bracket) = keys %$h;
				$c = $h->{$bracket};
			}
		}
		#area close
		push @areas , ([$o,$c],[$c2o{$c},$c2o{$o}]);
		# say "$o-$c $c2o{$c}-$c2o{$o}";
	}
	# for (sort {$$a[0] <=> $$b[0]} @areas){
	# 	say $$_[0].' '.$$_[1];
	# }

	return ($stk , $features, \@update , $write) if $#areas == -1;

	my ($annastart,$annastop,$ssPresentCount);
	for (sort {$$a[0] <=> $$b[0]} @areas){
		my ($sta, $sto) = ($$_[0],$$_[1]);

		$annastop = $sta if ! $annastop && $sta > $#ss/2 && $#areas > 0;
		$annastart = $sto+2 unless $annastop;

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

	my $hacasno = $rna_types=~/HACA-box/ ? 1 : 0;
	my @aca;
	my $i=$#cs;
	for(reverse @cs){
		if ($_=~/[a-zA-Z]/){
			unshift @aca, $i;
			$i--;
		}
		last if $#aca==12;
	}
	my @anna;
	my ($hp,$i,$bracket) = (0,0,0);
	for(@ss){ #count hairpins and get ananna locus between both hairpins
		$i++;
		unless($_=~/(\(|\))/){
			push @anna, $i-1 if $hp==1 && $bracket==0;
		} else {
			$bracket = $_ eq "(" ? $bracket+1 : $bracket-1;
			$hp++ if $bracket==1 && $_ eq "(";
		}
	}
	$hacasno = 0 unless $hp == 2;

	$minstructures = ($#areas +1)/2 unless $minstructures;
	for (keys %{$features}){
		my $f = $features->{$_};

		next if max($f->score,($f->get_tag_values('origscore'))[0]) >= 40;

		if ($#areas == 1){
			if (! exists $ssPresentCount->{$f->seq_id} || $ssPresentCount->{$f->seq_id} < 1){
				delete $features->{$_};
				$write = 1;
				$stk->remove_seq($stk->get_seq_by_id($f->seq_id));
				push @update , $f->seq_id.' '.$f->primary_tag.' S';
			}
		} else {
			if (! exists $ssPresentCount->{$f->seq_id} || $ssPresentCount->{$f->seq_id} < $minstructures){
				delete $features->{$_};
				$write = 1;
				$stk->remove_seq($stk->get_seq_by_id($f->seq_id));
				push @update , $f->seq_id.' '.$f->primary_tag.' S';
			} elsif ($hacasno){
				my $s = ($stk->get_seq_by_id($f->seq_id))->subseq($anna[0],$anna[-1]);
				my $s2 = ($stk->get_seq_by_id($f->seq_id))->subseq($aca[-11],$aca[-1]);
				$s=~s/[\W_]//g;
				$s2=~s/[\W_]//g;
				unless ( $s=~/a.{1,3}a.{1,4}a/i && $s2=~/a[^g]a.{2,8}$/i){
					delete $features->{$_};
					$write = 1;
					$stk->remove_seq($stk->get_seq_by_id($f->seq_id));
					push @update , $f->seq_id.' '.$f->primary_tag.' S';
				}
			}
		}
	}

	return ($stk , $features, \@update , $write);
}

sub sequence_filter {
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
			if ($consC/$allCons < 0.7 || ($f->primary_tag=~/_mir/i && $consC/$allCons < 0.8)){
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
	my ($self, $stk, $features, $constrains, $cs, $seedstk, $rna_types) = @_;
	$features = {map { $c++ => $_ } @{$features}} if ref($features) eq 'ARRAY';
	my $cdsno = $rna_types=~/CD-box/ ? 1 : 0;
	my @update;
	my $write;

	my @cs = split //,$cs;
	my ($seedss , $seedcs) = $self->get_ss_cs_from_file($seedstk);
	my @seedcs = split // , $seedcs;

	my %seedmap;
	my $j=0;
	for my $i (0..$#cs){
		$j++ while $i+$j<$#seedcs && $seedcs[$i+$j]=~/\W/;
		$seedmap{$i} = $i+$j;
	}

	my $seedaln = (Bio::AlignIO->new(-format  => 'stockholm', -file => $seedstk, -verbose => -1 ))->next_aln;
	for my $c (0..$#{$constrains}){
		my ($sta,$sto,$mm,$query) = @{$$constrains[$c]};
		$sta--;
		$sto--;
		# print "contrain: $query\n";
		# print @cs[$sta..$sto]; print "\t";
		# print @seedcs[$seedmap{$sta}..$seedmap{$sto}]; print "\n";

		my $hold=1;
		for $seedseq ($seedaln->each_seq) {
			my $subseq = $seedseq->subseq($seedmap{$sta}+1,$seedmap{$sto}+1);
			unless ($subseq=~/^\W/ && $subseq=~/\W$/){
				my @stkseq = split // , ($stk->get_seq_by_id($seedseq->id))->seq;
				my @seedseq = split //, $seedseq->seq;
				my %stkmap;
				my $j=0;
				for my $i (0..$#seedseq){
					next if $seedseq[$i]=~/\W/;
					$j++ while $i+$j<$#stkseq && $stkseq[$i+$j]=~/\W/;
					$stkmap{$i} = $i+$j;
				}
				# print $seedseq->seq."\n";
				# print $seedmap{$sta}."\t".$seedmap{$sto}."\t".$subseq."\n";
				# print @stkseq; print "\n";
				# print $stkmap{$seedmap{$sta}}."\t".$stkmap{$seedmap{$sto}}."\t";
				# print @stkseq[$stkmap{$seedmap{$sta}}..$stkmap{$seedmap{$sto}}]; print "\n";
				${$$constrains[$c]}[0] = $stkmap{$seedmap{$sta}};
				${$$constrains[$c]}[1] = $stkmap{$seedmap{$sto}};
				last;
			}		
		}
		${$$constrains[$c]}[3]='' unless $hold;
	}

	for my $k (keys %{$features}){
		my $f = $features->{$k};
		next if $f->score ne '.' && $f->score > 30 && $cdsno;

		my $stkseq = $stk->get_seq_by_id($f->seq_id);
		my @uga_ug;
		my @cu_ga;
		my $hold=1;
		for my $c (0..$#{$constrains}){
			my ($sta,$sto,$mm,$query) = @{$$constrains[$c]};
			next unless $query; #backmap unavailable
			#print $query."\t".$f->seq_id."\t".$stkseq->subseq($sta+1,$sto+1)."\n";

			$query = lc $query;
			$query=~s/\W//g;

			my ($seq) = $stk->get_seq_by_id($f->seq_id);
			my $seq = lc $stkseq->subseq($sta+1,$sto+1);
			$seq=~s/\W//g;
			my @seq = split //,$seq;
			
			my ($costs, @alnmap) = &gotoh($query,$seq);
			if ($costs*-1 > $mm){
				$hold=0;
				last;
			}

			@uga_ug = (defined $alnmap[3] && $seq[$alnmap[3]] ? $seq[$alnmap[3]] : '', defined $seq[4] && $seq[$alnmap[4]] ? $seq[$alnmap[4]] : '') if $c==0;
			@cu_ga = (defined $alnmap[0] && $seq[$alnmap[0]] ? $seq[$alnmap[0]] : '', defined $seq[1] && $seq[$alnmap[1]] ? $seq[$alnmap[1]] : '') if $c==1;
		}

		if ($hold && $cdsno && $#uga_ug>=1 && $#cu_ga>=1){
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
					$bpmm++ unless $cu_ga[1]=~/[cut]/;
				}
				case /[tu]/ {
					$bpmm++ unless $cu_ga[1]=~/[aut]/;
				}
				else {
					$bpmm++
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
					$bpmm++ unless $cu_ga[0]=~/[cut]/;
				}
				case /[tu]/ {
					$bpmm++ unless $cu_ga[0] eq 'a';
				}
				else {
					$bpmm++
				}
			}
			if ($hold){
				$hold=0 if $bpmm > 1;
				$hold=0 if $uga_ug[0]!~/[ut]/ && $cu_ga[1]!~/[ut]/;
			}
		}

		unless ($hold){
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
	my ($self,$stk,$tmpfile) = @_;

	my $out='';
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

sub cost {
	my ($q, $s) = @_;

	switch($q){
		case /[acgtu]/ {
			return $q eq $s ? 0 : -1;
		}
		case "r" {
			return $s=~/[ag]/ ? 0 : -1;
		}
		case "y" {
			return $s=~/[ctu]/ ? 0 : -1;
		}
		case "s" {
			return $s=~/[gc]/ ? 0 : -1;
		}
		case "w" {
			return $s=~/[atu]/ ? 0 : -1;
		}
		case "k" {
			return $s=~/[gtu]/ ? 0 : -1;
		}
		case "m" {
			return $s=~/[ac]/ ? 0 : -1;
		}
		case "b" {
			return $s=~/[cgtu]/ ? 0 : -1;
		}
		case "d" {
			return $s=~/[agtu]/ ? 0 : -1;
		}
		case "h" {
			return $s=~/[actu]/ ? 0 : -1;
		}
		case "v" {
			return $s=~/[acg]/ ? 0 : -1;
		}
		case "n" {
			return $s=~/[acgtu]/ ? 0 : -1;
		}
		else {}
	}
}

sub gotoh (){
	my ($q,$s) = @_;

	my @s = split //, $s;
	my @q = split //, $q;

	my $go=-1;
	my $ge=-1;

	my @M;
	my @H;
	my @V;

	for(my $i=0; $i<=scalar(@q); $i++){
		$M[$i][0]=$ge*($i-1)+$go;
		# $M[$i][0]=0;
		$H[$i][0]=-inf;
	}
	for(my $j=0; $j<=scalar(@s); $j++){
		$M[0][$j]=$ge*($j-1)+$go;
		# $M[0][$j]=0;
		$V[0][$j]=-inf;
	}
	$M[0][0]=0;

	my @maxposi;
	my @maxposj;
	my $max=-inf;
	my $c;
	for(my $i=1; $i<=scalar(@q); $i++){
		for(my $j=1; $j<=scalar(@s); $j++){
			$H[$i][$j]=max($H[$i][$j-1]+$ge,$M[$i][$j-1]+$go);
			$V[$i][$j]=max($V[$i-1][$j]+$ge,$M[$i-1][$j]+$go);
			$M[$i][$j]=max(max($M[$i-1][$j-1]+&cost($q[$i-1],$s[$j-1]),$H[$i][$j]),$V[$i][$j]);
		}
	}

	my $i=scalar(@q);
	my $j=scalar(@s);
	my $alS='';
	my $alQ='';
	my @map;
	$#map=$#q;
	while($i > 0 || $j > 0){
		if (($i > 0 && $j > 0) && $M[$i][$j]==$M[$i-1][$j-1]+&cost($q[$i-1],$s[$j-1])){
			$map[$i-1]=$j-1;
			$alS=$s[$j-1].$alS;
			$alQ=$q[$i-1].$alQ;
			$i--;
			$j--;
		} elsif ($j > 0 && $M[$i][$j]==$M[$i][$j-1]+$go) {
			$map[$i]=$j-1;
			$alS=$s[$j-1].$alS;
			$alQ="-".$alQ;
			$j--;
		} elsif ($i > 0 && $M[$i][$j]==$M[$i-1][$j]+$go){
			$map[$i-1]=$j;
			$alS="-".$alS;
			$alQ=$q[$i-1].$alQ;
			$i--;
		} elsif ($j > 0 && $M[$i][$j]==$H[$i][$j-1]+$ge){
			$map[$i]=$j-1;
			$alS=$s[$j-1].$alS;
			$alQ="-".$alQ;
			$j--;
		} else {
			$map[$i-1]=$j;
			$alS="-".$alS;
			$alQ=$q[$i-1].$alQ;
			$i--;
		}

	}
	# print "$alS\n";
	# print "$alQ\n";
	# print @map; print "\n";

	return ($M[scalar(@q)][scalar(@s)],@map);
}

1;
