package Bio::Gorap::Update;

use Archive::Extract;
use File::Spec::Functions;
use Bio::TreeIO;
use Bio::Gorap::Functions::STK;
use Bio::Gorap::Functions::CM;
use Bio::Gorap::DB::Taxonomy;
use File::Path qw(make_path rmtree);
use List::Util qw(min max);
use Try::Tiny;
use Symbol qw(gensym);
use IPC::Open3;
use Config::IniFiles;
use Switch;

sub dl_rfam {
	my ($self,$parameter,$taxdb) = @_;

	my $wgetpath = catdir($ENV{GORAP},'gorap','data','rfam_trees.tar.gz');
	my $url = 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed_tree.tar.gz';
	system("wget -O $wgetpath $url");
	my $extractpath = catdir($ENV{GORAP},'gorap','data','rfam_trees');
	make_path($extractpath);
	(Archive::Extract->new( archive => $wgetpath ))->extract( to => $extractpath );
	unlink $wgetpath;

	print "0% parsed\n";
	my @trees = glob(catfile($extractpath,'*'));
	my $percent = 1;
	my $step = ($#trees+1)/100;
	my $c=0;
	for (reverse @trees){
		$c++;
		if ($c > $percent * $step){
			print "$percent% parsed\n";
			$percent++;
		}
		open F,'<'.$_ or die $!;
		while(<F>){
			my @m = $_=~/[^_]+\/[^\[]+\[\d+\]/g;
			for (@m){
				my $acc = (split /\//,$_)[0];
				$_=~/\[(\d+)\]/;
				my $tmptaxid = $1;
				my $taxid;
				$taxid = $taxdb->accToTaxid->{$acc} if exists $taxdb->accToTaxid->{$acc};
				$taxid = $taxdb->getIDfromName($tmptaxid) unless $taxid;
				$taxid = $taxdb->getIDfromAccession($acc) unless $taxid;
				if ($taxid){
					$taxdb->accToTaxid->{$acc} = $taxid;
					my $name = $taxdb->getNameFromID($taxid);
					if ($name){
						$taxdb->taxidToName->{$taxid} = $name;
						$taxdb->nameToTaxid->{$name} = $taxid;
						$taxdb->accToName->{$acc} = $name;
					} else {
						$taxdb->accToName->{$acc} = $acc;
					}
				}
			}
		}
		close F;
	}
	rmtree($extractpath);
	print "100% parsed\n";

	open MAP , '>'.catfile($ENV{GORAP},'gorap','data','accSciTax.txt');
	print MAP $_."\t".$taxdb->accToName->{$_}."\t".$taxdb->accToTaxid->{$_}."\n" for(keys %{$taxdb->accToName});
	close MAP;

	$wgetpath=catdir($ENV{GORAP},'gorap','data','rfam.cm.gz');
	$url = 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz';
	system("wget --unlink -O $wgetpath $url");
	$extractpath = catfile($ENV{GORAP},'gorap','data','rfam.cm');
	(Archive::Extract->new( archive => $wgetpath ))->extract( to => $extractpath );
	unlink $wgetpath;

	open CM , '<'.$extractpath or die $!;
	my @tmp;
	my ($rf, $rna);
	while(<CM>){
		push @tmp, $_;
		if ($_=~/^INFERNAL/){
			next unless $rf;
			rmtree($_) for glob catdir($ENV{GORAP},'gorap','data','rfam',$rf.'*');
			$outpath=catdir($ENV{GORAP},'gorap','data','rfam',$rf.'_'.$rna);
			make_path($outpath);
			open CMOUT , '>'.catfile($outpath, $rf.'_'.$rna.'.cm') or die $!;
			pop @tmp;
			print CMOUT $_ for @tmp;
			close CMOUT;
			@tmp=();
			push @tmp, $_;
		} elsif ($_=~/NAME\s+(.+)/){
			$rna = $1;
		} elsif ($_=~/ACC\s+(.+)/){
			$rf = $1;
		}
	}
	rmtree($_) for glob catdir($ENV{GORAP},'gorap','data','rfam',$rf.'*');
	$outpath=catdir($ENV{GORAP},'gorap','data','rfam',$rf.'_'.$rna);
	make_path($outpath);
	open CMOUT , '>'.catfile($outpath, $rf.'_'.$rna.'.cm') or die $!;
	print CMOUT $_ for @tmp;
	close CMOUT;
	unlink $extractpath;

	$wgetpath = catdir($ENV{GORAP},'gorap','data','rfam.stk.gz');
	$url = 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz';
	system("wget -O $wgetpath $url");
	$extractpath = catfile($ENV{GORAP},'gorap','data','rfam.stk');
	(Archive::Extract->new( archive => $wgetpath ))->extract( to => $extractpath );
	unlink $wgetpath;

	print "0% parsed\n";
	$percent = 1;
	$step = $#trees / 10;
	$c=0;
	@prevlines=();
	my $map={};
	my $ss={};
	my $taxMapC={};
	open STK, '<'.$extractpath or die $!;
	while(<STK>){
		chomp $_;
		$_ =~ s/^\s+|\s+$//g;
		if($_=~/^\/\//){
			$c++;
			if ($c > $percent * $step){
				print "".($percent*10)."% parsed\n";
				$percent++;
			}

			my $maxlen = 0;
			for (keys %$map){
				my @line = split(/\s+/,join('',@{$map->{$_}}));
				$maxlen = max(length($line[1]),$maxlen);
			}
			$maxlen++;

			my $rfnr = (split(/\s+/,$prevlines[2]))[2];
			my $rna = (split(/\s+/,$prevlines[3]))[2];

			my $outpath = catdir($ENV{GORAP},'gorap','data','rfam',$rfnr.'_'.$rna);
			make_path($outpath);

			open STKOUT , '>'.catfile($outpath,$rfnr.'_'.$rna.'.stk') or die $!;
			print STKOUT $_."\n" for @prevlines;
			open FAOUT , '>'.catfile($outpath,$rfnr.'_'.$rna.'.fa') or die $!;
			for (keys %$map){
				my @line = split /\s+/ , join('',@{$map->{$_}});
				printf STKOUT ("%-".$maxlen."s %9s\n", $line[1] , $line[2]);
				$line[2] =~ s/\W//g;
				$line[2] = uc $line[2];
				print FAOUT '>'.$line[0]." ".$line[1]."\n".$line[2]."\n";
			}
			close FAOUT;
			my @line = split /\s+/ , join('',@{$ss->{'SS_cons'}});
			printf STKOUT ("%-".$maxlen."s %9s\n", $line[0].' '.$line[1], $line[2]);
			if (exists $ss->{'RF'}){
				@line = split(/\s+/,join('',@{$ss->{'RF'}}));
				printf STKOUT ("%-".$maxlen."s %9s\n", $line[0].' '.$line[1], $line[2]);
				print STKOUT "//\n";
				close STKOUT;
			} else {
				print STKOUT "//\n";
				close STKOUT;
				local *READER;
				my $pid = open3(gensym, \*READER, File::Spec->devnull , "cmbuild -F --v1p0 -O ".catfile($outpath,$rfnr.'_'.$rna.'.stk')." ".catfile($outpath,$rfnr.'_'.$rna.'.cm')." ".catfile($outpath,$rfnr.'_'.$rna.'.stk'));
				waitpid($pid, 0);
			}

			@prevlines=();
			$ss={};
			$map={};
			$taxMapC={};
		} elsif ($_=~/^#/ || $_=~/^\s*$/) {
			if ($_=~/^#=GC\s+(\S+)/){
				if (exists $ss->{$1}){
					push @{$ss->{$1}} , (split /\s+/ , $_)[2];
				} else {
					push @{$ss->{$1}} , $_;
				}
			} else {
				push @prevlines, $_;
			}
		} else {
			my @line = split(/\s+/,$_);
			if (exists $map->{$line[0]}){
				push @{$map->{$line[0]}} , $line[1];
			}else{
				my $acc = (split /\// , $line[0])[0];
				my $taxid;

				if (exists $taxdb->accToTaxid->{$acc}){
					$taxid = $taxdb->accToTaxid->{$acc};
				} else {
					$taxid = $taxdb->getIDfromAccession($acc) unless $taxid;
					$taxdb->accToTaxid->{$acc} = $taxid if $taxid;
				}
				if ($taxid){
					$taxMapC->{$taxid}++;
					my $name;
					if (exists $taxdb->accToName->{$acc}){
						$name = $taxdb->accToName->{$acc};
					} else {
						$name = $taxdb->getNameFromID($taxid);
					}
					if ($name){
						$taxdb->taxidToName->{$taxid} = $name;
						$taxdb->nameToTaxid->{$name} = $taxid;
						$taxdb->accToName->{$acc} = $name;
					} else {
						$taxdb->accToName->{$acc} = $acc;
					}
					push @{$map->{$line[0]}} , $taxid.'.'.$taxMapC->{$taxid}.' '.$taxdb->accToName->{$acc}.'.'.$taxMapC->{$taxid}.' '.$line[1];
				}
			}
		}
	}
	close STK;
	unlink $extractpath;
	print "100% parsed\n";

	open MAP , '>'.catfile($ENV{GORAP},'gorap','data','accSciTax.txt');
	print MAP $_."\t".$taxdb->accToName->{$_}."\t".$taxdb->accToTaxid->{$_}."\n" for(keys %{$taxdb->accToName});
	close MAP;
}

sub dl_silva {
	my ($self,$parameter,$taxdb) = @_;

	my $wgetpath = catfile($ENV{GORAP},'gorap','data','silvaList');
	my $url = 'ftp://ftp.arb-silva.de/living_tree/';
	my $version = `curl -s $url | awk \'{print \$NF}\' | sort -V | tail -1`;
	chomp $version;
	my $dir = `curl -s ftp://ftp.arb-silva.de/living_tree/$version/ | grep ^d | grep _SSU | awk \'{print \$NF}\'`;
	chomp $dir;

	$wgetpath = catfile($ENV{GORAP},'gorap','data','taxonomy','silva.newick');
	system("wget ftp://ftp.arb-silva.de/living_tree/$version/$dir/${dir}_tree.newick -O $wgetpath");

	print "0% parsed\n";
	my $tree = (Bio::TreeIO->new(-format => 'newick', -file => catfile($ENV{GORAP},'gorap','data','taxonomy','silva.newick') , -verbose => -1))->next_tree();
	my $percent = 1;
	my $step = $tree->number_nodes / 100;
	my $c=0;
	for my $node ( $tree->get_nodes ) {	
		$c++;
		if ($c > $percent * $step){
			print "$percent% parsed\n";
			$percent++;
		}

		next unless defined $node->id;
		$node->id(undef) unless $node->id;

		my $taxid;
		my ($name,$acc) = split /__+/,$node->id;
		for(($acc,$name)){
			next unless $_;
			$taxid = $taxdb->accToTaxid->{$_};
			$taxid = $taxdb->nameToTaxid->{$_} unless $taxid;
			$taxid = $taxdb->getIDfromAccession($_) unless $taxid;
			$taxid = $taxdb->getIDfromName($_) unless $taxid;
			print $_."\n" unless $taxid;
		}
		$node->id($taxid ? $taxid : undef);
	}
	print "100% parsed\n";
	
	Bio::TreeIO->new(-format => 'newick', -file => '>'.catfile($ENV{GORAP},'gorap','data','taxonomy','silva.newick') , -verbose => -1)->write_tree($tree);

	open MAP , '>'.catfile($ENV{GORAP},'gorap','data','accSciTax.txt');
	print MAP $_."\t".$taxdb->accToName->{$_}."\t".$taxdb->accToTaxid->{$_}."\n" for keys %{$taxdb->accToName};
	close MAP;
}

sub dl_silva_ambiguous {
	my ($self,$parameter,$taxdb) = @_;

	my $wgetpath = catfile($ENV{GORAP},'gorap','data','silvaList');
	my $url = 'ftp://ftp.arb-silva.de/current/Exports/taxonomy/';
	system("curl -s $url > $wgetpath");
	my $urlnewick;
	my $urltaxmap;
	open F , '<'.$wgetpath or die $!;
	while(<F>){
		my @l = split /\s+/, $_;
		$urlnewick = $url.$l[-1] if $l[-1] =~ /^tax_slv_ssu_.+\.tre$/;
		$urlacc = $url.$l[-1] if $l[-1] =~ /^tax_slv_ssu_.+\.acc_taxid$/;
		$urlacctax = $url.$l[-1] if $l[-1] =~ /^taxmap_embl_ssu_ref_.+\.txt.gz$/;
	}
	close F;
	unlink $wgetpath;

	my %acc2silva;
	my %silva2acc;
	$wgetpath = catfile($ENV{GORAP},'gorap','data','taxonomy','silva.acc');
	system("wget -O $wgetpath $urlacc");
	open F, '<'.$wgetpath or die $!;
	while(<F>){
		my @l = split /\s+/,$_;
		$acc2silva{(split /\./,$l[0])[0]} = $l[1];
		$silva2acc{$l[1]} = (split /\./,$l[0])[0];
	}
	close F;
	unlink $wgetpath;

	my %acc2tax;
	$wgetpath = catfile($ENV{GORAP},'gorap','data','taxonomy','silva.accTax.gz');
	my $extractpath = catfile($ENV{GORAP},'gorap','data','taxonomy','silva.accTax');
	system("wget -O $wgetpath $urlacctax");
	(Archive::Extract->new( archive => $wgetpath ))->extract( to => $extractpath );
	open F, '<'.$extractpath or die $!;
	while(<F>){
		my @l = split /\s+/,$_;
		$acc2tax{$l[0]} = $l[-1];
	}
	close F;
	unlink $wgetpath;
	unlink $extractpath;

	$wgetpath = catfile($ENV{GORAP},'gorap','data','taxonomy','silva.newick');
	system("wget -O $wgetpath $urlnewick");

	##!!!!not unambiguous:
	my $tree = (Bio::TreeIO->new(-format => 'newick', -file => $wgetpath , -verbose => -1))->next_tree();
	print "0% parsed\n";
	my $percent = 1;
	my $step = $tree->number_nodes / 10;
	my $c=0;
	for my $node ( $tree->get_nodes ) {
		$c++;
		if ($c > $percent * $step){
			print "".($percent*10)."% parsed\n";
			$percent++;
		}

		next unless defined $node->id;
		$node->id($acc2tax{$silva2acc{$node->id}});
	}
	print "100% parsed\n";
	Bio::TreeIO->new(-format => 'newick', -file => '>'.catfile($ENV{GORAP},'gorap','data','taxonomy','silva.newick') , -verbose => -1)->write_tree($tree);

	open MAP, '>'.catfile($ENV{GORAP},'gorap','data','silvaNcbi.txt') or die $!;
	print MAP $_."\t".$acc2silva{$_}."\t".$acc2tax{$_}."\n" for keys %acc2tax;
	close MAP;
}

sub dl_ncbi {
	my ($self,$parameter,$taxdb) = @_;

	print "Downloading NCBI taxonomy\n";
	my $wgetpath = catfile($ENV{GORAP},'gorap','data','taxdump.tar.gz');
	my $url = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz';
	system("wget --unlink -O $wgetpath $url");
	my $extractpath = catdir($ENV{GORAP},'gorap','data','taxonomy');
	(Archive::Extract->new( archive => $wgetpath ))->extract( to => $extractpath );
	unlink $wgetpath;

	my $taxdb = Bio::Gorap::DB::Taxonomy->new(
		parameter => $parameter
	);
	$taxdb->reindex();

	return $taxdb;
}

sub create_cfgs {
	my ($self,$parameter,$taxdb) = @_;

	make_path(catdir($ENV{GORAP},'gorap','config'));
	print "0% configured\n";
	my @stkfiles = glob(catfile($ENV{GORAP},'gorap','data','rfam','*','*.stk'));
	my $percent = 1;
	my $step = ($#stkfiles + 1) / 10;
	my $c=0;
	for my $stkfile (@stkfiles){
		$c++;
		if ($c > $percent * $step){
			print "".($percent*10)."% configured\n";
			$percent++;
		}

		my $bitscore = Bio::Gorap::Functions::STK->get_min_score($stkfile);
		my $rf_rna = Bio::Gorap::Functions::STK->get_rf_rna($stkfile);
		my ($ss , $cs) = Bio::Gorap::Functions::STK->get_ss_cs_from_file($stkfile);
		my ($rf,@rna) = split(/_/,$rf_rna);
		my $kingdoms;
		my @userdescription;
		my $types = '';
		my ($cfg) = glob catfile($ENV{GORAP},'gorap','config',$rf.'*');
		if (-e $cfg){
			$parameter->set_cfg($cfg);
			$kingdoms = $parameter->cfg->kingdoms;
			for (@{$parameter->cfg->constrains}){
				push @userdescription, $$_[-1] unless $userdescription[-1] eq $$_[-1];
			}
			$types = $parameter->cfg->types;
		}
		unlink $_ for glob catfile($ENV{GORAP},'gorap','config',$rf.'*');

		open CFG, '>'.catfile($ENV{GORAP},'gorap','config',$rf_rna.'.cfg') or die $!;
		print CFG '[family]'."\n";
		print CFG 'id='.$rf."\n";
		print CFG 'name='.join('_',@rna)."\n";
		print CFG "\n";
		print CFG '[cmd]'."\n";
		if ($types){
			$types = 'RNaseP' if $rf_rna=~/RNaseP/;
			$types = 'else' if $rf_rna=~/RF00002/;
			print CFG $self->get_cmd($types);
		} else {
			$types = 'rRNA' if $rf_rna =~/_rRNA/ && $rf_rna!~/RF00002/;
			$types = 'tRNA' if $rf_rna =~/\d_tRNA/;
			$types = 'RNaseP' if $rf_rna=~/RNaseP/;
			$types = 'CRISPR' if $rf_rna=~/CRISPR/;
			print CFG $self->get_cmd($types);
		}
		print CFG "\n";
		print CFG '[thresholds]'."\n";
		print CFG "evalue=1e-3\n";
		print CFG "bitscore=$bitscore\n";
		print CFG "\n";

		print CFG "[query]\n";
		print CFG 'fasta='.catfile('gorap','data','rfam',$rf_rna,$rf_rna.'.fa')."\n";
		print CFG 'stk='.catfile('gorap','data','rfam',$rf_rna,$rf_rna.'.stk')."\n";
		print CFG 'cm='.catfile('gorap','data','rfam',$rf_rna,$rf_rna.'.cm')."\n";
		print CFG "pseudogenes=all\n";

		my $cdsnorna = 0;
		if ($types){
			if ($types=~/snRNA/){
				if ($rf_rna=~/(F|f)ungi/ || $rf_rna=~/_Afu/ || $rf_rna=~/yeast/){
					print CFG  "kingdom=fungi\n";
				} elsif ($types=~/splicing/){
					print CFG "kingdom=euk\nkingdom=fungi\n";
				} else {
					print CFG  "kingdom=euk\nkingdom=fungi\nkingdom=arc\n";
				}
				$cdsnorna = 1 if $types =~ /CD-box/;
			} elsif($types=~/CRISPR/){
				print CFG "kingdom=arc\nkingdom=bac\n";
			} elsif($types=~/microRNA/){
				print CFG "kingdom=euk\n";
			} else {
				$types = '';
			}
		}
		unless($types){
			if ($rf_rna=~/(B|b)act/){
				print CFG "kingdom=bac\n";
			} elsif ($rf_rna=~/(E|e)uka/ || $rf_rna=~/(P|p)lant/ || $rf_rna=~/(M|m)etazo/ || $rf_rna=~/(P|p)rotoz/){
				print CFG $rf_rna=~/LSU/ ? "kingdom=euk\nkingdom=fungi\n" : "kingdom=euk\n";
			} elsif ($rf_rna=~/(A|a)rch/){
				print CFG "kingdom=arc\n";
			} elsif ($rf_rna=~/(F|f)ungi/ || $rf_rna=~/_Afu/ || $rf_rna=~/yeast/){
				print CFG "kingdom=fungi\n";
			} elsif ($rf_rna=~/_U[0-9](1|2|atac|_|$)/){
				print CFG "kingdom=euk\nkingdom=fungi\n";
			} elsif ($rf_rna=~/CRISPR/){
				print CFG "kingdom=arc\nkingdom=bac\n";
			} elsif ($rf_rna=~/_mir/i){
				print CFG "kingdom=euk\n";
			} elsif ($rf_rna=~/_ceN/ || $rf_rna=~/_DdR/ || $rf_rna=~/CD\d+$/ || $rf_rna=~/_sno_/ || $rf_rna=~/_SNOR/ || $rf_rna=~/_sno[A-Z]/ || $rf_rna=~/_sn\d/ || $rf_rna =~/_sn?o?s?n?o?[A-WYZ]+[a-z]?\d/){
				#not -.+- -> e.g. v-snoRNA-1 is viral:
				$cdsnorna = 1;
				$cdsnorna = 0 if $rf_rna=~/_snopsi/ || $rf_rna=~/_SNOR[A\d]/ || (($rf_rna=~/_sn\d/ || $rf_rna=~/_sno[A-Z]/) && $cs!~/UGA.{0,12}$/i && $cs=~/A[^G]A.{2,8}$/i);
				print CFG "kingdom=euk\nkingdom=fungi\nkingdom=arc\n";
			} elsif ($rf_rna=~/_6S/){
				print CFG "kingdom=arc\nkingdom=bac\n";
			} else {
				if (scalar(keys %$kingdoms) == 0 ){
					my $seqio = Bio::SeqIO->new( -format => 'Fasta' , -file => catfile($ENV{GORAP},'gorap','data','rfam',$rf_rna,$rf_rna.'.fa'), -verbose => -1);
					while(my $seqobj = $seqio->next_seq()){
						my @lineage = @{$taxdb->getLineageNodes((split /\./ ,$seqobj->id)[0])};
						next if $#lineage < 1;
						next if $lineage[0]->id == 28384 || $lineage[0]->id == 12908;
						$kingdoms->{'virus'} = 1 if $lineage[0]->id==10239 || $lineage[0]->id==12884;
						$kingdoms->{'euk'} = 1 if $lineage[1]->id==2759;
						$kingdoms->{'arc'} = 1 if $lineage[1]->id==2157;
						$kingdoms->{'bac'} = 1 if $lineage[1]->id==2;
						$kingdoms->{'fungi'} = 1 if $lineage[3] && $lineage[3]->id==4751;
					}
					$kingdoms->{'fungi'} = 1 if exists $kingdoms->{'euk'};
				}
				if (scalar keys %$kingdoms > 0){
					print CFG 'kingdom='.$_."\n" for sort keys %$kingdoms;
				} else {
					print CFG "kingdom=arc\nkingdom=bac\nkingdom=euk\nkingdom=fungi\nkingdom=virus\n";
				}
			}
		}
		print CFG "\n";
		print CFG '[constrains]'."\n";
        my @ss = split // , $ss;
        my @cs = split // , $cs;
        $ss='';
        $cs='';
        for my $i (0..$#cs){
          if ($cs[$i]=~/[a-zA-Z]/){
            $cs.=$cs[$i];
            $ss.=$ss[$i];
          }
        }
        $ss=~tr/<>{}[]/()()()/;
        $ss=~s/[^()]/./g;

        print CFG "structure=$ss\n";
		print CFG "conserved=$cs\n";
		if ($#userdescription != -1){
			print CFG 'constrain='.$_."\n" for @userdescription;
			print STDERR ":ERROR: manually re-check constrains for $rf_rna in ".catfile($ENV{GORAP},'gorap','config',$rf_rna.'.cfg')."\n" if length($cs) != length($userdescription[0]);
		} elsif($cdsnorna){
			print ":INFO: consider to add box constrains for $rf_rna in ".catfile($ENV{GORAP},'gorap','config',$rf_rna.'.cfg')."\n";
		}
		close CFG;
	}
	print "100% configures\n";
}

sub get_cmd {
	my ($self, $rna_type) = @_;

	my $options;
    switch ($rna_type) {
		case /rRNA/	{
			$options = "tool=barrnap\n";
			$options .= "parameter=--threads $cpus\n";
			$options .= "parameter=--kingdom $kingdom\n";
			$options .= "parameter=$genome\n";
			# $options = "tool=rnammer\n";
			# $options .= "parameter=-m ssu,lsu,tsu\n";
			# $options .= 'parameter=-S $kingdom'."\n";
			# $options .= 'parameter=-gff $output'."\n";
			# $options .= 'parameter=$genome'."\n";
		}
		case /tRNA/ {
			$options = "tool=tRNAscan-SE\n";
			$options .= "parameter=-q\n";
			$options .= "parameter=-b\n";
			$options .= "parameter=-Q\n";
			$options .= 'parameter=-$kingdom'."\n";
			$options .= 'parameter=$genome'."\n";
		}
		case /RNaseP/ {
			$options = "tool=Bcheck\n";
			$options .= 'parameter=-$kingdom'."\n";
			$options .= 'parameter=--ss'."\n";
			$options .= 'parameter=-o $output'."\n";
			$options .= 'parameter=$genome'."\n";
		}
		case /CRISPR/ {
			$options = "tool=crt\n";
			$options .= "parameter=-screen 1\n";
			$options .= "parameter=-minNR 4\n";
			$options .= "parameter=-minRL 10\n";
			$options .= "parameter=-maxRL 100\n";
			$options .= "parameter=-minSL 10\n";
			$options .= "parameter=-maxSL 100\n";
			$options .= "parameter=-searchWL 8\n";
			$options .= 'parameter=$genome'."\n";
		}
		else {
			$options = "tool=infernal\n";
			$options .= "tool=blast\n";
		}
    }

	return $options;
}

1;
