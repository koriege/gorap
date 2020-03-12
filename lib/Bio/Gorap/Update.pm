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
use File::Basename;

sub dl_rfam {
	my ($self,$parameter,$taxdb) = @_;

	my $wgetpath = catdir($ENV{GORAP},'db','data','rfam_trees.tar.gz');
	my $url = 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed_tree.tar.gz';
	system("wget -O $wgetpath $url");
	my $extractpath = catdir($ENV{GORAP},'db','data','rfam_trees');
	make_path($extractpath);
	(Archive::Extract->new( archive => $wgetpath ))->extract( to => $extractpath );
	unlink $wgetpath;

	print "0% parsed\n";
	my @trees = glob(catfile($extractpath,'*'));
	my $percent = 1;
	my $step = ($#trees+1)/100;
	my $c=0;
	for my $treefile (@trees){
		$c++;
		if ($c > $percent * $step){
			print "$percent% parsed\n";
			$percent++;
		}

		open F,"<$treefile" or die $!;
		while(<F>){
			chomp $_;
			@matches = $_=~/([^)(,]+\/[^)(,]+)/;
			for (@matches){
				my ($acc,$name) = split /\//;
				$acc =~ s/^[^_]*_//;
				$name =~ s/^[^_]*_//;
				$name =~ s/\[([^\]]+).*//;
				$name = $taxdb->flatname($name);
				$taxid = $1;
				unless($taxid=~/^\d+$/){ # sometimes rfam has newick format issues
					undef $taxid;
					$taxid = $taxdb->accToTaxid->{$acc} if $acc && exists $taxdb->accToTaxid->{$acc};
					$taxid = $taxdb->getIDfromAccession($acc) if $acc && ! $taxid;
					$taxid = $taxdb->getIDfromName($name) if $name && ! $taxid;
				}
				if ($taxid){
					$taxdb->accToTaxid->{$acc} = $taxid;
					$name = $taxdb->getNameFromID($taxid);
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

		# my $tree = (Bio::TreeIO->new(-format => 'newick', -file => $treefile , -verbose => -1))->next_tree(); # sometimes, rf00001 e.g., parser failes
		# for my $node ( $tree->get_nodes ) {	
		# 	next unless defined $node->is_Leaf;
		# 	$node->id=~/^[^_]+_([^\/]+)\/[^_]+_(\S+)/;
		# 	my $acc = $1;
		# 	my $name = taxdb->flatname($2);
		# 	my $taxid = $node->bootstrap;
		# }
	}
	rmtree($extractpath);
	print "100% parsed\n";

	open MAP , '>'.catfile($ENV{GORAP},'db','data','accSciTax.rfam');
	print MAP $_."\t".$taxdb->accToName->{$_}."\t".$taxdb->accToTaxid->{$_}."\n" for(keys %{$taxdb->accToName});
	close MAP;

	$wgetpath=catdir($ENV{GORAP},'db','data','rfam.cm.gz');
	$url = 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz';
	system("wget -O $wgetpath $url");
	$extractpath = catfile($ENV{GORAP},'db','data','rfam.cm');
	(Archive::Extract->new( archive => $wgetpath ))->extract( to => $extractpath );
	unlink $wgetpath;

	open CM , '<'.$extractpath or die $!;
	my @tmp;
	my ($rf, $rna);
	while(<CM>){
		push @tmp, $_;
		if ($_=~/^INFERNAL/){
			next unless $rf;
			rmtree($_) for glob catdir($ENV{GORAP},'db','data','rfam',$rf.'*');
			$outpath=catdir($ENV{GORAP},'db','data','rfam',$rf.'_'.$rna);
			make_path($outpath);
			open CMOUT , '>'.catfile($outpath, $rf.'_'.$rna.'.cm') or die $!;
			pop @tmp;
			print CMOUT $_ for @tmp;
			close CMOUT;
			@tmp=();
			push @tmp, $_;
		} elsif ($_=~/NAME\s+(.+)/){
			$rna = $taxdb->flatname($1);
		} elsif ($_=~/ACC\s+(.+)/){
			$rf = $1;
		}
	}
	rmtree($_) for glob catdir($ENV{GORAP},'db','data','rfam',$rf.'*');
	$outpath=catdir($ENV{GORAP},'db','data','rfam',$rf.'_'.$rna);
	make_path($outpath);
	open CMOUT , '>'.catfile($outpath, $rf.'_'.$rna.'.cm') or die $!;
	print CMOUT $_ for @tmp;
	close CMOUT;
	unlink $extractpath;


	$wgetpath = catdir($ENV{GORAP},'db','data','rfam.stk.gz');
	$url = 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz';
	system("wget -O $wgetpath $url");
	$extractpath = catfile($ENV{GORAP},'db','data','rfam.stk');
	(Archive::Extract->new( archive => $wgetpath ))->extract( to => $extractpath );
	unlink $wgetpath;

	print "0% parsed\n";
	$percent = 1;
	$step = ($#trees+1)/100;
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
				print "$percent% parsed\n";
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
			$rna = $taxdb->flatname($rna);

			my $outpath = catdir($ENV{GORAP},'db','data','rfam',$rfnr.'_'.$rna);
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
			$uid={};
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
					if ($acc=~/^URS0[^_]+_(\d+)$/){
						$taxid = $1;
					} else {
						$taxid = $taxdb->getIDfromAccession($acc);
					}
				}
				if ($taxid){
					$taxdb->accToTaxid->{$acc} = $taxid;
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
					$name = $taxdb->accToName->{$acc};
					$uid->{$name}++;
					push @{$map->{$line[0]}} , $taxid.'.'.$uid->{$name}.' '.$name.'.'.$uid->{$name}.' '.$line[1];
				}
			}
		}
	}
	close STK;

	unlink $extractpath;
	print "100% parsed\n";

	open MAP , '>'.catfile($ENV{GORAP},'db','data','accSciTax.rfam');
	print MAP $_."\t".$taxdb->accToName->{$_}."\t".$taxdb->accToTaxid->{$_}."\n" for(keys %{$taxdb->accToName});
	close MAP;
}

sub dl_silva {
	my ($self,$parameter,$taxdb) = @_;

	my $wgetpath = catfile($ENV{GORAP},'db','data','taxonomy','silva.newick.gz');
	my $url = 'ftp://ftp.arb-silva.de/current/Exports/taxonomy/*.tre.gz';
	system("wget --glob on $url -O $wgetpath");
	my $extractpath = catdir($ENV{GORAP},'db','data','taxonomy');
	(Archive::Extract->new( archive => $wgetpath ))->extract( to => $extractpath );
	unlink $wgetpath;
	
	return;

	$wgetpath = catfile($ENV{GORAP},'db','data','silvaList');
	$url = 'ftp://ftp.arb-silva.de/living_tree/';
	my $version = `curl -s $url | awk \'{print \$NF}\' | sort -V | tail -1`;
	chomp $version;
	my $file = `curl -s ftp://ftp.arb-silva.de/living_tree/$version/ | grep -F SSU_tree.newick | awk \'{print \$NF}\'`;
	chomp $file;

	$wgetpath = catfile($ENV{GORAP},'db','data','taxonomy','silva.newick');
	system("wget ftp://ftp.arb-silva.de/living_tree/$version/$file -O $wgetpath");

	print "0% parsed\n";
	my $tree = (Bio::TreeIO->new(-format => 'newick', -file => catfile($ENV{GORAP},'db','data','taxonomy','silva.newick') , -verbose => -1))->next_tree();
	my $percent = 1;
	my $step = $tree->number_nodes / 100;
	my $c=0;
	open ACC , '>'.catfile($ENV{GORAP},'db','data','accTax.silva');
	open SCI , '>'.catfile($ENV{GORAP},'db','data','sciTax.silva');
	for my $node ( $tree->get_nodes ) {	
		$c++;
		if ($c > $percent * $step){
			print "$percent% parsed\n";
			$percent++;
		}

		next unless defined $node->id;
		$node->id(undef) unless $node->id;

		my $taxid;
		my $nodeid = $node->id;
		$nodeid=~s/\W//g;
		my ($acc,$name) = split /_+/,$nodeid;
		if(! $name){
			$name = $acc;
			$acc = '';
		}
		$name=~s/\s+/_/g;
		if ($acc){
			$taxid = $taxdb->accToTaxid->{$acc};
			$taxid = $taxdb->getIDfromAccession($acc) unless $taxid;
			print ACC $acc."\t".$taxid."\n" if $taxid;
		} elsif ($name){
			$taxid = $taxdb->nameToTaxid->{$name};
			$taxid = $taxdb->getIDfromName($name) unless $taxid;
			print SCI $name."\t".$taxid."\n" if $taxid;
		}

		$node->id($taxid ? $taxid : undef);
	}
	close ACC;
	close SCI;
	print "100% parsed\n";

	Bio::TreeIO->new(-format => 'newick', -file => '>'.catfile($ENV{GORAP},'db','data','taxonomy','silva.newick') , -verbose => -1)->write_tree($tree);
}

sub dl_ncbi {
	my ($self,$parameter) = @_;

	print "Downloading NCBI taxonomy\n";
	my $wgetpath = catfile($ENV{GORAP},'db','data','taxdump.tar.gz');
	my $url = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz';
	system("wget -O $wgetpath $url");
	my $extractpath = catdir($ENV{GORAP},'db','data','taxonomy');
	(Archive::Extract->new( archive => $wgetpath ))->extract( to => $extractpath );
	unlink $wgetpath;

	my $taxdb = Bio::Gorap::DB::Taxonomy->new(
		parameter => $parameter
	);

	return $taxdb;
}

sub create_cfgs {
	my ($self,$parameter,$taxdb) = @_;

	make_path(catdir($ENV{GORAP},'db','config'));
	my %cfgs = map {my $f = basename($_); $f=~/(RF\d+)/; $1 => $_} glob catfile($ENV{GORAP},'db','config','*.cfg');
	print "0% configured\n";
	my @stkfiles = glob(catfile($ENV{GORAP},'db','data','rfam','*','*.stk'));
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
		my $cfg = $cfgs{$rf};
		if ($cfg){
			$parameter->set_cfg($cfg);
			$kingdoms = $parameter->cfg->kingdoms;
			for (@{$parameter->cfg->constrains}){
				push @userdescription, $$_[-1] unless $userdescription[-1] eq $$_[-1];
			}
			$types = $parameter->cfg->types;
			unlink $cfg;
			delete $cfgs{$rf};
		}

		open CFG, '>'.catfile($ENV{GORAP},'db','config',$rf_rna.'.cfg') or die $!;
		print CFG '[family]'."\n";
		print CFG 'id='.$rf."\n";
		print CFG 'name='.join('_',@rna)."\n";
		print CFG "\n";
		print CFG '[cmd]'."\n";
		if ($types){
			$types = 'RNaseP' if $rf_rna=~/RNaseP/;
			print CFG $self->get_cmd($types);
		} else {
			$types = 'rRNA' if $rf_rna =~/_rRNA/;
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
					my $seqio = Bio::SeqIO->new( -format => 'Fasta' , -file => catfile($ENV{GORAP},'db','data','rfam',$rf_rna,$rf_rna.'.fa'), -verbose => -1);
					while(my $seqobj = $seqio->next_seq()){
						my @lineage = @{$taxdb->getLineageNodes((split /\./ ,$seqobj->id)[0])};
						next if $#lineage < 1;
						next if $lineage[0] == 28384 || $lineage[0] == 12908;
						$kingdoms->{'virus'} = 1 if $lineage[0] == 10239 || $lineage[0] == 12884;
						$kingdoms->{'euk'} = 1 if $lineage[1] == 2759;
						$kingdoms->{'arc'} = 1 if $lineage[1] == 2157;
						$kingdoms->{'bac'} = 1 if $lineage[1] == 2;
						$kingdoms->{'fungi'} = 1 if $lineage[3] && $lineage[3] == 4751;
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
			print STDERR ":ERROR: manually re-check constrains for $rf_rna in ".catfile($ENV{GORAP},'db','config',$rf_rna.'.cfg')."\n" if length($cs) != length($userdescription[0]) && ! (length($userdescription[0])-length($cs) == 1 && $userdescription[0]=~/\|$/);
		} elsif($cdsnorna){
			print ":INFO: consider to add box constrains for $rf_rna in ".catfile($ENV{GORAP},'db','config',$rf_rna.'.cfg')."\n";
		}
		close CFG;
	}
	print "100% configured\n";
	unlink $cfgs{$_} for keys %cfgs;
}

sub get_cmd {
	my ($self, $rna_type) = @_;

	my $options;
	switch ($rna_type) {
		case /rRNA/	{
			$options = "tool=barrnap\n";
			$options .= 'parameter=--threads $cpus';
			$options .= ' --kingdom $kingdom';
			$options .= ' $genome'."\n";
			# $options = "tool=rnammer\n";
			# $options .= "parameter=-m ssu,lsu,tsu\n";
			# $options .= 'parameter=-S $kingdom'."\n";
			# $options .= 'parameter=-gff $output'."\n";
			# $options .= 'parameter=$genome'."\n";
		}
		case /tRNA/ {
			$options = "tool=tRNAscan-SE\n";
			$options .= "parameter=-L";
			$options .= " -q";
			$options .= " -Q";
			$options .= ' -$kingdom';
			$options .= ' $genome'."\n";
		}
		case /RNaseP/ {
			$options = "tool=Bcheck\n";
			$options .= 'parameter=-$kingdom';
			$options .= ' --ss';
			$options .= ' -o $output';
			$options .= ' $genome'."\n";
		}
		case /CRISPR/ {
			$options = "tool=crt\n";
			$options .= "parameter=-screen 1";
			$options .= " -minNR 4";
			$options .= " -minRL 10";
			$options .= " -maxRL 100";
			$options .= " -minSL 10";
			$options .= " -maxSL 100";
			$options .= " -searchWL 8";
			$options .= ' $genome'."\n";
		}
		else {
			$options = "tool=infernal\n";
			$options .= "tool=blast\n";
		}
	}

	return $options;
}

1;
