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

sub dl_rfam {
	my ($self,$parameter,$taxdb) = @_;

	my $wgetpath=catdir($ENV{GORAP},'data','rfam.cm.gz');
	my $url = 'ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.1_1.gz';
	system("wget --unlink -O $wgetpath $url");
	my $extractpath = catfile($ENV{GORAP},'data','rfam.cm');
	(Archive::Extract->new( archive => $wgetpath ))->extract( to => $extractpath );
	unlink $wgetpath;

	my $entryC = 0;
	open CM , '<'.$extractpath or die $!;
	my @prevlines=('','');
	my $cmstart=0;
	while(<CM>){	
		my $line = $_;		
		if ($line=~/^ACC\s+(.+)/){
			$entryC++;
			my $rf = $1;
			if ($prevlines[0] ne ''){
				close CMOUT;				
			}				
			my $rna = (split(/\s+/,$prevlines[1]))[1];

			rmtree($_) for glob catdir($ENV{GORAP},'data','rfam',$rf.'*');
			$outpath=catdir($ENV{GORAP},'data','rfam',$rf.'_'.$rna);	
			make_path($outpath);
			open CMOUT , '>'.catfile($outpath, $rf.'_'.$rna.'.cm') or die $!;
			print CMOUT $prevlines[0];
			print CMOUT $prevlines[1];
			print CMOUT $line;			 	
		} else {
			$cmstart=0 if $line=~/^INFERNAL/;
			$cmstart < 2 ? $prevlines[$cmstart++]=$line : print CMOUT $line; 		
		}
	}
	close CM;
	close CMOUT;
	unlink $extractpath;

	unless ($taxdb){
		$taxdb = Bio::Gorap::DB::Taxonomy->new(
			parameter => $parameter
		);
	}

	$wgetpath = catdir($ENV{GORAP},'data','rfam.stk.gz');
	$url = 'ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz';
	system("wget --unlink -O $wgetpath $url");
	$extractpath = catfile($ENV{GORAP},'data','rfam.stk');
	(Archive::Extract->new( archive => $wgetpath ))->extract( to => $extractpath );
	unlink $wgetpath;
	
	my $percent = 1;
	my $step = $entryC / 10;
	my $c=0;
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
			
			my $outpath = catdir($ENV{GORAP},'data','rfam',$rfnr.'_'.$rna);
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
				if (exists $taxdb->rfamToTaxid->{$acc}){
					my $tax = $taxdb->rfamToTaxid->{$acc};					
					$taxMapC->{$tax}++;
					push @{$map->{$line[0]}} , $tax.'.'.$taxMapC->{$tax}.' '.$taxdb->rfamToName->{$acc}.'.'.$taxMapC->{$tax}.' '.$line[1];
				} else {
					my $tax = $taxdb->getIDfromAccession($acc);
					if ($tax){
						$taxdb->rfamToTaxid->{$acc} = $tax;
						my $name = $taxdb->getNameFromID($tax);
						$taxdb->rfamToName->{$acc} = $name ? $name : $acc;
						$taxMapC->{$tax}++;
						push @{$map->{$line[0]}} , $tax.'.'.$taxMapC->{$tax}.' '.$taxdb->rfamToName->{$acc}.'.'.$taxMapC->{$tax}.' '.$line[1];
					} else {
						print "Skipping Rfam accession $acc - NCBI does not contain any related entry\n";
					}
				}
			}
		}
	}
	close STK;
	unlink $extractpath;

	open MAP , '>'.catfile($ENV{GORAP},'data','accSciTax.txt');
	print MAP $_."\t".$taxdb->rfamToName->{$_}."\t".$taxdb->rfamToTaxid->{$_}."\n" for(keys %{$taxdb->rfamToName}); 
	close MAP;

	#create all rfam queries related gorap configuration files
	#&create_cfgs($self,$parameter,$taxdb);
}

sub dl_silva {
	my ($self,$parameter,$taxdb) = @_;	

	my $wgetpath = catfile($ENV{GORAP},'data','silvaList');
	my $url = 'ftp://ftp.arb-silva.de/living_tree/';
	system("wget -q --unlink -O $wgetpath $url");		
	my $version = 0;
	my $link;
	open REL , '<'.$wgetpath or die $!;	
	while(<REL>){
		if ($_=~/>(.*LTP_release_)(\d+)\/</){
			if ($2 > $version){
				$link = $1.$2;
				$version=$2;
			}							
		}		
	}
	close REL;
	unlink $wgetpath;
		
	$wgetpath = catfile($ENV{GORAP},'data','taxonomy','silva16s.newick');
	for ( 'ftp://ftp.arb-silva.de/living_tree/'.$link.'/LTPs'.$version.'_SSU_tree.newick' ,
		'ftp://ftp.arb-silva.de/living_tree/'.$link.'/LTPs'.$version.'_SSU.tree.newick' ,
		'ftp://ftp.arb-silva.de/living_tree/'.$link.'/LTPs'.$version.'_SSU.ntree'){
		local *READER;
		my $pid = open3(gensym, \*READER, File::Spec->devnull , "wget -T 10 -q --unlink -O $wgetpath $_");
		if ($pid){
			waitpid($pid, 0);
			last;
		}
	}
	
	print "Parsing Silva tree - go and get some coffee ;)\n";

	my $tree = (Bio::TreeIO->new(-format => 'newick', -file => catfile($ENV{GORAP},'data','taxonomy','silva16s.newick') , -verbose => -1))->next_tree();			
	unless ($taxdb){
		$taxdb = Bio::Gorap::DB::Taxonomy->new(
			parameter => $parameter
		);
	}

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
				 
		my @tmp = split(/__+/,$node->id);
		if ($#tmp==0){			
			if (exists $taxdb->silvaToTaxid->{$tmp[0]}){
				$node->id($taxdb->silvaToTaxid->{$tmp[0]});		
			} else {
				my $id = $taxdb->getIDfromName((split(/_+/,$tmp[0]))[0]);			
				if ($id){
					$taxdb->silvaToTaxid->{$tmp[0]} = $id;
					$node->id($id);
				} else {
					$node->id(undef);
				}
			}
		} elsif ($#tmp == -1){
				$node->id(undef);			
		} else {		
			my $acc = $tmp[0]=~/^\D{0,3}\d+$/ ? $tmp[0] : $tmp[1];
			if (exists $taxdb->silvaToTaxid->{$acc}){
				$node->id($taxdb->silvaToTaxid->{$acc});		
			} else {
				my $id = $taxdb->getIDfromAccession($acc);
				if ($id){
					$taxdb->silvaToTaxid->{$acc} = $id;
					$node->id($id);
				} else {
					$id = $taxdb->getIDfromName($acc);
					if ($id){
						$taxdb->silvaToTaxid->{$acc} = $id;
						$node->id($id);
					} else {
						$node->id(undef);
					}
				}
			}	
		}			
	}
	
	Bio::TreeIO->new(-format => 'newick', -file => '>'.catfile($ENV{GORAP},'data','taxonomy','silva16s.newick') , -verbose => -1)->write_tree($tree);	

	open MAP, '>'.catfile($ENV{GORAP},'data','silvaNcbi.txt') or die $!;	
	print MAP $_."\t".$taxdb->silvaToTaxid->{$_}."\n" for keys %{$taxdb->silvaToTaxid};	
	close MAP;	
}

sub dl_ncbi {
	my ($self) = @_;

	print "Downloading NCBI taxonomy\n";
	my $wgetpath = catfile($ENV{GORAP},'data','taxdump.tar.gz');
	my $url = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz';
	system("wget --unlink -O $wgetpath $url");
	my $extractpath = catdir($ENV{GORAP},'data','taxonomy');
	(Archive::Extract->new( archive => $wgetpath ))->extract( to => $extractpath );	
	unlink $wgetpath;
}

sub create_cfgs {
	my ($self,$parameter,$taxdb) = @_;

	unless ($taxdb){
		$taxdb = Bio::Gorap::DB::Taxonomy->new(
			parameter => $parameter
		);
	}

	my @cmfiles = glob(catfile($ENV{GORAP},'data','rfam','*','*.cm'));
	my $percent = 1;
	my $step = ($#cmfiles + 1) / 10;
	my $c=0;
	for my $cmfile (@cmfiles){
		$c++;
		if ($c > $percent * $step){
			print "".($percent*10)."% configured\n";
			$percent++;
		}		

		my $bitscore = Bio::Gorap::Functions::CM->get_min_score($cmfile);
		my $rf_rna = Bio::Gorap::Functions::CM->get_rf_rna($cmfile);

		my @userdescription;
		my ($cfg) = glob catfile($ENV{GORAP},'parameter','config',(split(/_/,$rf_rna))[0].'*');
		if ($cfg){
			open F , '<'.$cfg or die $!;
			@userdescription = <F>;
			while($userdescription[-1]=~/^\s*$/){
				pop @userdescription;
			}
			close F;
			@userdescription = () if $userdescription[-1]=~/^#/;
		}

		unlink $_ for glob catfile($ENV{GORAP},'parameter','config',(split(/_/,$rf_rna))[0].'*');

		make_path(catdir($ENV{GORAP},'parameter','config'));
		open CFG, '>'.catfile($ENV{GORAP},'parameter','config',$rf_rna.'.cfg') or die $!;
		print CFG "#$rf_rna\n";
		print CFG "#system call of tool with linewise parameters\n";		
		if ($rf_rna =~/_rRNA/ && $rf_rna!~/RF00002/){
			print CFG "rnammer\n";
			print CFG "-m ssu,lsu,tsu\n";
			print CFG '-S $kingdom'."\n";
			print CFG '-gff $output'."\n";
			print CFG '$genome'."\n";
		} elsif($rf_rna=~/RNaseP/){
			print CFG "Bcheck\n";
			print CFG '-$kingdom'."\n";			
			print CFG '--ss'."\n";
			print CFG '-o $output'."\n";
			print CFG '$genome'."\n";		
		} elsif($rf_rna=~/\d_tRNA/){
			print CFG "tRNAscan-SE\n";
			print CFG "-q\n";
			print CFG "-b\n";
			print CFG "-Q\n";	
			print CFG '-$kingdom'."\n";
			print CFG '$genome'."\n";
		} elsif($rf_rna=~/CRISPR/){
			print CFG "crt\n";
			print CFG "-screen 1\n";
			print CFG "-minNR 4\n";
			print CFG "-minRL 10\n";
			print CFG "-maxRL 100\n";
			print CFG "-minSL 10\n";
			print CFG "-maxSL 100\n";
			print CFG "-searchWL 8\n";
			print CFG '$genome'."\n";
		} else {
			print CFG "infernal\n";
			print CFG "blast\n";
		}
		print CFG "#query fasta file\n";	
		print CFG catfile('data','rfam',$rf_rna,$rf_rna.'.fa')."\n";
		print CFG "#query stockholm file\n";
		print CFG catfile('data','rfam',$rf_rna,$rf_rna.'.stk')."\n";
		print CFG "#query covariance model\n";
		print CFG catfile('data','rfam',$rf_rna,$rf_rna.'.cm')."\n";
		print CFG "#e-value\n";
		print CFG "1e-3\n";
		print CFG "#bitscore\n";
		print CFG "$bitscore\n";
		print CFG "#maxcopies (pseudogenes): all/<INT>\n";
		print CFG "all\n";
		print CFG "#kingdoms: <bac/arc/euk/fungi/virus> comma seperated\n";
		

		if ($rf_rna=~/(B|b)act/){
			print CFG "bac\n";
		} elsif ($rf_rna=~/(E|e)uka/ || $rf_rna=~/(P|p)lant/ || $rf_rna=~/(M|m)etazo/ || $rf_rna=~/(P|p)rotoz/){
			print CFG $rf_rna=~/LSU/ ? "euk,fungi\n" : "euk\n";
		} elsif ($rf_rna=~/(A|a)rch/){
			print CFG "arc\n";
		} elsif ($rf_rna=~/(F|f)ungi/){
			print CFG "fungi\n";
		} elsif ($rf_rna=~/_U3$/){
			print CFG "euk\n";
		} elsif ($rf_rna=~/CRISPR/){
			print CFG "arc,bac\n";	
		} elsif ($rf_rna=~/_mir/ || $rf_rna=~/_MIR/){
			print CFG "euk\n";	
		} elsif ($rf_rna=~/_sn?o?(r|z|u|y)d?\d/i){
			print CFG "euk,fungi,arc\n";	
		} elsif ($rf_rna=~/_6S/){
			print CFG "arc,bac\n";	
		} else {				
			my $kingdoms;
			my $seqio = Bio::SeqIO->new( -format => 'Fasta' , -file => substr($cmfile,0,-2).'fa', -verbose => -1);
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
			if (scalar keys %$kingdoms > 0){
				print CFG join(',',sort keys %$kingdoms)."\n";	
			} else {
				print CFG "arc,bac,euk,fungi,virus\n";	
			}
		}
		my ($ss , $cs) = Bio::Gorap::Functions::STK->get_ss_cs_from_file(substr($cmfile,0,-2).'stk');
		my @ss = split // , $ss;
		my @cs = split // , $cs;
		my @csp;
		my @ssp; 
		for my $i (0..$#cs){
			if ( $cs[$i]=~/[a-zA-Z]/ || $ss[$i]=~/[<\[\(\{>\]\}\)]/){
				if ($#ssp>-1){
					if ( ($ss[$i]=~/[<\[\(\{]/ && $ssp[$#ssp]=~/[>\]\}\)]/) || ($ss[$i]=~/[>\]\}\)]/ && $ssp[$#ssp]=~/[<\[\(\{]/) ){
						push @csp , '|';
						push @ssp , '|';
					} elsif ($ss[$i]!~/[<\[\(\{>\]\}\)]/ && $ssp[$#ssp]=~/[<\[\(\{>\]\}\)]/){
						push @csp , '|';
						push @ssp , '|';
					} elsif ($ss[$i]=~/[<\[\(\{>\]\}\)]/ && $ssp[$#ssp]!~/[<\[\(\{>\]\}\)]/){
						push @csp , '|';
						push @ssp , '|';		
					}
				}	
				push @csp , $cs[$i];
				push @ssp , $ss[$i];
			}
		}
		if ($#userdescription > -1){
			my $usercs = $userdescription[-3];
			$usercs =~ s/\|//g;
			if ('#'.$cs eq $usercs ){
				print CFG $userdescription[-3];
				print CFG $userdescription[-2];
				print CFG $userdescription[-1];
			} else {
				print CFG '#'.join('',@csp)."\n";
				print CFG '#'.join('',@ssp)."\n\n";
			}
		} else {
			print CFG '#'.join('',@csp)."\n";
			print CFG '#'.join('',@ssp)."\n\n";	
		}
		close CFG;
	}
}

1;