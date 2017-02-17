#! /usr/bin/env bash
trap "trap - SIGTERM && kill -PIPE -- -$$" INT TERM SIGINT SIGTERM EXIT 
force=$([[ $1 == "force" ]] && echo 1 || echo '')

if [[ ! $GORAP ]]; then
	echo 'Setup $GORAP environment variable first.'
	exit 1
fi
mkdir -p $GORAP

echo '' > $GORAP/install.log
os='linux'
if [[ $OSTYPE == darwin* ]]; then 
	os='mac'
fi
bit=$(uname -m | awk '{if ($0 == "i686"){print "32"}else{print "64"}}')

progress (){
	x=0
	while [ 1 ]; do	
		x=$((x+1))
		if [ "$x" -eq "2" ]
			then
				x=0
				echo -en "\r/"
			else
				echo -en "\r\\"
		fi
		sleep 1
	done
}

getpaths (){
	case $toolid in
	'RAXML')
		tool='standard-RAxML-8.2.9'
		dlpath='https://github.com/stamatak/standard-RAxML/archive/v8.2.9.tar.gz'
	;;
	'MAFFT')
		tool='mafft-7.305-with-extensions'
		dlpath="http://mafft.cbrc.jp/alignment/software/$tool-src.tgz"
	;;
	'NEWICKTOPDF')
		tool='newicktopdf'
		dlpath="ftp://pbil.univ-lyon1.fr/pub/mol_phylogeny/njplot/$tool"
	;;
	'INFERNAL')
		tool='infernal-1.1.2'
		dlpath="http://eddylab.org/software/infernal/$tool.tar.gz"
	;;
	'INFERNAL1')
		tool='infernal-1.0'
		dlpath="http://eddylab.org/software/infernal/$tool.tar.gz"
	;;
	'GLIBC')
		tool='glibc-2.24'
		dlpath="https://ftp.gnu.org/gnu/libc/$tool.tar.gz"
	;;
	'RNABOB')
		tool='rnabob-2.2.1'
		dlpath="http://eddylab.org/software/rnabob/$tool.tar.gz"
	;;
	'BLAST')
		if [[ $bit -eq 32 ]]; then
			if [[ $os == 'mac' ]]; then
				dlpath='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/ncbi-blast-2.4.0+-universal-macosx.tar.gz'
			else
				dlpath='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-ia32-linux.tar.gz'
			fi
		else
			if [[ $os == 'mac' ]]; then
				dlpath=$(wget --spider --force-html -r -nd -np -l 1 -A *x64-macosx.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ 2>&1 | grep '^--.*tar.gz' | awk '{print $NF}')
			else
				dlpath=$(wget --spider --force-html -r -nd -np -l 1 -A *x64-linux.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ 2>&1 | grep '^--.*tar.gz' | awk '{print $NF}')
			fi
		fi
		tool=$(echo $dlpath | sed 's/%2B/+/g' | grep -Eo 'ncbi-blast-[0-9].*.tar.gz' | sed 's/.tar.gz$//')
	;;
	'TRNASCAN')
		tool='tRNAscan-SE-1.3.1'
		dlpath="http://www.rna.uni-jena.de/supplements/gorap/$tool-src.tar.gz"
	;;
	'HMMER')
		tool='hmmer-2.3.2'
		dlpath="http://eddylab.org/software/hmmer/2.3.2/$tool.tar.gz"
	;;
	'RNAMMER')
		tool='rnammer-1.2'
		dlpath="http://www.rna.uni-jena.de/supplements/gorap/$tool.tar.gz"		
	;;
	'BCHECK')
		tool='bcheck-0.6'
		# dlpath="http://rna.tbi.univie.ac.at/bcheck/Bcheck.tgz"
		dlpath="http://www.rna.uni-jena.de/supplements/gorap/$tool.tar.gz"
	;;
	'CRT')
		tool='crt-1.2'
		dlpath='http://www.room220.com/crt/CRT1.2-CLI.jar.zip'
	;;
	'SAMTOOLS')
		tool='samtools-0.1.19'
		dlpath='https://github.com/samtools/samtools/archive/0.1.19.tar.gz'
	;;
	'ZLIB')
		tool='zlib-1.2.11'
		dlpath="https://github.com/madler/zlib/archive/v1.2.11.tar.gz"
	;;
	'NCURSES')
		tool='ncurses-6.0'
		dlpath="https://ftp.gnu.org/pub/gnu/ncurses/$tool.tar.gz"
	;;
	*)
		exit 1
	;;
	esac
}

download (){
	echo "Installing $tool" | tee -a $GORAP/install.log
	progress &
	pid=$!
	if [[ $dlpath =~ (tgz|tar.gz)$ ]]; then
		wget -q -T 10 $dlpath -O $GORAP/$tool.tar.gz && tar -xzf $GORAP/$tool.tar.gz -C $GORAP && rm -f $GORAP/$tool.tar.gz
	else
		mkdir -p $GORAP/$tool/bin
		if [[ $dlpath =~ zip$ ]]; then
			wget -q -T 10 $dlpath -O $GORAP/$tool.zip && unzip -q -o -d $GORAP/$tool/bin $GORAP/$tool.zip && rm -f $GORAP/$tool.zip
		else 
			wget -q -T 10 $dlpath -O $GORAP/$tool/bin/$tool
		fi
	fi
	if [[ $? -gt 0 ]]; then
		kill $pid &> /dev/null
		wait $pid &> /dev/null
		echo
		echo "Downloading $tool failed - see $GORAP/install.log for details"
		exit 1
	fi
	sed -i'' "s/$toolid=.*/$toolid='$tool'/g" recompile_tool.sh
	bash recompile_tool.sh $tool &>> $GORAP/install.log
	excode=$?
	kill $pid &> /dev/null
	wait $pid &> /dev/null
	echo
	if [[ $excode -gt 0 ]]; then
		echo "Installing $tool failed - see $GORAP/install.log for details"
		exit 1
	fi
}

if [[ $bit -eq 32 ]] || [[ $os == 'mac' ]]; then
	echo 'Unexpected operating system detected. Software compilation necessary.'
	if [[ $bit -eq 32 ]]; then
		echo 'NOTE: Phylogeny reconstruction will not be available on your 32 bit OS'
	fi
	read -p 'Do you want to continue [y|n] ' in
	if [[ $in != 'y' ]] && [[ $in != 'Y' ]]; then
		echo 'Installation aborted...'
		exit
	fi
fi

if [[ $bit -eq 64 ]]; then
	toolid=RAXML
	getpaths
	ex=$([[ $(which raxml) ]] && echo $(which raxml) || echo $GORAP/$tool/bin/raxml)
	$ex -h &> /dev/null
	if [[ $? -gt 0 ]] || [[ $force ]]; then
		download
	fi

	toolid=MAFFT
	getpaths
	ex=$([[ $(which mafft) ]] && echo $(which mafft) || echo $GORAP/$tool/bin/mafft)
	if [[ $($ex -h 2>&1 | grep MAFFT | wc -l) -eq 0 ]] || [[ $force ]]; then
		download
	fi

	toolid=NEWICKTOPDF
	getpaths
	ex=$([[ $(which newicktopdf) ]] && echo $(which newicktopdf) || echo $GORAP/$tool/bin/newicktopdf)
	$ex -h &> /dev/null
	if [[ $? -gt 0 ]] || [[ $force ]]; then
		download
	fi
fi

if [[ $os != 'mac' ]]; then
	toolid=GLIBC
	getpaths
	if [[ $(which ldd) ]]; then
		gv=$(ldd --version | head -n 1 | awk '{split($NF,a,"."); print a[1]a[2]}')
	else 
		gv=0
	fi
	if [[ $gv -lt 214 ]]; then
		download
	fi
fi

toolid=INFERNAL
getpaths
ex=$([[ $(which cmsearch) ]] && echo $(which cmsearch) || echo $GORAP/$tool/bin/cmsearch)
if [[ ! $($ex -h 2>&1 | grep INFERNAL | awk '$3>=1.1{print 1}') ]] || [[ $force ]]; then
	download
fi

toolid=INFERNAL1
getpaths
ex=$([[ $(which cmsearch) ]] && echo $(which cmsearch) || echo $GORAP/$tool/bin/cmsearch)
if [[ ! $($ex -h 2>&1 | grep INFERNAL | awk '$3==1.0{print 1}') ]] || [[ $force ]]; then
	download
fi

toolid=RNABOB
getpaths
ex=$([[ $(which rnabob) ]] && echo $(which rnabob) || echo $GORAP/$tool/bin/rnabob)
$ex -h &> /dev/null
if [[ $? -gt 0 ]] || [[ $force ]]; then
	download
fi

toolid=BLAST
getpaths
ex=$([[ $(which blastn) ]] && echo $(which blastn) || echo $GORAP/$tool/bin/blastn)
if [[ $force ]] || [[ ! $($ex -h 2>&1 | grep DESCRIPTION -A 1 | tail -n 1 | awk '{split($NF,a,"."); if(a[2]>=2 && $NF~/\+$/){print 1}}') ]]; then
	download
fi

toolid=TRNASCAN
getpaths
ex=$([[ $(which tRNAscan-SE) ]] && echo $(which tRNAscan-SE) || echo $GORAP/$tool/bin/tRNAscan-SE)
export PERL5LIB=$GORAP/$tool:$PERL5LIB
$ex -h &> /dev/null
if [[ $? -gt 0 ]] || [[ $force ]]; then
	download
fi

toolid=HMMER
getpaths
ex=$([[ $(which hmmsearch) ]] && echo $(which hmmsearch) || echo $GORAP/$tool/bin/hmmsearch)
hmmer=$ex
if [[ ! $($ex -h | grep HMMER | awk '$3=="2.3.2" || $2=="2.3.2"{print 1}') ]] || [[ $force ]]; then
	download
	hmmer=$GORAP/$tool/bin/hmmsearch 
fi
toolid=RNAMMER
getpaths
ex=$([[ $(which rnammer) ]] && echo $(which rnammer) || echo $GORAP/$tool/bin/rnammer)
$ex -v &> /dev/null
if [[ $? -gt 0 ]] || [[ $force ]]; then
	download
	sed -i'' -e "s@HMMSEARCH_BINARY\s*=.*@HMMSEARCH_BINARY='$hmmer';@" $GORAP/$tool/bin/rnammer
fi

toolid=BCHECK
getpaths
ex=$([[ $(which Bcheck) ]] && echo $(which Bcheck) || echo $GORAP/$tool/bin/Bcheck)
$ex -h &> /dev/null
if [[ $? -gt 0 ]] || [[ $force ]]; then
	download
	# mkdir -p $GORAP/$tool/bin
	# mv $GORAP/$tool/* $GORAP/$tool/bin
	sed -i'' '/bob_version/,+3d' $GORAP/$tool/bin/Bcheck
fi

toolid=CRT
getpaths
download

if [[ ! -e /usr/include/zlib.h ]] || [[ $force ]]; then 
	toolid=ZLIB
	getpaths
	download
fi
if [[ ! -e /usr/include/ncurses.h ]] || [[ $force ]]; then 
	toolid=NCURSES
	getpaths
	download
fi
toolid=SAMTOOLS
getpaths
download

echo
echo 'All tools successfully installed'
