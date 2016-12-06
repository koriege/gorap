#!/bin/bash
if [[ ! $GORAP ]]; then
	echo 'Setup $GORAP environment variable first.'
	exit 1
fi
mkdir -p $GORAP

pwd=$PWD
rm -f install.log &> /dev/null
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

extract (){
	echo "Extracting $tool"
	# progress &
	# pid=$!
	tar -xzf $1 -C $GORAP
	if [[ $? -gt 0 ]]; then
		kill $pid &> /dev/null
		wait $pid &> /dev/null
		echo
		echo "Extracting $tool failed - see install.log for details"
		exit 1
	fi 
	# kill $pid &> /dev/null
	# wait $pid &> /dev/null
	# echo
}

recompile (){
	echo "Installing $tool"
	progress &
	pid=$!
	bash recompile_tool.sh $tool &>> $pwd/install.log
	if [[ $? -gt 0 ]]; then
		kill $pid &> /dev/null
		wait $pid &> /dev/null
		echo
		echo "Installing $tool failed - see install.log for details"
		exit 1
	fi 
	kill $pid &> /dev/null
	wait $pid &> /dev/null
	echo
}

download (){
	toolid=$tool
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
		dlpath=''
		if [[ ! -d $GORAP/$tool ]]; then
			echo
			echo "Downloading"	
			mkdir -p $GORAP/newicktopdf/bin
			echo "Extracting $tool"
			wget -q --show-progress -T 10 ftp://pbil.univ-lyon1.fr/pub/mol_phylogeny/njplot/newicktopdf -O $GORAP/newicktopdf/bin/newicktopdf
			echo "Installing $tool"
			echo "\\"
			chmod 755 $GORAP/newicktopdf/bin/newicktopdf
		fi
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
		dlpath=''
		tool=$(ls -d $GORAP/ncbi-blast*/ 2> /dev/null)
		if [[ $tool ]]; then
			tool=$(basename $tool)
		else
			echo
			echo "Downloading"
			rm -f ncbi-blast*.tar.gz
			if [[ $bit -eq 32 ]]; then
				if [[ $os == 'mac' ]]; then
					wget -q --show-progress -T 10 ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/ncbi-blast-2.4.0+-universal-macosx.tar.gz
				else
					wget -q --show-progress -T 10 ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-ia32-linux.tar.gz
				fi
			else
				if [[ $os == 'mac' ]]; then
					wget -q --show-progress -T 10 -r -nd -np -l 1 -A *x64-macosx.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
				else
					wget -q --show-progress -T 10 -r -nd -np -l 1 -A *x64-linux.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
				fi
			fi
			tool=$(ls ncbi-blast*.tar.gz | egrep -Eo 'ncbi-blast-[0-9]+\.[0-9]+\.[0-9]+\+')
			extract $(ls ncbi-blast*.tar.gz)
			rm -f ncbi-blast*.tar.gz
			echo "Installing $tool"
			echo "\\"
		fi
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
		tool='Bcheck-0.6'
		dlpath='http://rna.tbi.univie.ac.at/bcheck/Bcheck.tgz'
	;;
	'CRT')
		tool='crt-1.2'
		dlpath=''
		if [[ ! -d $GORAP/$tool ]]; then
			echo "Downloading"
			mkdir -p $GORAP/$tool/bin
			wget -q --show-progress -T 10 http://www.room220.com/crt/CRT1.2-CLI.jar.zip
			echo "Extracting $tool"
			unzip -q -d $GORAP/$tool/bin CRT1.2-CLI.jar.zip		
			rm -f CRT1.2-CLI.jar.zip
			echo "Installing $tool"
			echo "\\"
			mv $GORAP/$tool/bin/CRT1.2-CLI.jar $GORAP/$tool/bin/crt.jar
		fi
	;;
	'SAMTOOLS')
		tool='samtools-0.1.19'
		dlpath='https://github.com/samtools/samtools/archive/0.1.19.tar.gz'
	;;
	'ZLIB')
		tool='zlib-1.2.8'
		dlpath="http://zlib.net/$tool.tar.gz"
	;;
	'NCURSES')
		tool='ncurses-6.0'
		dlpath="https://ftp.gnu.org/pub/gnu/ncurses/$tool.tar.gz"
	;;
	*)
		echo
		echo "Downloading $tool failed - see install.log for details"
		exit 1
	;;
	esac
	if [[ $dlpath ]] && [[ ! -d $GORAP/$tool ]]; then
		echo "Downloading"
		wget -q --show-progress -T 10 $dlpath -O $tool.tar.gz
		if [[ $? -gt 0 ]]; then 
			echo "Download failed - see install.log for details"
			exit 1
		fi
		extract $tool.tar.gz
		rm -f $tool.tar.gz
	fi
	sed -i "s/$toolid=.*/$toolid='$tool'/g" recompile_tool.sh
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
	echo '64 bit operating system detected. Compiled software dependencies available.'
	read -p 'Do you want to download them instead of starting compilation [y|n]? ' in
	if [[ $in == 'y' ]] || [[ $in == 'Y' ]]; then
		echo
		echo 'Downloading software'
		wget -q --show-progress -T 10 www.rna.uni-jena.de/supplements/gorap/dependencies.tar.gz -O dependencies.tar.gz
		extract dependencies.tar.gz
	fi
fi

if [[ $bit -eq 64 ]]; then
	cd $pwd
	tool=RAXML
	ex=$(which raxml)
	if [[ $ex ]]; then
		raxml -h &> /dev/null
		if [ $? -gt 0 ]; then
			unset ex
		fi
	fi
	if [[ ! $ex ]]; then
		download
		$GORAP/$tool/bin/raxml -h &> /dev/null
		if [[ $? -gt 0 ]]; then				
			recompile
		fi
	fi

	cd $pwd
	tool=MAFFT
	ex=$(which mafft)
	if [[ $ex ]]; then
		if [[ ! $(mafft -h 2>&1 | grep MAFFT | wc -l) ]]; then
			unset ex
		fi
	fi
	if [[ ! $ex ]]; then
		download
		if [[ ! $($GORAP/$tool/bin/mafft -h 2>&1 | grep MAFFT) ]]; then
			recompile
		fi
	fi

	cd $pwd
	tool=NEWICKTOPDF
	ex=$(which newicktopdf)
	if [[ ! $ex ]]; then
		download
		$GORAP/$tool/bin/newicktopdf -h &> /dev/null
		if [[ $? -gt 0 ]]; then
			echo
			echo $tool' installation failed - see install.log for details'
			exit 1
		fi
	fi
fi

cd $pwd
tool=INFERNAL
ex=$(which cmsearch)
if [[ $ex ]]; then
	cmsearch -h &> /dev/null
	if [[ $? -gt 0 ]]; then
		unset ex
	else 
		ex=$(cmsearch -h | grep INFERNAL | awk '$3>=1.1{print 1}')
	fi
fi
if [[ ! $ex ]]; then
	$GORAP/$tool/bin/cmsearch -h &> /dev/null
	if [[ $? -gt 0 ]]; then
		recompile
	fi
fi
cd $pwd
ex=$(which esl-alimerge | wc -l)
if [[ $ex -gt 0 ]]; then
	esl-alimerge -h &> /dev/null
	if [[ $? -gt 0 ]]; then
		ex=0		
	fi			
fi
if [[ $ex -eq 0 ]]; then
	download
	$GORAP/$tool/bin/esl-alimerge -h &> /dev/null
	if [[ $? -gt 0 ]]; then
		tool='esl-alimerge'
		recompile
	fi
fi

cd $pwd
tool=INFERNAL1
ex=$(which cmsearch)
if [[ $ex ]]; then	
	cmsearch -h &> /dev/null
	if [[ $? -gt 0 ]]; then
		unset ex
	else
		ex=$(cmsearch -h | grep INFERNAL | awk '$3=="1.0"{print 1}')
	fi
fi
if [[ ! $ex ]]; then
	download
	$GORAP/$tool/bin/cmsearch -h &> /dev/null
	if [[ $? -gt 0 ]]; then		
		recompile
	fi
fi

cd $pwd
tool=RNABOB
ex=$(which rnabob)
if [[ $ex ]]; then
	rnabob -h &> /dev/null
	if [[ $? -gt 0 ]]; then
		unset ex
	fi
fi
if [[ ! $ex ]]; then
	download
	$GORAP/$tool/bin/rnabob -h &> /dev/null
	if [[ $? -gt 0 ]]; then		
		recompile
	fi
fi

cd $pwd
tool=BLAST
ex=$(which blastn) 
if [[ $ex ]]; then
	blastn -h &> /dev/null
	if [[ $? -gt 0 ]]; then
		unset ex
	else
		ex=$(blastn -h | grep DESCRIPTION -A 1 | tail -n 1 | awk '{split($NF,a,"."); if(a[2]>=2 && $NF~/\+$/){print 1}}')
	fi
fi
if [[ ! $ex ]]; then
	download
	$GORAP/$tool/bin/blastn -h &> /dev/null
	if [[ $? -gt 0 ]]; then
		echo
		echo $tool' installation failed - see install.log for details'
		exit 1
	fi
fi

cd $pwd
tool=TRNASCAN
ex=$(which tRNAscan-SE) 
if [[ $ex -gt 0 ]]; then
	tRNAscan-SE -h &> /dev/null
	if [[ $? -gt 0 ]]; then
		unset ex
	fi
fi
if [[ ! $ex ]]; then
	download
	export PERL5LIB=$GORAP/$tool:$PERL5LIB
	$GORAP/$tool/bin/tRNAscan-SE -h &> /dev/null
	if [[ $? -gt 0 ]]; then
		recompile
	fi
fi

cd $pwd
tool=HMMER
ex=$(which hmmsearch)
hmmer=$ex
if [[ $ex ]]; then
	hmmsearch -h &> /dev/null
	if [ $? -gt 0 ]; then
		unset ex
	else
		ex=$(hmmsearch -h | grep HMMER | awk '$3=="2.3.2" || $2=="2.3.2"{print 1}')
	fi
fi
if [[ ! $ex ]]; then
	download
	$GORAP/$tool/bin/hmmsearch -h &> /dev/null
	if [[ $? -gt 0 ]]; then
		recompile
	fi
	hmmer=$GORAP/$tool/bin/hmmsearch
fi

cd $pwd
tool=RNAMMER
ex=$(which rnammer)
if [[ $ex ]]; then
	rnammer -v &> /dev/null
	if [[ $? -gt 0 ]]; then
		unset ex
	fi
fi
if [[ ! $ex ]]; then
	download
	sed -iE "s@HMMSEARCH_BINARY\s*=.*@HMMSEARCH_BINARY='$hmmer';@" $GORAP/$tool/bin/rnammer
	$GORAP/$tool/bin/rnammer -v &> /dev/null
	if [[ $? -gt 0 ]]; then
		echo
		echo $tool' installation failed - see install.log for details'
		exit 1
	fi
fi

cd $pwd
tool=BCHECK
ex=$(which Bcheck)
if [[ $ex ]]; then
	Bcheck -h &> /dev/null
	if [[ $? -gt 0 ]]; then
		unset ex
	fi
fi
if [[ ! $ex ]]; then
	download
	mkdir -p $GORAP/$tool/bin
	mv $GORAP/$tool/* $GORAP/$tool/bin &> /dev/null
	sed -i '/bob_version/,+3d' $GORAP/$tool/bin/Bcheck
fi

cd $pwd
tool='CRT'
download

cd $pwd
tool=SAMTOOLS
download
samtool=$tool
if [[ ! -e /usr/include/zlib.h ]]; then 
	cd $pwd
	tool=ZLIB
	download
	recompile
fi
if [[ ! -e /usr/include/ncurses.h ]]; then 
	cd $pwd
	tool=NCURSES
	download
	recompile
fi
tool=$samtool
recompile

echo
echo 'All tools are successfully installed'
