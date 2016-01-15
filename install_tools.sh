#!/bin/bash
pwd=$PWD
bit=$(uname -m | awk '{if ($0 == "i686"){print "32-"}else{print ""}}')
if [[ $OSTYPE == darwin* ]]; then 
	bit='mac-'
fi
mkdir -p $GORAP
rm -f $pwd/install.log 2> $pwd/install.log > $pwd/install.log

progress(){
	x=0
	while [ 1 ]; do	
		x=$((x+1))
		if [ "$x" -eq "2" ]
			then
				x=0
				echo -en "\r/ for more information type: tail -f install.log"
			else
				echo -en "\r\\ for more information type: tail -f install.log"
		fi
		sleep 1
	done
}

download () {
	if [[ ! -d $GORAP/$tool ]] && [[ $ex -eq 0 ]]; then
		echo
		echo 'Downloading '$tool	
		echo		
		wget -T 10 -N www.rna.uni-jena.de/supplements/gorap/$bit$tool.tar.gz
		if [[ $? -gt 0 ]]; then
			echo
			echo $tool' download failed'
			exit 1
		fi 
		echo 'Extracting '$tool		
		tar -xzf $bit$tool.tar.gz -C $GORAP
		rm $bit$tool.tar.gz
	fi
}

recompile () {
	echo
	echo 'Installing '$tool
	echo
	progress &
	pid=$!
	bash recompile_tool.sh $tool 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then
		kill $pid 2>/dev/null > /dev/null
		wait $pid 2>/dev/null > /dev/null
		echo
		echo $tool' installation failed'
		exit 1
	fi 
	kill $pid 2>/dev/null > /dev/null
	wait $pid 2>/dev/null > /dev/null
	echo
}

if [[ $bit = '32-' ]]; then
	echo 'I am sorry, it is not possible to calculate phylogenies with RAxML on a 32 bit operating system'
	echo 'Skipping installation of RAxML and Mafft'
else
	echo 'Do you want to calculate phylogenies and'
	read -p 'therefor install RAxML and Mafft? [y|n] ' in
	bitbu=$bit
	bit=""
	if [[ $in = 'y' ]] || [[ $in = 'Y' ]]; then
		cd $pwd
		tool='RAxML-7.4.2'
		ex=$(which raxml | wc | awk '{print $1}')
		if [[ $ex -gt 0 ]]; then
			raxml -h 2>> $pwd/install.log >> $pwd/install.log
			if [ $? -gt 0 ]; then
				ex=0
			fi
		fi
		download
		if [[ $ex -eq 0 ]]; then
			$GORAP/$tool/bin/raxml -h 2>> $pwd/install.log >> $pwd/install.log
			if [[ $? -gt 0 ]]; then				
				recompile
			fi
		fi

		cd $pwd
		tool='mafft-7.017-with-extensions' 
		ex=$(which mafft | wc | awk '{print $1}')
		if [[ $ex -gt 0 ]]; then
			mafft --help 2>> $pwd/install.log >> $pwd/install.log
			if [[ $? -gt 0 ]]; then
				ex=0
			fi
		fi
		download
		if [[ $ex -eq 0 ]]; then
			$GORAP/$tool/bin/mafft -h 2>> $pwd/install.log >> $pwd/install.log
			if [[ $? -gt 0 ]]; then
				recompile
			fi
		fi
	fi
	bit=$bitbu
fi

cd $pwd
tool='infernal-1.1rc2'
ex=$(which cmsearch | wc | awk '{print $1}')
if [[ $ex -gt 0 ]]; then
	cmsearch -h 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then
		ex=0
	else 
		ex=$(cmsearch -h | grep INFERNAL | awk '{if($3>=1.1){print "1"}else{print "0"}}')
	fi
fi
download
if [[ $ex -eq 0 ]]; then
	$GORAP/$tool/bin/cmsearch -h 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then
		recompile
	fi
fi

cd $pwd
tool='infernal-1.1rc2'
ex=$(which esl-alimerge | wc | awk '{print $1}')
if [[ $ex -gt 0 ]]; then
	esl-alimerge -h 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then
		ex=0		
	fi			
fi
download
if [[ $ex -eq 0 ]]; then
	$GORAP/$tool/bin/esl-alimerge -h 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then
		tool='esl-alimerge'
		recompile
	fi
fi

cd $pwd
tool='infernal-1.0'
ex=$(which cmsearch | wc | awk '{print $1}')
if [[ $ex -gt 0 ]]; then	
	cmsearch -h 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then
		ex=0
	else
		ex=$(cmsearch -h | grep INFERNAL | awk '{if($3=="1.0"){print "1"}else{print "0"}}')
	fi
fi
download
if [[ $ex -eq 0 ]]; then
	$GORAP/$tool/bin/cmsearch -h 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then		
		recompile
	fi
fi

cd $pwd
tool='rnabob-2.2'
ex=0
ex=$(which rnabob | wc | awk '{print $1}')
if [[ $ex -gt 0 ]]; then
	rnabob -h 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then
		ex=0
	fi
fi
download
if [[ $ex -eq 0 ]]; then
	$GORAP/$tool/bin/rnabob -h 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then		
		recompile
	fi
fi

cd $pwd
tool='ncbi-blast-2.2.30+'
ex=$(which blastn | wc | awk '{print $1}') 
if [[ $ex -gt 0 ]]; then
	blastn -h 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then
		ex=0
	else
		ex=$(blastn -h | grep DESCRIPTION -A 1 | tail -n 1 | awk '{if ($3~/^2\.2\.2.+\+$/ || $3~/^2\.2\.3.+\+$/ || $3~/^2\.3\..+\+$/){print "1"}else{print "0"}}')		
	fi
fi
download
if [[ $ex -eq 0 ]]; then
	$GORAP/$tool/bin/blastn -h 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then		
		recompile
	fi
fi

cd $pwd
tool='tRNAscan-SE-1.3.1' 
ex=$(which tRNAscan-SE | wc | awk '{print $1}') 
if [[ $ex -gt 0 ]]; then
	tRNAscan-SE -h 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then
		ex=0
	fi
fi
download
if [[ $ex -eq 0 ]]; then
	$GORAP/$tool/bin/tRNAscan-SE -h 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then
		recompile
	fi
fi

cd $pwd
tool='hmmer-2.3.2'
ex=$(which hmmsearch | wc | awk '{print $1}')
if [[ $ex -gt 0 ]]; then
	hmmsearch -h 2>> $pwd/install.log >> $pwd/install.log
	if [ $? -gt 0 ]; then
		ex=0
	else
		ex=$(hmmsearch -h | grep HMMER | awk '{if($2=="2.3.2"){print "1"}else{print "0"}}')
	fi
fi
download
if [[ $ex -eq 0 ]]; then
	$GORAP/$tool/bin/hmmsearch -h 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then
		recompile
	fi
fi

cd $pwd
tool='rnammer-1.2'
ex=$(which rnammer | wc | awk '{print $1}')
if [[ $ex -gt 0 ]]; then
	rnammer -v 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then
		ex=0
	fi
fi
download

cd $pwd
tool='bcheck-0.6'
ex=$(which Bcheck | wc | awk '{print $1}')
if [[ $ex -gt 0 ]]; then
	Bcheck 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then
		ex=0
	fi
fi
download

cd $pwd
tool='crt-1.2'
ex=$(which CRT1.2-CLI.jar | wc | awk '{print $1}')
download

cd $pwd
tool='samtools-0.1.19'
ex=0
download
if [[ $ex -eq 0 ]]; then
	$GORAP/$tool/bin/samtools 2>> $pwd/install.log >> $pwd/install.log
	if [[ $? -gt 0 ]]; then
		recompile
	fi
fi

echo
echo 'All tools are successfully installed'