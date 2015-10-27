#!/bin/bash
pwd=`pwd`
bit=$(uname -m | awk '{if ($0 == "i686"){print "32-"}else{print ""}}')
if [[ $OSTYPE == darwin* ]]; then 
	bit='mac-'
fi
tools=$GORAP
mkdir -p  $tools
rm $pwd/install.log 2> $pwd/install.log > $pwd/install.log

progress(){
	x=0
	while [ 1 ]; do	
		x=$((x+1))
		if [ "$x" -eq "2" ]
			then
				x=0
				echo -en "\r/ for more information execute: tail -f install.log"
			else
				echo -en "\r\\ for more information execute: tail -f install.log"
		fi
		sleep 1
	done
}

download () {
	cd $tools
	if [ ! -d "$tool" -a "$ex" -eq "0" ]
		then			
			echo 'Downloading '$tool	
			echo ''
			wget -T 10 www.rna.uni-jena.de/supplements/gorap/$bit$tool.tar.gz
			echo 'Extracting '$bit$tool
			tar -xzf $bit$tool.tar.gz -C $tools	
			rm $bit$tool.tar.gz
			echo ''
		else
			echo 'Nothing to be done for '$tool
			echo ''
	fi
}

if [ "$bit" = "32-" ]
	then
		echo 'I am sorry, it is not possible to calculate SSU rRNA based phylogenies on a 32 bit operating system'
		echo 'Skipping installation of RAxML'
	else
		echo 'Do you want to calculate SSU rRNA based phylogenies and'
		read -p 'therefor install RAxML and Mafft? [y|n] ' in
		bit=""
		if [ "$in" = "y" -o "$in" = "Y" ]
			then
				tool='RAxML-7.4.2'
				ex=$(which raxml | wc | awk '{print $1}')
				if [ "$ex" -gt "0" ]; then
					raxml -h 2>> $pwd/install.log >> $pwd/install.log
					if [ $? -gt 0 ]; then
						ex=0
					fi
				fi
				download
				if [ "$ex" -eq "0" ]; then
					$tools/$tool/bin/raxml -h 2>> $pwd/install.log >> $pwd/install.log
					if [ $? -gt 0 ]; then
						cd $pwd
						bash recompile_tool.sh $tool
					fi
				fi

				tool='mafft-7.017-with-extensions' 
				ex=$(which mafft | wc | awk '{print $1}')
				if [ "$ex" -gt "0" ]; then
					mafft -h 2>> $pwd/install.log >> $pwd/install.log
					if [ $? -gt 0 ]; then
						ex=0
					fi
				fi
				download
				if [ "$ex" -eq "0" ]; then
					$tools/$tool/bin/mafft -h 2>> $pwd/install.log >> $pwd/install.log
					if [ $? -gt 0 ]; then
						cd $pwd
						bash recompile_tool.sh $tool
					fi
				fi
		fi
fi

tool='infernal-1.1rc2'
ex=$(which cmsearch | wc | awk '{print $1}')
if [ "$ex" -gt "0" ]
	then
		cmsearch -h 2>> $pwd/install.log >> $pwd/install.log
		if [ $? -gt 0 ]; then
			ex=0
		else 
			ex=$(cmsearch -h | grep INFERNAL | awk '{if($3>=1.1){print "1"}else{print "0"}}')
		fi
		download
	else
		download
fi
if [ "$ex" -eq "0" ]; then
	$tools/$tool/bin/cmsearch -h 2>> $pwd/install.log >> $pwd/install.log
	if [ $? -gt 0 ]; then
		cd $pwd
		bash recompile_tool.sh $tool
	fi
fi

infernal='infernal-1.1rc2'
ex=$(which esl-alimerge | wc | awk '{print $1}')
if [ "$ex" -gt "0" ]
	then
		esl-alimerge -h 2>> $pwd/install.log >> $pwd/install.log
		if [ $? -gt 0 ]; then
			ex=0		
			download				
			
			glibc=$(ldd --version | head -n 1 | awk '{split($NF,a,"."); print a[2]}')
			if [ "$glibc" -lt "14" ]; then
				tool='glibc-2.21'
				ex=0
				download

				echo 'Installing ESL-Alimerge'
				progress &
				pid=$!

				cd $tools/$tool
				mkdir -p build
				cd build
				../configure --prefix=$tools/$infernal/easel 2>> $pwd/install.log >> $pwd/install.log	
				make 2>> $pwd/install.log >> $pwd/install.log	
				make install 2>> $pwd/install.log >> $pwd/install.log	
				make clean 2>> $pwd/install.log >> $pwd/install.log	

				cd $tools/$infernal/easel
				cp -r include/* .
				
				./configure --prefix=$tools/$infernal 2>> $pwd/install.log >> $pwd/install.log	
				sed -i '/^CFLAGS/ c\CFLAGS = -O3 -fomit-frame-pointer -fstrict-aliasing -pthread -fPIC -I. -I..' Makefile							
				make 2>> $pwd/install.log >> $pwd/install.log	
				make install 2>> $pwd/install.log >> $pwd/install.log	
				make clean 2>> $pwd/install.log >> $pwd/install.log	

				kill $pid >/dev/null
			else
				echo 'Installing ESL-Alimerge'
				progress &
				pid=$!				

				cd $tools/$infernal/easel
				./configure --prefix=$tools/$infernal 2>> $pwd/install.log >> $pwd/install.log		
				make 2>> $pwd/install.log >> $pwd/install.log	
				make install 2>> $pwd/install.log >> $pwd/install.log	
				make clean 2>> $pwd/install.log >> $pwd/install.log

				kill $pid >/dev/null
			fi			
		fi		
	else
		download

		glibc=$(ldd --version | head -n 1 | awk '{split($NF,a,"."); print a[2]}')
		if [ "$glibc" -lt "14" ]; then
			tool='glibc-2.21'
			ex=0
			download

			echo 'Installing ESL-Alimerge'
			progress &
			pid=$!

			cd $tools/$tool
			mkdir -p build
			cd build
			../configure --prefix=$tools/$infernal/easel 2>> $pwd/install.log >> $pwd/install.log	
			make 2>> $pwd/install.log >> $pwd/install.log	
			make install 2>> $pwd/install.log >> $pwd/install.log	
			make clean 2>> $pwd/install.log >> $pwd/install.log	

			cd $tools/infernal-1.1rc2/easel
			cp -r include/* .
			
			./configure --prefix=$tools/$infernal 2>> $pwd/install.log >> $pwd/install.log
			sed -i '/^CFLAGS/ c\CFLAGS = -O3 -fomit-frame-pointer -fstrict-aliasing -pthread -fPIC -I. -I..' Makefile							
			make 2>> $pwd/install.log >> $pwd/install.log	
			make install 2>> $pwd/install.log >> $pwd/install.log	
			make clean 2>> $pwd/install.log >> $pwd/install.log	

			kill $pid >/dev/null
		else
			echo 'Installing ESL-Alimerge'
			progress &
			pid=$!

			cd $tools/$infernal/easel			
			./configure --prefix=$tools/$infernal 2>> $pwd/install.log >> $pwd/install.log
			make 2>> $pwd/install.log >> $pwd/install.log	
			make install 2>> $pwd/install.log >> $pwd/install.log	
			make clean 2>> $pwd/install.log >> $pwd/install.log

			kill $pid >/dev/null
		fi
fi

tool='infernal-1.0'
ex=$(which cmsearch | wc | awk '{print $1}')
if [ "$ex" -gt "0" ]
	then
		cmsearch -h 2>> $pwd/install.log >> $pwd/install.log
		if [ $? -gt 0 ]; then
			ex=0
		else
			ex=$(cmsearch -h | grep INFERNAL | awk '{if($3>=1 && $3<1.1){print "1"}else{print "0"}}')
		fi
		download
	else
		download
fi
if [ "$ex" -eq "0" ]; then
	$tools/$tool/bin/cmsearch -h 2>> $pwd/install.log >> $pwd/install.log
	if [ $? -gt 0 ]; then
		cd $pwd
		bash recompile_tool.sh $tool
	fi
fi

tool='rnabob-2.2'
ex=$(which rnabob | wc | awk '{print $1}')
if [ "$ex" -gt "0" ]; then
	rnabob -h 2>> $pwd/install.log >> $pwd/install.log
	if [ $? -gt 0 ]; then
		ex=0
	fi
fi
download
if [ "$ex" -eq "0" ]; then
	$tools/$tool/bin/rnabob -h 2>> $pwd/install.log >> $pwd/install.log
	if [ $? -gt 0 ]; then
		cd $pwd
		bash recompile_tool.sh $tool
	fi
fi

tool='hmmer-2.3.2'
ex=$(which hmmsearch | wc | awk '{print $1}')
if [ "$ex" -gt "0" ]
	then
		hmmsearch -h 2>> $pwd/install.log >> $pwd/install.log
		if [ $? -gt 0 ]; then
			ex=0
		else
			ex=$(hmmsearch -h | grep HMMER | awk '{if($2=="2.3.2"){print "1"}else{print "0"}}')
		fi
		download
	else
		download
fi
if [ "$ex" -eq "0" ]; then
	$tools/$tool/bin/hmmsearch -h 2>> $pwd/install.log >> $pwd/install.log
	if [ $? -gt 0 ]; then
		cd $pwd
		bash recompile_tool.sh $tool
	fi
fi

tool='ncbi-blast-2.2.27+'
ex=$(which blastn | wc | awk '{print $1}') 
if [ "$ex" -gt "0" ]
	then			
		blastn -h 2>> $pwd/install.log >> $pwd/install.log
		if [ $? -gt 0 ]; then
			ex=0
		else
			ex=$(blastn -h | grep DESCRIPTION -A 1 | tail -n 1 | awk '{print $3; if ($3~/^2\.2\.2.+\+$/){print "1"}else{print "0"}}')
		fi
fi
download
if [ "$ex" -eq "0" ]; then
	$tools/$tool/bin/blastn -h 2>> $pwd/install.log >> $pwd/install.log
	if [ $? -gt 0 ]; then
		cd $pwd
		bash recompile_tool.sh $tool
	fi
fi

tool='rnammer-1.2'
ex=$(which rnammer | wc | awk '{print $1}')
if [ "$ex" -gt "0" ]; then
	rnammer -v 2>> $pwd/install.log >> $pwd/install.log
	if [ $? -gt 0 ]; then
		ex=0
	fi
fi
download
if [ "$ex" -eq "0" ]; then
	$tools/$tool/bin/rnammer -v 2>> $pwd/install.log >> $pwd/install.log
	if [ $? -gt 0 ]; then
		cd $pwd
		bash recompile_tool.sh $tool
	fi
fi

tool='bcheck-0.6'
ex=$(which Bcheck | wc | awk '{print $1}')
if [ "$ex" -gt "0" ]; then
	Bcheck 2>> $pwd/install.log >> $pwd/install.log
	if [ $? -gt 0 ]; then
		ex=0
	fi
fi
download
if [ "$ex" -eq "0" ]; then
	$tools/$tool/bin/Bcheck 2>> $pwd/install.log >> $pwd/install.log
	if [ $? -gt 0 ]; then
		cd $pwd
		bash recompile_tool.sh $tool
	fi
fi

tool='crt-1.2'
ex=$(which CRT1.2-CLI.jar | wc | awk '{print $1}')
download
		
samtools='samtools-0.1.19'
ex=0
download
if [ -d "$tools/$samtools" ]
	then
		
		zlib='zlib-1.2.8'
		ex=0
		download						

		ncurses='ncurses-5.9'
		ex=0
		download
		
		echo 'Installing '$samtools
		
		progress &
		pid=$!

		cd $tools/$zlib
		./configure --prefix=$tools/$samtools 2>> $pwd/install.log >> $pwd/install.log	
		make 2>> $pwd/install.log >> $pwd/install.log
		make install 2>> $pwd/install.log >> $pwd/install.log
		make clean 2>> $pwd/install.log >> $pwd/install.log
		
		cd $tools/$ncurses
		./configure --prefix=$tools/$samtools 2>> $pwd/install.log >> $pwd/install.log
		make 2>> $pwd/install.log >> $pwd/install.log
		make install 2>> $pwd/install.log >> $pwd/install.log
		make clean 2>> $pwd/install.log >> $pwd/install.log

		cd $tools/$samtools
		cp -r include/* .
        cp -r lib/* .        
        #cp -r $tools/$zlib/lib/* .		
		make 2>> $pwd/install.log >> $pwd/install.log
		                       
		mkdir -p bin 
		cp samtools bin/samtools
		
		kill $pid >/dev/null
fi

tool='tRNAscan-SE-1.3.1' 
ex=$(which tRNAscan-SE | wc | awk '{print $1}') 
if [ "$ex" -gt "0" ]; then
	tRNAscan-SE -h 2>> $pwd/install.log >> $pwd/install.log
	if [ $? -gt 0 ]; then
		ex=0
	fi
fi
download
if [ "$ex" -eq "0" ]
	then
		echo 'Installing '$tool
		progress &
		pid=$!
				
		cd $tools/tRNAscan-SE-1.3.1
		make 2>> $pwd/install.log >> $pwd/install.log
		make install 2>> $pwd/install.log >> $pwd/install.log
		make clean 2>> $pwd/install.log >> $pwd/install.log
		cd $pwd
		
		kill $pid >/dev/null
fi