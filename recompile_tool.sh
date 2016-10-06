#!/bin/bash
if [[ ! $GORAP ]]; then
	echo 'Setup $GORAP environment variable first.'
	exit 1
fi

retool=$1
pwd=$PWD

if [[ ! $retool ]]; then
	echo 'Do you forgot to add a tool name as parameter? Or'
	read -p 'do you want to recompile all tools [y|n] ' in
	if [[ $in = 'y' ]] || [[ $in = 'Y' ]]; then		
		retool='all'
	else 
		read -p 'Please enter the correct tool name as used in '$GORAP'/<toolname-version> or hit return to exit' retool
		if [[ "$retool" = '' ]]; then
			echo 'Recompilation aborted...'
			exit
		fi
	fi
fi

cd $pwd
RAXML='standard-RAxML-8.2.9'
tool=$RAXML
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then	
	if [[ -d $GORAP/$tool ]]; then
		if [[ $OSTYPE == darwin* ]]; then 
			cpu=$(sysctl -n machdep.cpu.brand_string | awk '{if( $0~/(i5-|i7-)/ || ($0~/Opteron/ && ($NF~/^X/ || $NF > 3000 && $NF < 8000)) ){print 1}}')
		else
			cpu=$(grep -m 1 ^model\ name /proc/cpuinfo | awk '{if( $0~/(i5-|i7-)/ || ($0~/Opteron/ && ($NF~/^X/ || $NF > 3000 && $NF < 8000)) ){print 1}}')
		fi	

		if [[ ! $cpu ]]; then
			cd $GORAP/$tool
			make clean -f Makefile.SSE3.PTHREADS.gcc
			make -f Makefile.SSE3.PTHREADS.gcc
			if [[ $? -gt 0 ]]; then				
				exit 1
			fi 
			mkdir -p bin 
			mv raxmlHPC-PTHREADS-SSE3 bin/raxml 		
			make clean -f Makefile.SSE3.PTHREADS.gcc			
		else
			cd $GORAP/$tool 
			make clean -f Makefile.AVX.PTHREADS.gcc 
			make -f Makefile.AVX.PTHREADS.gcc 
			if [[ $? -gt 0 ]]; then				
				exit 1
			fi
			mkdir -p bin 
			mv raxmlHPC-PTHREADS-AVX bin/raxml
			make clean -f Makefile.AVX.PTHREADS.gcc
		fi
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

cd $pwd
MAFFT='mafft-7.305-with-extensions'
tool=$MAFFT
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d $GORAP/$tool ]]; then
		cd $GORAP/$tool/
		mkdir bin
		cd $GORAP/$tool/core
		sed -iE "s@PREFIX\s*=.*@PREFIX=$GORAP/$tool@" Makefile
		make clean
		make
		if [[ $? -gt 0 ]]; then
			exit 1
		fi
		make install
		make clean
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

cd $pwd
INFERNAL='infernal-1.1.2'
tool=$INFERNAL
infernal=$tool
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d $GORAP/$tool ]]; then
		cd $GORAP/$tool
		./configure --prefix=$PWD
		make clean
		make
		if [[ $? -gt 0 ]]; then				
			exit 1
		fi
		make install
		make clean
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

cd $pwd
tool='esl-alimerge'
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d $GORAP/$infernal ]]; then
		cd $GORAP/$infernal/easel
		./configure --prefix=$GORAP/$infernal
		if [[ ! $(egrep -E 'CFLAGS\s*=' Makefile | grep 'I\.\/include') ]]; then
			sed -iE -r 's/(CFLAGS\s*=.+)/\1 -I\.\/include/' Makefile
		fi
		make clean
		make
		if [[ $? -gt 0 ]]; then
			
			if [[ $(which ldd) ]]; then
				gv=$(ldd --version | head -n 1 | awk '{split($NF,a,"."); print a[1]a[2]}')
			else 
				gv=0
			fi
			if [[ $gv -lt 214 ]]; then
				cd $pwd
				tool='glibc-2.24'
				if [[ ! -d $GORAP/$tool ]]; then
					wget -q --show-progress -T 10 https://ftp.gnu.org/gnu/libc/$tool.tar.gz -O $tool.tar.gz
					tar -xzf $tool.tar.gz -C $GORAP
				fi
				
				cd $GORAP/$tool
				mkdir -p build
				cd build
				../configure --prefix=$GORAP/$infernal/easel
				make clean
				make
				if [[ $? -gt 0 ]]; then				
					exit 1
				fi
				make install
				make clean

				cd $GORAP/$infernal/easel
				if [[ ! $(egrep -E 'CFLAGS\s*=' Makefile | grep 'I\.\/include') ]]; then
					sed -iE -r 's/(CFLAGS\s*=.+)/\1 -I\.\/include/' Makefile
				fi
				make clean
				make
				if [[ $? -gt 0 ]]; then				
					exit 1
				fi
			fi
		fi
		make install 
		make clean	
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

cd $pwd
INFERNAL1='infernal-1.0'
tool=$INFERNAL1
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d $GORAP/$tool ]]; then		
		cd $GORAP/$tool
		./configure --prefix=$PWD
		make clean
		make
		if [[ $? -gt 0 ]]; then				
			exit 1
		fi
		make install
		make clean
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

cd $pwd
RNABOB='rnabob-2.2.1'
tool=$RNABOB
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d $GORAP/$tool ]]; then			
		cd $GORAP/$tool
		mkdir -p bin
		mkdir -p man/man1
		sed -i '/^BINDIR/ c\BINDIR = $(CURDIR)/bin' Makefile
		sed -i '/^MANDIR/ c\MANDIR = $(CURDIR)/man' Makefile
		make clean
		make
		if [[ $? -gt 0 ]]; then				
			exit 1
		fi
		make install
		make clean
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

cd $pwd
TRNASCAN='tRNAscan-SE-1.3.1'
tool=$TRNASCAN
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d "$GORAP/$tool" ]]; then
		cd $GORAP/$tool
		mkdir -p bin
		mkdir -p man/man1
		mkdir -p lib/tRNAscan-SE
		sed -i '/^BINDIR/ c\BINDIR = $(CURDIR)/bin' Makefile 
		sed -i '/^MANDIR/ c\MANDIR = $(CURDIR)/man' Makefile
		sed -i '/^LIBDIR/ c\LIBDIR = $(CURDIR)/lib/tRNAscan-SE' Makefile
		make clean
		make
		if [[ $? -gt 0 ]]; then				
			exit 1
		fi
		make install
		make clean
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

cd $pwd
HMMER='hmmer-2.3.2'
tool=$HMMER
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then	
	if [[ -d $GORAP/$tool ]]; then
		cd $GORAP/$tool
		make clean
		./configure --enable-threads --enable-lfs --prefix=$PWD
		make
		if [[ $? -gt 0 ]]; then				
			exit 1
		fi
		make install
		make clean
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

cd $pwd
SAMTOOLS='samtools-0.1.19'
tool=$SAMTOOLS
samtool=$tool

ZLIB='zlib-1.2.8'
tool=$ZLIB
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d "$GORAP/$tool" ]]; then
		cd $GORAP/$tool 
		make clean
		./configure --prefix=$GORAP/$samtool
		make
		if [[ $? -gt 0 ]]; then				
			exit 1
		fi
		make install
		make clean
	else 
		echo 'Please run install_tools.sh first, then try again'		
	fi
fi

NCURSES='ncurses-6.0'
tool=$NCURSES
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d "$GORAP/$tool" ]]; then
		cd $GORAP/$tool 
		make clean
		./configure --prefix=$GORAP/$samtool
		make
		if [[ $? -gt 0 ]]; then
			exit 1
		fi
		make install
		make clean
		cp $GORAP/$samtool/include/ncurses/* $GORAP/$samtool/include/
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

tool=$samtool
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d "$GORAP/$tool" ]]; then
		cd $GORAP/$tool
		if [[ ! $(egrep -E 'CFLAGS\s*=' Makefile | grep 'I\.\/include') ]]; then
			sed -i '/^CFLAGS/ c\CFLAGS = -g -Wall -O2 -fPIC -I./include -L./lib #-m64 #-arch ppc' Makefile
		fi
		make clean
		make 
		if [[ $? -gt 0 ]]; then
			exit 1
		fi 
		mkdir -p bin 
		mv samtools bin/samtools
		mv libbam.a libbam.x
		make clean
		mv libbam.x libbam.a
	else 
		echo 'Please run install_tools.sh first, then try again'		
	fi
fi
