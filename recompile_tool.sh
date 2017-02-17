#! /usr/bin/env bash

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
			make clean -f Makefile.SSE3.PTHREADS.gcc; make -f Makefile.SSE3.PTHREADS.gcc && make clean -f Makefile.SSE3.PTHREADS.gcc
			if [[ $? -gt 0 ]]; then				
				exit 1
			fi 
			mkdir -p bin 
			mv raxmlHPC-PTHREADS-SSE3 bin/raxml
		else
			cd $GORAP/$tool 
			make clean -f Makefile.AVX.PTHREADS.gcc; make -f Makefile.AVX.PTHREADS.gcc && make clean -f Makefile.AVX.PTHREADS.gcc
			if [[ $? -gt 0 ]]; then				
				exit 1
			fi
			mkdir -p bin 
			mv raxmlHPC-PTHREADS-AVX bin/raxml
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
		cd $GORAP/$tool/core
		make clean;	sed -i'''' -r "s@PREFIX\s*=.*@PREFIX=$GORAP/$tool@" Makefile && make PREFIX=$GORAP/$tool && make install && make clean
		if [[ $? -gt 0 ]]; then
			exit 1
		fi
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

cd $pwd
NEWICKTOPDF='newicktopdf'
tool=$NEWICKTOPDF
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d $GORAP/$tool ]]; then
		chmod 755 $GORAP/$tool/bin/$tool
		if [[ $? -gt 0 ]]; then
			exit 1
		fi
	else
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

cd $pwd
GLIBC='glibc-2.24'
tool=$GLIBC
glibc=$tool
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d $GORAP/$tool ]]; then
		cd $GORAP/$tool
		rm -rf built; mkdir -p built
		cd built
		../configure --prefix=$PWD && make && make install && make clean
		if [[ $? -gt 0 ]]; then
			exit 1
		fi
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

cd $pwd
INFERNAL='infernal-1.1.2'
tool=$INFERNAL
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d $GORAP/$tool ]]; then
		cd $GORAP/$tool
		make clean; ./configure --prefix=$PWD && make && make install && make clean
		if [[ $? -gt 0 ]]; then				
			exit 1
		fi

		cd $GORAP/$infernal/easel
		if [[ $(which ldd) ]]; then
			gv=$(ldd --version | head -n 1 | awk '{split($NF,a,"."); print a[1]a[2]}')
		else 
			gv=0
		fi
		if [[ $gv -lt 214 ]]; then
			make clean; ./configure --prefix=$GORAP/$tool CFLAGS="-I$GORAP/$glibc/built/include" LDFLAGS="-I$GORAP/$glibc/built/lib" && sed -i'' -r 's/^\s*(C|LD)FLAGS\s*=\s*/override \1FLAGS += /' Makefile && make CFLAGS="-I$GORAP/$glibc/built/include" LDFLAGS="-I$GORAP/$glibc/built/lib" && make install && make clean
		else 
			make clean; ./configure --prefix=$GORAP/$tool && make && make install && make clean
		fi
		if [[ $? -gt 0 ]]; then
			exit 1
		fi
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
		make clean; ./configure --prefix=$PWD && make && make install && make clean
		if [[ $? -gt 0 ]]; then				
			exit 1
		fi
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

cd $pwd
CRT='crt-1.2'
tool=$CRT
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d $GORAP/$tool ]]; then	
		cp $GORAP/$tool/bin/CRT1.2-CLI.jar $GORAP/$tool/bin/crt.jar
		if [[ $? -gt 0 ]]; then
			exit 1
		fi
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
		sed -i'' '/^BINDIR/ c\BINDIR = $(CURDIR)/bin' Makefile
		sed -i'' '/^MANDIR/ c\MANDIR = $(CURDIR)/man' Makefile
		make clean; make && make install && make clean
		if [[ $? -gt 0 ]]; then				
			exit 1
		fi
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

cd $pwd
TRNASCAN='tRNAscan-SE-1.3.1'
tool=$TRNASCAN
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d $GORAP/$tool ]]; then
		cd $GORAP/$tool
		mkdir -p bin
		mkdir -p man/man1
		mkdir -p lib/tRNAscan-SE
		sed -i'' '/^BINDIR/ c\BINDIR = $(CURDIR)/bin' Makefile 
		sed -i'' '/^MANDIR/ c\MANDIR = $(CURDIR)/man' Makefile
		sed -i'' '/^LIBDIR/ c\LIBDIR = $(CURDIR)/lib/tRNAscan-SE' Makefile
		make clean; make && make install && make clean
		if [[ $? -gt 0 ]]; then				
			exit 1
		fi
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
		make clean;	./configure --enable-threads --enable-lfs --prefix=$PWD && make && make install && make clean
		if [[ $? -gt 0 ]]; then				
			exit 1
		fi
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

cd $pwd
ZLIB='zlib-1.2.8'
tool=$ZLIB
zlib=$tool
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d "$GORAP/$tool" ]]; then
		cd $GORAP/$tool 
		rm -rf built; mkdir -p built
		make clean; ./configure --prefix=$PWD/built && make && make install && make clean
		if [[ $? -gt 0 ]]; then				
			exit 1
		fi
	else 
		echo 'Please run install_tools.sh first, then try again'		
	fi
fi

cd $pwd
NCURSES='ncurses-6.0'
tool=$NCURSES
ncurses=$tool
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d "$GORAP/$tool" ]]; then
		cd $GORAP/$tool
		rm -rf built; mkdir -p built
		make clean; ./configure --prefix=$PWD/built && make && make install && make clean
		if [[ $? -gt 0 ]]; then
			exit 1
		fi
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi

cd $pwd
SAMTOOLS='samtools-0.1.19'
tool=$SAMTOOLS
if [[ $tool == $retool ]] || [[ $retool == 'all' ]]; then
	if [[ -d "$GORAP/$tool" ]]; then
		cd $GORAP/$tool
		cfl='-fPIC '
		ldfl=''
		if [[ ! -e /usr/include/zlib.h ]]; then 
			cfl="-I$GORAP/$zlib/built/include "
			ldfl="-L$GORAP/$zlib/built/lib "
		fi
		if [[ ! -e /usr/include/ncurses.h ]]; then 
			cfl="$cfl-I$GORAP/$ncurses/built/include -I$GORAP/$ncurses/built/include/ncurses"
			ldfl="$ldfl-L$GORAP/$ncurses/built/lib"
		fi
		make clean; sed -i'' -r 's/^\s*CFLAGS\s*=\s*/override CFLAGS += /' Makefile && sed -i'' 's/-lcurses/-lncurses/' Makefile && make CFLAGS="$cfl $ldfl"
		if [[ $? -gt 0 ]]; then
			exit 1
		fi 
		mkdir -p bin 
		cp samtools bin/samtools
	else 
		echo 'Please run install_tools.sh first, then try again'
	fi
fi
