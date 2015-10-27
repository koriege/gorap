#!/bin/bash
retool=$1

if [ "$retool" = "" ]
	then
	
		echo 'Do you forgot to add a tool name as parameter? Or'
		read -p 'do you want to recompile all tools [y|n] ' in
		if [ "$in" = "y" -o "$in" = "Y" ]
			then
				retool='all'
			else 
				read -p 'Please enter the correct tool name as used in <~/.gorap/[toolname-version]> or hit return to exit' retool
				if [ "$retool" = "" ]
					then
						echo 'Recompilation aborted...'
						exit
				fi
		fi
fi

pwd=`pwd`
tools=$GORAP
mkdir -p  $tools

download () {
	if [ ! -d "$tool" -a "$ex" -eq "0" ]
		then
			echo 'Downloading '$tool	
			wget -T 10 -N www.rna.uni-jena.de/supplements/gorap/$tool.tar.gz
			echo 'Extracting '$tool
			tar -xzf $tool.tar.gz		
			rm $tool.tar.gz
	fi
}

tool='mafft-7.017-with-extensions'
if [ "$tool" = "$retool" ] || [ "$retool" = "all" ]
	then	
		if [ -d "$tools/$tool" ]
			then
				cd $tools/$tool/core
				make 
				make install
				make clean
			else 
				echo 'please run install_tools.sh first, then try again'
		fi
fi

tool='infernal-1.1rc2'
if [ "$tool" = "$retool" ] || [ "$retool" = "all" ]
	then	
		if [ -d "$tools/$tool" ]
			then
				cd $tools/$tool
				./configure --prefix=`pwd`
				make
				make install
				make clean

				glibc=$(ldd --version | head -n 1 | awk '{split($NF,a,"."); print a[2]}')
				if [ "$glibc" -lt "14" ]; then
					tool='glibc-2.21'
					ex=0
					download
					
					cd $tools/$tool
					mkdir -p build
					cd build
					../configure --prefix=$tools/infernal-1.1rc2/easel
					make
					make install
					make clean

					cd $tools/infernal-1.1rc2/easel
					cp -r include/* .
					
					./configure --prefix=$tools/infernal-1.1rc2
					sed -i '/^CFLAGS/ c\CFLAGS = -O3 -fomit-frame-pointer -fstrict-aliasing -pthread -fPIC -I. -I..' Makefile							
					make
					make install
					make clean
				else
					cd $tools/infernal-1.1rc2/easel
					./configure --prefix=$tools/infernal-1.1rc2
					make 
					make install 
					make clean
				fi

			else 
				echo 'please run install_tools.sh first, then try again'
		fi
fi

tool='infernal-1.0'
if [ "$tool" = "$retool" ] || [ "$retool" = "all" ]
	then	
		if [ -d "$tools/$tool" ]
			then
				cd $tools/$tool
				./configure --prefix=`pwd` 
				make
				make install
				make clean
			else 
				echo 'please run install_tools.sh first, then try again'
		fi
fi

tool='rnabob-2.2'
if [ "$tool" = "$retool" ] || [ "$retool" = "all" ]
	then	
		if [ -d "$tools/$tool" ]
			then
				cd $tools/$tool
				make
				make install
				make clean
			else 
				echo 'please run install_tools.sh first, then try again'
		fi
fi

tool='ncbi-blast-2.2.27+'
if [ "$tool" = "$retool" ] || [ "$retool" = "all" ]
	then	
		if [ -d "$tools/$tool" ]
			then
				cd $tools/$tool
				./configure --without-debug --with-strip --with-mt --with-build-root=`pwd`
				cd $tools/$tool/src
				make
				make all_r
				make clean
			else 
				echo 'please run install_tools.sh first, then try again'
		fi
fi

tool='tRNAscan-SE-1.3.1' 
if [ "$tool" = "$retool" ] || [ "$retool" = "all" ]
	then	
		if [ -d "$tools/$tool" ]
			then
				cd $tools/$tool
				make
				make install
				make clean
			else 
				echo 'please run install_tools.sh first, then try again'
		fi
fi

tool='hmmer-2.3.2'
if [ "$tool" = "$retool" ] || [ "$retool" = "all" ]
	then	
		if [ -d "$tools/$tool" ]
			then
				cd $tools/$tool
				./configure --enable-threads --enable-lfs --prefix=`pwd` 
				make
				make install
				make clean
			else 
				echo 'please run install_tools.sh first, then try again'
		fi
fi

samtool='samtools-0.1.19'
if [ "$tool" = "$retool" ] || [ "$retool" = "all" ]
	then	
		if [ -d "$tools/$samtool" ]
			then
				tool='zlib-1.2.8'		
				ex=0
				download				
				
				cd $tools/$tool 
				./configure --prefix=$tools/$samtool				
				make
				make install
				make clean
				
				tool='ncurses-5.9'
				ex=0
				download
				
				cd $tools/$tool 
				./configure --prefix=$tools/$samtool
				make
				make install
				make clean

				cd $tools/$samtool
				make clean
				make
				
				cp -r include/* .
				cp -r lib/* .
								
				mkdir -p bin 
				cp samtools bin/samtools
			else 
				echo 'please run install_tools.sh first, then try again'
		fi
fi

tool='RAxML-7.4.2'
if [ "$tool" = "$retool" ] || [ "$retool" = "all" ]
	then	
		if [ -d "$tools/$tool" ]
			then
				os=`uname -s`
				if [ "$os" == "Darwin" ]
					then
						cpu=$(sysctl -n machdep.cpu.brand_string | egrep '(i5|i7)' | wc | awk '{print $1}')		
					else
						cpu=$(grep -m 1 ^model\ name /proc/cpuinfo | egrep '(i5|i7)' | wc | awk '{print $1}')
				fi	

				if [ "$cpu" -eq "0" ]
					then
						cd $tools/$tool
						make -f Makefile.SSE3.PTHREADS.gcc 
						mkdir -p bin 
						mv raxmlHPC-PTHREADS-SSE3 bin/raxml 						
					else
						cd $tools/$tool 
						make -f Makefile.AVX.PTHREADS.gcc 
						mkdir -p bin 
						mv raxmlHPC-PTHREADS-AVX bin/raxml						
				fi
			else 
				echo 'please run install_tools.sh first, then try again'
		fi
fi
