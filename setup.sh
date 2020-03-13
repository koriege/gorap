#! /usr/bin/env bash
# (c) Konstantin Riege
trap 'die' INT TERM
trap 'sleep 1; kill -PIPE $(pstree -p $$ | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+") &> /dev/null' EXIT

die() {
	[[ $* ]] && echo ":ERROR: $*" || echo ":ERROR: failed"
	exit 1
}

progressbar() {
	local mod=0
	while true; do
		((++mod))
		case $mod in
			[15]) echo -en "\r|";;
			[26]) echo -en "\r/";;
			[37]) echo -en "\r-";;
			4) echo -en "\r\\";;
			8) echo -en "\r\\"; mod=0;;
		esac
		sleep 0.2
	done

	return 0
}

progresslog() {
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage:
			-v [0|1|2] | verbosity level
			-o <file>  | path to
		EOF
		return 0
	}
	local OPTIND arg mandatory log verbosity
	while getopts 'v:o:' arg; do
		case $arg in
			v)	((++mandatory)); verbosity=$OPTARG;;
			o)	((++mandatory)); log="$OPTARG";;
			*)	_usage;	return 1;;
		esac
	done
	[[ $mandatory -lt 2 ]] && { _usage; return 1; }

	case $verbosity in
		0)	progressbar &
			{ tail -f $log 2>&1 | grep -E --line-buffered '^\s*(:INFO:|:ERROR:|:BENCHMARK:|:WARNING:)'; } &
			;;
		1)	progressbar &
			{ tail -f $log 2>&1 | grep -E --line-buffered '^\s*(:INFO:|:CMD:|:ERROR:|:BENCHMARK:|:WARNING:)'; } &
			;;
		2)	{ tail -f $log 2>&1; } &;;
		*)	_usage;	return 1;;
	esac
	
	return 0
}

############### GLOBAL VARS ###############

export SRC=$(readlink -e $(dirname $0))
VERSION=$(perl -e 'use File::Spec::Functions; BEGIN{ unshift @INC,"$ENV{SRC}/lib"}; use Bio::Gorap::Gorap; print Bio::Gorap::Gorap->VERSION')
THREADS=$(cat /proc/cpuinfo | grep -cF processor)
VERBOSITY=0
INSTALL=false
UPGRADE=false
INSDIR=''
declare -a TOOLS
for t in $(grep -oE '^install_[^[:space:]\(]+' $0); do
	t=$(cut -d '_' -f 2 <<< $t)
	TOOLS+=($t)
done

############### USAGE ###############

usage() {
	cat <<- EOF
		DESCRIPTION
		Setup Gorap, the genome wide ncRNA anntation pipeline

		VERSION
		$VERSION

		SYNOPSIS
		$(basename $0) -d [path]

		OPTIONS
		-h | --help            # print this message
		-d | --dir [path]      # installation path
		-t | --threads [num]   # threads to use for comilation - default: $THREADS
		-v | --verbosity [num] # verbosity level (0|1|2)
		-u | --update          # update Gorap data bases

		REFERENCES
		(c) Konstantin Riege
		konstantin{.}riege{a}leibniz-fli{.}de
	EOF

	exit 0
}

checkopt() {
	local arg=false
	case $1 in
		-h | --h | -help | --help) usage;;
		-u | --u | -update | --update) TOOLS=('db');;
		-t | --t | -threads | --threads) arg=true; THREADS=$2;;
		-v | --v | -verbosity | --verbosity) arg=true; VERBOSITY=$2;;
		-l | --l | -list | --list) arg=true; TOOLS=(); mapfile -d ',' -t <<< $2; for t in "${MAPFILE[@]}"; do TOOLS+=($(printf '%s' $t)); done;;
		-d | --d | -dir | --dir) arg=true; INSDIR=$2; [[ $(mkdir -p $INSDIR &> /dev/null; echo $?) -gt 0 ]] && die 'check your installation path';;
		-*) die "illegal option $1";;
		*) die "illegal option $2";;
	esac
	$arg && {
		[[ ! $2 ]] && die "argument missing for option $1"
		[[ "$2" =~ ^- ]] && die "illegal argument $2 for option $1"
		return 0
	} || {
		[[ $2 ]] && [[ ! "$2" =~ ^- ]] && die "illegal argument $2 for option $1"
		return 0
	}
}

############### FUNCTIONS ###############

makeclean() {
	rm -rf built
	[[ $1 ]] && make clean -f $1 &> /dev/null || make clean &> /dev/null
	return 0
}

install_conda() {
	local url version
	{	rm -rf $INSDIR/conda && \
		mkdir -p $INSDIR/conda && \
		cd $INSDIR/conda && \
		url='https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh' && \
		wget -q --show-progress --progress=bar:force --waitretry 1 --tries 5 --retry-connrefused -N -O miniconda.sh $url && \
		bash miniconda.sh -b -f -p $PWD && \
		return 0
	} || return 1
}

install_conda-env() {
	# use python 2 env due to bcheck
	{	source $INSDIR/conda/bin/activate base && \
		conda env remove -y -n gorap
		conda create -y -n gorap python=2 && \
		conda install -y -n gorap --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
			gcc_linux-64 readline make automake xz zlib bzip2 pigz pbzip2 unzip ncurses htslib \
			perl perl-threaded perl-dbi perl-app-cpanminus perl-bioperl perl-bio-eutilities perl-moose perl-bio-db-sam perl-postscript \
			perl-archive-extract perl-list-moreutils perl-try-tiny perl-math-round perl-hash-merge perl-test-more perl-extutils-makemaker \
			perl-file-temp raxml openjdk mafft trnascan-se=2.0.0 hmmer2 infernal barrnap blast rnabob && \
		conda clean -y -a && \
		return 0
	} || return 1
}

install_njplot() {
	{	source $INSDIR/conda/bin/activate gorap && \
		rm -rf $INSDIR/njplot/2.3 && \
		mkdir -p $INSDIR/njplot/2.3/bin && \
		cd $INSDIR/njplot/2.3 && \
		url='ftp://pbil.univ-lyon1.fr/pub/mol_phylogeny/njplot/newicktopdf' && \
		wget -q --show-progress --progress=bar:force --waitretry 1 --tries 5 --retry-connrefused -N -O bin/newicktopdf $url && \
		chmod 755 bin/newicktopdf
		return 0
	} || return 1
}

install_crt() {
	{	source $INSDIR/conda/bin/activate gorap && \
		rm -rf $INSDIR/crt/1.2/ && \
		mkdir -p $INSDIR/crt/1.2/bin && \
		cd $INSDIR/crt && \
		url='http://www.room220.com/crt/CRT1.2-CLI.jar.zip' && \
		wget -q --show-progress --progress=bar:force --waitretry 1 --tries 5 --retry-connrefused -N -O crt.zip $url && \
		unzip -q -o -d 1.2 crt.zip && \
		rm -f crt.zip && \
		touch 1.2/bin/crt && \
		chmod 755 1.2/bin/crt
	} || return 1

	cd $INSDIR/crt/1.2/bin
	cat <<- 'EOF' > crt
		#!/usr/bin/env bash
		exec java -cp $GORAP/crt/1.2/CRT1.2-CLI.jar crt $*
	EOF
	return 0
}

install_bcheck() {
	{	source $INSDIR/conda/bin/activate gorap && \
		rm -rf $INSDIR/bcheck/0.6 && \
		mkdir -p $INSDIR/bcheck/0.6 && \
		cd $INSDIR/bcheck/0.6 && \
		url='http://rna.tbi.univie.ac.at/bcheck/Bcheck.tgz' && \
		wget -q --show-progress --progress=bar:force --waitretry 1 --tries 5 --retry-connrefused -N -O bcheck.tar.gz $url && \
		tar -xzf bcheck.tar.gz && \
		rm -f bcheck.tar.gz
		mv Bcheck-0.6 bin && \
		sed -i "s/if bob_version != '2.1'/if bob_version != '2.2.1'/" bin/Bcheck && \
		return 0
	} || return 1
}

install_bcheck-infernal() {
	{	source $INSDIR/conda/bin/activate gorap && \
		rm -rf $INSDIR/infernal/1.0 && \
		mkdir -p $INSDIR/infernal && \
		cd $INSDIR/infernal && \
		url='http://eddylab.org/software/infernal/infernal-1.0.tar.gz' && \
		wget -q --show-progress --progress=bar:force --waitretry 1 --tries 5 --retry-connrefused -N -O infernal1.tar.gz $url && \
		tar -xzf infernal1.tar.gz && \
		rm -f infernal1.tar.gz && \
		mv infernal-1.0 1.0 && \
		cd 1.0 && \
		./configure --prefix=$PWD && \
		makeclean && \
		make -j $THREADS && \
		make install && \
		return 0
	} || return 1
}

install_jquery(){
	{	source $INSDIR/conda/bin/activate gorap && \
		rm -rf $INSDIR/jquery && \
		mkdir -p $INSDIR/jquery && \
		cd $INSDIR/jquery && \
		url='https://mottie.github.io/tablesorter/js/jquery.tablesorter.js' && \
		wget -q --show-progress --progress=bar:force --waitretry 1 --tries 5 --retry-connrefused -N -O jquery.tablesorter.js $url && \
		url=$(curl -s https://jquery.com/download/ | grep -m 1 -F 'Download the compressed' | grep -oE 'http[^"]+') && \
		wget -q --show-progress --progress=bar:force --waitretry 1 --tries 5 --retry-connrefused -N -O jquery.js $url && \
		return 0
	} || return 1
}

install_gorap(){
	{	source $INSDIR/conda/bin/activate gorap && \
		rm -rf $INSDIR/db $INSDIR/bin && \
		mkdir -p $INSDIR/db $INSDIR/bin && \
		cd $SRC && \
		tar -xzf $(ls -v data-*.tar.gz | tail -1) -C $INSDIR/db && \
		cpanm --reinstall . && \
		makeclean && \
		rm -f Makefile.old && \
		touch $INSDIR/bin/gorap && \
		chmod 755 $INSDIR/bin/gorap
	} || return 1

	for f in $(find $INSDIR/conda/envs/gorap/lib/ -type f -name 'SimpleAlign.pm'); do
		chmod 644 $f
		sed -i '/contains no residues/d' $f
		sed -i '/Sequence excluded/d' $f
		sed -i '/Replacing one sequence/d' $f
	done

	cd $INSDIR/bin
	cat <<- 'EOF' >> gorap
		#!/usr/bin/env bash
		[[ -z $GORAP ]] && {
		    echo "Export GORAP environment variable pointing towards installation directory and try again!"
		    exit 1
		}
		source $GORAP/conda/bin/activate gorap 
		Gorap.pl $*
	EOF
}

############### MAIN ###############

[[ $# -eq 0 ]] && usage
[[ $# -eq 1 ]] && [[ ! $1 =~ ^- ]] && die "illegal option $1"
for i in $(seq 1 $#); do
	if [[ ${!i} =~ ^- ]]; then
		j=$((i+1))
		checkopt "${!i}" "${!j}" || die
	else 
		((++i))
	fi
done

[[ $INSDIR ]] || die "mandatory parameter -d missing"
INSDIR=$(readlink -e $INSDIR)

LOG=$INSDIR/install.log

echo ":INFO: installation started. please be patient." > $LOG || die "cannot access $LOG"
progresslog -v $VERBOSITY -o $LOG

for i in "${TOOLS[@]}"; do
	echo ":INFO: installing $i" >> $LOG
	install_$i >> $LOG 2> >(tee -a $LOG >&2) || die 
done

echo ":INFO: success" >> $LOG
exit 0
