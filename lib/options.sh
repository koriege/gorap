#! /usr/bin/env bash
# (c) Konstantin Riege

function options::usage() {
	commander::print {COMMANDER[0]}<<- EOF
		DESCRIPTION
		GORAP

		VERSION
		$VERSION
		utilizing bashbone $BASHBONE_VERSION

		SYNOPSIS
		gorap.sh -g genome.fa

		OPTIONS
		-h       | --help                     : prints this message
		-v       | --verbosity [value]        : set level of verbosity. default: 0
		                                        0 - get simple status updates
		                                        1 - get status updates and commands
		                                        2 - get full output
		-d       | --download                 : download and setup latest databases to $DATADIR (takes a while!)
		-o       | --out [path]               : output directory. default: $OUTDIR
		-l       | --log [path]               : output directory. default: $OUTDIR/run.log
		-tmp     | --tmp                      : temporary directory. default: $TMPDIR/rippchen.XXXXXXXXXX
		-k       | --keep                     : keep temporary and unnecessary files
		-t       | --threads [value]          : number of threads. default: $THREADS

		-g       | --genome [path]            : genome fasta input. without, only preprocessing is performed (see dlgenome.sh)
		-gtf     | --gtf [path]               : annotation gtf input. default: [-g].gtf (see dlgenome.sh)

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	exit 1
}

function options::checkopt (){
	local arg=false
	case $1 in
		-h        | --help) options::usage || exit 0;;

		-tmp      | --tmp) arg=true; TMPDIR="$2";;
		-k        | --keep) CLEANUP=false;;
		-d        | --download) DOWNLOAD=true;;
		-v        | --verbosity) arg=true; VERBOSITY=$2;;
		-t        | --threads) arg=true; THREADS=$2;;

		-g        | --genome) arg=true; GENOME="$2";;
		-gtf      | --gtf) arg=true; GTF="$2";;
		-o        | --out) arg=true; OUTDIR="$2";;
		-l        | --log) arg=true; LOG="$2";;

		-*) commander::printerr "illegal option $1"; return 1;;
		*) commander::printerr "illegal option $2"; return 1;;
	esac

	$arg && {
		[[ ! $2 ]] && commander::printerr "argument missing for option $1" && return 1
		[[ "$2" =~ ^- ]] && commander::printerr "illegal argument $2 for option $1" && return 1
		return 0
	} || {
		[[ $2 ]] && [[ ! "$2" =~ ^- ]] && commander::printerr "illegal argument $2 for option $1" && return 1
		return 0
	}
}
