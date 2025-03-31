#! /usr/bin/env bash
# (c) Konstantin Riege

function annotate::_align(){
	function _usage() {
		commander::print {COMMANDER[0]}<<- EOF
			annotate::infernal runs infernal
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-g <genome>   | path to
			-f <cfgfile>  | path to
			-b <bed>      | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory threads genome tmpdir="${TMPDIR:-/tmp}" cfg bed execute=true maxmemory
	declare -n _cmds1_align
	while getopts '1:t:g:f:b:M:' arg; do
		case $arg in
			1) execute=false; _cmds1_align=$OPTARG;;
			t) ((++mandatory)); threads=$OPTARG;;
			g) ((++mandatory)); genome="$OPTARG";;
			f) ((++mandatory)); cfg="$OPTARG";;
			b) ((++mandatory)); bed="$OPTARG";;
			M) maxmemory="$OPTARG";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage

	if $execute; then
		unset _cmds1_align
		declare -a _cmds1_align
	fi
	local instances memory
	read -r instances memory < <(configure::memory_by_instances -i 1 -M "$maxmemory")

	# no need to convert T -> U
	commander::makecmd -a _cmds1_align -s '|' -c {COMMANDER[0]}<<-CMD {COMMANDER[1]}<<-'CMD' {COMMANDER[2]}<<-CMD {COMMANDER[3]}<<-CMD
		bedtools getfasta -bedOut -s -fi "$genome" -bed "$bed"
	CMD
		perl -lane 'print ">GORAP.$F[3]"; print $F[-1]=~tr/acgt/ACGT/r'
	CMD
		tee -i "${bed%.*}.fa"
	CMD
		cmalign --mxsize $memory --cpu $threads --noprob --sfile "${bed%.*}.scores" --mapali "${cfg%.*}.stk" -o "${bed%.*}.stk" "${cfg%.*}.cm" -
	CMD
	# esl-alimerge --rna -o "${bed%.*}.full.stk" "${cfg%.*}.stk" "${bed%.*}.stk"

	if $execute; then
		commander::runcmd -c infernal -v -b -t 1 -a _cmds1_align
	fi
}

function annotate::infernal(){
	function _usage() {
		commander::print {COMMANDER[0]}<<- EOF
			annotate::infernal runs infernal
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-g <genome>   | path to
			-o <outdir>   | path to
			-r <cfgfiles> | array of
			-M <maxmem>   | amount of
			-F            | force infernal
		EOF
		return 1
	}

	local OPTIND arg mandatory threads genome outdir tmpdir="${TMPDIR:-/tmp}" skip=false force=false maxmemory
	declare -n _cfgfiles_infernal
	while getopts 'S:s:t:g:o:r:M:F' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			g) ((++mandatory)); genome="$OPTARG";;
			o) ((++mandatory)); outdir="$OPTARG/infernal"; mkdir -p "$outdir";;
			r) ((++mandatory)); _cfgfiles_infernal=$OPTARG;;
			M) maxmemory=$OPTARG;;
			F) force=true;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage

	declare -a cmds1 cmds2

	# TODO
	local f v params="-E 0.01"
	[[ $cutoff ]] && params+=" -T $cutoff" || params+=" --cut_nc"
	# TODO length filter !!! sto-sta > cmstat/2 bzw in output header is "Query:  HgcF  [CLEN=168]"
	# TODO user filter

	for f in "${_cfgfiles_infernal[@]}"; do
		# declare variables CFG_id CFG_bitscore ...
		source <(sed -nE '/^\w/{s/^/declare -a CFG_/;s/=.*/=\(\)/p}' "$f")
		# assign variables
		source <(sed -nE '/^\w/{s/^/CFG_/;s/=/+=\("/;s/$/"\)/p}' "$f")

		! $force && [[ "${CFG_tool[*],,}" != *infernal* ]] && continue
		commander::makecmd -a cmds1 -s '|' -c {COMMANDER[0]}<<-CMD {COMMANDER[1]}<<-'CMD' {COMMANDER[2]}<<-CMD {COMMANDER[3]}<<-CMD {COMMANDER[4]}<<-CMD
			cmsearch
				--noali
				--cpu $threads
				"${f%.*}.cm"
				"$genome"
		CMD
			perl -M'List::Util qw(min max)' -lane 'print STDERR; next if /(^#|^\s*$)/ || $#F < 11; print join "\t",($F[5],min($F[6],$F[7])-1,max($F[6],$F[7]),".",$F[3],$F[8])'
		CMD
			helper::sort -t $threads -M "$maxmemory"
		CMD
			bedtools merge
				-i -
				-s
				-c 4,5,6
				-o distinct,max,distinct
		CMD
			awk -v OFS='\t' -v n=$CFG_name '{\$4="$CFG_name"":"NR; print}' > "$outdir/$CFG_id.full.bed"; annotate::_align -t $threads -g "$genome" -f "$f" -b "$outdir/$CFG_id.full.bed"
		CMD
	done

	if $skip; then
		commander::printcmd -a cmds1
	else
		commander::runcmd -c infernal -v -b -t 1 -a cmds1
		commander::runcmd -c infernal -v -b -t 1 -a cmds2
	fi

	return 0
}
