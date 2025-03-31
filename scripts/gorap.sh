#! /usr/bin/env bash
# (c) Konstantin Riege

source "$(dirname "$(dirname "$(readlink -e "$0")")")/activate.sh" -i $BASHBONE -l true -c false -r true -x cleanup -a "$@" || exit 1

cleanup() {
	[[ -e "$LOG" ]] && {
		echo "date: $(date)" | tee -ia "$LOG"
		[[ $1 -eq 0 ]] && echo "success" | tee -ia "$LOG" || echo "failed" | tee -ia "$LOG"
	}
	[[ -e "$CLEANUP_TMPDIR" ]] && {
		find -L "$CLEANUP_TMPDIR" -type f -name "cleanup.*" -exec rm -f "{}" \; &> /dev/null || true
		find -L "$CLEANUP_TMPDIR" -depth -type d -name "cleanup.*" -exec rm -rf "{}" \; &> /dev/null || true
	}
	${CLEANUP:=true} && {
		[[ -e "$CLEANUP_TMPDIR" ]] && {
			find -L "$CLEANUP_TMPDIR" -type f -exec rm -f "{}" \; &> /dev/null || true
			find -L "$CLEANUP_TMPDIR" -depth -type d -exec rm -rf "{}" \; &> /dev/null || true
			rm -rf "$CLEANUP_TMPDIR"
		}
	}
	return 0
}

VERSION=$version
CMD="$(basename "$0") $*"
THREADS=$(grep -cF processor /proc/cpuinfo)
VERBOSITY=0
OUTDIR="$PWD/results"
TMPDIR="${TMPDIR:-$OUTDIR}"

BASHBONE_ERROR="parameterization issue"
options::parse "$@"
bashbone -c

[[ ! $LOG ]] && LOG="$OUTDIR/run.log"
BASHBONE_ERROR="cannot access $LOG"
mkdir -p "$(dirname "$LOG")"
commander::printinfo "gorap $VERSION utilizing bashbone $BASHBONE_VERSION started with command: $CMD" | tee -i "$LOG"
LOG="$(realpath -se "$LOG")"

BASHBONE_ERROR="cannot access $OUTDIR"
mkdir -p "$OUTDIR"
OUTDIR="$(realpath -se "$OUTDIR")"

BASHBONE_ERROR="cannot access $TMPDIR"
if [[ $PREVIOUSTMPDIR ]]; then
	TMPDIR="$(realpath -se "$PREVIOUSTMPDIR")"
else
	mkdir -p "$TMPDIR"
	TMPDIR="$(realpath -se "$TMPDIR")"
	TMPDIR="$(mktemp -d -p "$TMPDIR" rippchen.XXXXXXXXXX)"
fi
CLEANUP_TMPDIR="$TMPDIR"


# progress::log -v $VERBOSITY -o "$LOG" -f pipeline::annotate


cfgs=(data/rfam/RF00003/RF00003.cfg data/rfam/RF00004/RF00004.cfg)

annotate::get_thresholds

annotate::infernal -t $THREADS -o $OUTDIR -g test/reference.fa -r cfgs

