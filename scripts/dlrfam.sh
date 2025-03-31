#! /usr/bin/env bash
# (c) Konstantin Riege

source "$(dirname "$0")/../bashbone_lite.sh" -a "$@" || exit 1

############################################

usage(){
	cat <<- EOF
		DESCRIPTION
		$(basename $0) downloads most recent Rfam and NCBI taxonomy

		VERSION
		1.0.0

		INPUT OPTIONS
		-h | --help               : prints this message
		-t | --threads [value]    : threads - predicted default: $threads
		-o | --out [path]         : output directory - default: $PWD

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF

	return 1
}

alias dl="_bashbone_wrapper dl"
function dl(){
	mkdir -p $PWD/data/

	wget -O - "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/Rfam.fa.gz" | gzip -kcd | grep -oE '^>[^/]+' | sed -E 's/^>//;s/\.[0-9]+$//' | psort -u > $PWD/data/acc

	wget -O - "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/rfamseq.txt.gz" | gzip -kcd | cut -f 1,4 | sed -E 's/\.[0-9]+\t/\t/' | grep -Fw -f $PWD/data/acc | psort -k2,2n > $PWD/data/acc2tax
	wget -O - "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/taxonomy.txt.gz" | gzip -kcd | awk -F '\t' -v OFS='\t' '{print $1,$(NF-1)}' | grep -Fw -f <(cut -f 2 $PWD/data/acc2tax) | sed -E 's/\t[[({]([[:alnum:]]+)[]})]/\t\1/; s/([^[:alnum:]]).*\1$//; s/[^[:alnum:]]*[[({].*[]})]$//; s/[^[:alnum:]]+$//; s/[^[:alnum:]_-]+/_/2g; s/_+/_/g' | psort -k1,1n > $PWD/data/tax2sci
	# match: 1234 [Name]_name
	# match: 1234 Name_name_(foo)
	# match: 1234 Name_name_'foo'
	# match: 1234 Name_name_'
	# match: 1234 Name_name_sp.
	awk -v OFS='\t' -v f=$PWD/data/tax2sci 'BEGIN{getline < f; tax2=$1; sci=$2}{acc=$1;tax1=$2; while(tax2<tax1){getline < f; tax2=$1; sci=$2} print acc,tax1,sci}' $PWD/data/acc2tax | psort -k1,1V > $PWD/data/acc2tax2sci


	wget -O - "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/full_region.txt.gz" | gzip -kcd | awk -v OFS='\t' '{print $5,$1,$2}' | sed -E 's/\.[0-9]+$//' | grep -Fw -f $PWD/data/acc | psort -k2,3V -k1,1g | uniq -f 1 | awk -v OFS='\t' '{print $3,$2,$1}' | psort -k1,2V > $PWD/data/acc2fam2bs

# TODO some sci names empty RF03878
	awk -v OFS='\t' -v f=$PWD/data/acc2tax2sci 'BEGIN{getline < f; acc2=$1; tax=$2; sci=$3}{acc1=$1; fam=$2; bs=$3; while(acc1!=acc2){getline < f; acc2=$1; tax=$2; sci=$3} print acc1,tax,sci,fam,bs}' $PWD/data/acc2fam2bs | psort -k4,4V -k5,5g -k1,1V > $PWD/data/rfam.tsv

	# tsv
	# attention incomplete! e.g. RF00041 is missing, because there are no full alignments. all sequences are seeds. ftp does not provide seed bitscores
	mkdir -p $PWD/data/rfam
	perl -salne 'if ($f ne $F[3]){close OUT; open OUT,">","$o/$F[3].tsv" or die $!; $f=$F[3]}; print OUT $_' -- -o="$PWD/data/rfam" $PWD/data/rfam.tsv

	# CM
	mkdir -p $PWD/data/rfam
	wget -O - "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.tar.gz" | gzip -kcd | tar -x -C $PWD/data/rfam

	# ncbi + silva
	mkdir -p $PWD/data/taxonomy
	wget -O - "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" | gzip -kcd | tar -x -C $PWD/data/taxonomy
	# index
	perl -M'Bio::DB::Taxonomy' -se 'Bio::DB::Taxonomy->new(-source => "flatfile",-nodesfile => $x, -namesfile => $y, -directory => $d, -force => 1, -verbose => -1)' -- -x=$PWD/data/taxonomy/nodes.dmp -y=$PWD/data/taxonomy/names.dmp -d=$PWD/data/taxonomy
	wget -O - "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/tax_slv_ssu_138.2.tre.gz" | gzip -kcd > $PWD/data/taxonomy/silva.newick

	# seed stk and fasta
	# may use Bio::Index::Stockholm and $id->fetch_aln() instead of Bio::AlignIO::Stockholm
	# to alter ids: don't do $s->id("new") this is not a valid setter here!
	# => use $aln->displayname("old","new")
	# stockholm peculiarity: id is ABC.1/124-456 but $s->id or $s->display_name gives ABC only which is not a valid id to use in $aln->displayname()
	# in contrast to write_aln $s->id cannot be changed via $aln->set_displayname_flat/set_displayname_normal
	# => use $s->get_nse() to get unique seqname/start-end formatted key to be used as "old"
	# see Bio::LocatableSeq <- Bio::Align::AlignI->next_aln->Bio::SimpleAlign->each_seq->Bio::PrimarySeq-> (in other cases: Bio::Seq or Bio::SeqFeatureI)
	# further note: aln->id == name , $aln->accession = RF0000



	wget -O - "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz" | gzip -kcd | perl -M'Bio::AlignIO' -sE 'BEGIN{open F,"<",$f or die $!; while(<F>){chomp; @F=split; $acc2tax{$F[0]}=$F[1]; $acc2sci{$F[0]}=$F[2]}} $aio=Bio::AlignIO->new(-format => "stockholm", -fh => \*STDIN, -verbose => -1); while($aln = $aio->next_aln){ open FA,">","$o/".$aln->accession.".fa" or die $!; open STK,">","$o/".$aln->accession.".stk" or die $!; $str=""; open $STKFD,">",\$str; %n=(); $aln->set_displayname_flat; for $s ($aln->each_seq){$newid=$acc2sci{$s->id}.".".$acc2tax{$s->id}; $n{$newid}++; $no=$n{$newid}; $newid.=".$no"; say FA ">$newid ".$s->id."\n".uc(($s->seq=~s/[^acgtnACGTNU]//gr)=~tr/Uu/Tt/r); $aln->displayname($s->get_nse,"$newid")} (Bio::AlignIO->new(-format  => "stockholm", -fh => $STKFD, -verbose => -1))->write_aln($aln); @stk = split/\n/,$str; for (@stk){say STK unless /^#=GS/}; close FA; close STK; close $STKFD; }' -- -f=/u/people/kriege/workspace/bash/gorap/data/acc2tax2sci -o=/u/people/kriege/workspace/bash/gorap/data/rfam




# TODO ids -> sciname.taxid.repl and repl continue from full.fa to add seeds there
	# add from full fa those sequences listed in tsv to seed fasta
	mkdir -p $PWD/data/rfam
	for f in $(curl -s https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/ | grep -oE 'RF[0-9]+\.fa\.gz' | psort -Vu); do
		wget -O - "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/$f" | gzip -kcd | sed -E 's@^>(.+)\.[0-9]+/.*@\1@;t;s/(.*)/\U\1/' | paste - - | psort -u | perl -slane 'BEGIN{open F,"<","$f" or die $!; while(<F>){@F=split; $acc2sci{$F[0]}=$F[2]; $acc2tax{$F[0]}=$F[1];} close F;} $n{$acc2tax{$F[0]}}++; $no=$n{$acc2tax{$F[0]}}; print ">$acc2tax{$F[0]}.$no $acc2sci{$F[0]} $F[0]\n$F[1]" if exists $acc2sci{$F[0]}' -- -f="$PWD/data/acc2tax2sci" > "$PWD/data/rfam/$(basename $f .fa.gz).full.fa"
	done

	# grep what is missing in full tsv (should be seed only), run cmaling to add scores
	for f in data/rfam/*.full.fa; do
		cut -f 1 "${f%.full.fa}.tsv" 2> /dev/null | grep -vFw -f - <(cat "$f" "${f%.full.fa}.fa" | paste - -) | sed 's/\t/\n/;s/ /./g' > tmp; [[ -s tmp ]] || continue; cmalign -g --mxsize 10000 --cpu 32 --noprob -o /dev/null "${f%.full.fa}.cm" tmp | perl -slane 'next if /^#/; @o=split/\./,$F[1]; $acc=pop @o; $tax=shift @o; shift @o; $sci=join".",@o; print join"\t",($F[6],$acc,$tax,$sci,$i)' -- -i="$(basename "$f" .full.fa)" | psort -k2,3V -k1,1g | uniq -f 1 | awk -v OFS='\t' '{print $2,$3,$4,$5,$1}' > ${f%.full.fa}.tsv
		rm -f tmp
	done

	# finally add family name and types to tsv
	for f in data/rfam/*.tsv; do
		types="$(grep -E '^#=GF (ID|TP)' "${f%.*}.stk" | sed -nE 's/^\S+\s+\S+\s+//;s/\s+//g;s/;$//;1h;1!H;${x;s/\n/;/gp}')"
		awk -v n="$types" -v OFS='\t' '{print $0,n}' "$f" > tmp
		mv tmp "$f"
	done
}

update_cfg(){
	# TODO mach fuer alle constrains!!!
	head -1 RF00016.tsv | awk '{print "accession="$4; print "type="$NF}' > new.cfg
	if [[ -s /u/people/kriege/workspace/bash/gorap/test/RF00016.toupdate.cfg ]]; then
		grep -qF CD-box new.cfg && { grep -m 1 -qF '|' /u/people/kriege/workspace/bash/gorap/test/RF00016.toupdate.cfg || echo "consider CD-box constrains"; }
		grep -qF miRNA new.cfg && { grep -m 1 -qF '|' /u/people/kriege/workspace/bash/gorap/test/RF00016.toupdate.cfg || echo "consider miRNA stem constrains"; }
		grep conserved /u/people/kriege/workspace/bash/gorap/test/RF00016.toupdate.cfg | sed 's/^/>/;s/=/\n/' | cmalign --mxsize 20000 --cpu 32 --noprob --mapali RF00016.stk RF00016.cm - | tail -n 4 | head -3 | awk '/^#/{print sprintf("%-9s",$2),$3}!/^#/{print sprintf("%-9s",$1),toupper($2)}' | cat - <(tail -n +3 /u/people/kriege/workspace/bash/gorap/test/RF00016.toupdate.cfg | sed -E 's/( |=)(\S+)/ \U\2/') | perl -M'List::Util qw(max)' -slanE 'print $_ if $debug; push @n,$F[0]; push @l,[(split//,$F[-1])]; $len=max($len,length($F[-1])-1); END{ print "--" if $debug; $i=0; for(0..$len){ if($l[4]->[$i] ne $l[0]->[$i]){ splice $l[$_]->@*,$i,0,"." for 3..$#l} $i++} if($debug){ sprintf("%-9s",$n[$_])." ".join("",$l[$_]->@*) for 0..$#l; print "--"} $n[1]=$n[3]; $n[2]=$n[4]; if($debug){ print $n[$_]." ".join("",$l[$_]->@*) for (1,2,5..$#l); print "--"} for my $i (0..$len){ if((join "",map { $l[$_]->[$i] } (1,2,5..$#l)) =~ /^\.+$/){ delete $l[$_]->@[$i] for (1,2,5..$#l)}} print $n[$_]."=".join("",$l[$_]->@*) for (1,2,5..$#l) }' -- -debug=0 >> new.cfg
	else
		grep -qF CD-box new.cfg && echo "consider CD-box constrains"
		grep -qF miRNA new.cfg && echo "consider miRNA stem constrains"
		tail -n 3 RF00016.stk | head -2 | perl -M'List::Util qw(max)' -lanE 'push @l,[(split//,$F[-1])]; $len=max($len,length($F[-1])-1); END{ push @l,[(".")x$len]; @n=("conserved","structure","constrain"); for my $i (0..$len){ if((join "",map { $l[$_]->[$i] } (0..2)) =~ /^\.+$/){ delete $l[$_]->@[$i] for (0..2)}} print $n[$_]."=".join("",$l[$_]->@*) for (1,0,2) }' >> new.cfg
	fi
}
