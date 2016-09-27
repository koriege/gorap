#! /bin/bash

if [[ $1 == 'on' ]]; then
	if [[ -e $GORAPLIB/lib/perl5/Bio/SimpleAlign.pm ]]; then
		chmod 644 $GORAPLIB/lib/perl5/Bio/AlignIO/stockholm.pm
		chmod 644 $GORAPLIB/lib/perl5/Bio/SimpleAlign.pm

		cp $GORAPLIB/lib/perl5/Bio/SimpleAlign.pm $GORAP/SimpleAlign.pm_backup
		cp $GORAPLIB/lib/perl5/Bio/AlignIO/stockholm.pm $GORAP/stockholm.pm_backup		
		sed -i 's/%-4s%-22s%-3s%s/%-4s%-100s%-3s%s/g' $GORAPLIB/lib/perl5/Bio/AlignIO/stockholm.pm		
		sed -i '/contains no residues/d' $GORAPLIB/lib/perl5/Bio/SimpleAlign.pm
		sed -i '/Sequence excluded/d' $GORAPLIB/lib/perl5/Bio/SimpleAlign.pm
		sed -i '/Replacing one sequence/d' $GORAPLIB/lib/perl5/Bio/SimpleAlign.pm
	fi
else
	if [[ -e $GORAP/SimpleAlign.pm_backup ]] && [[ -e $GORAPLIB/lib/perl5/Bio/SimpleAlign.pm ]]; then
		cp $GORAP/SimpleAlign.pm_backup $GORAPLIB/lib/perl5/Bio/SimpleAlign.pm
		cp $GORAP/stockholm.pm_backup $GORAPLIB/lib/perl5/Bio/AlignIO/stockholm.pm
	fi
fi
