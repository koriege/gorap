#! /bin/bash
sed -i 's/%-4s%-22s%-3s%s/%-4s%-100s%-3s%s/g' $GORAPLIB/lib/perl5/Bio/AlignIO/stockholm.pm
sed -i '/contains no residues/d' $GORAPLIB/lib/perl5/Bio/SimpleAlign.pm
sed -i '/Sequence excluded/d' $GORAPLIB/lib/perl5/Bio/SimpleAlign.pm
sed -i '/Replacing one sequence/d' $GORAPLIB/lib/perl5/Bio/SimpleAlign.pm