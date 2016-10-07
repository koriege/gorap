package Bio::Gorap::Evaluation::HTML;

use Moose;
use File::Spec::Functions;
use File::Basename;

sub create {
	my ($self,$parameter,$gffdb,$stkdb,$rf_rnas,$filename) = @_;
	my $data_pos = tell DATA;
	my $tablehtml=0;
	open INDEX , '>'.catfile($parameter->output,'index.html') or die $!;
	open HTML , '>'.catfile($parameter->output,'html',$filename.'.html') or die $!;
	while (<DATA>){		
		if ($_ =~ /DOCTYPE/){
			$tablehtml++;			
		}
		if ($tablehtml==1){
			if ($_=~/option value/){
				print INDEX $_;
				for (glob catfile($parameter->output,'html','*.html')){							
					my $file = basename $_;
					next if $file eq $filename.'.html';
					print INDEX '<option value="'.catfile('html',$file).'">'.substr($file,0,-5).'</option>'."\n";
				}	
				print INDEX '<option value="'.catfile('html',$filename.'.html').'">'.$filename.'</option>'."\n";		
			} else {
				print INDEX $_;	
			}
		} else {
			if ($_=~/Parameter/){
				print HTML $_;
				if ($parameter->has_file){
					open F , '<'.$parameter->file or die $!;
					while(<F>){
						chomp $_;
						print HTML $_.'<br>'."\n";
					}
					close F;
					print HTML '<br>'."\n";
					print HTML 'Gorap.pl -file '.$parameter->file.' ';
					print HTML '-ss ' if $parameter->strandspec;
					print HTML '-notax ' if $parameter->notax;
					print HTML '-sort ' if $parameter->sort;
					print HTML '-nofi ' if $parameter->nofilter;
					print HTML '-nobl ' if $parameter->noblast;
					print HTML '-noc' unless $parameter->check_overlaps;
					print HTML "\n";
				} else {	
					print HTML "Gorap.pl -i ".join(',\\<br>'."\n",@{$parameter->genomes}).'\\<br>'."\n";
					print HTML "-a ".join(',',@{$parameter->abbreviations}).'\\<br>'."\n";
					print HTML "-o ".$parameter->output.'\\<br>'."\n";
					print HTML '-c '.$parameter->threads.'\\<br>'."\n";
					print HTML '-k '.join(',',keys %{$parameter->kingdoms}).'\\<br>'."\n";
					if ($#{$parameter->queries}>-1){				
						print HTML '-q '.$parameter->querystring.'\\<br>'."\n" if $parameter->querystring;
					} else {
						print HTML '-q 0\\<br>'."\n";
					}
					print HTML '-minl '.$parameter->denovolength.'\\<br>'."\n" unless $parameter->denovolength == 50;
					print HTML '-minh '.$parameter->denovoheigth.'\\<br>'."\n" unless $parameter->denovoheigth == 1000;
					print HTML '-r '.$parameter->rank.'\\<br>'."\n" if $parameter->has_rank;
					print HTML '-s '.$parameter->species.'\\<br>'."\n" if $parameter->has_species;
					print HTML '-og '.join('\\<br>'."\n",@{$parameter->outgroups}).'\\<br>'."\n" if $parameter->has_outgroups;
					print HTML '-oga '.join(',',@{$parameter->ogabbreviations}).'\\<br>'."\n" if $#{$parameter->ogabbreviations} > -1;				
					if ($parameter->has_gffs){
						my $gffs='';
						$gffs .= join(',\\<br>',@{$_}).',\\<br>' for @{$parameter->{'gffs'}};
						$gffs = substr($gffs,0,-6);
						print HTML '-g '.$gffs.'\\<br>'."\n";
					}
					if ($parameter->has_bams){
						my $bams='';
						$bams .= join(',\\<br>',@{$_}).',\\<br>' for @{$parameter->{'bams'}};
						$bams = substr($bams,0,-6);
						print HTML '-b '.$bams.'\\<br>'."\n";
					}
					print HTML '-strand\\<br>'."\n" if $parameter->strandspec;
					print HTML '-notax\\<br>'."\n" if $parameter->taxonomy;
					print HTML '-sort\\<br>'."\n" if $parameter->sort;
					print HTML '-nofi\\<br>'."\n" if $parameter->nofilter;
					print HTML '-nobl\\<br>'."\n" if $parameter->noblast;
					print HTML '-noo\\<br>'."\n" unless $parameter->check_overlaps;
					print HTML '-t '.$parameter->{'tmp'}.'<br>'."\n";
				}
			} elsif ($_=~/Used data/){
				print HTML $_;
				print HTML '<table class="tablesorter">'."\n";
				print HTML '<thead>'."\n";
				print HTML '<tr>'."\n";
				print HTML '<th>Genome</th>'."\t".'<th>File</th>'."\n";
				print HTML '</tr>'."\n";
				print HTML '</thead>'."\n";
				print HTML '<tbody>'."\n";
				for (0..$#{$parameter->genomes}){
					print HTML '<tr>'."\n";					
					print HTML '<td>'.${$parameter->abbreviations}[$_].'</td>'."\t";
					print HTML '<td><a href="'.${$parameter->genomes}[$_].'">'.basename(${$parameter->genomes}[$_]).'</a></td>'."\n";
					print HTML '</tr>'."\n";	
				}
				for (0..$#{$parameter->outgroups}){
					print HTML '<tr>'."\n";					
					print HTML '<td>'.${$parameter->ogabbreviations}[$_].'</td>'."\t";
					print HTML '<td><a href="'.${$parameter->outgroups}[$_].'">'.basename(${$parameter->outgroups}[$_]).'</a></td>'."\n";
					print HTML '</tr>'."\n";	
				}			
				print HTML '</tbody>'."\n";				
				print HTML '</table>'."\n";
			} elsif($_=~/ncRNA annotation/) {
				print HTML $_;
				print HTML '<table class="tablesorter">'."\n";
				print HTML '<thead>'."\n";
				print HTML '<tr>'."\n";
				print HTML '<th>Genome</th>'."\t".'<th>GFF</th>'."\t".'<th>FASTA</th>'."\n";
				print HTML '</tr>'."\n";
				print HTML '</thead>'."\n";
				print HTML '<tbody>'."\n";
				for (0..$#{$parameter->genomes}){
					print HTML '<tr>'."\n";					
					print HTML '<td>'.${$parameter->abbreviations}[$_].'</td>'."\t";
					print HTML '<td><a href="../annotations/'.${$parameter->abbreviations}[$_].'.final.orig.gff">final</a> / <a href="../annotations/'.${$parameter->abbreviations}[$_].'.orig.gff">unreliable</a></td>';
					print HTML '<td><a href="../annotations/'.${$parameter->abbreviations}[$_].'.final.orig.fa">final</a> / <a href="../annotations/'.${$parameter->abbreviations}[$_].'.orig.fa">unreliable</a></td>'."\n";
					print HTML '</tr>'."\n";	
				}
				for (0..$#{$parameter->outgroups}){
					print HTML '<tr>'."\n";					
					print HTML '<td>'.${$parameter->ogabbreviations}[$_].'</td>'."\t";
					print HTML '<td><a href="../annotations/'.${$parameter->ogabbreviations}[$_].'.final.orig.gff">final</a> / <a href="../annotations/'.${$parameter->ogabbreviations}[$_].'.orig.gff">unreliable</a></td>';
					print HTML '<td><a href="../annotations/'.${$parameter->ogabbreviations}[$_].'.final.orig.fa">final</a> / <a href="../annotations/'.${$parameter->ogabbreviations}[$_].'.orig.fa">unreliable</a></td>'."\n";
					print HTML '</tr>'."\n";	
				}
				print HTML '</tbody>'."\n";				
				print HTML '</table>'."\n";
			} elsif($_=~/ncRNA alignments/){
				print HTML $_;
				print HTML '<table class="tablesorter">'."\n";
				print HTML '<thead>'."\n";
				print HTML '<tr>'."\n";
				print HTML '<th>ncRNA</th>'."\t".'<th>Rfam Accession</th>'."\t".'<th>STK</th>';
				for (@{$parameter->abbreviations},@{$parameter->ogabbreviations}){
					print HTML "\t".'<th>'.$_.'</th>';
				}
				print HTML "\n";				
				print HTML '</tr>'."\n";
				print HTML '</thead>'."\n";
				print HTML '<tbody>'."\n";
				for my $rf_rna (keys %{$rf_rnas}){
					my $anno;
					my @counts;
					for (@{$parameter->abbreviations},@{$parameter->ogabbreviations}){
						push @counts , ($#{$gffdb->get_features($rf_rna,[$_],'!')}+1);
						$anno = 1 if $counts[-1] > 0;
					}
					next unless $anno;
					my ($rf , @rna) = split /_/ , $rf_rna;
					print HTML '<tr>'."\n";
					print HTML '<td>'.join('_',@rna).'</td>'."\t";
					print HTML '<td>'.$rf.'</td>'."\t";
					my $path = $stkdb->idToPath->{$rf_rna};
					print HTML $path ? '<td><a href="../alignments/'.$path.'">STK</a></td>' : '<td>NA</tp>';
					for (@counts){												
						print HTML "\t".'<td>'.$_.'</td>';
					}
					print HTML "\n";
					print HTML '</tr>'."\n";
				}				
				print HTML '</tbody>'."\n";				
				print HTML '</table>'."\n";
			} else {
				print HTML $_;
			}
		}
	}
	my $outdir=catdir($parameter->output,'phylogeny-'.$parameter->label);
	if ($parameter->has_outgroups && (-e catfile($outdir,'SSU.mafft.eps') || 
	-e catfile($outdir,'RNome.stk.eps') ||
	-e catfile($outdir,'coreRNome.stk.eps') ||
	-e catfile($outdir,'coreRNome.mafft.eps') )){		
		print HTML '<a name="phylo"></a>'."\n";
		if (-e catfile($outdir,'SSU.mafft.eps')){
			print HTML '<div class="staticbox" id="phylo">'."\n";
			print HTML '<h3>Phylogeny SSU</h3>'."\n";
			print HTML 'Tree based on all SSU rRNA predictions, built from Mafft alignment.<br><br>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/SSU.mafft.eps">EPS</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/SSU.mafft.pdf">PDF</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/SSU.mafft.tree">Newick</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/SSU.mafft">Alignment</a><br><br>'."\n";			
			print HTML '<embed src="../phylogeny-'.$parameter->label.'/SSU.mafft.pdf#view=FitH" width="400" height="400" type="application/pdf">'."\n";
			print HTML '</div>'."\n";
		}
		if (-e catfile($outdir,'RNome.stk.eps')){
			print HTML '<div class="staticbox" id="phylo">'."\n";
			print HTML '<h3>Phylogeny RNome</h3>'."\n";
			print HTML 'Tree based on all ncRNA predictions except rRNAs and tRNAs, built from Stockholm super-alignment.<br><br>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/RNome.stk.eps">EPS</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/RNome.stk.pdf">PDF</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/RNome.stk.tree">Newick</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/RNome.stkfa">Alignment</a><br><br>'."\n";
			print HTML '<embed src="../phylogeny-'.$parameter->label.'/RNome.stk.pdf#view=FitH" width="400" height="400" type="application/pdf">'."\n";
			print HTML '</div>'."\n";
		}
		if (-e catfile($outdir,'core50RNome.stk.eps')){
			print HTML '<div class="staticbox" id="phylo">'."\n";
			print HTML '<h3>Phylogeny core50RNome</h3>'."\n";
			print HTML 'Tree based on ncRNA predictions, present in 50% of given species and except rRNAs and tRNAs, built from Stockholm super-alignment.<br><br>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/core50RNome.stk.eps">EPS</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/core50RNome.stk.pdf">PDF</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/core50RNome.stk.tree">Newick</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/core50RNome.stkfa">Alignment</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/INFO_core50RNome.txt">Included ncRNAs</a><br><br>'."\n";
			print HTML '<embed src="../phylogeny-'.$parameter->label.'/core50RNome.stk.pdf#view=FitH" width="400" height="400" type="application/pdf">'."\n";
			print HTML '</div>'."\n";
		}
		if (-e catfile($outdir,'coreRNome.stk.eps')){
			print HTML '<div class="staticbox" id="phylo">'."\n";
			print HTML '<h3>Phylogeny coreRNome (Stockholm)</h3>'."\n";
			print HTML 'Tree based on ncRNAs present in all species, built from concatenated Stockholm alignments.<br><br>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/coreRNome.stk.eps">EPS</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/coreRNome.stk.pdf">PDF</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/coreRNome.stk.tree">Newick</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/coreRNome.stkfa">Alignment</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/INFO_coreRNome.txt">Included ncRNAs</a><br><br>'."\n";
			print HTML '<embed src="../phylogeny-'.$parameter->label.'/coreRNome.stk.pdf#view=FitH" width="400" height="400" type="application/pdf">'."\n";
			print HTML '</div>'."\n";
		}
		if (-e catfile($outdir,'coreRNome.mafft.eps')){
			print HTML '<div class="staticbox" id="phylo">'."\n";
			print HTML '<h3>Phylogeny coreRNome (Mafft)</h3>'."\n";
			print HTML 'Tree based on ncRNAs present in all species, built by Mafft from concatenated FASTA sequences.<br><br>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/coreRNome.mafft.eps">EPS</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/coreRNome.mafft.pdf">PDF</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/coreRNome.mafft.tree">Newick</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/coreRNome.mafft">Alignment</a>'."\n";
			print HTML '<a href="../phylogeny-'.$parameter->label.'/INFO_coreRNome.txt">Included ncRNAs</a><br><br>'."\n";
			print HTML '<embed src="../phylogeny-'.$parameter->label.'/coreRNome.mafft.pdf#view=FitH" width="400" height="400" type="application/pdf">'."\n";
			print HTML '</div>'."\n";
		}
	}
	print HTML "\n";
	print HTML '</div>'."\n";
	print HTML '</body>'."\n";
	print HTML '</html>'."\n";
	close HTML;
	close INDEX;
	seek DATA, $data_pos, 0;
}

1;

__DATA__

<!DOCTYPE html>
<html>

<head>
	<title>GORAP</title>

	<style> 
		button {
			width: 80px
		}
		#header {
			background-color:black;
			color:white;
			text-align:center;
			padding:1px;
		}
		#nav { 
			float: left;
			line-height:30px;
			background-color:black;
			color:white;
			width:100px;
			padding:5px;	      
		}
		html {
			display: table;
			margin: auto;
			font-family:arial;
			font-size: 8pt;
			text-align: center;
			overflow: auto;
		}
	</style>

	<script src="http://www.rna.uni-jena.de/supplements/src/jquery/jquery-latest.js"></script> 

	<script>
		function getDocHeight() {
			return Math.max(
				Math.max(document.body.scrollHeight, document.documentElement.scrollHeight),
				Math.max(document.body.offsetHeight, document.documentElement.offsetHeight),
				Math.max(document.body.clientHeight, document.documentElement.clientHeight)
			);
		}
		function getDocWidth() {
			return Math.max(
				Math.max(document.body.scrollWidth, document.documentElement.scrollWidth),
				Math.max(document.body.offsetWidth, document.documentElement.offsetWidth),
				Math.max(document.body.clientWidth, document.documentElement.clientWidth)
			);
		}
		function onResize() {
			location.reload();
			location.href = location.href;
		}
		function select() {	
			var element = document.getElementById("select");
			var value = element.options[element.selectedIndex].value;
			if (value != "select") {
				document.getElementById("frame").src = value;
			} else {
				document.getElementById("frame").src = "";
			}
			document.getElementById("frame").onload = function() {
				var x = document.getElementById("frame");
				var y = (x.contentDocument || x.contentWindow);					
				if (y.getElementById("phylo")) {
					document.getElementById("bphylo").style.display="inline";
				} else {
					document.getElementById("bphylo").style.display="none";
				}
				if (y.getElementById("param")) {
					document.getElementById("bparam").style.display="inline";
					document.getElementById("bdb").style.display="inline";
					document.getElementById("banno").style.display="inline";
					document.getElementById("baln").style.display="inline";
				} else {
					document.getElementById("bparam").style.display="none";
					document.getElementById("bdb").style.display="none";
					document.getElementById("banno").style.display="none";
					document.getElementById("baln").style.display="none";
				}
			}
		}
		function resize(){
			document.getElementById("nav").setAttribute("style","height:" + (getDocHeight()-130) + "px");
			document.getElementById("frame").height=getDocHeight()-130;
			document.getElementById("frame").width=getDocWidth()-130;
		}
		function jump(target){
			var element = document.getElementById("select");
			var value = element.options[element.selectedIndex].value;
			if (value != "select") {
				document.getElementById("frame").src = value+target;
				document.getElementById("bparam").style.display="inline";
				document.getElementById("bdb").style.display="inline";
				document.getElementById("banno").style.display="inline";
				document.getElementById("baln").style.display="inline";
				document.getElementById("bphylo").style.display="inline";
			} else {
				document.getElementById("bparam").style.display="none";
				document.getElementById("bdb").style.display="none";
				document.getElementById("banno").style.display="none";
				document.getElementById("baln").style.display="none";
				document.getElementById("bphylo").style.display="none";
			}
		}
	</script>
</head>

<body onResize="onResize()">
<div id="header">
	<!-- reloading (F5) and/or resizing the window may solve displaying issues -->
	<h2>GORAP results</h2>
	<select id="select" onChange="select()">
		<option value="select">Select</option>
	</select>
	<br>
	<br>
	<a href="https://github.com/rna-hta-jena/gorap/blob/master/README" style="text-decoration: none; font-size: 12pt; color: rgb(0,255,0)"><b>>>Manual<<</b></a>
</div>

<div id="nav">
	<button type="button" id="bparam" onclick="jump('#param')" style="display:none;">Parameter</button>
	<br>
	<button type="button" id="bdb" onclick="jump('#db')" style="display:none;">Data</button>
	<br>
	<button type="button" id="banno" onclick="jump('#anno')" style="display:none;">Annotation</button>
	<br>
	<button type="button" id="baln" onclick="jump('#aln')" style="display:none;">Alignments</button>
	<br>
	<button type="button" id="bphylo" onclick="jump('#phylo')" style="display:none;">Phylogeny</button>
</div>

<iframe id="frame" frameborder="0" onLoad="resize()"></iframe>

</body>

</html>

<!DOCTYPE html>
<html>

<head>
	<style> 
		.staticbox {
			float: left;
			margin: 10px;
			padding: 10px;
			width: 400px;
			overflow: auto;
		}
		.box {
			float: left;
			margin: 10px;
			padding: 10px;
			overflow: auto;
		}
		html {
			display: table;
			margin: auto;
			font-family:arial;
			font-size: 8pt;
			overflow: auto;
		}

		table.tablesorter {
			font-family:arial;
			background-color: #CDCDCD;
			margin:10px 0pt 10px;
			font-size: 8pt;
			width: 100%;
			text-align: left;
		}
		table.tablesorter thead tr th, table.tablesorter tfoot tr th {
			background-color: #e6EEEE;
			border: 1px solid #FFF;
			font-size: 8pt;
			padding: 4px;
		}
		table.tablesorter thead tr .header {
			background-image: url(http://tablesorter.com/themes/blue/bg.gif);
			background-repeat: no-repeat;
			background-position: center right;
			cursor: pointer;
			padding-right: 20px;
		}
		table.tablesorter tbody td {
			color: #3D3D3D;
			padding: 4px;
			background-color: #FFF;
			vertical-align: top;
		}
		table.tablesorter thead tr .headerSortUp {
			background-image: url(http://tablesorter.com/themes/blue/asc.gif);
		}
		table.tablesorter thead tr .headerSortDown {
			background-image: url(http://tablesorter.com/themes/blue/desc.gif);
		}
		table.tablesorter thead tr .headerSortDown, table.tablesorter thead tr .headerSortUp {
			background-color: #8dbdd8;
		}
	</style>

	<script src="http://www.rna.uni-jena.de/supplements/src/jquery/jquery-latest.js"></script> 
	<script src="http://www.rna.uni-jena.de/supplements/src/jquery/jquery.tablesorter.js"></script>
    
    <script>
		$(document).ready(function() {
			$("table").tablesorter({
				sortList: [[1,0]] 
			});
		});
	</script>
</head>

<body>
<div>

<div class="staticbox" id="param">
	<a name="param"></a>
	<h3>Parameter</h3>
</div>

<div class="staticbox">
	<a name="db"></a>
	<h3>Used data</h3>
</div>

<div class="staticbox">
	<a name="anno"></a>
	<h3>ncRNA annotation</h3>
</div>

<div class="box">
	<a name="aln"></a>
	<h3>ncRNA alignments</h3>
</div>
