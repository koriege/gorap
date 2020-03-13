# GORAP -  Genomewide ncRNA Annotation Pipeline

**Gorap** screens genomic sequences for all non-coding RNAs present in the **Rfam database** using either a generalized strategy which includes multiple filters or utilizing specialized software. Furthermore, Gorap provides ncRNA based reconstruction of phylogenetic trees and is able to perform *de novo* predictions including TPM calculations from RNA-Seq experiments. RNA family specific screening options like command lines, threshold and constrains can be easily setupped.

**Gorap** is shipped with a setup routine which installs the Perl library ``Bio::Gorap``, all necessary third-party software and databases without any requirements. ``Bio::Gorap`` is a distribution of Perl modules to provide software wrappers and functions for efficient Fasta, GFF, Stockholm alignment file and taxonomic tree manipulation.


## Outline

1. [Used software and databases](#used-software-and-databases)
2. [Installation](#installation)
3. [Update Gorap](#update-gorap)
4. [Run Gorap](#run-gorap)
5. [Results](#results)
6. [Configure Gorap](#configure-gorap)
7. [Contact](#contact)


## Used software and databases

| Tool        | Source                                              |
| ----------- | --------------------------------------------------- |
| Infernal    | <http://eddylab.org/infernal>                       |
| Blast       | <https://blast.ncbi.nlm.nih.gov>                    |
| tRNAscan-SE | <http://lowelab.ucsc.edu/tRNAscan-SE>               |
| Barrnap     | <https://github.com/tseemann/barrnap>               |
| CRT         | <http://room220.com/crt>                            |
| RAxML       | <https://cme.h-its.org/exelixis/web/software/raxml> |
| htslib      | <https://htslib.org>                                |
| BioPerl     | <https://bioperl.org>                               |

| Database      | Source                              |
| ------------- | ----------------------------------- |
| Rfam DB       | <http://rfam.xfam.org>              |
| NCBI Taxonomy | <https://ncbi.nlm.nih.gov/taxonomy> |
| Silva DB      | <https://arb-silva.de>              |


## Installation

Set the ``GORAP`` environment variable, assigned to its installation directory
```
export GORAP=/path/to/install/dir
```
Optionally, Store the ``GORAP`` environment variable permanently to your bashrc
```
echo "export GORAP=$GORAP" >> ~/.bashrc
```
Download Gorap
```
git clone --recursive https://github.com/koriege/gorap.git
git checkout $(git describe --tags)
```
Enter the source directory
```
cd gorap
```
Install ``Bio::Gorap``, Rfam and Silva DB plus third-party tools
```
./setup.sh -d $GORAP
```
Download latest NCBI Taxonomy
```
$GORAP/bin/gorap -update ncbi
```


## Update Gorap

Download Gorap
```
git clone --recursive https://github.com/koriege/gorap.git
git checkout $(git describe --tags)
```
Enter the source directory
```
cd gorap
```
Update Gorap and replace databases plus configurations from enclosed data archive
```
./setup.sh -u -d $GORAP
```
Optionally, download latest (bleeding edge) databases. **This requires manual config file adaptions!**
```
$GORAP/bin/gorap -update all
```


## Run Gorap

### Run Gorap using command line parameters
Screen *E. coli* genome for 6S (Rfam entry 13), and RNaseP (Rfam entries 9-11) on 8 threads given taxonomic information

- **Hint 1**: Multiple runs can be saved into the same output directory
- **Hint 2**: Restart a run with same parameters plus ``-skip`` option fastforwards to downstream analysis
- **Hint 3**: Restart a run with same output directory and label plus ``-refresh`` option will update annotation files (GFF,FASTA) according to present alignments and sequences

```
$GORAP/bin/gorap -q 13 -c 8 -k bac -s 562 -r 543 -i $GORAP/db/example/ecoli.fa -a ecoli -l ecoli_ncrna -o results
$GORAP/bin/gorap -q 9:11 -c 8 -k bac -s 562 -r 543 -i $GORAP/db/example/ecoli.fa -a ecoli -l ecoli_ncrna -o results
$GORAP/bin/gorap -i $GORAP/db/example/ecoli.fa -a ecoli -l ecoli_ncrna -o results -refresh
```

### Run Gorap using a parameter file

- Hint 1: See template parameter.txt file
- Hint 2: Command line parameters priorize parameter file settings

```
$GORAP/bin/gorap -file parameter.txt
$GORAP/bin/gorap -file parameter.txt -skip -sort -l ecoli_sorted
```

### Run Gorap and perform phylogenetic tree reconstruction

**Example for annotation and SSU/RNome based phylogeny reconstruction in one step**

- **Hint 1**: A given outgroup genome will be screened only for ncRNAs predicted in genome files given by -i
- **Hint 2**: For phylogeny reconstruction on curated alignments see below

```
$GORAP/bin/gorap -i <FASTA>,<FASTA>,<FASTA>,<FASTA> -k bac -og <FASTA>
```

**Example for annotation and SSU/RNome based phylogeny reconstruction on curated alignments**

- **Hint 1**: Using labels is highly recommended to update initially created HTML output
- **Hint 2**: Manual curation of alignments can be done after step 1

```
$GORAP/bin/gorap -i <FASTA>,<FASTA>,<FASTA>,<FASTA> -k bac -l mylabel
$GORAP/bin/gorap -i <FASTA>,<FASTA>,<FASTA>,<FASTA> -k bac -l mylabel -refresh
$GORAP/bin/gorap -i <FASTA>,<FASTA>,<FASTA> -og <FASTA> -k bac -l mylabel -skip
```


## Results


- **index.html** file  
	inspect results using a webbrowser

- **annotations** directory  
	**\*.orig** - GFF and FASTA files with original sequence ids  
	**\*.passed** - GFF and FASTA  files with predictions that passed all filters  
	**\*** - GFF and FASTA files with predictions which failed at least one filter  

	GFF filter tags - sorted in order of application

	- **L** - Length filter  
		Applied if a predicted gene is shorter than 40% of the seed consensus sequence.  
	- **B** - Bitscore filter  
		Applied if the score of a prediction gene is below the minimum of scores of taxonomically related species. If no information is passed to Gorap or genes of related species can be found in the seed alignemnt, Gorap uses thee Rfam suggested noise cutoff.  
	- **S** - Structure filter  
		Applied if a predicted gene supports less than 50% conserved seed hairpins.  
	- **P** - Primary sequence filter
		Applied if a predicted gene either supports less than 70% of nucleotides which are strongly conserved (>=90%) within the seed alignment. User defined constrains replace the default conservation based filter strategy.  
	   	**Exception**: C/D-Box snoRNA familiy constrains must cover C and D box motifs which are subsequently screened for a proper kink-turn motif.  
	- **C** - Copy number filter  
		Applied according to the number of allowed pseudogenes as setupped in the RNA family specific configuration files. The primary genes are seleted by their Bitscores.  
	- **O** - Overlap filter  
		Applied if multiple ncRNA predictions intersect each other with more than 50%. All but the highest scored ncRNAs will be filtered.  
	- **!** - All filters passed  
	- **X** - User filter  
		Introduced upon executing Gorap with ``-refresh`` paramter to flag manually deleted sequences.

- **alignments** directory  
	Alignments in Stockholm format which hold all predicted genes (**!** - tagged) plus Rfam seed sequences

- **phylogeny** directory  
	phylogenetic trees in newick and PDF format  
	- **RNome** - all ncRNA predictions (except rRNAs and tRNAs), using a super-matrix approach
	- **core50RNome** - ncRNA predictions (except rRNAs and tRNAs), present in >=50% of given input FASTA files, using a super-matrix approach  
	- **coreRNomeSTK** - ncRNA predictions (except rRNAs and tRNAs) present in all input FASTA files, using a super-matrix approach  
	- **\*.stk** - super-matrix approach contrained by secondary structure based on horizontally concatenated Stockhol alignments
	- **\*.mafft** - super-matrix approach from aligned and horizontally concatenated ncRNA seqeunces

- **meta** directory  
	Stockholm alignments for each Gorap filter applied


## Configure Gorap


Configuration files are located in the ``$GORAP/db/config`` directory and enable RNA family specific definitions of command line(s), thresholds and constrains.

### Custom constrains

Define consensus sequence dependend mismatch constrains as shown below.

**Excpetion**: C/D-Box snoRNA familiy constrains must cover C and D box motifs by **four** region constrains.
Gorap comes with pre-defined constrains for many C/D-Box snoRNAs.

Example:
```
[constrains]
conserved=uuGCAAUGAUGuUAagAAUUUCUUcacCUGAAuuaaaCcuUGAaGuucAAAaauCGAGCUUUUUAACaCUGAGCaaa
constrain=.....|..1...|......................................................|.1..|....
constrain=......|0.|...........................................................|0.|....
```

Use constrains to analogously set the number of allowed mismatches (``..|.<number>.|..``) in a specific region of **non**-snoRNAs.
The given STK consensus can be modified according to IUPAC definition, but must not change in length.

### Custom tools

**Requirement**: The output format of yet unsupported tools must be GFF3

The command line can be defined using this list of available placeholders

- $genome - internally substituted by a fasta file path
- $kingdom - internally substituted by bac, arc, euk, fungi or virus
- $cpus - internally substituted by number of defined worker threads
- $output - internally substituted by a temporary file path  
	without $output , STDOUT will be parsed

Example:
```
[cmd]
tool=barrnap
parameter=--threads $cpus --kingdom $kingdom $genome
```

### Custom RNA families

Create a new configuration file in the ``$GORAP/db/config`` directory with a file name that follows this naming scheme.

``<RFxxxxx_description.cfg>``

To use the default screening methods Infernal and/or Blast, setup the following tool definitions in your config file.

```
[cmd]
tool=Infernal
tool=Blast
```

According to the chosen naming scheme, place the folling files into a new directory under ``$GORAP/db/data/rfam/<RFxxxxx_description>``

- A Stockholm alignment file ``<RFxxxxx_description.stk>``
- A covariance model ``<RFxxxxx_description.cm>``
- A FASTA file ``<RFxxxxx_description.fa>``


## Contact


Suggestions? Bugs? Troubles? Please let me know!

<konstantin.riege@leibniz-fli.de>

See also <https://www.rna.uni-jena.de>