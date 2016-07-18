use ExtUtils::MakeMaker;
use strict;
use warnings;

my @scripts = qw (Gorap.pl Gosno.pl);

WriteMakefile(
	'NAME'	=>	'Bio::Gorap',
	'AUTHOR'	=>	q{Konstantin Riege},
	'VERSION_FROM'	=>	'lib/Bio/Gorap/Gorap.pm',
	'ABSTRACT_FROM'	=>	'lib/Bio/Gorap/Gorap.pm',
	'LICENSE'	=>	'perl',
	'PREREQ_PM'	=>	{
		'Getopt::Long'	=>	'2.0',
		'Bio::Root::Version'	=>	'1.00690001',
		'Bio::DB::Sam'	=>	'1.39',
		'Bio::Index::Fasta'	=>	'0',
		'Bio::SimpleAlign'	=>	'0',
		'Bio::DB::SeqFeature::Store'	=>	'0',
		'Bio::AlignIO'	=>	'0',
		'Bio::DB::Taxonomy'	=>	'0',
		'Bio::Tree::Tree'	=>	'0',
		'Bio::TreeIO'	=>	'0',
		'Bio::DB::EUtilities'	=>	'0',
		'Cwd'	=>	'0',
		'File::Spec::Functions'	=>	'0',
		'File::Basename'	=>	'0',
		'File::Path'	=> '0',
		'Bio::SeqIO'	=> '0',
		'List::Util'	=>	'0',
		'List::MoreUtils' => '0',
		'Try::Tiny'	=>	'0',
		'IPC::Cmd'	=>	'0',
		'IPC::Run'	=>	'0',
		'IPC::Open3'	=>	'0',
		'Symbol'	=>	'0',
		'Switch'	=>	'0',
		'Pod::Usage'	=>	'0',
		'IO::Select'	=>	'0',
		'IO::Pipe'	=>	'0',
		'Moose::Role'	=>	'0',		
		'POSIX'	=>	'0',
		'Moose'	=> '0',
		'Encode'	=> '0',
		'Tree::Simple'	=> 	'0',
		'List::MoreUtils'	=>	'0',
		'Bio::Tree::Draw::Cladogram' => '0',
		'Math::Round'	=>	'0',
		'PostScript::TextBlock' => '0',
		'Sys::MemInfo' => '0',
		'File::Path' => '0',
		'Config::IniFiles' => '0',
		'Hash::Merge' => 0
	},
	'MIN_PERL_VERSION'	=>	'5.010',
	'EXE_FILES'	=>	[ map { "$_" } @scripts ],
	'CONFIGURE_REQUIRES'	=>	{
		'ExtUtils::MakeMaker'	=>	0
	},
	'BUILD_REQUIRES'	=>	{
		'Test::More' => 0
	},
	dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
	clean => { FILES => 'Bio-Gorap-*' },
);
