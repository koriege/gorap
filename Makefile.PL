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
		'List::MoreUtils' => '0',
		'Try::Tiny'	=>	'0',
		'Switch'	=>	'0',
		'Math::Round'	=>	'0',
		'Sys::MemInfo' => '0',
		'PostScript::TextBlock' => '0',
		'Config::IniFiles' => '0',
		'Hash::Merge' => '0',
		'Moose'	=> '0',
		'Bio::Root::Version'	=>	'0',
		'Bio::DB::Sam'	=>	'0',
		'Bio::DB::EUtilities'	=>	'0',
		'Tree::Simple'	=> 	'0',
		'Archive::Extract' => '0'
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
