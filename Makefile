# This Makefile is for the Bio::Gorap extension to perl.
#
# It was generated automatically by MakeMaker version
# 6.66 (Revision: 66600) from the contents of
# Makefile.PL. Don't edit this file, edit Makefile.PL instead.
#
#       ANY CHANGES MADE HERE WILL BE LOST!
#
#   MakeMaker ARGV: ()
#

#   MakeMaker Parameters:

#     ABSTRACT_FROM => q[lib/Bio/Gorap/Gorap.pm]
#     AUTHOR => [q[Konstantin Riege]]
#     BUILD_REQUIRES => { Test::More=>q[0] }
#     CONFIGURE_REQUIRES => { ExtUtils::MakeMaker=>q[0] }
#     EXE_FILES => [q[Gorap.pl], q[Gorap_noTaxonomy.pl], q[Evaluation.pl]]
#     LICENSE => q[perl]
#     MIN_PERL_VERSION => q[5.010]
#     NAME => q[Bio::Gorap]
#     PREREQ_PM => { Bio::DB::EUtilities=>q[0], Encode=>q[0], IO::Select=>q[0], Bio::Root::Version=>q[1.00690001], Bio::SimpleAlign=>q[0], IPC::Cmd=>q[0], List::Util=>q[0], Cwd=>q[0], Bio::Tree::Draw::Cladogram=>q[0], Getopt::Long=>q[2.0], Test::More=>q[0], File::Path=>q[0], Bio::TreeIO=>q[0], Bio::DB::Sam=>q[1.39], Pod::Usage=>q[0], List::MoreUtils=>q[0], Bio::DB::SeqFeature::Store=>q[0], Bio::DB::Taxonomy=>q[0], Try::Tiny=>q[0], IPC::Run=>q[0], IPC::Open3=>q[0], Bio::Tree::Tree=>q[0], Bio::Index::Fasta=>q[0], File::Spec::Functions=>q[0], Moose::Role=>q[0], Switch=>q[0], IO::Pipe=>q[0], Moose=>q[0], Tree::Simple=>q[0], Symbol=>q[0], Bio::SeqIO=>q[0], Bio::AlignIO=>q[0], File::Basename=>q[0], Math::Round=>q[0], POSIX=>q[0] }
#     TEST_REQUIRES => {  }
#     VERSION_FROM => q[lib/Bio/Gorap/Gorap.pm]
#     clean => { FILES=>q[Bio-Gorap-*] }
#     dist => { COMPRESS=>q[gzip -9f], SUFFIX=>q[gz] }

# --- MakeMaker post_initialize section:


# --- MakeMaker const_config section:

# These definitions are from config.sh (via /usr/lib/perl/5.18/Config.pm).
# They may have been overridden via Makefile.PL or on the command line.
AR = ar
CC = cc
CCCDLFLAGS = -fPIC
CCDLFLAGS = -Wl,-E
DLEXT = so
DLSRC = dl_dlopen.xs
EXE_EXT = 
FULL_AR = /usr/bin/ar
LD = cc
LDDLFLAGS = -shared -L/usr/local/lib -fstack-protector
LDFLAGS =  -fstack-protector -L/usr/local/lib
LIBC = 
LIB_EXT = .a
OBJ_EXT = .o
OSNAME = linux
OSVERS = 3.2.0-58-generic
RANLIB = :
SITELIBEXP = /usr/local/share/perl/5.18.2
SITEARCHEXP = /usr/local/lib/perl/5.18.2
SO = so
VENDORARCHEXP = /usr/lib/perl5
VENDORLIBEXP = /usr/share/perl5


# --- MakeMaker constants section:
AR_STATIC_ARGS = cr
DIRFILESEP = /
DFSEP = $(DIRFILESEP)
NAME = Bio::Gorap
NAME_SYM = Bio_Gorap
VERSION = v2.0
VERSION_MACRO = VERSION
VERSION_SYM = v2_0
DEFINE_VERSION = -D$(VERSION_MACRO)=\"$(VERSION)\"
XS_VERSION = v2.0
XS_VERSION_MACRO = XS_VERSION
XS_DEFINE_VERSION = -D$(XS_VERSION_MACRO)=\"$(XS_VERSION)\"
INST_ARCHLIB = blib/arch
INST_SCRIPT = blib/script
INST_BIN = blib/bin
INST_LIB = blib/lib
INST_MAN1DIR = blib/man1
INST_MAN3DIR = blib/man3
MAN1EXT = 1p
MAN3EXT = 3pm
INSTALLDIRS = site
INSTALL_BASE = /home/koriege/perl5
DESTDIR = 
PREFIX = $(INSTALL_BASE)
INSTALLPRIVLIB = $(INSTALL_BASE)/lib/perl5
DESTINSTALLPRIVLIB = $(DESTDIR)$(INSTALLPRIVLIB)
INSTALLSITELIB = $(INSTALL_BASE)/lib/perl5
DESTINSTALLSITELIB = $(DESTDIR)$(INSTALLSITELIB)
INSTALLVENDORLIB = $(INSTALL_BASE)/lib/perl5
DESTINSTALLVENDORLIB = $(DESTDIR)$(INSTALLVENDORLIB)
INSTALLARCHLIB = $(INSTALL_BASE)/lib/perl5/x86_64-linux-gnu-thread-multi
DESTINSTALLARCHLIB = $(DESTDIR)$(INSTALLARCHLIB)
INSTALLSITEARCH = $(INSTALL_BASE)/lib/perl5/x86_64-linux-gnu-thread-multi
DESTINSTALLSITEARCH = $(DESTDIR)$(INSTALLSITEARCH)
INSTALLVENDORARCH = $(INSTALL_BASE)/lib/perl5/x86_64-linux-gnu-thread-multi
DESTINSTALLVENDORARCH = $(DESTDIR)$(INSTALLVENDORARCH)
INSTALLBIN = $(INSTALL_BASE)/bin
DESTINSTALLBIN = $(DESTDIR)$(INSTALLBIN)
INSTALLSITEBIN = $(INSTALL_BASE)/bin
DESTINSTALLSITEBIN = $(DESTDIR)$(INSTALLSITEBIN)
INSTALLVENDORBIN = $(INSTALL_BASE)/bin
DESTINSTALLVENDORBIN = $(DESTDIR)$(INSTALLVENDORBIN)
INSTALLSCRIPT = $(INSTALL_BASE)/bin
DESTINSTALLSCRIPT = $(DESTDIR)$(INSTALLSCRIPT)
INSTALLSITESCRIPT = $(INSTALL_BASE)/bin
DESTINSTALLSITESCRIPT = $(DESTDIR)$(INSTALLSITESCRIPT)
INSTALLVENDORSCRIPT = $(INSTALL_BASE)/bin
DESTINSTALLVENDORSCRIPT = $(DESTDIR)$(INSTALLVENDORSCRIPT)
INSTALLMAN1DIR = $(INSTALL_BASE)/man/man1
DESTINSTALLMAN1DIR = $(DESTDIR)$(INSTALLMAN1DIR)
INSTALLSITEMAN1DIR = $(INSTALL_BASE)/man/man1
DESTINSTALLSITEMAN1DIR = $(DESTDIR)$(INSTALLSITEMAN1DIR)
INSTALLVENDORMAN1DIR = $(INSTALL_BASE)/man/man1
DESTINSTALLVENDORMAN1DIR = $(DESTDIR)$(INSTALLVENDORMAN1DIR)
INSTALLMAN3DIR = $(INSTALL_BASE)/man/man3
DESTINSTALLMAN3DIR = $(DESTDIR)$(INSTALLMAN3DIR)
INSTALLSITEMAN3DIR = $(INSTALL_BASE)/man/man3
DESTINSTALLSITEMAN3DIR = $(DESTDIR)$(INSTALLSITEMAN3DIR)
INSTALLVENDORMAN3DIR = $(INSTALL_BASE)/man/man3
DESTINSTALLVENDORMAN3DIR = $(DESTDIR)$(INSTALLVENDORMAN3DIR)
PERL_LIB = /usr/share/perl/5.18
PERL_ARCHLIB = /usr/lib/perl/5.18
LIBPERL_A = libperl.a
FIRST_MAKEFILE = Makefile
MAKEFILE_OLD = Makefile.old
MAKE_APERL_FILE = Makefile.aperl
PERLMAINCC = $(CC)
PERL_INC = /usr/lib/perl/5.18/CORE
PERL = /usr/bin/perl
FULLPERL = /usr/bin/perl
ABSPERL = $(PERL)
PERLRUN = $(PERL)
FULLPERLRUN = $(FULLPERL)
ABSPERLRUN = $(ABSPERL)
PERLRUNINST = $(PERLRUN) "-I$(INST_ARCHLIB)" "-I$(INST_LIB)"
FULLPERLRUNINST = $(FULLPERLRUN) "-I$(INST_ARCHLIB)" "-I$(INST_LIB)"
ABSPERLRUNINST = $(ABSPERLRUN) "-I$(INST_ARCHLIB)" "-I$(INST_LIB)"
PERL_CORE = 0
PERM_DIR = 755
PERM_RW = 644
PERM_RWX = 755

MAKEMAKER   = /usr/share/perl/5.18/ExtUtils/MakeMaker.pm
MM_VERSION  = 6.66
MM_REVISION = 66600

# FULLEXT = Pathname for extension directory (eg Foo/Bar/Oracle).
# BASEEXT = Basename part of FULLEXT. May be just equal FULLEXT. (eg Oracle)
# PARENT_NAME = NAME without BASEEXT and no trailing :: (eg Foo::Bar)
# DLBASE  = Basename part of dynamic library. May be just equal BASEEXT.
MAKE = make
FULLEXT = Bio/Gorap
BASEEXT = Gorap
PARENT_NAME = Bio
DLBASE = $(BASEEXT)
VERSION_FROM = lib/Bio/Gorap/Gorap.pm
OBJECT = 
LDFROM = $(OBJECT)
LINKTYPE = dynamic
BOOTDEP = 

# Handy lists of source code files:
XS_FILES = 
C_FILES  = 
O_FILES  = 
H_FILES  = 
MAN1PODS = Gorap.pl \
	Gorap_noTaxonomy.pl
MAN3PODS = Gorap.pl \
	Gorap_noTaxonomy.pl \
	RefreshOutput.pl \
	lib/Bio/Gorap/Gorap.pm

# Where is the Config information that we are using/depend on
CONFIGDEP = $(PERL_ARCHLIB)$(DFSEP)Config.pm $(PERL_INC)$(DFSEP)config.h

# Where to build things
INST_LIBDIR      = $(INST_LIB)/Bio
INST_ARCHLIBDIR  = $(INST_ARCHLIB)/Bio

INST_AUTODIR     = $(INST_LIB)/auto/$(FULLEXT)
INST_ARCHAUTODIR = $(INST_ARCHLIB)/auto/$(FULLEXT)

INST_STATIC      = 
INST_DYNAMIC     = 
INST_BOOT        = 

# Extra linker info
EXPORT_LIST        = 
PERL_ARCHIVE       = 
PERL_ARCHIVE_AFTER = 


TO_INST_PM = Evaluation.pl \
	Gorap.pl \
	Gorap_noTaxonomy.pl \
	RefreshOutput.pl \
	lib/Bio/Gorap/CFG.pm \
	lib/Bio/Gorap/DB/BAM.pm \
	lib/Bio/Gorap/DB/Fasta.pm \
	lib/Bio/Gorap/DB/GFF.pm \
	lib/Bio/Gorap/DB/STK.pm \
	lib/Bio/Gorap/DB/Taxonomy.pm \
	lib/Bio/Gorap/Evaluation/Clustering.pm \
	lib/Bio/Gorap/Evaluation/HTML.pm \
	lib/Bio/Gorap/Evaluation/Statistics.pm \
	lib/Bio/Gorap/Functions/CM.pm \
	lib/Bio/Gorap/Functions/STK.pm \
	lib/Bio/Gorap/Functions/ToolParser.pm \
	lib/Bio/Gorap/Gorap.pm \
	lib/Bio/Gorap/Parameter.pm \
	lib/Bio/Gorap/ThrListener.pm \
	lib/Bio/Gorap/Tool/Bcheck.pm \
	lib/Bio/Gorap/Tool/Blast.pm \
	lib/Bio/Gorap/Tool/Crt.pm \
	lib/Bio/Gorap/Tool/Default.pm \
	lib/Bio/Gorap/Tool/Infernal.pm \
	lib/Bio/Gorap/Tool/Rnammer.pm \
	lib/Bio/Gorap/Tool/Trnascanse.pm \
	lib/Bio/Gorap/ToolI.pm \
	lib/Bio/Gorap/Update.pm

PM_TO_BLIB = Gorap.pl \
	$(INST_LIB)/Bio/Gorap.pl \
	lib/Bio/Gorap/CFG.pm \
	blib/lib/Bio/Gorap/CFG.pm \
	lib/Bio/Gorap/DB/STK.pm \
	blib/lib/Bio/Gorap/DB/STK.pm \
	lib/Bio/Gorap/DB/GFF.pm \
	blib/lib/Bio/Gorap/DB/GFF.pm \
	lib/Bio/Gorap/DB/Fasta.pm \
	blib/lib/Bio/Gorap/DB/Fasta.pm \
	lib/Bio/Gorap/Update.pm \
	blib/lib/Bio/Gorap/Update.pm \
	lib/Bio/Gorap/Tool/Default.pm \
	blib/lib/Bio/Gorap/Tool/Default.pm \
	lib/Bio/Gorap/Gorap.pm \
	blib/lib/Bio/Gorap/Gorap.pm \
	lib/Bio/Gorap/ThrListener.pm \
	blib/lib/Bio/Gorap/ThrListener.pm \
	lib/Bio/Gorap/DB/Taxonomy.pm \
	blib/lib/Bio/Gorap/DB/Taxonomy.pm \
	lib/Bio/Gorap/Tool/Infernal.pm \
	blib/lib/Bio/Gorap/Tool/Infernal.pm \
	lib/Bio/Gorap/Tool/Bcheck.pm \
	blib/lib/Bio/Gorap/Tool/Bcheck.pm \
	lib/Bio/Gorap/Evaluation/Clustering.pm \
	blib/lib/Bio/Gorap/Evaluation/Clustering.pm \
	lib/Bio/Gorap/ToolI.pm \
	blib/lib/Bio/Gorap/ToolI.pm \
	RefreshOutput.pl \
	$(INST_LIB)/Bio/RefreshOutput.pl \
	lib/Bio/Gorap/Tool/Trnascanse.pm \
	blib/lib/Bio/Gorap/Tool/Trnascanse.pm \
	lib/Bio/Gorap/Functions/ToolParser.pm \
	blib/lib/Bio/Gorap/Functions/ToolParser.pm \
	lib/Bio/Gorap/Functions/STK.pm \
	blib/lib/Bio/Gorap/Functions/STK.pm \
	lib/Bio/Gorap/Evaluation/HTML.pm \
	blib/lib/Bio/Gorap/Evaluation/HTML.pm \
	lib/Bio/Gorap/Tool/Rnammer.pm \
	blib/lib/Bio/Gorap/Tool/Rnammer.pm \
	Gorap_noTaxonomy.pl \
	$(INST_LIB)/Bio/Gorap_noTaxonomy.pl \
	lib/Bio/Gorap/Tool/Blast.pm \
	blib/lib/Bio/Gorap/Tool/Blast.pm \
	lib/Bio/Gorap/DB/BAM.pm \
	blib/lib/Bio/Gorap/DB/BAM.pm \
	lib/Bio/Gorap/Tool/Crt.pm \
	blib/lib/Bio/Gorap/Tool/Crt.pm \
	lib/Bio/Gorap/Parameter.pm \
	blib/lib/Bio/Gorap/Parameter.pm \
	lib/Bio/Gorap/Functions/CM.pm \
	blib/lib/Bio/Gorap/Functions/CM.pm \
	Evaluation.pl \
	$(INST_LIB)/Bio/Evaluation.pl \
	lib/Bio/Gorap/Evaluation/Statistics.pm \
	blib/lib/Bio/Gorap/Evaluation/Statistics.pm


# --- MakeMaker platform_constants section:
MM_Unix_VERSION = 6.66
PERL_MALLOC_DEF = -DPERL_EXTMALLOC_DEF -Dmalloc=Perl_malloc -Dfree=Perl_mfree -Drealloc=Perl_realloc -Dcalloc=Perl_calloc


# --- MakeMaker tool_autosplit section:
# Usage: $(AUTOSPLITFILE) FileToSplit AutoDirToSplitInto
AUTOSPLITFILE = $(ABSPERLRUN)  -e 'use AutoSplit;  autosplit($$$$ARGV[0], $$$$ARGV[1], 0, 1, 1)' --



# --- MakeMaker tool_xsubpp section:


# --- MakeMaker tools_other section:
SHELL = /bin/sh
CHMOD = chmod
CP = cp
MV = mv
NOOP = $(TRUE)
NOECHO = @
RM_F = rm -f
RM_RF = rm -rf
TEST_F = test -f
TOUCH = touch
UMASK_NULL = umask 0
DEV_NULL = > /dev/null 2>&1
MKPATH = $(ABSPERLRUN) -MExtUtils::Command -e 'mkpath' --
EQUALIZE_TIMESTAMP = $(ABSPERLRUN) -MExtUtils::Command -e 'eqtime' --
FALSE = false
TRUE = true
ECHO = echo
ECHO_N = echo -n
UNINST = 0
VERBINST = 0
MOD_INSTALL = $(ABSPERLRUN) -MExtUtils::Install -e 'install([ from_to => {@ARGV}, verbose => '\''$(VERBINST)'\'', uninstall_shadows => '\''$(UNINST)'\'', dir_mode => '\''$(PERM_DIR)'\'' ]);' --
DOC_INSTALL = $(ABSPERLRUN) -MExtUtils::Command::MM -e 'perllocal_install' --
UNINSTALL = $(ABSPERLRUN) -MExtUtils::Command::MM -e 'uninstall' --
WARN_IF_OLD_PACKLIST = $(ABSPERLRUN) -MExtUtils::Command::MM -e 'warn_if_old_packlist' --
MACROSTART = 
MACROEND = 
USEMAKEFILE = -f
FIXIN = $(ABSPERLRUN) -MExtUtils::MY -e 'MY->fixin(shift)' --


# --- MakeMaker makemakerdflt section:
makemakerdflt : all
	$(NOECHO) $(NOOP)


# --- MakeMaker dist section:
TAR = tar
TARFLAGS = cvf
ZIP = zip
ZIPFLAGS = -r
COMPRESS = gzip -9f
SUFFIX = gz
SHAR = shar
PREOP = $(NOECHO) $(NOOP)
POSTOP = $(NOECHO) $(NOOP)
TO_UNIX = $(NOECHO) $(NOOP)
CI = ci -u
RCS_LABEL = rcs -Nv$(VERSION_SYM): -q
DIST_CP = best
DIST_DEFAULT = tardist
DISTNAME = Bio-Gorap
DISTVNAME = Bio-Gorap-v2.0


# --- MakeMaker macro section:


# --- MakeMaker depend section:


# --- MakeMaker cflags section:


# --- MakeMaker const_loadlibs section:


# --- MakeMaker const_cccmd section:


# --- MakeMaker post_constants section:


# --- MakeMaker pasthru section:

PASTHRU = LIBPERL_A="$(LIBPERL_A)"\
	LINKTYPE="$(LINKTYPE)"\
	LD="$(LD)"\
	PREFIX="$(PREFIX)"\
	INSTALL_BASE="$(INSTALL_BASE)"


# --- MakeMaker special_targets section:
.SUFFIXES : .xs .c .C .cpp .i .s .cxx .cc $(OBJ_EXT)

.PHONY: all config static dynamic test linkext manifest blibdirs clean realclean disttest distdir



# --- MakeMaker c_o section:


# --- MakeMaker xs_c section:


# --- MakeMaker xs_o section:


# --- MakeMaker top_targets section:
all :: pure_all manifypods
	$(NOECHO) $(NOOP)


pure_all :: config pm_to_blib subdirs linkext
	$(NOECHO) $(NOOP)

subdirs :: $(MYEXTLIB)
	$(NOECHO) $(NOOP)

config :: $(FIRST_MAKEFILE) blibdirs
	$(NOECHO) $(NOOP)

help :
	perldoc ExtUtils::MakeMaker


# --- MakeMaker blibdirs section:
blibdirs : $(INST_LIBDIR)$(DFSEP).exists $(INST_ARCHLIB)$(DFSEP).exists $(INST_AUTODIR)$(DFSEP).exists $(INST_ARCHAUTODIR)$(DFSEP).exists $(INST_BIN)$(DFSEP).exists $(INST_SCRIPT)$(DFSEP).exists $(INST_MAN1DIR)$(DFSEP).exists $(INST_MAN3DIR)$(DFSEP).exists
	$(NOECHO) $(NOOP)

# Backwards compat with 6.18 through 6.25
blibdirs.ts : blibdirs
	$(NOECHO) $(NOOP)

$(INST_LIBDIR)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_LIBDIR)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_LIBDIR)
	$(NOECHO) $(TOUCH) $(INST_LIBDIR)$(DFSEP).exists

$(INST_ARCHLIB)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_ARCHLIB)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_ARCHLIB)
	$(NOECHO) $(TOUCH) $(INST_ARCHLIB)$(DFSEP).exists

$(INST_AUTODIR)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_AUTODIR)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_AUTODIR)
	$(NOECHO) $(TOUCH) $(INST_AUTODIR)$(DFSEP).exists

$(INST_ARCHAUTODIR)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_ARCHAUTODIR)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_ARCHAUTODIR)
	$(NOECHO) $(TOUCH) $(INST_ARCHAUTODIR)$(DFSEP).exists

$(INST_BIN)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_BIN)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_BIN)
	$(NOECHO) $(TOUCH) $(INST_BIN)$(DFSEP).exists

$(INST_SCRIPT)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_SCRIPT)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_SCRIPT)
	$(NOECHO) $(TOUCH) $(INST_SCRIPT)$(DFSEP).exists

$(INST_MAN1DIR)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_MAN1DIR)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_MAN1DIR)
	$(NOECHO) $(TOUCH) $(INST_MAN1DIR)$(DFSEP).exists

$(INST_MAN3DIR)$(DFSEP).exists :: Makefile.PL
	$(NOECHO) $(MKPATH) $(INST_MAN3DIR)
	$(NOECHO) $(CHMOD) $(PERM_DIR) $(INST_MAN3DIR)
	$(NOECHO) $(TOUCH) $(INST_MAN3DIR)$(DFSEP).exists



# --- MakeMaker linkext section:

linkext :: $(LINKTYPE)
	$(NOECHO) $(NOOP)


# --- MakeMaker dlsyms section:


# --- MakeMaker dynamic section:

dynamic :: $(FIRST_MAKEFILE) $(INST_DYNAMIC) $(INST_BOOT)
	$(NOECHO) $(NOOP)


# --- MakeMaker dynamic_bs section:

BOOTSTRAP =


# --- MakeMaker dynamic_lib section:


# --- MakeMaker static section:

## $(INST_PM) has been moved to the all: target.
## It remains here for awhile to allow for old usage: "make static"
static :: $(FIRST_MAKEFILE) $(INST_STATIC)
	$(NOECHO) $(NOOP)


# --- MakeMaker static_lib section:


# --- MakeMaker manifypods section:

POD2MAN_EXE = $(PERLRUN) "-MExtUtils::Command::MM" -e pod2man "--"
POD2MAN = $(POD2MAN_EXE)


manifypods : pure_all  \
	Gorap.pl \
	Gorap_noTaxonomy.pl \
	Gorap_noTaxonomy.pl \
	Gorap.pl \
	lib/Bio/Gorap/Gorap.pm \
	RefreshOutput.pl
	$(NOECHO) $(POD2MAN) --section=$(MAN1EXT) --perm_rw=$(PERM_RW) \
	  Gorap.pl $(INST_MAN1DIR)/Gorap.pl.$(MAN1EXT) \
	  Gorap_noTaxonomy.pl $(INST_MAN1DIR)/Gorap_noTaxonomy.pl.$(MAN1EXT) 
	$(NOECHO) $(POD2MAN) --section=$(MAN3EXT) --perm_rw=$(PERM_RW) \
	  Gorap_noTaxonomy.pl $(INST_MAN3DIR)/Bio::Gorap_noTaxonomy.$(MAN3EXT) \
	  Gorap.pl $(INST_MAN3DIR)/Bio::Gorap.$(MAN3EXT) \
	  lib/Bio/Gorap/Gorap.pm $(INST_MAN3DIR)/Bio::Gorap::Gorap.$(MAN3EXT) \
	  RefreshOutput.pl $(INST_MAN3DIR)/Bio::RefreshOutput.$(MAN3EXT) 




# --- MakeMaker processPL section:


# --- MakeMaker installbin section:

EXE_FILES = Gorap.pl Gorap_noTaxonomy.pl Evaluation.pl

pure_all :: $(INST_SCRIPT)/Gorap_noTaxonomy.pl $(INST_SCRIPT)/Gorap.pl $(INST_SCRIPT)/Evaluation.pl
	$(NOECHO) $(NOOP)

realclean ::
	$(RM_F) \
	  $(INST_SCRIPT)/Gorap_noTaxonomy.pl $(INST_SCRIPT)/Gorap.pl \
	  $(INST_SCRIPT)/Evaluation.pl 

$(INST_SCRIPT)/Gorap_noTaxonomy.pl : Gorap_noTaxonomy.pl $(FIRST_MAKEFILE) $(INST_SCRIPT)$(DFSEP).exists $(INST_BIN)$(DFSEP).exists
	$(NOECHO) $(RM_F) $(INST_SCRIPT)/Gorap_noTaxonomy.pl
	$(CP) Gorap_noTaxonomy.pl $(INST_SCRIPT)/Gorap_noTaxonomy.pl
	$(FIXIN) $(INST_SCRIPT)/Gorap_noTaxonomy.pl
	-$(NOECHO) $(CHMOD) $(PERM_RWX) $(INST_SCRIPT)/Gorap_noTaxonomy.pl

$(INST_SCRIPT)/Gorap.pl : Gorap.pl $(FIRST_MAKEFILE) $(INST_SCRIPT)$(DFSEP).exists $(INST_BIN)$(DFSEP).exists
	$(NOECHO) $(RM_F) $(INST_SCRIPT)/Gorap.pl
	$(CP) Gorap.pl $(INST_SCRIPT)/Gorap.pl
	$(FIXIN) $(INST_SCRIPT)/Gorap.pl
	-$(NOECHO) $(CHMOD) $(PERM_RWX) $(INST_SCRIPT)/Gorap.pl

$(INST_SCRIPT)/Evaluation.pl : Evaluation.pl $(FIRST_MAKEFILE) $(INST_SCRIPT)$(DFSEP).exists $(INST_BIN)$(DFSEP).exists
	$(NOECHO) $(RM_F) $(INST_SCRIPT)/Evaluation.pl
	$(CP) Evaluation.pl $(INST_SCRIPT)/Evaluation.pl
	$(FIXIN) $(INST_SCRIPT)/Evaluation.pl
	-$(NOECHO) $(CHMOD) $(PERM_RWX) $(INST_SCRIPT)/Evaluation.pl



# --- MakeMaker subdirs section:

# none

# --- MakeMaker clean_subdirs section:
clean_subdirs :
	$(NOECHO) $(NOOP)


# --- MakeMaker clean section:

# Delete temporary files but do not touch installed files. We don't delete
# the Makefile here so a later make realclean still has a makefile to use.

clean :: clean_subdirs
	- $(RM_F) \
	  lib$(BASEEXT).def perlmain.c \
	  core.[0-9] $(BASEEXT).def \
	  core.[0-9][0-9][0-9] core.[0-9][0-9][0-9][0-9] \
	  perl.exe *perl.core \
	  core pm_to_blib \
	  $(MAKE_APERL_FILE) core.[0-9][0-9] \
	  tmon.out $(BASEEXT).exp \
	  $(BASEEXT).x $(BASEEXT).bso \
	  so_locations core.*perl.*.? \
	  perl$(EXE_EXT) pm_to_blib.ts \
	  core.[0-9][0-9][0-9][0-9][0-9] *$(OBJ_EXT) \
	  $(INST_ARCHAUTODIR)/extralibs.ld $(BOOTSTRAP) \
	  perl mon.out \
	  $(INST_ARCHAUTODIR)/extralibs.all MYMETA.json \
	  MYMETA.yml blibdirs.ts \
	  *$(LIB_EXT) 
	- $(RM_RF) \
	  blib Bio-Gorap-* 
	- $(MV) $(FIRST_MAKEFILE) $(MAKEFILE_OLD) $(DEV_NULL)


# --- MakeMaker realclean_subdirs section:
realclean_subdirs :
	$(NOECHO) $(NOOP)


# --- MakeMaker realclean section:
# Delete temporary files (via clean) and also delete dist files
realclean purge ::  clean realclean_subdirs
	- $(RM_F) \
	  $(MAKEFILE_OLD) $(FIRST_MAKEFILE) 
	- $(RM_RF) \
	  $(DISTVNAME) 


# --- MakeMaker metafile section:
metafile : create_distdir
	$(NOECHO) $(ECHO) Generating META.yml
	$(NOECHO) $(ECHO) '---' > META_new.yml
	$(NOECHO) $(ECHO) 'abstract: '\''A Perl distribution for genome wide ncRNA annotation based on various software '\''' >> META_new.yml
	$(NOECHO) $(ECHO) 'author:' >> META_new.yml
	$(NOECHO) $(ECHO) '  - '\''Konstantin Riege'\''' >> META_new.yml
	$(NOECHO) $(ECHO) 'build_requires:' >> META_new.yml
	$(NOECHO) $(ECHO) '  Test::More: 0' >> META_new.yml
	$(NOECHO) $(ECHO) 'configure_requires:' >> META_new.yml
	$(NOECHO) $(ECHO) '  ExtUtils::MakeMaker: 0' >> META_new.yml
	$(NOECHO) $(ECHO) 'dynamic_config: 1' >> META_new.yml
	$(NOECHO) $(ECHO) 'generated_by: '\''ExtUtils::MakeMaker version 6.66, CPAN::Meta::Converter version 2.120921'\''' >> META_new.yml
	$(NOECHO) $(ECHO) 'license: perl' >> META_new.yml
	$(NOECHO) $(ECHO) 'meta-spec:' >> META_new.yml
	$(NOECHO) $(ECHO) '  url: http://module-build.sourceforge.net/META-spec-v1.4.html' >> META_new.yml
	$(NOECHO) $(ECHO) '  version: 1.4' >> META_new.yml
	$(NOECHO) $(ECHO) 'name: Bio-Gorap' >> META_new.yml
	$(NOECHO) $(ECHO) 'no_index:' >> META_new.yml
	$(NOECHO) $(ECHO) '  directory:' >> META_new.yml
	$(NOECHO) $(ECHO) '    - t' >> META_new.yml
	$(NOECHO) $(ECHO) '    - inc' >> META_new.yml
	$(NOECHO) $(ECHO) 'requires:' >> META_new.yml
	$(NOECHO) $(ECHO) '  Bio::AlignIO: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Bio::DB::EUtilities: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Bio::DB::Sam: 1.39' >> META_new.yml
	$(NOECHO) $(ECHO) '  Bio::DB::SeqFeature::Store: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Bio::DB::Taxonomy: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Bio::Index::Fasta: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Bio::Root::Version: 1.00690001' >> META_new.yml
	$(NOECHO) $(ECHO) '  Bio::SeqIO: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Bio::SimpleAlign: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Bio::Tree::Draw::Cladogram: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Bio::Tree::Tree: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Bio::TreeIO: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Cwd: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Encode: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  File::Basename: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  File::Path: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  File::Spec::Functions: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Getopt::Long: 2.0' >> META_new.yml
	$(NOECHO) $(ECHO) '  IO::Pipe: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  IO::Select: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  IPC::Cmd: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  IPC::Open3: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  IPC::Run: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  List::MoreUtils: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  List::Util: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Math::Round: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Moose: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Moose::Role: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  POSIX: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Pod::Usage: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Switch: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Symbol: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Tree::Simple: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  Try::Tiny: 0' >> META_new.yml
	$(NOECHO) $(ECHO) '  perl: 5.010' >> META_new.yml
	$(NOECHO) $(ECHO) 'version: v2.0.0' >> META_new.yml
	-$(NOECHO) $(MV) META_new.yml $(DISTVNAME)/META.yml
	$(NOECHO) $(ECHO) Generating META.json
	$(NOECHO) $(ECHO) '{' > META_new.json
	$(NOECHO) $(ECHO) '   "abstract" : "A Perl distribution for genome wide ncRNA annotation based on various software ",' >> META_new.json
	$(NOECHO) $(ECHO) '   "author" : [' >> META_new.json
	$(NOECHO) $(ECHO) '      "Konstantin Riege"' >> META_new.json
	$(NOECHO) $(ECHO) '   ],' >> META_new.json
	$(NOECHO) $(ECHO) '   "dynamic_config" : 1,' >> META_new.json
	$(NOECHO) $(ECHO) '   "generated_by" : "ExtUtils::MakeMaker version 6.66, CPAN::Meta::Converter version 2.120921",' >> META_new.json
	$(NOECHO) $(ECHO) '   "license" : [' >> META_new.json
	$(NOECHO) $(ECHO) '      "perl_5"' >> META_new.json
	$(NOECHO) $(ECHO) '   ],' >> META_new.json
	$(NOECHO) $(ECHO) '   "meta-spec" : {' >> META_new.json
	$(NOECHO) $(ECHO) '      "url" : "http://search.cpan.org/perldoc?CPAN::Meta::Spec",' >> META_new.json
	$(NOECHO) $(ECHO) '      "version" : "2"' >> META_new.json
	$(NOECHO) $(ECHO) '   },' >> META_new.json
	$(NOECHO) $(ECHO) '   "name" : "Bio-Gorap",' >> META_new.json
	$(NOECHO) $(ECHO) '   "no_index" : {' >> META_new.json
	$(NOECHO) $(ECHO) '      "directory" : [' >> META_new.json
	$(NOECHO) $(ECHO) '         "t",' >> META_new.json
	$(NOECHO) $(ECHO) '         "inc"' >> META_new.json
	$(NOECHO) $(ECHO) '      ]' >> META_new.json
	$(NOECHO) $(ECHO) '   },' >> META_new.json
	$(NOECHO) $(ECHO) '   "prereqs" : {' >> META_new.json
	$(NOECHO) $(ECHO) '      "build" : {' >> META_new.json
	$(NOECHO) $(ECHO) '         "requires" : {' >> META_new.json
	$(NOECHO) $(ECHO) '            "Test::More" : "0"' >> META_new.json
	$(NOECHO) $(ECHO) '         }' >> META_new.json
	$(NOECHO) $(ECHO) '      },' >> META_new.json
	$(NOECHO) $(ECHO) '      "configure" : {' >> META_new.json
	$(NOECHO) $(ECHO) '         "requires" : {' >> META_new.json
	$(NOECHO) $(ECHO) '            "ExtUtils::MakeMaker" : "0"' >> META_new.json
	$(NOECHO) $(ECHO) '         }' >> META_new.json
	$(NOECHO) $(ECHO) '      },' >> META_new.json
	$(NOECHO) $(ECHO) '      "runtime" : {' >> META_new.json
	$(NOECHO) $(ECHO) '         "requires" : {' >> META_new.json
	$(NOECHO) $(ECHO) '            "Bio::AlignIO" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Bio::DB::EUtilities" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Bio::DB::Sam" : "1.39",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Bio::DB::SeqFeature::Store" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Bio::DB::Taxonomy" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Bio::Index::Fasta" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Bio::Root::Version" : "1.00690001",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Bio::SeqIO" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Bio::SimpleAlign" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Bio::Tree::Draw::Cladogram" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Bio::Tree::Tree" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Bio::TreeIO" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Cwd" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Encode" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "File::Basename" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "File::Path" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "File::Spec::Functions" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Getopt::Long" : "2.0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "IO::Pipe" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "IO::Select" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "IPC::Cmd" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "IPC::Open3" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "IPC::Run" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "List::MoreUtils" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "List::Util" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Math::Round" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Moose" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Moose::Role" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "POSIX" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Pod::Usage" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Switch" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Symbol" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Tree::Simple" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "Try::Tiny" : "0",' >> META_new.json
	$(NOECHO) $(ECHO) '            "perl" : "5.010"' >> META_new.json
	$(NOECHO) $(ECHO) '         }' >> META_new.json
	$(NOECHO) $(ECHO) '      }' >> META_new.json
	$(NOECHO) $(ECHO) '   },' >> META_new.json
	$(NOECHO) $(ECHO) '   "release_status" : "stable",' >> META_new.json
	$(NOECHO) $(ECHO) '   "version" : "v2.0.0"' >> META_new.json
	$(NOECHO) $(ECHO) '}' >> META_new.json
	-$(NOECHO) $(MV) META_new.json $(DISTVNAME)/META.json


# --- MakeMaker signature section:
signature :
	cpansign -s


# --- MakeMaker dist_basics section:
distclean :: realclean distcheck
	$(NOECHO) $(NOOP)

distcheck :
	$(PERLRUN) "-MExtUtils::Manifest=fullcheck" -e fullcheck

skipcheck :
	$(PERLRUN) "-MExtUtils::Manifest=skipcheck" -e skipcheck

manifest :
	$(PERLRUN) "-MExtUtils::Manifest=mkmanifest" -e mkmanifest

veryclean : realclean
	$(RM_F) *~ */*~ *.orig */*.orig *.bak */*.bak *.old */*.old 



# --- MakeMaker dist_core section:

dist : $(DIST_DEFAULT) $(FIRST_MAKEFILE)
	$(NOECHO) $(ABSPERLRUN) -l -e 'print '\''Warning: Makefile possibly out of date with $(VERSION_FROM)'\''' \
	  -e '    if -e '\''$(VERSION_FROM)'\'' and -M '\''$(VERSION_FROM)'\'' < -M '\''$(FIRST_MAKEFILE)'\'';' --

tardist : $(DISTVNAME).tar$(SUFFIX)
	$(NOECHO) $(NOOP)

uutardist : $(DISTVNAME).tar$(SUFFIX)
	uuencode $(DISTVNAME).tar$(SUFFIX) $(DISTVNAME).tar$(SUFFIX) > $(DISTVNAME).tar$(SUFFIX)_uu

$(DISTVNAME).tar$(SUFFIX) : distdir
	$(PREOP)
	$(TO_UNIX)
	$(TAR) $(TARFLAGS) $(DISTVNAME).tar $(DISTVNAME)
	$(RM_RF) $(DISTVNAME)
	$(COMPRESS) $(DISTVNAME).tar
	$(POSTOP)

zipdist : $(DISTVNAME).zip
	$(NOECHO) $(NOOP)

$(DISTVNAME).zip : distdir
	$(PREOP)
	$(ZIP) $(ZIPFLAGS) $(DISTVNAME).zip $(DISTVNAME)
	$(RM_RF) $(DISTVNAME)
	$(POSTOP)

shdist : distdir
	$(PREOP)
	$(SHAR) $(DISTVNAME) > $(DISTVNAME).shar
	$(RM_RF) $(DISTVNAME)
	$(POSTOP)


# --- MakeMaker distdir section:
create_distdir :
	$(RM_RF) $(DISTVNAME)
	$(PERLRUN) "-MExtUtils::Manifest=manicopy,maniread" \
		-e "manicopy(maniread(),'$(DISTVNAME)', '$(DIST_CP)');"

distdir : create_distdir distmeta 
	$(NOECHO) $(NOOP)



# --- MakeMaker dist_test section:
disttest : distdir
	cd $(DISTVNAME) && $(ABSPERLRUN) Makefile.PL 
	cd $(DISTVNAME) && $(MAKE) $(PASTHRU)
	cd $(DISTVNAME) && $(MAKE) test $(PASTHRU)



# --- MakeMaker dist_ci section:

ci :
	$(PERLRUN) "-MExtUtils::Manifest=maniread" \
	  -e "@all = keys %{ maniread() };" \
	  -e "print(qq{Executing $(CI) @all\n}); system(qq{$(CI) @all});" \
	  -e "print(qq{Executing $(RCS_LABEL) ...\n}); system(qq{$(RCS_LABEL) @all});"


# --- MakeMaker distmeta section:
distmeta : create_distdir metafile
	$(NOECHO) cd $(DISTVNAME) && $(ABSPERLRUN) -MExtUtils::Manifest=maniadd -e 'exit unless -e q{META.yml};' \
	  -e 'eval { maniadd({q{META.yml} => q{Module YAML meta-data (added by MakeMaker)}}) }' \
	  -e '    or print "Could not add META.yml to MANIFEST: $$$${'\''@'\''}\n"' --
	$(NOECHO) cd $(DISTVNAME) && $(ABSPERLRUN) -MExtUtils::Manifest=maniadd -e 'exit unless -f q{META.json};' \
	  -e 'eval { maniadd({q{META.json} => q{Module JSON meta-data (added by MakeMaker)}}) }' \
	  -e '    or print "Could not add META.json to MANIFEST: $$$${'\''@'\''}\n"' --



# --- MakeMaker distsignature section:
distsignature : create_distdir
	$(NOECHO) cd $(DISTVNAME) && $(ABSPERLRUN) -MExtUtils::Manifest=maniadd -e 'eval { maniadd({q{SIGNATURE} => q{Public-key signature (added by MakeMaker)}}) } ' \
	  -e '    or print "Could not add SIGNATURE to MANIFEST: $$$${'\''@'\''}\n"' --
	$(NOECHO) cd $(DISTVNAME) && $(TOUCH) SIGNATURE
	cd $(DISTVNAME) && cpansign -s



# --- MakeMaker install section:

install :: pure_install doc_install
	$(NOECHO) $(NOOP)

install_perl :: pure_perl_install doc_perl_install
	$(NOECHO) $(NOOP)

install_site :: pure_site_install doc_site_install
	$(NOECHO) $(NOOP)

install_vendor :: pure_vendor_install doc_vendor_install
	$(NOECHO) $(NOOP)

pure_install :: pure_$(INSTALLDIRS)_install
	$(NOECHO) $(NOOP)

doc_install :: doc_$(INSTALLDIRS)_install
	$(NOECHO) $(NOOP)

pure__install : pure_site_install
	$(NOECHO) $(ECHO) INSTALLDIRS not defined, defaulting to INSTALLDIRS=site

doc__install : doc_site_install
	$(NOECHO) $(ECHO) INSTALLDIRS not defined, defaulting to INSTALLDIRS=site

pure_perl_install :: all
	$(NOECHO) umask 022; $(MOD_INSTALL) \
		$(INST_LIB) $(DESTINSTALLPRIVLIB) \
		$(INST_ARCHLIB) $(DESTINSTALLARCHLIB) \
		$(INST_BIN) $(DESTINSTALLBIN) \
		$(INST_SCRIPT) $(DESTINSTALLSCRIPT) \
		$(INST_MAN1DIR) $(DESTINSTALLMAN1DIR) \
		$(INST_MAN3DIR) $(DESTINSTALLMAN3DIR)
	$(NOECHO) $(WARN_IF_OLD_PACKLIST) \
		$(SITEARCHEXP)/auto/$(FULLEXT)


pure_site_install :: all
	$(NOECHO) umask 02; $(MOD_INSTALL) \
		read $(SITEARCHEXP)/auto/$(FULLEXT)/.packlist \
		write $(DESTINSTALLSITEARCH)/auto/$(FULLEXT)/.packlist \
		$(INST_LIB) $(DESTINSTALLSITELIB) \
		$(INST_ARCHLIB) $(DESTINSTALLSITEARCH) \
		$(INST_BIN) $(DESTINSTALLSITEBIN) \
		$(INST_SCRIPT) $(DESTINSTALLSITESCRIPT) \
		$(INST_MAN1DIR) $(DESTINSTALLSITEMAN1DIR) \
		$(INST_MAN3DIR) $(DESTINSTALLSITEMAN3DIR)
	$(NOECHO) $(WARN_IF_OLD_PACKLIST) \
		$(PERL_ARCHLIB)/auto/$(FULLEXT)

pure_vendor_install :: all
	$(NOECHO) umask 022; $(MOD_INSTALL) \
		$(INST_LIB) $(DESTINSTALLVENDORLIB) \
		$(INST_ARCHLIB) $(DESTINSTALLVENDORARCH) \
		$(INST_BIN) $(DESTINSTALLVENDORBIN) \
		$(INST_SCRIPT) $(DESTINSTALLVENDORSCRIPT) \
		$(INST_MAN1DIR) $(DESTINSTALLVENDORMAN1DIR) \
		$(INST_MAN3DIR) $(DESTINSTALLVENDORMAN3DIR)

doc_perl_install :: all

doc_site_install :: all
	$(NOECHO) $(ECHO) Appending installation info to $(DESTINSTALLSITEARCH)/perllocal.pod
	-$(NOECHO) umask 02; $(MKPATH) $(DESTINSTALLSITEARCH)
	-$(NOECHO) umask 02; $(DOC_INSTALL) \
		"Module" "$(NAME)" \
		"installed into" "$(INSTALLSITELIB)" \
		LINKTYPE "$(LINKTYPE)" \
		VERSION "$(VERSION)" \
		EXE_FILES "$(EXE_FILES)" \
		>> $(DESTINSTALLSITEARCH)/perllocal.pod

doc_vendor_install :: all


uninstall :: uninstall_from_$(INSTALLDIRS)dirs
	$(NOECHO) $(NOOP)

uninstall_from_perldirs ::

uninstall_from_sitedirs ::
	$(NOECHO) $(UNINSTALL) $(SITEARCHEXP)/auto/$(FULLEXT)/.packlist

uninstall_from_vendordirs ::



# --- MakeMaker force section:
# Phony target to force checking subdirectories.
FORCE :
	$(NOECHO) $(NOOP)


# --- MakeMaker perldepend section:


# --- MakeMaker makefile section:
# We take a very conservative approach here, but it's worth it.
# We move Makefile to Makefile.old here to avoid gnu make looping.
$(FIRST_MAKEFILE) : Makefile.PL $(CONFIGDEP)
	$(NOECHO) $(ECHO) "Makefile out-of-date with respect to $?"
	$(NOECHO) $(ECHO) "Cleaning current config before rebuilding Makefile..."
	-$(NOECHO) $(RM_F) $(MAKEFILE_OLD)
	-$(NOECHO) $(MV)   $(FIRST_MAKEFILE) $(MAKEFILE_OLD)
	- $(MAKE) $(USEMAKEFILE) $(MAKEFILE_OLD) clean $(DEV_NULL)
	$(PERLRUN) Makefile.PL 
	$(NOECHO) $(ECHO) "==> Your Makefile has been rebuilt. <=="
	$(NOECHO) $(ECHO) "==> Please rerun the $(MAKE) command.  <=="
	$(FALSE)



# --- MakeMaker staticmake section:

# --- MakeMaker makeaperl section ---
MAP_TARGET    = perl
FULLPERL      = /usr/bin/perl

$(MAP_TARGET) :: static $(MAKE_APERL_FILE)
	$(MAKE) $(USEMAKEFILE) $(MAKE_APERL_FILE) $@

$(MAKE_APERL_FILE) : $(FIRST_MAKEFILE) pm_to_blib
	$(NOECHO) $(ECHO) Writing \"$(MAKE_APERL_FILE)\" for this $(MAP_TARGET)
	$(NOECHO) $(PERLRUNINST) \
		Makefile.PL DIR= \
		MAKEFILE=$(MAKE_APERL_FILE) LINKTYPE=static \
		MAKEAPERL=1 NORECURS=1 CCCDLFLAGS=


# --- MakeMaker test section:

TEST_VERBOSE=0
TEST_TYPE=test_$(LINKTYPE)
TEST_FILE = test.pl
TEST_FILES = 
TESTDB_SW = -d

testdb :: testdb_$(LINKTYPE)

test :: $(TEST_TYPE) subdirs-test

subdirs-test ::
	$(NOECHO) $(NOOP)

	$(NOECHO) $(ECHO) 'No tests defined for $(NAME) extension.'

test_dynamic :: pure_all

testdb_dynamic :: pure_all
	PERL_DL_NONLAZY=1 $(FULLPERLRUN) $(TESTDB_SW) "-I$(INST_LIB)" "-I$(INST_ARCHLIB)" $(TEST_FILE)

test_ : test_dynamic

test_static :: test_dynamic
testdb_static :: testdb_dynamic


# --- MakeMaker ppd section:
# Creates a PPD (Perl Package Description) for a binary distribution.
ppd :
	$(NOECHO) $(ECHO) '<SOFTPKG NAME="$(DISTNAME)" VERSION="$(VERSION)">' > $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '    <ABSTRACT>A Perl distribution for genome wide ncRNA annotation based on various software </ABSTRACT>' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '    <AUTHOR>Konstantin Riege</AUTHOR>' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '    <IMPLEMENTATION>' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <PERLCORE VERSION="5,010,0,0" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Bio::AlignIO" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Bio::DB::EUtilities" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE VERSION="1.39" NAME="Bio::DB::Sam" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Bio::DB::SeqFeature::Store" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Bio::DB::Taxonomy" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Bio::Index::Fasta" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE VERSION="1.00690001" NAME="Bio::Root::Version" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Bio::SeqIO" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Bio::SimpleAlign" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Bio::Tree::Draw::Cladogram" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Bio::Tree::Tree" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Bio::TreeIO" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Cwd::" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Encode::" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="File::Basename" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="File::Path" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="File::Spec::Functions" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE VERSION="2" NAME="Getopt::Long" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="IO::Pipe" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="IO::Select" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="IPC::Cmd" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="IPC::Open3" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="IPC::Run" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="List::MoreUtils" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="List::Util" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Math::Round" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Moose::" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Moose::Role" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="POSIX::" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Pod::Usage" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Switch::" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Symbol::" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Tree::Simple" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <REQUIRE NAME="Try::Tiny" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <ARCHITECTURE NAME="x86_64-linux-gnu-thread-multi-5.18" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '        <CODEBASE HREF="" />' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '    </IMPLEMENTATION>' >> $(DISTNAME).ppd
	$(NOECHO) $(ECHO) '</SOFTPKG>' >> $(DISTNAME).ppd


# --- MakeMaker pm_to_blib section:

pm_to_blib : $(FIRST_MAKEFILE) $(TO_INST_PM)
	$(NOECHO) $(ABSPERLRUN) -MExtUtils::Install -e 'pm_to_blib({@ARGV}, '\''$(INST_LIB)/auto'\'', q[$(PM_FILTER)], '\''$(PERM_DIR)'\'')' -- \
	  Gorap.pl $(INST_LIB)/Bio/Gorap.pl \
	  lib/Bio/Gorap/CFG.pm blib/lib/Bio/Gorap/CFG.pm \
	  lib/Bio/Gorap/DB/STK.pm blib/lib/Bio/Gorap/DB/STK.pm \
	  lib/Bio/Gorap/DB/GFF.pm blib/lib/Bio/Gorap/DB/GFF.pm \
	  lib/Bio/Gorap/DB/Fasta.pm blib/lib/Bio/Gorap/DB/Fasta.pm \
	  lib/Bio/Gorap/Update.pm blib/lib/Bio/Gorap/Update.pm \
	  lib/Bio/Gorap/Tool/Default.pm blib/lib/Bio/Gorap/Tool/Default.pm \
	  lib/Bio/Gorap/Gorap.pm blib/lib/Bio/Gorap/Gorap.pm \
	  lib/Bio/Gorap/ThrListener.pm blib/lib/Bio/Gorap/ThrListener.pm \
	  lib/Bio/Gorap/DB/Taxonomy.pm blib/lib/Bio/Gorap/DB/Taxonomy.pm \
	  lib/Bio/Gorap/Tool/Infernal.pm blib/lib/Bio/Gorap/Tool/Infernal.pm \
	  lib/Bio/Gorap/Tool/Bcheck.pm blib/lib/Bio/Gorap/Tool/Bcheck.pm \
	  lib/Bio/Gorap/Evaluation/Clustering.pm blib/lib/Bio/Gorap/Evaluation/Clustering.pm \
	  lib/Bio/Gorap/ToolI.pm blib/lib/Bio/Gorap/ToolI.pm \
	  RefreshOutput.pl $(INST_LIB)/Bio/RefreshOutput.pl \
	  lib/Bio/Gorap/Tool/Trnascanse.pm blib/lib/Bio/Gorap/Tool/Trnascanse.pm \
	  lib/Bio/Gorap/Functions/ToolParser.pm blib/lib/Bio/Gorap/Functions/ToolParser.pm \
	  lib/Bio/Gorap/Functions/STK.pm blib/lib/Bio/Gorap/Functions/STK.pm \
	  lib/Bio/Gorap/Evaluation/HTML.pm blib/lib/Bio/Gorap/Evaluation/HTML.pm \
	  lib/Bio/Gorap/Tool/Rnammer.pm blib/lib/Bio/Gorap/Tool/Rnammer.pm \
	  Gorap_noTaxonomy.pl $(INST_LIB)/Bio/Gorap_noTaxonomy.pl \
	  lib/Bio/Gorap/Tool/Blast.pm blib/lib/Bio/Gorap/Tool/Blast.pm \
	  lib/Bio/Gorap/DB/BAM.pm blib/lib/Bio/Gorap/DB/BAM.pm \
	  lib/Bio/Gorap/Tool/Crt.pm blib/lib/Bio/Gorap/Tool/Crt.pm \
	  lib/Bio/Gorap/Parameter.pm blib/lib/Bio/Gorap/Parameter.pm \
	  lib/Bio/Gorap/Functions/CM.pm blib/lib/Bio/Gorap/Functions/CM.pm \
	  Evaluation.pl $(INST_LIB)/Bio/Evaluation.pl \
	  lib/Bio/Gorap/Evaluation/Statistics.pm blib/lib/Bio/Gorap/Evaluation/Statistics.pm 
	$(NOECHO) $(TOUCH) pm_to_blib


# --- MakeMaker selfdocument section:


# --- MakeMaker postamble section:


# End.
