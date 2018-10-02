#-------------------------------------------------------------------------------
#
# directories indicating where the source files are, where the object files
# should be created and where to save the resulting executable file

SRCSDIR := ../sources
OBJSDIR := ./obj
DESTDIR := .

#-------------------------------------------------------------------------------
#
# include configuration
#
$(info )
ifeq ($(wildcard make.config),make.config)
$(info Using customized file 'make.config'.)
$(info )
include make.config
else
$(warning Could not find customized file 'make.config'!)
$(info )
$(info File 'make.default' is an example with all available compilation time   \
       options.  You can modify it, but it is better to copy it to a new file  \
       'make.config' and then customize it.  This will also make disappear     \
       the following 15 second pause every time you compile.)
$(info )
$(shell sleep 15)
include make.default
endif

#-------------------------------------------------------------------------------
#
# host name
#
HOST := $(shell hostname)

include hosts/default
ifeq ($(wildcard hosts/$(HOST)),hosts/$(HOST))
$(info Using customized compiler setup from 'hosts/$(HOST)' file.)
$(info )
include hosts/$(HOST)
else
$(info Using default compiler setup from file 'hosts/default'.                 \
       Do not modify it!)
$(info Instead, copy this file to 'hosts/$(HOST)' and customize compilers and  \
       compilation options there.  This will also get rid of the following     \
       15 second pause every time you compile.)
$(info )
$(shell sleep 15)
endif

#-------------------------------------------------------------------------------
#
# check flag conditions
#

#-------------------------------------------------------------------------------
#
# pass additional flags to the compiler
#
# compiler
#
FFLAGS  += -D${COMPILER}

# precision
#
FFLAGS  += -DPREC=${PREC}

# add module path to compiler options
#
ifeq ($(COMPILER),GNU)
FFLAGS  += -J $(OBJSDIR)
else
FFLAGS  += -module $(OBJSDIR)
endif

#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .o

#-------------------------------------------------------------------------------

name     = particles

sources := $(wildcard $(SRCSDIR)/*.F90)
objects := $(sources:$(SRCSDIR)/%.F90=$(OBJSDIR)/%.o)
modules := $(sources:$(SRCSDIR)/%.F90=$(OBJSDIR)/%.mod)

all: $(name).x

$(name).x: $(objects) | $(DESTDIR)
	$(LD) $(LDFLAGS) $(objects) $(LIBS) -o $(DESTDIR)/$(name).x

$(OBJSDIR)/%.o : $(SRCSDIR)/%.F90
	$(FC) -c $(FFLAGS) $< -o $@

$(objects): | $(OBJSDIR)

$(OBJSDIR):
	mkdir -p $(OBJSDIR)

$(DESTDIR):
	mkdir -p $(DESTDIR)

makedeps: $(sources)
	./mkdeps.sh $(SRCSDIR) $(OBJSDIR) > makedeps

.PHONY: clean

clean:
	rm -f $(objects) $(modules) $(DESTDIR)/$(name).x makedeps
	if [ -d $(OBJSDIR) ]; then rmdir $(OBJSDIR); fi

#-------------------------------------------------------------------------------
#
# object file dependencies (generated using mkdeps.sh script)
#
include makedeps

#-------------------------------------------------------------------------------
