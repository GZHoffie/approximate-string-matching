FOLDER_BIN=bin
FOLDER_BUILD=build

UNAME=$(shell uname)

CC=gcc
CPP=g++

LD_FLAGS=-lm
CC_FLAGS=-Wall -g
ifeq ($(UNAME), Linux)
  LD_FLAGS+=-lrt 
endif

AR=ar
AR_FLAGS=-rsc

SUBDIRS=benchmark/util
       
LIB_WFA=$(FOLDER_BUILD)/libwfa.a
LIB_WFA_SO=$(FOLDER_BUILD)/libwfa.so


all: CC_FLAGS+=-O3 -fPIC
all: MODE=all
all: setup
all: $(SUBDIRS)



setup:
	@mkdir -p $(FOLDER_BIN) $(FOLDER_BUILD)


clean:
	rm -rf $(FOLDER_BIN) $(FOLDER_BUILD)
	
###############################################################################
# Subdir rule
###############################################################################
export
$(SUBDIRS):
	$(MAKE) --directory=$@ all
	
tools: $(SUBDIRS)
	$(MAKE) --directory=$@ $(MODE)

.PHONY: $(SUBDIRS) tools;