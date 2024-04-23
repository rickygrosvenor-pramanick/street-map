# debug or release
CONF=debug

PROG := ssmap
CFLAGS := -Wall -std=gnu99
LOADLIBS := -lm

ifeq ($(CONF),debug)
CFLAGS += -g -O0 -ggdb
else ifeq ($(CONF),release)
CFLAGS += -O3
else
$(error CONF must be either debug or release)
endif

SOURCES := $(wildcard *.c)
OBJECTS := $(SOURCES:.c=.o)

all: depend $(PROG)

$(PROG): $(OBJECTS)
	$(CC) -o $@ $(CFLAGS) $^ $(LOADLIBS)

.PHONY: clean zip
clean:
	rm -f *.o depend.mk $(PROG) *.exe *.stackdump *~

zip: clean
	tar cvf ../a2-$(notdir $(shell pwd)).tar * 

depend:
	$(CC) -MM $(SOURCES) > depend.mk

ifeq (depend.mk,$(wildcard depend.mk))
include depend.mk
endif