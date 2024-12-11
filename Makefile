CC=gcc
CFLAGS=-c -Wall
LDFLAGS= -lrt -lpthread -lm
SOURCES= composite.c composite_calcs.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLES=composite 

all: $(SOURCES) $(EXECUTABLES)

# pull in dependency info for *existing* .o files
-include $(OBJECTS:.o=.d)

composite: composite.o composite_calcs.o
	$(CC) composite.o composite_calcs.o $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@
	$(CC) -MM $< > $*.d

clean:
	rm -rf $(OBJECTS) $(EXECUTABLES) *.d