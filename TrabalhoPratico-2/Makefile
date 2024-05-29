CC = g++
CFLAGS = -Wall -c
LFLAGS = -Wall
SRCDIR = src
OBJDIR = obj
BINDIR = bin

all: $(BINDIR)/main

$(BINDIR)/main: $(OBJDIR)/main.o
	$(CC) $(LFLAGS) $(OBJDIR)/main.o -o $(BINDIR)/main

$(OBJDIR)/main.o: $(SRCDIR)/main.cpp
	$(CC) $(CFLAGS) $(SRCDIR)/main.cpp -o $(OBJDIR)/main.o
	$(CC) $(CFLAGS) $(SRCDIR)/convexHull.cpp -o $(OBJDIR)/convexHull.o

clean:
	rm -rf $(BINDIR)/* $(OBJDIR)/*

run: all
	./$(BINDIR)/main $(word 3, $(MAKECMDGOALS))

%:
	@: