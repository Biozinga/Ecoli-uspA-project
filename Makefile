# Makefile

CC = gcc
CFLAGS = -Wall -Wextra -std=c11
INCLUDES = -Iinclude
SRCDIR = src
OBJDIR = obj
BINDIR = bin

SOURCES = $(wildcard $(SRCDIR)/*.c)
OBJECTS = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SOURCES))
EXEC = $(BINDIR)/projet_bioinfo

all: $(EXEC)

$(EXEC): $(OBJECTS)
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $^
	@echo "Compilation terminée. Exécutable : $(EXEC)"

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	rm -rf $(OBJDIR) $(BINDIR)
	@echo "Nettoyage terminé."
