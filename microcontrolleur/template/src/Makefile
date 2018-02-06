
CC=g++
CFLAGS=-W -Wall -ansi -pedantic
LDFLAGS=-lstdc++
EXEC=main
# SRC= $(wildcard *.c)
SRC= main.cpp Interface.cpp Verification.cpp Module.cpp Watcher.cpp
OBJ= $(SRC:.c=.o)
INCL=../include

all: $(EXEC)

main: $(OBJ)
	@$(CC) -I$(INCL) -o $@ $^ $(LDFLAGS)

%.o: %.c
	@$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
	@rm -rf *.o

mrproper: clean
	@rm -rf $(EXEC)
