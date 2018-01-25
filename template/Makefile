
CC=g++
CFLAGS=-W -Wall -ansi -pedantic
LDFLAGS=-lstdc++
EXEC=main
DIRSRC=./src
SRC= $(wildcard $(DIRSRC)/*.cpp)
INCL=./include
OBJ= $(SRC:.c=.o)

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
