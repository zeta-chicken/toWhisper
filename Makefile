PROGRAM = toWhisper
OBJS = vector.o wave.o toWhisper.o
CC = gcc
CFLAGS = -Wall -O3

.SUFFIXES: .o .c

$(PROGRAM): $(OBJS)
	$(CC) -o $(PROGRAM) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) $(PROGRAM) $(OBJS)
