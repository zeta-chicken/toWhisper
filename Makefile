PROGRAM = toWhisper
OBJS = vector.o wave.o toWhisper.o
CC = gcc
CFLAGS = -Wall -O3 -std=gnu11
LDFLAGS=-lm

.SUFFIXES: .o .c

$(PROGRAM): $(OBJS)
	$(CC) -o $(PROGRAM) $^ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) $(PROGRAM) $(OBJS) $(LDFLAGS)
