CC=gcc
CFLAGS=
LDFLAGS=
LIBS= -lm
EXE=exe
NPROCS?=1

SRCS= question1.c utest.c
OBJS= $(patsubst %.c, %.o, $(SRCS))

%.o: %.cpp
	$(CC) -c $< $(CFLAGS)

$(EXE): main.o $(OBJS)
	$(CC) -o $@ main.o $(OBJS) $(LDFLAGS) $(LIBS)

all: $(EXE)

clean:
	$(RM) $(EXE) *.o

run:
	@if [ "$(NPROCS)" = "1" ]; then ./$(EXE); \
	  else mpirun -n $(NPROCS) ./$(EXE); fi
