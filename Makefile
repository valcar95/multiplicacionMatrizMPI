CC = mpicc
OPTS = -Wall

DEBUG ?= 1
ifeq ($(DEBUG), 1)
    OPTS += -g -O0
else
    OPTS += -O3
endif

all: mat_vec_mult_mpi.out take_time.out

mat_vec_mult_mpi.out: mat_vec_mult_mpi.o
        $(CC) -o $@ $< $(LIBS)

take_time.out: take_time.o
        $(CC) -o $@ $< $(LIBS)

%.o: %.c
        $(CC) $(OPTS) -c $< -o $@ $(LIBS)

clean:
        rm -f *.o *.out *~

run:
        mpirun --oversubscribe -np 4 ./mat_vec_mult_mpi.out
        mpirun --oversubscribe -np 4 ./take_time.out

run-dist:
        mpirun --oversubscribe -np 16 -host localhost,worker1,worker2,worker3 -mca pml ob1 -mca btl tcp,self -mca btl_tcp_if_include br0,lo ./take_time.out