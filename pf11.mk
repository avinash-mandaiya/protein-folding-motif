PROGRAM = pf11
FILES.c = pro_fold11.c motifproj.c 
FILES.h = pro_fold11.h 
FILES.o = ${FILES.c:.c=.o}
	
CC      = gcc
SFLAGS  = -std=c11
OFLAGS  = -O
CFLAGS  = ${SFLAGS} ${OFLAGS} 
LDFLAGS = -lm

all: ${PROGRAM} 

${PROGRAM}: ${FILES.o}
	${CC} -o $@ ${CFLAGS} ${FILES.o} ${LDFLAGS}

pro_fold8.o: ${FILES.h}
motifproj.o: ${FILES.h}

ProID        = 5DMA
fracSP       = 1.0
fracOS       = 1.0
totHBB       = 150
totHBS       = 50
CHstep       = 10
id           = t1
beta         = 0.5
maxiter      = 20000
iterstride   = 100
stoperr      = 0.0005
epsilon      = 0.01
seed         = 0

run: ${PROGRAM}
	$(addprefix ./,${PROGRAM}) ${ProID} ${fracSP} ${fracOS} ${totHBB} ${totHBS} ${CHstep} ${id} ${beta} ${maxiter} ${iterstride} ${stoperr} ${epsilon} ${seed}

.PHONY: run 
