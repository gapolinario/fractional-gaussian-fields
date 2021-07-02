WARN = -Wmissing-prototypes -Wall #-Winline
OPTI = -O3 -finline-functions -fomit-frame-pointer -DNDEBUG \
-fno-strict-aliasing --param max-inline-insns-single=1800
STD = -std=c99
CC = gcc
CCFLAGS = $(OPTI) $(WARN) $(STD)

fgf1: FGF1D.c
	$(CC) -o $@ FGF1D.c -lfftw3 -lm

fgf2: FGF2D.c
	$(CC) -o $@ FGF2D.c -lfftw3 -lm

fgf3: FGF3D.c
	$(CC) -o $@ FGF3D.c -lfftw3 -lm
