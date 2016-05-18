AR     = ar rcv
RANLIB = ranlib 

HEADERS = f2c.h blas.h blasp.h
FILES = dnrm2.o dasum.o daxpy.o dcopy.o ddot.o dscal.o dsymv.o dgemv.o dtrsv.o\
xerbla.o lsame.o dpotf2.o dgemm.o dsyrk.o dsysv.o dsytrf.o\
dsytrs.o ieeeck.o dlasyf.o dsytf2.o dlasyf.o daxpy.o dtrsv.o\
ieeeck.o dsyr.o dswap.o scopy.o dger.o idamax.o ilaenv.o s_cmp.o s_copy.o

CFLAGS = $(OPTFLAGS) 
FFLAGS = $(OPTFLAGS)

blas: $(FILES) $(HEADERS)
	$(AR) f2c.a $(FILES)  
	$(RANLIB) f2c.a

clean:
	- rm -f *.o
	- rm -f *.a
	- rm -f *~

.c.o:
	$(CC) $(CFLAGS) -c $*.c

#.f.o:   $(FC) $(FFLAGS) -c $*.f
