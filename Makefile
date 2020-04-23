
# Define implicit build rules
.SUFFIXES :
%.o: %.F90
	$(MPIF90) -c $< $(FFLAGS) 

%.o: %.f90
	$(MPIF90) -c $< $(FFLAGS) 


# Include platform-specific parameters
include make.inc


OBJS = \
constants.o \
control.o \
decide_input.o \
error.o \
es_tools.o \
flinal_tools.o \
io_ascii.o \
io_espresso.o \
io_global.o \
io_w90.o \
io_wfc.o \
mp_global.o \
op_tools.o \
read_tools.o \
sc_readin.o \
sc_tools.o \
version.o \
write_kpt_info.o


# Build rules

all: sc.x

sc.x: shiftcurrent.o $(OBJS)
	$(LD) -o $@ $^ -liotk -lfftw3 $(IFLAGS) $(LIBS)

test: test.x

test.x: test.o $(OBJS)
	$(LD) $^ -o -liotk -lfftw3 $@ $(IFLAGS) $(LIBS) 


.PHONY : clean
clean:
	$(RM) *.o *.mod *.optrpt *.x

include make.depend

