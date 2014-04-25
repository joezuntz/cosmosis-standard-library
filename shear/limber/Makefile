include ${COSMOSIS_DIR}/config/compilers.mk
include ${COSMOSIS_DIR}/config/subdirs.mk


CFLAGS +=  -std=c99
LDFLAGS += -lgsl -lgslcblas  -lm

CFLAGS += -fPIC



AR = ar rc

OBJS = limber.o interp2d.o
LIBNAME = liblimber.a
SHLIBNAME = interface.so
EXEC=test_limber

default:lib  $(EXEC)
all: $(EXEC)

OBJS1=$(OBJS)

$(EXEC): lib test_limber.o
	$(CC) $(CFLAGS) -o $@ $(OBJS1)  test_limber.o $(LDFLAGS) -L. -llimber


lib: $(OBJS)
	$(AR) $(LIBNAME) $(OBJS)


.PHONY : clean

clean: 
	rm -rf *.o $(LIBNAME) $(LIBNAME).dSYM $(EXEC) $(EXEC).dSYM $(SHLIBNAME) $(SHLIBNAME).dSYM
