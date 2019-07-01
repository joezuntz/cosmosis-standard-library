include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_CFLAGS = -std=c99 -Wall -Wextra -Werror -Wno-unknown-pragmas -pedantic
USER_CFLAGS += -Ofast  -march=native -mfpmath=sse
USER_LDFLAGS = -lm

.PHONY: all clean

all: wigner_d.so

clean:
	rm -f wigner_d.so

wigner_d.so: wigner_d.c
	$(CC) -shared $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)
