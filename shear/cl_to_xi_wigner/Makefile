include ${COSMOSIS_SRC_DIR}/config/compilers.mk


USER_CFLAGS = -std=c99 -Wall -Wextra -Werror -Wno-unknown-pragmas -pedantic
USER_CFLAGS += -Ofast
USER_LDFLAGS = -lm


all: tessore_wigner.so

clean:
	rm -f tessore_wigner.so

tessore_wigner.so: transform.c wigner_d.c
	$(CC) -shared $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)


.PHONY: all clean
