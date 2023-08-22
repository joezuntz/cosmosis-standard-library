#include "errors.h"
#include <gsl/gsl_errno.h>
#include <stdio.h>
gsl_error_handler_t * original_handler = NULL;

void replacement_handler (const char * reason, const char * file, int line, int gsl_errno){
    fprintf(stderr, "GSL error: %s\n", reason);
    fprintf(stderr, "Location: %s:%d\n", file, line);
    fprintf(stderr, "GSL error code: %d\n", gsl_errno);
    fprintf(stderr, "error: %s\n", gsl_strerror(gsl_errno));

}

void enable_gsl_error_handling() {
    original_handler = gsl_set_error_handler(&replacement_handler);

}

void disable_gsl_error_handling() {
    gsl_set_error_handler(original_handler);
}
