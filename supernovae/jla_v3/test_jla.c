void * setup(void * options);

#include <dlfcn.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "datablock/c_datablock.h"

int main(void) {
        /* Load the library */
        void* lib = dlopen("./jla.so", RTLD_NOW | RTLD_LOCAL);
	if (!lib) {
          printf("%s\n", dlerror());
	}
        assert(lib);

        /* Get a callable pointer to the setup function */
        void* fcn = dlsym(lib, "setup");
        assert(fcn);
        void*(*psetup)(void*) = (void*(*)(void*))(fcn);

        /* Create a datablock with the right entries */
        c_datablock* config = make_c_datablock();
        assert(config);
        c_datablock_put_double(config, "module_options", "scriptmcut", 10.0);
        c_datablock_put_string(config, "module_options", "data_dir", "data");
        c_datablock_put_string(config, "module_options", "data_file", "jla_lcparams.txt");
        c_datablock_put_string(config, "module_options", "mag_covmat_file", "jla_v0_covmatrix.dat");
        c_datablock_put_string(config, "module_options", "stretch_covmat_file", "jla_va_covmatrix.dat");
        c_datablock_put_string(config, "module_options", "colour_covmat_file", "jla_vb_covmatrix.dat");
        c_datablock_put_string(config, "module_options", "mag_stretch_covmat_file", "jla_v0a_covmatrix.dat");
        c_datablock_put_string(config, "module_options", "mag_colour_covmat_file", "jla_v0b_covmatrix.dat");
        c_datablock_put_string(config, "module_options", "stretch_colour_covmat_file", "jla_vab_covmatrix.dat");

        /* Test that our function behaves correctly */
        void* calculator = (psetup)(config);
        assert(calculator);
        free(calculator);

        /* Clean up */
        dlclose(lib);           
        return 0;
}
