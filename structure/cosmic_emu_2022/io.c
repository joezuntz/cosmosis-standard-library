#include "stdio.h"
#include "string.h"
#include "emu.h"


void dump_0d(FILE * f, double val){
    int c = sizeof(double);
    fwrite(&val, c, 1, f);
}

void dump_1d(FILE * f, int n, double val[n]){
    int c = sizeof(double);
    fwrite(val, c, n, f);
}
void dump_2d(FILE * f, int m, int n, double val[m][n]){
    int c = sizeof(double);
    for (int i=0; i<m; i++){
        fwrite(val[i], c, n, f);
    }
}

double read_0d(FILE * f){
    double val;
    int c = sizeof(double);
    fread(&val, c, 1, f);
    return val;
}

void read_1d(FILE * f, int n, double val[n]){
    int c = sizeof(double);
    fread(val, c, n, f);
}
void read_2d(FILE * f, int m, int n, double val[m][n]){
    int c = sizeof(double);
    for (int i=0; i<m; i++){
        fread(val[i], c, n, f);
    }
}




void dump_params(const char * filename, emu_params * p){
    FILE * f = fopen(filename, "w");
    dump_2d(f, 111, 8, p->x);
    dump_1d(f, 8, p->xmin);
    dump_1d(f, 8, p->xrange);
    dump_1d(f, 8, p->xmax);
    dump_1d(f, 8, p->z);
    dump_1d(f, 8, p->z_asc);
    dump_1d(f, 351, p->mode);
    dump_2d(f, 2808, 45, p->K);
    dump_2d(f, 45, 111, p->w);
    dump_1d(f, 2808, p->mean);
    dump_0d(f, p->sd);
    dump_2d(f, 45, 8, p->beta);
    dump_1d(f, 45, p->lamws);
    dump_1d(f, 45, p->lamz);
}

void read_params(const char * filename, emu_params * p){
    FILE * f = fopen(filename, "r");
    read_2d(f, 111, 8, p->x);
    read_1d(f, 8, p->xmin);
    read_1d(f, 8, p->xrange);
    read_1d(f, 8, p->xmax);
    read_1d(f, 8, p->z);
    read_1d(f, 8, p->z_asc);
    read_1d(f, 351, p->mode);
    read_2d(f, 2808, 45, p->K);
    read_2d(f, 45, 111, p->w);
    read_1d(f, 2808, p->mean);
    p->sd = read_0d(f);
    read_2d(f, 45, 8, p->beta);
    read_1d(f, 45, p->lamws);
    read_1d(f, 45, p->lamz);
    fclose(f);
}


// int main(void){

//     emu_params p;


//     memcpy(p.x, x, sizeof x);
//     memcpy(p.xmin, xmin, sizeof xmin);
//     memcpy(p.xrange, xrange, sizeof xrange);
//     memcpy(p.xmax, xmax, sizeof xmax);
//     memcpy(p.z, z, sizeof z);
//     memcpy(p.z_asc, z_asc, sizeof z_asc);
//     memcpy(p.mode, mode, sizeof mode);
//     memcpy(p.K, K, sizeof K);
//     memcpy(p.w, w, sizeof w);
//     memcpy(p.mean, mean, sizeof mean);
//     memcpy(p.beta, beta, sizeof beta);
//     memcpy(p.lamws, lamws, sizeof lamws);
//     memcpy(p.lamz, lamz, sizeof lamz);

//     p.sd = sd;

//     dump_params("P_tot.dat", &p);

//     return 0;
// }
