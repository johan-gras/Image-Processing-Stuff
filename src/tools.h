#ifndef _TOOLS_H_
#define _TOOLS_H_

static const char PATH_SEQUENCE[] = "Sequences\\Lbox\\ppm\\lbox%03d.ppm";

long nrl,max_y,ncl,max_x;

byte** rgb8_to_byte(rgb8** img, long max_x, long max_y);

int** byte_to_int(byte** img);

byte** diff_img(byte **img1, byte **img2, long max_x, long max_y, long seuil);

void load_batch_ppm(int a, int b, char* path_source, byte ***tab);

void load_batch_rgb(int a, int b, const char* path_source, rgb8 ***tab);

void load_batch_grey(int a, int b, char* path_source, byte ***tab);

void save_batch(int a, int b, char* path_source, byte ***tab);

void save_batch_ppm(int a, int b, char* path_source, rgb8 ***tab);

void free_batch(int size, byte ***tab);

void free_batch_rgb(int size, rgb8 ***tab);

void save_caracterisation(int img, int nb_region, int* size, int barycentre[][2], int ecart_type[][3], double *axe, int* avgGrey, int avgColor[][3], int histGrey[][256], int range[][4]);

void load_caracterisation(int img, int *nb_region, int* size, int **barycentre, int **ecart_type, double *axe, int* avgGrey, int **avgColor, int **histGrey, int **range);

#endif