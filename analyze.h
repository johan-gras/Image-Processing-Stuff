#ifndef _ANALYZE_H_
#define _ANALYZE_H_

byte** etiquetage_intuitif(byte** img);

void etiquetage_batch(int size, byte*** batch, byte*** batch_result);

void caracterisation(byte** img, rgb8** rgb, byte** grey, int* nb_region, int* size, int barycentre[][2], int ecart_type[][3], double *axe, int* avgGrey, int avgColor[][3], int histGrey[][256], int range[][4]);

void generate_json();

#endif