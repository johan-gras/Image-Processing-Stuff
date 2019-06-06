#include<stdio.h>
#include<stdlib.h>
#include <math.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "tools.h"
#include "morpho.h"


void change_eti(int **m, long max_x, long max_y, int a, int b){
    long i, j;
    for (j=0 ; j<=max_y ; j++) {
        for (i = 0; i <= max_x; i++) {
            if (m[j][i] == a)
                m[j][i] = b;
        }
    }
}


byte** etiquetage_intuitif(byte** img){
    int** m=imatrix(0, max_y, 0, max_x);
    byte** m_final=bmatrix(0, max_y, 0, max_x);
    long i, j;
    int a, b, eti = 1, region = 1;
    byte correspondance[255];

    for (i=0 ; i<255 ; i++)
        correspondance[i] = 255;

    for (j=0 ; j<=max_y ; j++){
        for(i=0 ; i<=max_x ; i++)
            m[j][i] = 0;
    }

    for (j=0 ; j<=max_y ; j++) {
        for (i = 0; i <= max_x; i++){
            if (img[j][i] == 255){
                a = i-1>=0 ? m[j][i-1] : 0;
                b = j-1>=0 ? m[j-1][i] : 0;

                if (a==0 && b==0){
                    m[j][i] = eti;
                    eti++;
                }
                else if (a==b || a==0)
                    m[j][i] = b;
                else if (b==0)
                    m[j][i] = a;
                else{
                    m[j][i] = b;
                    change_eti(m, max_x, max_y, a, b);
                }
            }
        }
    }

    for (j=0 ; j<=max_y ; j++){
        for(i=0 ; i<=max_x ; i++){
            if (m[j][i] == 0)
                m_final[j][i] = 255;
            else {
                if (correspondance[m[j][i]] == 255){
                    correspondance[m[j][i]] = region;
                    m_final[j][i] = region;
                    region++;
                } else
                    m_final[j][i] = correspondance[m[j][i]];
            }
        }
    }

    free_imatrix(m,0,max_y,0,max_x);
    return m_final;
}


void etiquetage_batch(int size, byte*** batch, byte*** batch_result){
    int i;

    for (i=0 ; i<size ; i++){
        batch_result[i] = etiquetage_intuitif(batch[i]);
    }
}


void caracterisation(byte** img, rgb8** rgb, byte** grey, int* nb_region, int* size, int barycentre[][2], int ecart_type[][3], double *axe, int* avgGrey, int avgColor[][3], int histGrey[][256], int range[][4]){
    int i, j, z, region;

    for (i=0 ; i<255 ; i++){
        range[i][0] = 1000000;
        range[i][2] = 1000000;
    }

    for (j=0 ; j<=max_y ; j++) {
        for (i = 0; i <= max_x; i++){
            region = img[j][i];
            if (region != 255)
            {
                // Nb Region
                if (region > *nb_region)
                    *nb_region = region;

                // Size calculation
                size[region-1]++;

                // Barycentre sum
                barycentre[region-1][0] += i;
                barycentre[region-1][1] += j;

                // Range limit calculation
                if (range[region-1][0] > i)
                    range[region-1][0] = i;
                else if (range[region-1][1] < i)
                    range[region-1][1] = i;
                if (range[region-1][2] > j)
                    range[region-1][2] = j;
                else if (range[region-1][3] < j)
                    range[region-1][3] = j;
            }
        }
    }

    for (z=0 ; z<*nb_region ; z++){
        // Barycentre avg
        barycentre[z][0] /= size[z];
        barycentre[z][1] /= size[z];

        for (j=0 ; j<=max_y ; j++) {
            for (i = 0; i <= max_x; i++) {
                if (img[j][i] == z+1){
                    // Ecart Type Calculation
                    ecart_type[z][0] += (i - barycentre[z][0]) * (i - barycentre[z][0]);
                    ecart_type[z][1] += (j - barycentre[z][1]) * (j - barycentre[z][1]);
                    ecart_type[z][2] += (i - barycentre[z][0]) * (j - barycentre[z][1]);

                    // Grey and Color Calculation
                    avgGrey[z] += grey[j][i];
                    avgColor[z][0] += rgb[j][i].r;
                    avgColor[z][1] += rgb[j][i].g;
                    avgColor[z][2] += rgb[j][i].b;
                    histGrey[z][grey[j][i]]++;
                }
            }
        }
        // Ecart Type Calculation
        ecart_type[z][0] /= size[z];
        ecart_type[z][1] /= size[z];
        ecart_type[z][2] /= size[z];

        // Grey and Color Calculation
        avgGrey[z] /= size[z];
        avgColor[z][0] /= size[z];
        avgColor[z][1] /= size[z];
        avgColor[z][2] /= size[z];

        // Axe calculation
        axe[z] = (180.0 / M_PI) * atan( (double)(2 * ecart_type[z][2]) / (ecart_type[z][0] - ecart_type[z][1]) );
    }
}

void generate_json(){
    int a = 1, b = 870;
    byte** batch[b-a+1];
    byte** clean_batch[b-a+1];
    byte** etiq_batch[b-a+1];
    byte** batch_grey[b-a+1];
    rgb8** batch_rgb[b-a+1];

    // Median
    byte** median = median_temp(a, b);
    SavePGM_bmatrix(median,nrl,max_y,ncl,max_x,"median.pgm");

    // Batches
    move_binary_ref(35, a, b, median, batch);
    clean_binary(b-a+1, batch, clean_batch);
    etiquetage_batch(b-a+1, clean_batch, etiq_batch);
    save_batch(a, b, "Result_Move\\r%03d.pgm", batch);
    save_batch(a, b, "Result_Clean\\r%03d.pgm", clean_batch);
    save_batch(a, b, "Result_Etiq\\r%03d.pgm", etiq_batch);
    load_batch_grey(a, b, "Sequences\\Lbox\\ppm\\lbox%03d.ppm", batch_grey);
    load_batch_rgb(a, b, "Sequences\\Lbox\\ppm\\lbox%03d.ppm", batch_rgb);

    // Caracterisation
    for (long i=a; i<=b ; i++){
        int nb = 0;
        int size[255] = {0};
        int barycentre[255][2] = {0};
        int ecart_type[255][3] = {0};
        double axe[255] = {0};
        int avgGrey[255] = {0};
        int avgColor[255][3] = {0};
        int histGrey[255][256] = {0};
        int range[255][4] = {0};

        caracterisation(etiq_batch[i-a], batch_rgb[i-a], batch_grey[i-a], &nb, size, barycentre, ecart_type, axe, avgGrey, avgColor, histGrey, range);
        save_caracterisation(i, nb, size, barycentre, ecart_type, axe, avgGrey, avgColor, histGrey, range);
    }

    // Free
    free_bmatrix(median,nrl,max_y,ncl,max_x);
    free_batch(b-a+1, batch);
    free_batch(b-a+1, clean_batch);
    free_batch(b-a+1, etiq_batch);
    free_batch(b-a+1, batch_grey);
    free_batch_rgb(b-a+1, batch_rgb);

    /*
     * printf("Nb :%d       size : %d\n", nb, size[0]);
    printf("Bary x : %d          Bary y : %d\n", barycentre[0][0], barycentre[0][1]);
    printf("Ecart 1 : %d        2: %d        3: %d\n", ecart_type[0][0], ecart_type[0][1], ecart_type[0][2]);
    printf("Axe : %f\n", axe[0]);
    printf("AvgGrey : %d\n", avgGrey[0]);
    printf("AvgRGB : %d  %d  %d\n", avgColor[0][0], avgColor[0][1], avgColor[0][2]);
     */
}
