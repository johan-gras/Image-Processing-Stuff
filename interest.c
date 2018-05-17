#include<stdio.h>
#include<stdlib.h>
#include <math.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "tools.h"

#define LAMBDA 0.004

int sobelX[3][3] = {{-1, 0, 1},
                       {-2, 0, 2},
                       {-1, 0, 1}};
int sobelX_normalizer = 4;

int sobelY[3][3] = {{-1, -2, -1},
                       {0,  0,  0},
                       {1,  2,  1}};
int sobelY_normalizer = 4;

int gaussien[3][3] = {{1,2,1},
                        {2,-4,2},
                        {1,2,1}};
int gaussien_normalizer = 12;


int** convolution(int** img, int mask[3][3], int mask_normalizer){
    int** m_final=imatrix(0, max_y, 0, max_x);

    for (int j = 0 ; j<=max_y ; j++) {
        for (int i = 0; i <= max_x; i++) {
            m_final[j][i] = 0;
        }
    }

    for (int j = 1 ; j<max_y ; j++){
        for (int i = 1 ; i<max_x ; i++){
            int sum = 0;

            for (int j2=0 ; j2<=2 ; j2++){
                for (int i2=0 ; i2<=2 ; i2++){
                    sum += img[j + j2 - 1][i + i2 - 1] * mask[j2][i2];
                }
            }
            m_final[j][i] = sum / mask_normalizer;
        }
    }

    return m_final;
}

int** dot(int** img, int** img2){
    int** m_final=imatrix(0, max_y, 0, max_x);

    for (int j = 0 ; j<=max_y ; j++){
        for (int i = 0 ; i<=max_x ; i++){
            m_final[j][i] = img[j][i] * img2[j][i];
        }
    }

    return m_final;
}

void harris(){
    int a = 1, b = 50;
    byte** batch[b-a+1];
    int** img;
    int **Ix, **Iy, **IxIy, **Ix_square, **Iy_square, **Ix_gaussian, **Iy_gaussian, **IxIy_gaussian, **harris;

    // Load batch
    load_batch_ppm(a, b, "Sequences\\Lbox\\ppm\\lbox%03d.ppm", batch);
    img = byte_to_int(batch[49]);

    // Calcul
    Ix = convolution(img, sobelX, sobelX_normalizer);
    Iy = convolution(img, sobelY, sobelY_normalizer);

    IxIy = dot(Ix, Iy);
    Ix_square = dot(Ix, Ix);
    Iy_square = dot(Iy, Iy);

    Ix_gaussian = convolution(Ix_square, gaussien, gaussien_normalizer);
    Iy_gaussian = convolution(Iy_square, gaussien, gaussien_normalizer);
    IxIy_gaussian = convolution(IxIy, gaussien, gaussien_normalizer);

    // Harris method
    harris = imatrix(0, max_y, 0, max_x);
    byte** test = bmatrix(0, max_y, 0, max_x);

    for (int j = 1 ; j<max_y ; j++) {
        for (int i = 1; i < max_x; i++) {
            harris[j][i] = Ix_gaussian[j][i] * Iy_gaussian[j][i] - IxIy_gaussian[j][i] - LAMBDA * pow(Ix_gaussian[j][i] + Iy_gaussian[j][i], 2);
            if (harris[j][i] > 1000000)
                test[j][i] = 255;
            else
                test[j][i] = 0;
        }
    }

    // Painting
    SavePGM_bmatrix(test,nrl,max_y,ncl,max_x,"test.ppm");

    // Free
    free_batch(b-a+1, batch);
    free_imatrix(img, 0, max_y, 0, max_x);
    free_imatrix(Ix, 0, max_y, 0, max_x);
    free_imatrix(Iy, 0, max_y, 0, max_x);
    free_imatrix(IxIy, 0, max_y, 0, max_x);
    free_imatrix(Ix_square, 0, max_y, 0, max_x);
    free_imatrix(Iy_square, 0, max_y, 0, max_x);
    free_imatrix(Ix_gaussian, 0, max_y, 0, max_x);
    free_imatrix(Iy_gaussian, 0, max_y, 0, max_x);
    free_imatrix(IxIy_gaussian, 0, max_y, 0, max_x);
    free_imatrix(harris, 0, max_y, 0, max_x);
}
