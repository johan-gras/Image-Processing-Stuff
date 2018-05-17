#include<stdio.h>
#include<stdlib.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "tools.h"


byte** morpho(byte** img, long max_x, long max_y, byte (*fun)(byte**, long, long)){
    long i, j;
    byte** m=bmatrix(0, max_y, 0, max_x);

    for(j=0 ; j<max_y ; j+=max_y-1) {
        for (i = 0; i < max_x; i++)
            m[j][i] = img[j][i];
    }

    for(i=0 ; i<max_x ; i+=max_x-1) {
        for (j = 1; j < max_y; j++)
            m[j][i] = img[j][i];
    }

    for(j=1 ; j<max_y-1 ; j++){
        for(i=1 ; i<max_x-1 ; i++)
            m[j][i] = (*fun)(img, i, j);
    }

    return m;
}


byte erosion_op(byte** img, long i, long j){
    return img[j][i]*img[j-1][i]*img[j+1][i]*img[j][i-1]*img[j][i+1];
}


byte dilatation_op(byte** img, long i, long j){
    return (img[j][i] || img[j-1][i] || img[j+1][i] || img[j][i-1] || img[j][i+1])* 255;
}


byte** composition(byte** img, long max_x, long max_y, long n, byte (*fun1)(byte**, long, long), byte (*fun2)(byte**, long, long)){
    long i;
    byte** m= morpho(img, max_x, max_y, fun1);
    byte** m_t;

    for (i=1 ; i<n ; i++){
        m_t = morpho(m, max_x, max_y, fun1);
        free_bmatrix(m,0,max_y,0,max_x);
        m=m_t;
    }

    for (i=0 ; i<n ; i++){
        m_t = morpho(m, max_x, max_y, fun2);
        free_bmatrix(m,0,max_y,0,max_x);
        m=m_t;
    }

    return m;
}


byte** ouverture(byte** img, long max_x, long max_y, long n){
    return composition(img, max_x, max_y, n, &erosion_op, &dilatation_op);
}


byte** fermeture(byte** img, long max_x, long max_y, long n){
    return composition(img, max_x, max_y, n, &dilatation_op, &erosion_op);
}


void move_binary(long seuil, int a, int b, byte*** tab2){
    byte** tab1[b-a+1];
    long i;

    load_batch_ppm(a, b, "Sequences\\Lbox\\ppm\\lbox%03d.ppm", tab1);

    for (i=a ; i<=b-1 ; i++){
        tab2[i-a] = diff_img(tab1[i-a], tab1[i-a+1], max_x, max_y, seuil);
    }

    save_batch(a, b-1, "Result\\r%03d.pgm", tab2);
    free_batch(b-a+1, tab1);
}


void move_binary_ref(long seuil, int a, int b, byte** ref, byte*** tab2){
    byte** tab1[b-a+1];
    long i;

    load_batch_ppm(a, b, "Sequences\\Lbox\\ppm\\lbox%03d.ppm", tab1);

    for (i=a ; i<=b ; i++){
        tab2[i-a] = diff_img(tab1[i-a], ref, max_x, max_y, seuil);
    }

    save_batch(a, b-1, "Result\\r%03d.pgm", tab2);
    free_batch(b-a+1, tab1);
}


void clean_binary(int size, byte ***tab, byte ***tab_r){
    byte** b_1;
    long i;

    for (i=0 ; i<size ; i++){
        b_1 = fermeture(tab[i], max_x, max_y, 5);
        tab_r[i] = ouverture(b_1, max_x, max_y, 4);
        free_bmatrix(b_1,nrl,max_y,ncl,max_x);
    }

    save_batch(1, size, "Result_Clean\\r%03d.pgm", tab_r);
}


byte** moyenne_temp(int a, int b){
    byte** tab[b-a];
    long i, j, time;

    //Init matrice
    load_batch_ppm(a, b, "Sequences\\Lbox\\ppm\\lbox%03d.ppm", tab);
    byte** moyenne = bmatrix(nrl, max_y, ncl, max_x);
    int **m = imatrix(nrl, max_y, ncl, max_x);
    for(j=0 ; j<=max_y ; j++) {
        for (i = 0; i <= max_x; i++)
            m[j][i] = 0;
    }

    //Somme
    for (time=a ; time<=b ; time++){

        for(j=0 ; j<=max_y ; j++){
            for(i=0 ; i<=max_x ; i++) {
                m[j][i] += tab[time - a][j][i];
            }
        }
    }

    //Normalisation
    for(j=0 ; j<=max_y ; j++) {
        for (i = 0; i <= max_x; i++) {
            moyenne[j][i] = m[j][i] / (b-a);
        }
    }

    //Freeing
    free_batch(b-a, tab);
    free_imatrix(m,nrl,max_y,ncl,max_x);

    return moyenne;
}

byte** median_temp(int a, int b){
    byte** tab[b-a];
    long i, j, l, time, v;

    //Init matrice
    load_batch_ppm(a, b, "Sequences\\Lbox\\ppm\\lbox%03d.ppm", tab);
    byte** median = bmatrix(nrl, max_y, ncl, max_x);;
    int ***hist = i3tensor(nrl, max_y, ncl, max_x, 0, 255);
    for(l=0; l<256 ; l++){
        for(j=0 ; j<=max_y ; j++) {
            for (i = 0; i <= max_x; i++)
                hist[j][i][l] = 0;
        }
    }

    //Histo
    for (time=a ; time<=b ; time++){

        for(j=0 ; j<=max_y ; j++){
            for(i=0 ; i<=max_x ; i++) {
                hist[j][i][ tab[time - a][j][i] ]++;
            }
        }
    }

    //Find median values
    for(j=0 ; j<=max_y ; j++){
        for(i=0 ; i<=max_x ; i++) {
            v = 0;
            for (l=0 ; (v<(b-a)/2)  ; l++){
                v+= hist[j][i][l];
            }
            median[j][i] = l-1;
        }
    }

    //Freeing
    free_batch(b-a, tab);
    free_i3tensor(hist, nrl, max_y, ncl, max_x, 0, 255);

    return median;
}
