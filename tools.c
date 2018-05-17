#include<stdio.h>
#include<stdlib.h>

#include "def.h"
#include "nrio.h"
#include "nralloc.h"
#include "tools.h"
#include "analyze.h"
#include "morpho.h"
#include "cJSON.h"


byte** rgb8_to_byte(rgb8** img, long max_x, long max_y){
    byte** m=bmatrix(0, max_y, 0, max_x);
    long i, j;

    for(j=0 ; j<=max_y ; j++){
        for(i=0 ; i<=max_x ; i++)
            m[j][i] = (img[j][i].g + img[j][i].g + img[j][i].r) / 3;
    }

    return m;
}


int** byte_to_int(byte** img){
    int** m=imatrix(0, max_y, 0, max_x);
    long i, j;

    for(j=0 ; j<=max_y ; j++){
        for(i=0 ; i<=max_x ; i++)
            m[j][i] = img[j][i];
    }

    return m;
}


byte** diff_img(byte **img1, byte **img2, long max_x, long max_y, long seuil){
    byte** m=bmatrix(0, max_y, 0, max_x);
    long i, j;

    for(j=0 ; j<=max_y ; j++){
        for(i=0 ; i<=max_x ; i++)
            m[j][i] = (abs(img1[j][i] - img2[j][i]) > seuil) * 255;
    }

    return m;
}


void load_batch_ppm(int a, int b, char* path_source, byte ***tab){
    rgb8** rgb;
    long i;
    char path[256];

    for (i=a ; i<=b ; i++){
        sprintf(path, path_source, i);
        rgb=LoadPPM_rgb8matrix(path,&nrl,&max_y,&ncl,&max_x);
        tab[i-a] = rgb8_to_byte(rgb, max_x, max_y);
        free_rgb8matrix(rgb,nrl,max_y,ncl,max_x);
    }
}

void load_batch_rgb(int a, int b, const char* path_source, rgb8 ***tab){
    long i;
    char path[256];

    for (i=a ; i<=b ; i++){
        sprintf(path, path_source, i);
        tab[i-a]=LoadPPM_rgb8matrix(path,&nrl,&max_y,&ncl,&max_x);
    }
}


void load_batch_grey(int a, int b, char* path_source, byte ***tab){
    rgb8** rgb;
    long i;
    char path[256];

    for (i=a ; i<=b ; i++){
        sprintf(path, path_source, i);
        rgb=LoadPPM_rgb8matrix(path,&nrl,&max_y,&ncl,&max_x);
        tab[i-a] = rgb8_to_byte(rgb, max_x, max_y);
        free_rgb8matrix(rgb,nrl,max_y,ncl,max_x);
    }
}


void save_batch(int a, int b, char* path_source, byte ***tab){
    long i;
    char path[256];

    for (i=a; i<=b ; i++){
        sprintf(path, path_source, i);
        SavePGM_bmatrix(tab[i-a],nrl,max_y,ncl,max_x,path);
    }
}


void save_batch_ppm(int a, int b, char* path_source, rgb8 ***tab){
    long i;
    char path[256];

    for (i=a; i<=b ; i++){
        sprintf(path, path_source, i);
        SavePPM_rgb8matrix(tab[i-a],nrl,max_y,ncl,max_x,path);
    }
}


void free_batch(int size, byte ***tab){
    long i;

    for (i=0; i<size ; i++){
        free_bmatrix(tab[i],nrl,max_y,ncl,max_x);
    }
}


void free_batch_rgb(int size, rgb8 ***tab){
    long i;

    for (i=0; i<size ; i++){
        free_rgb8matrix(tab[i],nrl,max_y,ncl,max_x);
    }
}


void save_caracterisation(int img, int nb_region, int* size, int barycentre[][2], int ecart_type[][3], double *axe, int* avgGrey, int avgColor[][3], int histGrey[][256], int range[][4]){
    cJSON *regions = NULL;
    cJSON *hist = NULL;
    cJSON *monitor = cJSON_CreateObject();
    char *string = NULL;
    int i, j;
    char path[256];

    if (cJSON_AddNumberToObject(monitor, "nb_region", nb_region) == NULL) exit(-1);

    regions = cJSON_AddArrayToObject(monitor, "regions");
    if (regions == NULL) exit(-1);

    for (i=0 ; i<nb_region ; i++){
        cJSON *region = cJSON_CreateObject();

        if (cJSON_AddNumberToObject(region, "size", size[i]) == NULL) exit(-1);
        if (cJSON_AddNumberToObject(region, "barycentre_x", barycentre[i][0]) == NULL) exit(-1);
        if (cJSON_AddNumberToObject(region, "barycentre_y", barycentre[i][1]) == NULL) exit(-1);

        if (cJSON_AddNumberToObject(region, "range_x_min", range[i][0]) == NULL) exit(-1);
        if (cJSON_AddNumberToObject(region, "range_x_max", range[i][1]) == NULL) exit(-1);
        if (cJSON_AddNumberToObject(region, "range_y_min", range[i][2]) == NULL) exit(-1);
        if (cJSON_AddNumberToObject(region, "range_y_max", range[i][3]) == NULL) exit(-1);

        if (cJSON_AddNumberToObject(region, "ecart_type_x2", ecart_type[i][0]) == NULL) exit(-1);
        if (cJSON_AddNumberToObject(region, "ecart_type_y2", ecart_type[i][1]) == NULL) exit(-1);
        if (cJSON_AddNumberToObject(region, "ecart_type_xy", ecart_type[i][2]) == NULL) exit(-1);

        if (cJSON_AddNumberToObject(region, "axe", axe[i]) == NULL) exit(-1);

        if (cJSON_AddNumberToObject(region, "avg_grey", avgGrey[i]) == NULL) exit(-1);
        if (cJSON_AddNumberToObject(region, "avg_R", avgColor[i][0]) == NULL) exit(-1);
        if (cJSON_AddNumberToObject(region, "avg_G", avgColor[i][1]) == NULL) exit(-1);
        if (cJSON_AddNumberToObject(region, "avg_B", avgColor[i][2]) == NULL) exit(-1);

        hist = cJSON_AddArrayToObject(region, "Histogram");
        if (hist == NULL) exit(-1);

        for (j=0 ; j<256 ; j++){
            cJSON_AddItemToArray(hist, cJSON_CreateNumber(histGrey[i][j]));
        }

        cJSON_AddItemToArray(regions, region);
    }

    string = cJSON_Print(monitor);
    if (string == NULL) {
        fprintf(stderr, "Failed to print monitor.\n");
    }

    cJSON_Delete(monitor);

    sprintf(path, "Data\\img%03d.json", img);
    FILE *fp = fopen(path, "w");
    if (fp != NULL)
    {
        fputs(string, fp);
        fclose(fp);
    }
}


void load_caracterisation(int img, int *nb_region, int* size, int **barycentre, int **ecart_type, double *axe, int* avgGrey, int **avgColor, int **histGrey, int **range){
    char path[256];
    char *fcontent = NULL;
    sprintf(path, "Data\\img%03d.json", img);

    FILE *f = fopen(path, "r");
    if(f) {
        fseek(f, 0, SEEK_END);
        int fsize = ftell(f);
        rewind(f);

        fcontent = (char*) malloc(sizeof(char) * fsize);
        fread(fcontent, 1, fsize, f);

        fclose(f);
    }

    fclose(f);
    cJSON *json = cJSON_Parse(fcontent);
    if (json == NULL)
        printf("Error parsing\n");

    const cJSON *nb = cJSON_GetObjectItemCaseSensitive(json, "nb_region");
    *nb_region = nb->valueint;

    const cJSON *regions = cJSON_GetObjectItemCaseSensitive(json, "regions");
    const cJSON *region = NULL;
    int r=0;
    cJSON_ArrayForEach(region, regions)
    {
        cJSON *size_j = cJSON_GetObjectItemCaseSensitive(region, "size");
        cJSON *barycentre_x = cJSON_GetObjectItemCaseSensitive(region, "barycentre_x");
        cJSON *barycentre_y = cJSON_GetObjectItemCaseSensitive(region, "barycentre_y");
        cJSON *range_x_min = cJSON_GetObjectItemCaseSensitive(region, "range_x_min");
        cJSON *range_x_max = cJSON_GetObjectItemCaseSensitive(region, "range_x_max");
        cJSON *range_y_min = cJSON_GetObjectItemCaseSensitive(region, "range_y_min");
        cJSON *range_y_max = cJSON_GetObjectItemCaseSensitive(region, "range_y_max");
        cJSON *ecart_type_x2 = cJSON_GetObjectItemCaseSensitive(region, "ecart_type_x2");
        cJSON *ecart_type_y2 = cJSON_GetObjectItemCaseSensitive(region, "ecart_type_y2");
        cJSON *ecart_type_xy = cJSON_GetObjectItemCaseSensitive(region, "ecart_type_xy");
        cJSON *axe_j = cJSON_GetObjectItemCaseSensitive(region, "axe");
        cJSON *avg_grey_json = cJSON_GetObjectItemCaseSensitive(region, "avg_grey");
        cJSON *avg_r = cJSON_GetObjectItemCaseSensitive(region, "avg_R");
        cJSON *avg_g = cJSON_GetObjectItemCaseSensitive(region, "avg_G");
        cJSON *avg_b = cJSON_GetObjectItemCaseSensitive(region, "avg_B");

        size[r] = size_j->valueint;
        barycentre[r][0] = barycentre_x->valueint;
        barycentre[r][1] = barycentre_y->valueint;
        ecart_type[r][0] = ecart_type_x2->valueint;
        ecart_type[r][1] = ecart_type_y2->valueint;
        ecart_type[r][2] = ecart_type_xy->valueint;
        axe[r] = axe_j->valueint;
        avgGrey[r] = avg_grey_json->valueint;
        avgColor[r][0] = avg_r->valueint;
        avgColor[r][1] = avg_g->valueint;
        avgColor[r][2] = avg_b->valueint;
        range[r][0] = range_x_min->valueint;
        range[r][1] = range_x_max->valueint;
        range[r][2] = range_y_min->valueint;
        range[r][3] = range_y_max->valueint;

        r++;
    }

    cJSON_Delete(json);
}
