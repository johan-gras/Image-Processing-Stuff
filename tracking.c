#include <stdio.h>
#include <stdlib.h>
#include "def.h"
#include "morpho.h"
#include "tools.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "analyze.h"
#include "interest.h"
#include "tracking.h"

#define A 396
#define B 800

void focus(rgb8 **img, int *barycentre, int range[4], int lost){
    // Barycentre
    img[barycentre[1]][barycentre[0]].r = 255;

    // Vertical
    for (int x=range[0] ; x<=range[1] ; x++){
        if (lost){
            img[range[2]][x].b = 255;
            img[range[3]][x].b = 255;
        }
        else{
            img[range[2]][x].r = 255;
            img[range[3]][x].r = 255;
        }
    }
    // Horizontal
    for (int y=range[2] ; y<=range[3] ; y++){
        if (lost){
            img[y][range[0]].b = 255;
            img[y][range[1]].b = 255;
        }
        else{
            img[y][range[0]].r = 255;
            img[y][range[1]].r = 255;
        }
    }
}


int range_check(int *range_1, int *range_2){
    if (range_1[1] >= range_2[0] && range_1[0] <= range_2[1] && range_1[3] >= range_2[2] && range_1[2] <= range_2[3])
        return 1;
    return 0;
}


void tracking(){
    // Batch
    rgb8** batch[B-A+1];
    load_batch_rgb(A, B, PATH_SEQUENCE, batch);
    // Caracteristique
    int *nb_regions = (int*)calloc(B-A+1, sizeof(int));
    int **size = imatrix(0, B-A+1, 0, 255);
    int ***barycentre = i3tensor(0, B-A+1, 0, 255, 0, 2);
    int ***ecart_type = i3tensor(0, B-A+1, 0, 255, 0, 3);
    double **axe = dmatrix(0, B-A+1, 0, 255);
    int **avgGrey = imatrix(0, B-A+1, 0, 255);
    int ***avgColor = i3tensor(0, B-A+1, 0, 255, 0, 3);
    int **histGrey = imatrix(0, 255, 0, 256);
    int ***range = i3tensor(0, B-A+1, 0, 255, 0, 4);
    // Region tracking
    int region_track[B-A+1] = {0};
    region_track[0] = 0;
    int last_seen = 0;
    int nb_region_conflict = 0;
    int region_conflict[255] = {0};
    // Mean caracteristique
    int sumGrey = 0;
    int sumColor[3] = {0};
    int n = 0;

    // First image
    load_caracterisation(A, &nb_regions[0], size[0], barycentre[0], ecart_type[0], axe[0], avgGrey[0], avgColor[0], histGrey, range[0]);
    focus(batch[0], barycentre[0][region_track[0]], range[0][region_track[0]], 0);

    // Others images
    for (long i=A+1; i<=B ; i++){
        // Load new features
        load_caracterisation(i, &nb_regions[i - A], size[i - A], barycentre[i - A], ecart_type[i - A], axe[i - A], avgGrey[i - A], avgColor[i - A], histGrey, range[i - A]);

        // Select the more likely region based on features
        int region;
        for (region=0 ; region<nb_regions[i-A] ; region++){
            int barycentre_x = barycentre[i-A][region][0];
            int barycentre_y = barycentre[i-A][region][1];
            int *old_range = range[last_seen][region_track[last_seen]];

            // Barycentre within the old range
            if (barycentre_x >= old_range[0] && barycentre_x <= old_range[1] && barycentre_y >= old_range[2] && barycentre_y <= old_range[3]){
                region_conflict[nb_region_conflict] = region;
                nb_region_conflict++;
            }

            // Tracking lost
            if (last_seen != i-A-1){
                // Actual range within the old range
                if (range_check(range[i-A][region], old_range)){
                    region_conflict[nb_region_conflict] = region;
                    nb_region_conflict++;
                }
            }
        }

        // Region not found
        if (nb_region_conflict == 0) {
            region_track[i - A] = -1;
            focus(batch[i - A], barycentre[last_seen][region_track[last_seen]], range[last_seen][region_track[last_seen]], 1);
        }

        // Region found
        else {
            // Take a decision between conflicts
            if (nb_region_conflict > 1){
                int min = 10000000;
                for (int index_region=0 ; index_region<nb_region_conflict ; index_region++){
                    int temp = abs((sumColor[0] / n) - avgColor[i - A][region_conflict[index_region]][0]) +
                            abs((sumColor[1] / n) - avgColor[i - A][region_conflict[index_region]][1]) +
                            abs((sumColor[2] / n) - avgColor[i - A][region_conflict[index_region]][2]);
                    if (temp < min){
                        min = temp;
                        region = region_conflict[index_region];
                    }
                }
            } else
                region = region_conflict[0];

            sumGrey += avgGrey[i - A][region];
            sumColor[0] += avgColor[i - A][region][0];
            sumColor[1] += avgColor[i - A][region][1];
            sumColor[2] += avgColor[i - A][region][2];
            n++;
            region_track[i - A] = region;
            last_seen = i-A;
            focus(batch[i - A], barycentre[i - A][region_track[i-A]], range[i - A][region_track[i-A]], 0);
        }

        nb_region_conflict = 0;
    }

    save_batch_ppm(A, B, "Tracking\\tracking%03d.ppm", batch);

    // Free
    free_batch_rgb(B-A+1, batch);
    free(nb_regions);
    free_imatrix(size, 0, B-A+1, 0, 255);
    free_i3tensor(barycentre, 0, B-A+1, 0, 255, 0, 2);
    free_i3tensor(ecart_type, 0, B-A+1, 0, 255, 0, 3);
    free_dmatrix(axe, 0, B-A+1, 0, 255);
    free_imatrix(avgGrey, 0, B-A+1, 0, 255);
    free_i3tensor(avgColor, 0, B-A+1, 0, 255, 0, 3);
    free_imatrix(histGrey, 0, B-A+1, 0, 255);
    free_i3tensor(range, 0, B-A+1, 0, 255, 0, 4);
}
