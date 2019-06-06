#ifndef _INTEREST_H_
#define _INTEREST_H_

int** convolution(int** img, int mask[3][3], int mask_normalizer);

int** dot(int** img, int** img2);

void harris();

#endif