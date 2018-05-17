#ifndef _MORPHO_H_
#define _MORPHO_H_

byte** morpho(byte** img, long max_x, long max_y, byte (*fun)(byte**, long, long));

byte erosion_op(byte** img, long i, long j);

byte dilatation_op(byte** img, long i, long j);

byte** composition(byte** img, long max_x, long max_y, long n, byte (*fun1)(byte**, long, long), byte (*fun2)(byte**, long, long));

byte** ouverture(byte** img, long max_x, long max_y, long n);

byte** fermeture(byte** img, long max_x, long max_y, long n);

void move_binary(long seuil, int a, int b, byte*** tab2);

void move_binary_ref(long seuil, int a, int b, byte** ref, byte*** tab2);

void clean_binary(int size, byte ***tab, byte ***tab_r);

byte** moyenne_temp(int a, int b);

byte** median_temp(int a, int b);

#endif