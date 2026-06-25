#ifndef SITO_H
 #define SITO_H
 typedef int ta[17][17];
 typedef int tu[17];
#endif

int connected(int N, ta A);
int integral(int N, ta A);
void BMKdecode(char * BUFFOR, int *N, ta A);
