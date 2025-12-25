
#include <stdio.h>
#include <stdlib.h>



int main(void) {
    int *a = NULL;
    unsigned int n = 1024;
    int i = 1024;
    a = (int *)malloc(n * sizeof(int));
    while (1) {
    i--;
    a[i] = 1000000;
    if (i < 0) break;
    }
    free(a);
    a = NULL;
    return 0;
}
