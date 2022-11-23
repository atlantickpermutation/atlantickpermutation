#include <stdio.h>
// Constants Function for w
int computeComplecity(int p, int l, int L, int s)
{
    int d = 5, delta = 4, D = 3;
    return s * (10 * d - 1 + 6 * p + delta * (l + 2 * p - 2)) + s * (14 * p + 4 * d - 2 + 2 * D * (p + L)) + 3 * s;
}

void main()
{

    printf("\n\np\tl\tL\ts\tz=2s\tsize\tOps_Nbr\n");
    printf("--------------------------------------------------\n");

    //printf("------- Complexity 192 ----------\n");
    printf("%d\t%d\t%d\t%d\t%4d\t%4d\t%d\n", 3, 2, 1, 7, 14, 192, computeComplecity(3, 2, 1, 7));
    //printf("-----------------------------------\n\n");

    //printf("------- Complexity 256 ----------\n");
    printf("%d\t%d\t%d\t%d\t%4d\t%4d\t%d\n", 4, 2, 1, 7, 14, 256, computeComplecity(4, 2, 1, 7));
    //printf("-----------------------------------\n\n");

    //printf("------- Complexity 320 ----------\n");
    printf("%d\t%d\t%d\t%d\t%4d\t%4d\t%d\n", 5, 2, 1, 7, 14, 320, computeComplecity(5, 2, 1, 7));
    //printf("-----------------------------------\n\n");

    //printf("------- Complexity 384 ----------\n");
    printf("%d\t%d\t%d\t%d\t%4d\t%4d\t%d\n", 6, 2, 1, 7, 14, 384, computeComplecity(6, 2, 1, 7));
    //printf("-----------------------------------\n\n");

    //printf("------- Complexity 512 ----------\n");
    printf("%d\t%d\t%d\t%d\t%4d\t%4d\t%d\n", 8, 2, 1, 7, 14, 512, computeComplecity(8, 2, 1, 7));
    //printf("-----------------------------------\n\n");

    //printf("------- Complexity 768 ----------\n");
    printf("%d\t%d\t%d\t%d\t%4d\t%4d\t%d\n", 12, 4, 1, 9, 18, 768, computeComplecity(12, 4, 1, 9));
    //printf("-----------------------------------\n\n");

    //printf("------- Complexity 1024 ----------\n");
    printf("%d\t%d\t%d\t%d\t%4d\t%4d\t%d\n", 16, 6, 1, 9, 18, 1024, computeComplecity(16, 6, 1, 9));
    //printf("-----------------------------------\n\n");

    //printf("------- Complexity 1600 ----------\n");
    printf("%d\t%d\t%d\t%d\t%4d\t%4d\t%d\n", 25, 6, 2, 9, 18, 1600, computeComplecity(25, 6, 2, 9));
    //printf("-----------------------------------\n\n");

    //printf("------- Complexity 2048 ----------\n");
    printf("%d\t%d\t%d\t%d\t%4d\t%4d\t%d\n", 32, 8, 2, 11, 22, 2048, computeComplecity(32, 6, 2, 11));
    //printf("-----------------------------------\n\n");

    //printf("------- Complexity 4096 ----------\n");
    printf("%d\t%d\t%d\t%d\t%4d\t%4d\t%d\n", 64, 8, 4, 11, 22, 4096, computeComplecity(64, 8, 4, 11));
    //printf("-----------------------------------\n\n");
}