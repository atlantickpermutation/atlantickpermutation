#include <stdio.h>
#include <stdint.h>

// Function for extended Euclidean Algorithm
int gcdExtended(int a, int b, int *x, int *y);

// Function to find modulo inverse of a
int modInverse(int A, int M)
{
    int x, y;
    int g = gcdExtended(A, M, &x, &y);
    if (g != 1)
        return 0;
    else
    {

        // m is added to handle negative x
        int res = (x % M + M) % M;
        return res;
    }
}

// Function for extended Euclidean Algorithm
int gcdExtended(int a, int b, int *x, int *y)
{

    // Base Case
    if (a == 0)
    {
        *x = 0, *y = 1;
        return b;
    }

    // To store results of recursive call
    int x1, y1;
    int gcd = gcdExtended(b % a, a, &x1, &y1);

    // Update x and y using results of recursive
    // call
    *x = y1 - (b / a) * x1;
    *y = x1;

    return gcd;
}

///////////////////////////////
// Shift constants a, b, d, e and h
//////////////////////////////

// Constants Function for a
uint32_t Constants_a(uint32_t constants_a[])
{
    for (int k = 0; k <= 15; k++)
    {
        constants_a[k] = 1 + ((1 + 15 * modInverse(k + 1, 17)) % 16);
    }
}

// Constants Function for b
uint32_t Constants_b(uint32_t constants_b[])
{
    for (int k = 0; k <= 15; k++)
    {
        constants_b[k] = 1 + ((2 + 13 * modInverse(2 * (k + 1), 17)) % 16);
    }
}

// Constants Function for e
uint32_t Constants_e(uint32_t constants_e[])
{
    for (int k = 0; k <= 15; k++)
    {
        constants_e[k] = 1 + ((3 + 11 * modInverse(3 * (k + 1), 17)) % 16);
    }
}

// Constants Function for h
uint32_t Constants_h(uint32_t constants_h[])
{
    for (int k = 0; k <= 15; k++)
    {
        constants_h[k] = 1 + ((4 + 9 * modInverse(4 * (k + 1), 17)) % 16);
    }
}

// Constants Function for v
uint32_t Constants_v(uint32_t constants_v[])
{
    for (int k = 0; k <= 15; k++)
    {
        constants_v[k] = 1 + ((5 + 7 * modInverse(5 * (k + 1), 17)) % 16);
    }
}

// Constants Function for w
uint32_t Constants_w(uint32_t constants_w[])
{
    for (int k = 0; k <= 15; k++)
    {
        constants_w[k] = 1 + ((6 + 5 * modInverse(6 * (k + 1), 17)) % 16);
    }
}

void main()
{
    int constants_a[16], constants_b[16], constants_e[16], constants_h[16], constants_v[16], constants_w[16];

    printf("// ------- Contants a ----------\n");
    Constants_a(constants_a);
    printf("int a[] = {");
    for (int k = 0; k <= 15; k++)
    {
        printf("%d, ", constants_a[k]);
    }
    printf("};");
    printf("\n\n");
    //----------------------------------------

    printf("// ------- Contants b ----------\n");
    Constants_b(constants_b);
    printf("int b[] = {");
    for (int k = 0; k <= 15; k++)
    {
        printf("%d, ", constants_b[k]);
    }
    printf("};");
    printf("\n\n");
    //----------------------------------------

    printf("// ------- Contants e ----------\n");
    Constants_e(constants_e);
    printf("int e[] = {");
    for (int k = 0; k <= 15; k++)
    {
        printf("%d, ", constants_e[k]);
    }
    printf("};");
    printf("\n\n");
    //----------------------------------------

    printf("// ------- Contants h ----------\n");
    Constants_h(constants_h);
    printf("int h[] = {");
    for (int k = 0; k <= 15; k++)
    {
        printf("%d, ", constants_h[k]);
    }
    printf("};");
    printf("\n\n");
    //----------------------------------------

    printf("// ------- Contants v ----------\n");
    Constants_v(constants_v);
    printf("int v[] = {");
    for (int k = 0; k <= 15; k++)
    {
        printf("%d, ", constants_v[k]);
    }
    printf("};");
    printf("\n\n");
    //----------------------------------------

    printf("// ------- Contants w ----------\n");
    Constants_w(constants_w);
    printf("int w[] = {");
    for (int k = 0; k <= 15; k++)
    {
        printf("%d, ", constants_w[k]);
    }
    printf("};");
    printf("\n\n");
    //----------------------------------------
}