/** Atlantick Permutation and its auxiliary functions
 * Implementors: Atlantick Team
 * License: CC0, attribution kindly requested. Blame taken too,
 * but not liability.
 */

////////////////////////////////////////////////////////////////////////////////////////////
//  This implementation is for q=1 (One SuperBlock), 32-bits words and a block size of 192 //
////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdint.h>

// Round Constantes (for 26 rounds at most)
uint32_t Rc[] = {0x100cd809, 0x8cc76981, 0xbbfdf360, 0xef5210ad, 0x14deb5c6, 0x035e2f35, 0x1b5e8771, 0x527b79bf, 0x55b5b64c, 0x463314df, 0x3aee3d1a, 0x05221a65, 0x9a550a2d, 0xee3ed3bb, 0x75e0a6f4, 0x1d12d04e, 0x7f895ac6, 0x86a82676, 0x49855833, 0xf96a8950, 0x9b15d051, 0x004381ad, 0xc58fe08b, 0xa40763e1, 0x29254a93, 0xe4cf79bb};

// Left and Right Rotation constants

// ------- Contants a ----------
int a[] = {1, 9, 12, 5, 11, 15, 13, 3, 16, 6, 4, 8, 14, 7, 10, 2};

// ------- Contants b ----------
int b[] = {8, 12, 10, 6, 15, 5, 2, 3, 16, 1, 14, 4, 13, 9, 7, 11};

// ------- Contants e ----------
int e[] = {6, 5, 10, 2, 12, 15, 3, 11, 8, 16, 4, 7, 1, 9, 14, 13};

// ------- Contants h ----------
int h[] = {10, 12, 15, 5, 11, 2, 3, 13, 6, 16, 1, 8, 14, 4, 7, 9};

// ------- Contants v ----------
int v[] = {7, 10, 14, 16, 15, 2, 13, 11, 8, 6, 1, 4, 3, 5, 9, 12};

// ------- Contants w ----------
int w[] = {6, 9, 12, 16, 11, 4, 2, 14, 5, 1, 15, 8, 3, 7, 10, 13};

/*************************************************
* Names:        Mapd, MapDelta, MapD
*
* Description: Elementary Diffusion Functions.
*
* Arguments:   - uint32_t x: a 32 bits word
*              - int v: the right shift parameter
*              - int w: the left shift parameter
*
* Returns a 32 bits word (success)
**************************************************/
uint32_t Mapd(uint32_t x, int v, int w)
{
    return (~x >> v) ^ (x << w) ^ x;
}

uint32_t MapDelta(uint32_t x, int v, int w)
{
    return (x >> v) ^ (x << w) ^ x;
}

uint32_t MapD(uint32_t x, int v, int w)
{
    return (x >> v) ^ (x << w);
}

/*************************************************
* Names:        AtlantickLinearLayer
*
* Description: Atlantick Linear.
*
* Arguments:   - uint32_t State: a block (point to an array) 
*                       of 2*p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns a pointer to the output block of length 2*p 32 bits
*         words(success)
**************************************************/
uint32_t *AtlantickLinearLayer(uint32_t *State, int p)
{
    uint32_t C, l = 2, temp[l + 1];

    ///////////////////////////////
    State[2 * p - 1] = State[2 * p - 1] ^ Mapd(State[0], v[0], w[0]);
    ///////////////////////////////

    ////////////// MIXA in dimension 2p-1
    temp[0] = Mapd(State[2 * p - 1], v[1], w[1]);
    temp[1] = Mapd(State[1], v[2], w[2]);
    temp[2] = State[2];
    for (int i = 3; i <= 2 * p - 2; i++)
    {
        temp[2] ^= State[i];
    }

    temp[0] = temp[0] ^ temp[1] ^ temp[2];
    temp[0] = Mapd(temp[0], v[3], w[3]);

    State[0] = State[0] ^ temp[0];
    ///////////////////////////

    ////////////// MIXA in dimension 2
    temp[0] = Mapd(State[0], v[4], w[4]);
    temp[1] = Mapd(State[2 * p - 1], v[5], w[5]);

    C = temp[0] ^ temp[1];
    for (int j = 0; j < l; j++)
    {
        temp[j] = MapDelta(C, v[j % 16], w[j % 16]);
    }

    // Sequential implementation (Compute in parallel if wanted)
    int j, u;
    for (int k = 1; k <= 2 * p - 2; k++)
    {
        j = k % l;
        u = (k - j) / l;

        C = MapDelta(temp[j], e[u % 16], a[u % 15]);
        State[k] = State[k] ^ C;
    }
    ///////////////////////////

    ////////////// MIXA in dimension 2p-1
    temp[0] = Mapd(State[2 * p - 2], v[6], w[6]);
    temp[1] = Mapd(State[0], v[7], w[7]);
    temp[2] = State[2];
    for (int i = 3; i <= 2 * p - 2; i++)
    {
        temp[2] ^= State[i];
    }

    temp[0] = temp[0] ^ temp[1] ^ temp[2];
    temp[0] = Mapd(temp[0], v[8], w[8]);

    State[2 * p - 1] = State[2 * p - 1] ^ temp[0];
    ///////////////////////////

    ///////////////////////////
    State[0] = State[0] ^ Mapd(State[2 * p - 1], v[9], w[9]);
    ///////////////////////////

    return State;
}

/*************************************************
* Names:        InvAtlantickLinearLayer
*
* Description: Inverse Atlantick Linear.
*
* Arguments:   - uint32_t State: a block (point to an array) 
*                       of 2*p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns a pointer to the output block of length 2*p 32 bits
*         words(success)
**************************************************/
uint32_t *InvAtlantickLinearLayer(uint32_t *State, int p)
{
    uint32_t C, l = 2, temp[l + 1];

    ///////////////////////////
    State[0] = State[0] ^ Mapd(State[2 * p - 1], v[9], w[9]);
    ///////////////////////////

    ////////////// MIXA in dimension 2p-1
    temp[0] = Mapd(State[2 * p - 2], v[6], w[6]);
    temp[1] = Mapd(State[0], v[7], w[7]);
    temp[2] = State[2];
    for (int i = 3; i <= 2 * p - 2; i++)
    {
        temp[2] ^= State[i];
    }

    temp[0] = temp[0] ^ temp[1] ^ temp[2];
    temp[0] = Mapd(temp[0], v[8], w[8]);

    State[2 * p - 1] = State[2 * p - 1] ^ temp[0];
    ///////////////////////////

    ////////////// MIXA in dimension 2
    temp[0] = Mapd(State[0], v[4], w[4]);
    temp[1] = Mapd(State[2 * p - 1], v[5], w[5]);

    C = temp[0] ^ temp[1];
    for (int j = 0; j < l; j++)
    {
        temp[j] = MapDelta(C, v[j % 16], w[j % 16]);
    }

    // Sequential implementation (Compute in parallel if wanted)
    int j, u;
    for (int k = 2 * p - 2; k >= 1; k--)
    {
        j = k % l;
        u = (k - j) / l;

        C = MapDelta(temp[j], e[u % 16], a[u % 15]);
        State[k] = State[k] ^ C;
    }
    ///////////////////////////

    ////////////// MIXA in dimension 2p-1
    temp[0] = Mapd(State[2 * p - 1], v[1], w[1]);
    temp[1] = Mapd(State[1], v[2], w[2]);
    temp[2] = State[2];
    for (int i = 3; i <= 2 * p - 2; i++)
    {
        temp[2] ^= State[i];
    }

    temp[0] = temp[0] ^ temp[1] ^ temp[2];
    temp[0] = Mapd(temp[0], v[3], w[3]);

    State[0] = State[0] ^ temp[0];
    ///////////////////////////

    ///////////////////////////////
    State[2 * p - 1] = State[2 * p - 1] ^ Mapd(State[0], v[0], w[0]);
    ///////////////////////////////

    return State;
}

/*************************************************
* Names:        NonLinear
*
* Description: Elementary non linear map.
*
* Arguments:   - uint32_t tempB: a block (point to an array) 
*                       of p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns 0 (success)
**************************************************/
int NonLinear(uint32_t *tempB, int p)
{
    uint32_t tmp1 = tempB[0], tmp2 = tempB[1];
    for (int i = 0; i <= p - 3; i++)
    {
        tempB[i] = (~tempB[i + 1]) & tempB[i + 2];
    }
    tempB[p - 2] = (~tempB[p - 1]) & tmp1;
    tempB[p - 1] = (~tmp1) & tmp2;
    return 0;
}

/*************************************************
* Names:        MIXASTAR
*
* Description: MIXA* map.
*
* Arguments:   - uint32_t tempR: a block (point to an array) 
*                       of p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns a pointer to the output block of length p 32 bits
*         words(success)
**************************************************/
uint32_t *MIXASTAR(uint32_t *tempR, int p)
{
    uint32_t C, temp[3], tempB[p], *tempo;
    ////////////// MIXA STAR
    temp[0] = Mapd(tempR[p - 1], v[10], w[10]);
    temp[1] = Mapd(tempR[0], v[11], w[11]);

    C = temp[0] ^ temp[1];
    temp[0] = MapD(C, v[0], w[0]);

    // Sequential implementation (Parallel implementation possible)
    for (int k = 0; k <= p - 1; k++)
    {
        tempB[k] = tempR[k] ^ MapD(temp[0], b[k % 16], h[k % 15]);
    }

    tempo = &tempB[0];
    return tempo;
}

/*************************************************
* Names:        MIXB
*
* Description: MIXB map.
*
* Arguments:   - uint32_t tempR: a block (point to an array) 
*                       of p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns a pointer to the output block of length p 32 bits
*         words(success)
**************************************************/
uint32_t *MIXB(uint32_t *tempR, int p)
{
    uint32_t *tempB;

    //////////////////// MIXA*
    tempB = MIXASTAR(tempR, p);
    //////////////////////////

    ////////////// NonLinear
    NonLinear(tempB, p);
    //////////////////////////

    return tempB;
}

/*************************************************
* Names:        SumBlock
*
* Description: Sum of two 32 bits words Block of size p.
*
* Arguments:   - uint32_t input1: a block (point to an array) 
*                       of p 32 bits words
*              - uint32_t input2: a block (point to an array) 
*                       of p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns 0 (success)
**************************************************/
int SumBlock(uint32_t input1[], uint32_t input2[], int p)
{
    for (int i = 0; i < p; i++)
    {
        input1[i] = input1[i] ^ input2[i];
    }
    return 0;
}

/*************************************************
* Names:        RotBlock
*
* Description: RotBlock map.
*
* Arguments:   - uint32_t tempB: a block (point to an array) 
*                       of p 32 bits words
*              - int p: a parameter that defines the
*                       block size             
*
* Returns a pointer to the output block of length p 32 bits
*         words(success)
**************************************************/
uint32_t *RotBlock(uint32_t *tempB, int p)
{
    uint32_t tmp[p], *output;
    tmp[p - 1] = tempB[0];
    for (int i = 0; i <= p - 2; i++)
        tmp[i] = tempB[i + 1];

    output = &tmp[0];
    return output;
}

/*************************************************
* Names:        AtlantickSPbox
*
* Description: Atlantick SPbox.
*
* Arguments:   - uint32_t State: a block (point to an array) 
*                       of 2*p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns a pointer to the output block of length 2*p 32 bits
*         words(success)
**************************************************/
uint32_t *AtlantickSPbox(uint32_t *State, int p)
{
    uint32_t tempL[p], tempR[p];
    // Sequential computation (Compute in parallel if wanted)

    for (int i = 0; i <= p - 1; i++)
    {
        tempL[i] = State[i];
        tempR[i] = State[p + i];
    }

    SumBlock(tempR, RotBlock(tempL, p), p);
    SumBlock(tempL, MIXB(tempR, p), p);
    SumBlock(tempR, MIXB(tempL, p), p);
    SumBlock(tempL, RotBlock(tempR, p), p);

    for (int i = 0; i <= p - 1; i++)
    {
        State[i] = tempL[i];
        State[p + i] = tempR[i];
    }
    return State;
}

/*************************************************
* Names:        InvAtlantickSPbox
*
* Description: Inverse Atlantick SPbox.
*
* Arguments:   - uint32_t State: a block (point to an array) 
*                       of 2*p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns a pointer to the output block of length 2*p 32 bits
*         words(success)
**************************************************/
uint32_t *InvAtlantickSPbox(uint32_t *State, int p)
{
    uint32_t tempL[p], tempR[p];
    // Sequential computation (Compute in parallel if wanted)
    for (int i = 0; i <= p - 1; i++)
    {
        tempL[i] = State[i];
        tempR[i] = State[p + i];
    }

    SumBlock(tempL, RotBlock(tempR, p), p);
    SumBlock(tempR, MIXB(tempL, p), p);
    SumBlock(tempL, MIXB(tempR, p), p);
    SumBlock(tempR, RotBlock(tempL, p), p);

    for (int i = 0; i <= p - 1; i++)
    {
        State[i] = tempL[i];
        State[p + i] = tempR[i];
    }
    return State;
}

/*************************************************
* Names:        AtlantickRoundFunction
*
* Description: Atlantick Round Permutation.
*
* Arguments:   - uint32_t M: a block (point to an array) 
*                       of 2*p 32 bits words
*              - int p: a parameter that defines the
*                       block size  
*              - int r: the current round
*
* Returns a pointer to the output block of length 2*p 32 bits
*         words(success)
**************************************************/
uint32_t *AtlantickRoundFunction(uint32_t *M, int p, int r)
{
    M = AtlantickLinearLayer(M, p);
    M[0] = M[0] ^ Rc[2 * r];
    M[p] = M[p] ^ (uint32_t)p ^ Rc[2 * r + 1];
    M = AtlantickSPbox(M, p);
    return M;
}

/*************************************************
* Names:        AtlantickPermutation
*
* Description: Atlantick Permutation.
*
* Arguments:   - uint32_t M: a block (point to an array) 
*                       of 2*p 32 bits words
*              - int p: a parameter that defines the
*                       block size
*              - int s: the number of rounds    
*
* Returns a pointer to the output block of length 2*p 32 bits
*         words(success)
**************************************************/
uint32_t *AtlantickPermutation(uint32_t *M, int p, int s)
{
    for (int r = 0; r <= s - 1; r++)
    {
        M = AtlantickRoundFunction(M, p, r);
    }
    return M;
}

/*************************************************
* Names:        InvAtlantickRoundFunction
*
* Description: Inverse Atlantick Round Permutation.
*
* Arguments:   - uint32_t M: a block (point to an array) 
*                       of 2*p 32 bits words
*              - int p: a parameter that defines the
*                       block size
*              - int r: the current round
*
* Returns a pointer to the output block of length 2*p 32 bits
*         words(success)
**************************************************/
uint32_t *InvAtlantickRoundFunction(uint32_t *M, int p, int r)
{
    M = InvAtlantickSPbox(M, p);
    M[p] = M[p] ^ (uint32_t)p ^ Rc[2 * r + 1];
    M[0] = M[0] ^ Rc[2 * r];
    M = InvAtlantickLinearLayer(M, p);
    return M;
}

/*************************************************
* Names:        InvAtlantickPermutation
*
* Description: Inverse Atlantick Permutation.
*
* Arguments:   - uint32_t M: a block (point to an array) 
*                       of 2*p 32 bits words
*              - int p: a parameter that defines the
*                       block size
*              - int s: the number of rounds    
*
* Returns a pointer to the output block of length 2*p 32 bits
*         words(success)
**************************************************/
uint32_t *InvAtlantickPermutation(uint32_t *M, int p, int s)
{
    for (int r = s - 1; r >= 0; r--)
    {
        M = InvAtlantickRoundFunction(M, p, r);
    }
    return M;
}
