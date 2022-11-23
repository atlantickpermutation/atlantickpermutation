########################################################################
# Atlantick Permutation and its auxiliary functions                    #
# Implementors: Atlantick Team                                         #
# License: CC0, attribution kindly requested. Blame taken too,         #        
# but not liability.                                                   #
########################################################################

########################################################################
#   This implementation is for q=1 (One SuperBlock) and 32-bit words   #
########################################################################

import numpy as np

# Round Constantes (for 26 rounds at most)
Rc = [
    0x100cd809, 0x8cc76981, 0xbbfdf360, 0xef5210ad, 0x14deb5c6, 0x035e2f35, 0x1b5e8771, 
    0x527b79bf, 0x55b5b64c, 0x463314df, 0x3aee3d1a, 0x05221a65, 0x9a550a2d, 0xee3ed3bb, 
    0x75e0a6f4, 0x1d12d04e, 0x7f895ac6, 0x86a82676, 0x49855833, 0xf96a8950, 0x9b15d051, 
    0x004381ad, 0xc58fe08b, 0xa40763e1, 0x29254a93, 0xe4cf79bb]

# Left and Right Rotation constants

# ------- Contants a ----------
a = [1, 9, 12, 5, 11, 15, 13, 3, 16, 6, 4, 8, 14, 7, 10, 2]

# ------- Contants b ----------
b = [8, 12, 10, 6, 15, 5, 2, 3, 16, 1, 14, 4, 13, 9, 7, 11]

# ------- Contants e ----------
e = [6, 5, 10, 2, 12, 15, 3, 11, 8, 16, 4, 7, 1, 9, 14, 13]

# ------- Contants h ----------
h = [10, 12, 15, 5, 11, 2, 3, 13, 6, 16, 1, 8, 14, 4, 7, 9]

# ------- Contants v ----------
v = [7, 10, 14, 16, 15, 2, 13, 11, 8, 6, 1, 4, 3, 5, 9, 12]

# ------- Contants w ----------
w = [6, 9, 12, 16, 11, 4, 2, 14, 5, 1, 15, 8, 3, 7, 10, 13]


def Mapd(x, v, w):
    return ((0xFFFFFFFF ^ x) >> v) ^ ((x << w) % (1 << 32)) ^ x


def MapDelta(x, v, w):
    return (x >> v) ^ ((x << w) % (1 << 32)) ^ x


def MapD(x, v, w):
    return (x >> v) ^ ((x << w) % (1 << 32))


#### Atlantick Linear
def AtlantickLinearLayer(State, p):
    l = 8; temp = np.zeros(l + 2, dtype=np.longlong)

    ###############/
    State[2 * p - 1] = State[2 * p - 1] ^ Mapd(State[0], v[0], w[0])
    ###############/

    ####### MIXA in dimension 2p-1
    temp[0] = Mapd(State[2 * p - 1], v[1], w[1])
    temp[1] = Mapd(State[1], v[2], w[2])
    temp[2] = State[2]
    for i in range(3, 2 * p - 1):
        temp[2] ^= State[i]

    temp[0] = temp[0] ^ temp[1] ^ temp[2]
    temp[0] = Mapd(temp[0], v[3], w[3])

    State[0] = State[0] ^ temp[0]
    #############/

    ####### MIXA in dimension 2
    temp[0] = Mapd(State[0], v[4], w[4])
    temp[1] = Mapd(State[2 * p - 1], v[5], w[5])

    C = temp[0] ^ temp[1]
    for j in range(0, l):
        temp[j] = MapDelta(C, v[j % 16], w[j % 16])

    # Sequential implementation (Compute in parallel if wanted)
    for k in range(1, 2 * p - 1):
        j = k % l
        u = (k - j) // l

        C = MapDelta(temp[j], e[u % 16], a[u % 15])
        State[k] = State[k] ^ C
    
    #############/

    ####### MIXA in dimension 2p-1
    temp[0] = Mapd(State[2 * p - 2], v[6], w[6])
    temp[1] = Mapd(State[0], v[7], w[7])
    temp[2] = State[2]
    for i in range(3, 2 * p - 1):
        temp[2] ^= State[i]
    
    temp[0] = temp[0] ^ temp[1] ^ temp[2]
    temp[0] = Mapd(temp[0], v[8], w[8])

    State[2 * p - 1] = State[2 * p - 1] ^ temp[0]
    #############/

    #############/
    State[0] = State[0] ^ Mapd(State[2 * p - 1], v[9], w[9])
    #############/

    return State


def InvAtlantickLinearLayer(State, p):
    l = 8; temp = np.zeros(l + 2, dtype=np.longlong)

    #############/
    State[0] = State[0] ^ Mapd(State[2 * p - 1], v[9], w[9])
    #############/

    ####### MIXA in dimension 2p-1
    temp[0] = Mapd(State[2 * p - 2], v[6], w[6])
    temp[1] = Mapd(State[0], v[7], w[7])
    temp[2] = State[2]
    for i in range(3, 2 * p - 1):
        temp[2] ^= State[i]

    temp[0] = temp[0] ^ temp[1] ^ temp[2]
    temp[0] = Mapd(temp[0], v[8], w[8])

    State[2 * p - 1] = State[2 * p - 1] ^ temp[0]
    #############/

    ####### MIXA in dimension 2
    temp[0] = Mapd(State[0], v[4], w[4])
    temp[1] = Mapd(State[2 * p - 1], v[5], w[5])

    C = temp[0] ^ temp[1]
    for j in range(0, l):
        temp[j] = MapDelta(C, v[j % 16], w[j % 16])

    # Sequential implementation (Compute in parallel if wanted)
    for k in range(2 * p - 2, 0, -1):
        j = k % l
        u = (k - j) // l

        C = MapDelta(temp[j], e[u % 16], a[u % 15])
        State[k] = State[k] ^ C
    #############/

    ####### MIXA in dimension 2p-1
    temp[0] = Mapd(State[2 * p - 1], v[1], w[1])
    temp[1] = Mapd(State[1], v[2], w[2])
    temp[2] = State[2]
    for i in range(3, 2 * p - 1):
        temp[2] ^= State[i]

    temp[0] = temp[0] ^ temp[1] ^ temp[2]
    temp[0] = Mapd(temp[0], v[3], w[3])

    State[0] = State[0] ^ temp[0]
    #############/

    ###############/
    State[2 * p - 1] = State[2 * p - 1] ^ Mapd(State[0], v[0], w[0])
    ###############/

    return State


######/ NonLinear
def NonLinear(tempB, p):
    tmp1 = tempB[0]; tmp2 = tempB[1]
    for i in range(p - 2):
        tempB[i] = (~tempB[i + 1]) & tempB[i + 2]

    tempB[p - 2] = (~tempB[p - 1]) & tmp1
    tempB[p - 1] = (~tmp1) & tmp2

######/ MIXA STAR
def MIXASTAR(tempR, p):
    L = 4
    temp = np.zeros(L, dtype=np.longlong)
    tempB = np.zeros(p, dtype=np.longlong)
    ####### MIXA STAR
    temp[0] = Mapd(tempR[p - 1], v[10], w[10])
    temp[1] = Mapd(tempR[0], v[11], w[11])

    temp[2] = tempR[1]
    for i in range(2, p - 1):
        temp[2] ^= tempR[i]

    C = temp[0] ^ temp[1] ^ temp[2]
    for j in range(L):
        temp[j] = MapD(C, v[j % 16], w[j % 16])

    # Sequential implementation (Parallel implementation possible)
    for k in range(p):
        j = k % L
        u = (k - j) // L
        tempB[k] = tempR[k] ^ MapD(temp[j], b[u % 16], h[u % 15])
    
    return tempB


######/MIXB
def MIXB(tempR, p):
    ####### First MIXA*
    tempB = MIXASTAR(tempR, p)
    #############

    ####### NonLinear
    NonLinear(tempB, p)

    return tempB


def SumBlock(input1: np.array, input2: np.array, p):
    for i in range(p):
        input1[i] = input1[i] ^ input2[i]


def RotBlock(tempB: np.array, p):
    tmp = np.zeros(p, dtype=np.longlong)
    tmp[p - 1] = tempB[0]
    for i in range(p - 1):
        tmp[i] = tempB[i + 1]

    return tmp


def AtlantickSPbox(State, p):
    tempL = np.zeros(p, dtype=np.longlong)
    tempR = np.zeros(p, dtype=np.longlong)
    # Sequential computation (Compute in parallel if wanted)

    for i in range(p):
        tempL[i] = State[i]
        tempR[i] = State[p + i]
    

    SumBlock(tempR, RotBlock(tempL, p), p)
    SumBlock(tempL, MIXB(tempR, p), p)
    SumBlock(tempR, MIXB(tempL, p), p)
    SumBlock(tempL, RotBlock(tempR, p), p)

    for i in range(p):
        State[i] = tempL[i]
        State[p + i] = tempR[i]
    
    return State


def InvAtlantickSPbox(State, p):
    tempL = np.zeros(p, dtype=np.longlong)
    tempR = np.zeros(p, dtype=np.longlong)
    # Sequential computation (Compute in parallel if wanted)
    for i in range(p):
        tempL[i] = State[i]
        tempR[i] = State[p + i]

    SumBlock(tempL, RotBlock(tempR, p), p)
    SumBlock(tempR, MIXB(tempL, p), p)
    SumBlock(tempL, MIXB(tempR, p), p)
    SumBlock(tempR, RotBlock(tempL, p), p)

    for i in range(p):
        State[i] = tempL[i]
        State[p + i] = tempR[i]
    
    return State

def AtlantickRoundFunction(M, p, r):
    M = AtlantickLinearLayer(M, p)
    M[0] = M[0] ^ Rc[2 * r]
    M[p] = M[p] ^ p ^ Rc[2 * r + 1]
    M = AtlantickSPbox(M, p)
    return M

def AtlantickPermutation(M, p, s):
    for r in range(s):
        M = AtlantickRoundFunction(M, p, r)        

    return M


def InvAtlantickRoundFunction(M, p, r):
    M = InvAtlantickSPbox(M, p)
    M[p] = M[p] ^ p ^ Rc[2 * r + 1]
    M[0] = M[0] ^ Rc[2 * r]
    M = InvAtlantickLinearLayer(M, p)
    return M

def InvAtlantickPermutation(M, p, s):
    for r in range(s - 1, -1, -1):
        M = InvAtlantickRoundFunction(M, p, r)

    return M


if __name__ == '__main__':
    p = 64; s = 11
    M = [
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
    ]

    S = AtlantickPermutation(M, p, s)
    S = [v % (1 << 32) for v in S]
    print([hex(v) for v in S])
    S = InvAtlantickPermutation(S, p, s)
    print([hex(v) for v in S])

    # S = AtlantickLinearLayer(M, p)
    # print(S)
    # S = InvAtlantickLinearLayer(S, p)
    # print(S)

    # S = AtlantickSPbox(M, p)
    # print(S)
    # S = InvAtlantickSPbox(S, p)

