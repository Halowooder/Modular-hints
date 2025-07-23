from sage.all import *
from random import randint
from tqdm import tqdm
import numpy


G = GF(3329)
M = matrix(G, 128, 128, immutable=False)
for i in range(128):
    for j in range(128):
        M[i, j] = G(17) ** ((2 * i + 1) * j)
Zero = matrix(G, 128, 128)
M2 = block_matrix([[M, Zero], [Zero, M]])

P1 = matrix(G, 256, 256)
for i in range(128):
    P1[2 * i, i] = 1

for j in range(128):
    P1[2 * j + 1, 128 + j] = 1

M3 = P1 * M2 / P1

save(M3, "kyber_nttmatrix.sobj")


G = GF(3329)
_s = vector(G, 256)
for i in range(len(_s)):
    r = 0
    for _ in range(2):
        r += randint(0,1) - randint(0,1)
    _s[i] = G(r)

M = load("kyber_nttmatrix.sobj")
M = matrix(G, M)
hat_s = M * _s
z = vector(G, 128)
for i in range(128):
    z[i] = -hat_s[2 * i] / hat_s[2 * i + 1]

B = matrix(G, 128, 256)
for i in range(128):
    B[i, 2 * i] = 1
    B[i, 2 * i + 1] = z[i]
assert list(B * hat_s).count(0) == 128

BB = B * M
N = B.dimensions()[1]
assert list(BB * _s).count(0) == 128

for delta in tqdm(range(1)):
    Bbar = deepcopy(BB)
    for i in range(delta):
        P = identity_matrix(G, 128, 128)
        for j in range(i, 128):
            P[j, i] = Bbar[j, i]
            P[j, j] = -Bbar[i, i]
        Bbar = P * Bbar
    Bbar = Bbar[delta:, delta:]

    s = deepcopy(_s[delta:])
    assert list(Bbar * s).count(0) == 128 - delta

    B1 = Bbar.matrix_from_columns([2 * i for i in range(Bbar.ncols() // 2)])
    B2 = Bbar.matrix_from_columns([2 * i + 1 for i in range(Bbar.ncols() // 2)])
    s1 = [s.change_ring(ZZ)[2 * i] for i in range(Bbar.ncols() // 2)]
    s2 = [s.change_ring(ZZ)[2 * i + 1] for i in range(Bbar.ncols() // 2)]

    # B1 = Bbar[:, :-128 + delta]
    # B2 = Bbar[:, -128 + delta:]
    # s1 = list(s[:-128 + delta].change_ring(ZZ))
    # s2 = list(s[-128 + delta:].change_ring(ZZ))
    # C = (B2 ** -1 * B1).change_ring(ZZ)
    # assert sum([abs(_) for _ in list((C * vector(s1) + vector(s2)) % 3329)]) == 0

    C = (B1**-1 * B2).change_ring(ZZ)
    assert sum([abs(_) for _ in list((C * vector(s2) + vector(s1)) % 3329)]) == 0

    L1 = identity_matrix(ZZ, (C.T).nrows(), (C.T).nrows())
    L2 = matrix(ZZ, (C.T).ncols(), (C.T).nrows(), 0)
    L3 = identity_matrix(ZZ, (C.T).ncols(), (C.T).ncols())
    L = block_matrix(ZZ, [[L1, C.T], [L2, (3329 * L3)]], subdivide=False)

    for i in range(len(s1)):
        if s1[i] > (3329 - 1) // 2:
            s1[i] = s1[i] - 3329
        s1[i] = s1[i] * -1

    for i in range(len(s2)):
        if s2[i] > (3329 - 1) // 2:
            s2[i] = s2[i] - 3329

    deltab = lambda beta: ((pi * beta) ** (1 / beta) * (beta / (2 * pi * exp(1)))) ** (1 / (2 * beta - 2))
    r0 = vector(s1 + s2).norm().n()
    print(r0)
    gh = (sqrt(L.nrows() / 2 / pi / exp(1)) * (3329) ** ((128 - delta) / (256 - delta))).n()

    for bb in range(L.nrows(), 40, -1):
        if sqrt(bb) * numpy.std(vector(s1 + s2)) < deltab(bb) ** (2 * bb - L.nrows() - 1) * (3329) ** ((128 - delta) / (256 - delta)):
            beta = bb

    print("delta:{}, sigma:{}, r0:{},gh:{}, beta:{}".format(delta, numpy.std(vector(s1 + s2)), r0, gh, beta))
    with open("kyber_data/Matrix_{}".format(L.dimensions()[0]), "w") as f:
        f.write(str(L))

    with open("kyber_data/shortvector_{}".format(L.dimensions()[0]), "w") as f:
        f.write(str(s1 + s2))
