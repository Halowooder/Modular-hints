from sage.all import *
from tqdm import tqdm
import numpy


G = GF(7681)
s = vector(G, 768)
for i in range(len(s)):
    r = 0
    for _ in range(2):
        r += randint(0,1) - randint(0,1)
    s[i] = G(r * 3)
s[-1] = s[-1] + 1
M = load("nttru_nttmatrix.sobj")
M = matrix(G, M)
hat_s = M * s
z1, z2 = vector(G, 256), vector(G, 256)
for i in range(256):
    z1[i] = -hat_s[3 * i] / hat_s[3 * i + 2]
    z2[i] = -hat_s[3 * i + 1] / hat_s[3 * i + 2]
B = matrix(G, 512, 768)
for i in range(256):
    B[2 * i, 3 * i] = 1
    B[2 * i, 3 * i + 2] = z1[i]
    B[2 * i + 1, 3 * i + 1] = 1
    B[2 * i + 1, 3 * i + 2] = z2[i]

assert list(B * hat_s).count(0) == 512

Bbar = B * M
N = B.dimensions()[1]
assert list(Bbar * s).count(0) == 512

s = [int(_) for _ in list(s)]

for i in range(len(s)):
    if s[i] > 7680 / 2:
        s[i] -= 7681

s = (vector(s[:-1]) / 3).change_ring(G)

P = identity_matrix(G, 512, 512)
for j in range(512):
    P[j, -1] = Bbar[j, -1]
    P[j, j] = -Bbar[-1, -1]
Bbar = P * Bbar
Bbar = Bbar[:-1, :-1]
assert list(Bbar * s).count(0) == 511



delta = 265

for i in tqdm(range(delta)):
    P = identity_matrix(G, 511, 511)
    for j in range(i, 511):
        P[j, i] = Bbar[j, i]
        P[j, j] = -Bbar[i, i]
    Bbar = P * Bbar
Bbar = Bbar[delta:, delta:]
s = s[delta:]
assert list(Bbar * s).count(0) == 511 - delta

B1 = Bbar[:, : -(N * 2 // 3 - 1) + delta]
B2 = Bbar[:, -(N * 2 // 3 - 1) + delta :]
s1 = list(s[: -(N * 2 // 3 - 1) + delta].change_ring(ZZ))
s2 = list(s[-(N * 2 // 3 - 1) + delta :].change_ring(ZZ))
C = (B2**-1 * B1).change_ring(ZZ)
assert sum([abs(_) for _ in list((C * vector(s1) + vector(s2)) % 7681)]) == 0


L = block_matrix(
    ZZ,
    [
        [identity_matrix(ZZ, (C.T).nrows(), (C.T).nrows()), C.T],
        [
            matrix(ZZ, (C.T).ncols(), (C.T).nrows(), 0),
            (7681 * matrix(ZZ, (C.T).ncols(), (C.T).ncols(), 1)),
        ],
    ],
    subdivide=False,
)

for i in range(len(s1)):
    if s1[i] > (7681 - 1) // 2:
        s1[i] = s1[i] - 7681
    s1[i] = s1[i] * -1

for i in range(len(s2)):
    if s2[i] > (7681 - 1) // 2:
        s2[i] = s2[i] - 7681

deltab = lambda beta: ((pi * beta) ** (1 / beta) * (beta / (2 * pi * exp(1)))) ** (
    1 / (2 * beta - 2)
)
r0 = vector(s1 + s2).norm().n()
gh = (sqrt(L.nrows() / 2 / pi / exp(1)) * (7681) ** ((511 - delta) / (767 - delta))).n()

for bb in range(L.nrows(), 40, -1):
    if sqrt(bb) * numpy.std(vector(s1 + s2)) < deltab(bb) ** (
        2 * bb - L.nrows() - 1
    ) * (7681) ** ((511 - delta) / (767 - delta)):
        beta = bb
print(numpy.std(vector(s1 + s2)))

with open("Nttru_data/Matrix_{}".format(L.dimensions()[0]),"w") as f:
    f.write(str(L))
with open("Nttru_data/shortvector_{}".format(L.dimensions()[0]),"w") as f:
    f.write(str(s1 + s2))


print(delta, beta, r0, gh)