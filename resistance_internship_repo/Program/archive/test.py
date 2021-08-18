import numpy as np

N1 = np.array([[-1,0,0,0],[1,-1,-1,0],[0,1,0,-1],[0,0,1,1]])
V1 = np.array([2,1,1,1])

print(N1)
print(V1)
print(N1 * V1)
print(np.sum(N1 * V1))


N2 = np.array([[-1,0,0],[1,-1,-1],[0,1,-1],[0,0,1]])
V2 = np.array([2,1,1])

print("\n")
print(N2)
print(V2)
print(N2 * V2)
print(np.sum(N2 * V2))


outf = np.array(N2*V2)
outf[outf >0] = 0
inf = np.array(N2*V2)
inf[inf <0] = 0

print(outf)
print(inf)
