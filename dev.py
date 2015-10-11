import GPd
import scipy as sp

k = GPd.gen_mat32_k_d([1.1,0.2,0.3])

print k([0.1,0.2],[0.12,0.25],[[0],[1]],[[1]])