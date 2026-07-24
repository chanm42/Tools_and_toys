import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.random as npr
import math
from numpy.linalg import eig
# Following imports so I can visualize on my linux machine
import os
os.environ["QT_LOGGING_RULES"] = "qt.qpa.wayland.textinput=false"

def symmetrize_matrix(M):
    MT = np.transpose(M)
    SM = 0.5 * (M + MT)
    return SM

def normalize(v):
    l = len(v)

    v_sum = 0.0
    for i in range(0,l):
        v_i    = abs(v[i])
        v_sum += v_i * v_i

    v_n = v / np.sqrt(v_sum)
    return v_n

def power_iterate_v1(amat,guess_vec,num_iter,exact_vec):
    #normalize the guess vector
    guess_vec_norm = normalize(guess_vec)

    v_k      = None
    prev_vec = guess_vec_norm
    for k in range(num_iter):
        v_k = np.dot(amat,prev_vec)

        #error out on zero vector
        v_k_norm = normalize(v_k)
        print("Guess #",k,"\t",v_k_norm)

        assert math.sqrt(np.dot(v_k_norm,v_k_norm)) != 0.0, "Power iteration returned the zero vector"

        plt.quiver(0,0,v_k_norm[0],v_k_norm[1],label='v_k',angles='xy',scale_units='xy',scale=1,color='red')
        plt.quiver(0,0,exact_vec[0],exact_vec[1],label='exact',angles='xy',scale_units='xy',scale=1,color='blue')
        plt.xlim(-1.2,1.2)
        plt.ylim(-1.2,1.2)
        plt.legend()
        plt.show()
        prev_vec = v_k_norm

    return v_k_norm

def test1A():
    mat_size  = 2
    num_iter  = 5 
    amat      = np.eye(mat_size,mat_size)
    guess_vec = np.zeros(mat_size)
    for i in range(0,mat_size):
        guess_vec[i] = npr.rand()
        for j in range(0,mat_size):
            amat[i,j] = npr.rand()

    amat        = symmetrize_matrix(amat)
    evals,evecs = eig(amat)
    idx         = np.argmax(evals)
    exact_vec   = evecs[:,idx]
    approx_vec  = power_iterate_v1(amat,guess_vec,num_iter,exact_vec)
    print("Exact eigenvector: ",exact_vec,"\n")
    print("Iterated eigenvector: ",approx_vec)

test1A()
