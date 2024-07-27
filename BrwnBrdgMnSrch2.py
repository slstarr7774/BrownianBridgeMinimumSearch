
import numpy as np
import matplotlib.pyplot as plt
import imageio

def BBFI(n,x):
    x_pr=np.zeros(2*n+1)
    for k_ctr in range(n+1):
        x_pr[2*k_ctr]=x[k_ctr]
    for k_ctr in range(n):
        x_pr[2*k_ctr+1] = x_pr[2*k_ctr]/2 + x_pr[2*k_ctr+2]/2 + 1/np.sqrt(2*n)*np.random.standard_normal()
    return(x_pr)

def Init(d):
    x = np.array([0,0])
    t = np.array(range(2**d+1))/2**d
    for r in range(d):
        x_pr = BBFI(2**r,x)
        x = x_pr
    return(np.array([t,x]))

def ScndCrtfct(d,K,B):
    mn_val = B[K]
    m_lst = np.ones(2**d)*mn_val
    for j_ctr in range(2**d):
        k_ctr=j_ctr+1
        if k_ctr<2**(d-2)+1 or k_ctr>2**(d-2)+2**(d-1):
            x = B[k_ctr-1]
            y = B[k_ctr]
            m = (x+y)/2-np.sqrt((x-y)**2/4-np.log(np.random.uniform())/2**(d+1))
            m_lst[k_ctr-1] = m
    return(m_lst)

N=4
fig, axs = plt.subplots(4,figsize=(12, 12))
d=7
output = Init(d)
t = output[0]
x = output[1]
K = np.argmin(x)
t_star = np.array([t[K]])
x_pr = np.concatenate((x[K:2**d+1:1],x[1:K+1:1]),axis=None)
B_tilde = np.concatenate((x_pr[2**(d-1):2**d+1:1],x_pr[1:2**(d-1)+1:1]),axis=None)-x[K]
B_hat = np.sqrt(2)*B_tilde[2**(d-2):2**(d-2)+2**(d-1)+1:1]
m_lst = ScndCrtfct(d,2**(d-1),B_tilde)
axs[0].plot(t,np.sqrt(2)*B_tilde,color='black')
for k in range(2**d):
    axs[0].plot(np.array([t[k],t[k+1]]),np.sqrt(2)*np.array([m_lst[k],m_lst[k]]),color='grey')
for n_ctr in range(N-1):
    B = BBFI(2**(d-1),B_hat)
    K_min = np.argmin(B)
    if K_min<2**(d-1)-2**(d-3) or K_min>2**(d-1)+2**(d-3):
        print(quit)
        quit()
    ScndCrtfct(d,K_min,B)
    t_star_pr = np.concatenate((t_star,np.array(t[K_min])),axis=None)
    t_star=t_star_pr
    B_hat = np.sqrt(2)*(B[K_min-2**(d-2):K_min+2**(d-2)+1:1]-B[K_min])
    m_lst = ScndCrtfct(d,K_min,B)
    axs[n_ctr+1].plot(t,np.sqrt(2)*(B-B[K_min]),color='black')
    for k in range(2**d):
        axs[n_ctr+1].plot(np.array([t[k],t[k+1]]),np.sqrt(2)*np.array([m_lst[k],m_lst[k]]),color='grey')
