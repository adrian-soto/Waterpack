import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

Afile='trash.dat'
A = np.loadtxt(Afile)
A = A.astype(int)

G=nx.from_numpy_matrix(A)

plt.matshow(A)
#plt.show()

NCC=nx.number_connected_components(G)
print "# of connected components =", NCC
nx.draw(G, with_labels=True)


# Prepare subplots
fig = plt.figure()
axD=fig.add_subplot(221)
axA=fig.add_subplot(222)
axL=fig.add_subplot(223)
axC=fig.add_subplot(224)


edges=np.arange(0,max(nx.degree(G).values())+2)
axD.hist(sorted(nx.degree(G).values()), edges)
axD.set_xlabel('Node degree')
axD.set_ylabel('Frequency')
axD.title.set_text('Node degree distribution of A')


A=nx.to_numpy_matrix(G)
axA.plot(np.sort(np.linalg.eigvals( A ) ) )
axA.set_xlim(0,np.shape(A)[0])
axA.set_xlabel('eigenvector index')
axA.set_ylabel('eigenvalue')
axA.title.set_text('Eigenvalue spectrum of A')


L=nx.normalized_laplacian_matrix(G)
axL.plot(np.sort(np.linalg.eigvals( L.todense() ) ) )
axL.set_xlim(0,np.shape(L)[0])
axL.set_xlabel('eigenvector')
axL.set_ylabel('eigenvalue')
axL.title.set_text('Eigenvalue spectrum of L_n')
    #plt.show()

C=nx.degree_centrality(G)
axC.plot(np.sort(C.values()) )
axC.set_xlim(0,np.shape(L)[0])
axC.set_xlabel('node index')
axC.set_ylabel('centrality')
axC.title.set_text('Degree centrality spectrum')

plt.tight_layout()
#plt.show()



A1 =A.astype(int)
A2 =np.dot(A1 , A1)
A3 =np.dot(A2 , A1)
A4 =np.dot(A3 , A1)
A5 =np.dot(A4 , A1)
A6 =np.dot(A5 , A1)
A7 =np.dot(A6 , A1)
A8 =np.dot(A7 , A1)
A9 =np.dot(A8 , A1)
A10=np.dot(A9 , A1)

Nmol=A.shape
Nmol=Nmol[0]

print Nmol
for i in range(0, Nmol):
    print i, A1[i,i], A2[i,i], A3[i,i], A4[i,i], A5[i,i], A6[i,i], A7[i,i], A8[i,i], A9[i,i], A10[i,i]
    


