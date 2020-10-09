import numpy as np
import time
from datetime import date
import matplotlib.pyplot as plt


from numpy import *
import pickle
from sklearn.cluster import KMeans
from sklearn import neighbors
from sklearn.tree import DecisionTreeRegressor
import matplotlib
from getCode import *

#stop
matplotlib.rcParams.update({'font.size': 14})
zKu1L=[]
zKa1L=[]
zsfcL=[]
recL=pickle.load(open('belowZc.pklz','rb'))
recL[:,0]=(recL[:,0]-28)/10.
a=np.nonzero(recL[:,-1]<200)
recL=recL[a[0],:]
#recL[:,2:3]=0
#recL[:,5]=0
n_neighbors=30
knn = neighbors.KNeighborsRegressor(n_neighbors, weights='distance')
#knn= DecisionTreeRegressor(max_depth=3)
from sklearn.ensemble import RandomForestRegressor
#sklearn.linear_model.Ridge
from sklearn.linear_model import Ridge

knn =RandomForestRegressor(n_estimators=10,random_state=0)
#knn =DecisionTreeRegressor()#n_estimators=15,random_state=0)
clf = Ridge(alpha=1.0)

nx=recL.shape[0]
for i in range(10):
    r=np.random.random(nx)
    a=np.nonzero(r<0.75)

    x_train=recL[a[0],0:5]
    y_train=(recL[a[0],5]+1e-9)

    b=np.nonzero(r>0.75)

    x_valid=recL[b[0],0:5]
    y_valid=(recL[b[0],5]+1e-9)

    knn.fit(x_train,y_train)
    clf.fit(x_train,log(y_train))
    yp=knn.predict(x_valid)
    ypr=clf.predict(x_valid)
    print(np.corrcoef((yp),(y_valid)))
    #print(np.corrcoef(exp(ypr),(y_valid)))

fout=open("estimators.c","w")
for i in range(10):
    s=get_code("getSfcRain%2.2i"%i,knn.estimators_[i],'x',x_train.shape[1])
    fout.write(s)
s=get_codeEns("getSfcRain","x",10)
fout.write(s)
fout.close()
