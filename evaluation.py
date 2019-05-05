import numpy as np
def rmse(true,pred):
    if len(true) != len(pred):
        print("true and pred do not have same length")
        exit
    sum = 0
    for i in range(len(true)):
        sum = np.sum((true[i]-pred[i])*(true[i]-pred[i]))
    return np.sqrt(sum/len(true))

def average(pred):
    sum = np.zeros(pred[0].shape)
    for i in range(len(pred)):
        sum = sum+pred[i]
    return sum/len(pred)
def pixel_to_mm(data,factor):
    data = np.concatenate((data[:,0:1]*0.4,data[:,1:2]*factor,data[:,2:3]*0.4,data[:,3:4]*factor,data[:,4:5]*0.4,data[:,5:6]*factor,data[:,6:7]*0.4,data[:,7:8]*factor),axis=1)
    return data
def dx_theta_phi(data):
    dx = np.sqrt((data[:,0:1]-data[:,2:3])*(data[:,0:1]-data[:,2:3])+
                 (data[:,4:5]-data[:,6:7])*(data[:,4:5]-data[:,6:7])+
                 (data[:,1:2]-data[:,3:4])*(data[:,1:2]-data[:,3:4]))
    theta = (data[:,0:1]-data[:,2:3])/dx
    theta = np.degrees(np.arccos(theta))
    phi = (data[:,4:5]-data[:,6:7])/np.maximum((data[:,1:2]-data[:,3:4]),0.000001)
    phi = np.degrees(np.arctan(phi))
    return dx,theta,phi
