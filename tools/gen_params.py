

import numpy as np
from astropy.constants import G,c
from astropy.cosmology import WMAP9 as cosmo
import scipy.stats as st
import scipy.interpolate
import tarfile
import scipy
from astropy.io import fits
import astropy.table as tb
from numpy import linalg as LA
import matplotlib.pyplot as plt
from astropy.io import fits



#Read Cl,l
number_of_spectrum=4

l_txt=open("single_cosmology/test_output_single/shear_cl/ell.txt")
l_list=l_txt.readlines()[1:]
l_list=[float(l) for l in l_list]
l_len=len(l_list)
Cl=np.zeros(shape=(number_of_spectrum,number_of_spectrum,l_len))
for row in range(number_of_spectrum):
    for colume in range(row+1):
        cl_cell_txt=open("single_cosmology/test_output_single/shear_cl/bin_"+str(row+1)+"_"+str(colume+1)+".txt")
        cl_cell_list=cl_cell_txt.readlines()[1:]
        cl_cell_list=[float(cl) for cl in cl_cell_list]
        Cl[row][colume]=cl_cell_list
for row in range(number_of_spectrum):
    for colume in range(row+1,number_of_spectrum):
        Cl[row][colume]=Cl[colume][row]
Cl_rot=Cl
Cl=np.swapaxes(Cl,0,2)
Cl=np.swapaxes(Cl,1,2)

#put in Nl

Nl=np.zeros(shape=(l_len,number_of_spectrum,number_of_spectrum))
sigma_list=[0.25,0.28,0.26,0.27]
n_eff_list=[1.47,1.46,1.5,0.73]
for l in range(l_len):
    for index in range(number_of_spectrum):
        #Nl[l][index][index]=sigma_list[index]**2/(n_eff_list[index]*3600*57**2*np.sqrt(2*l_list[l]+1))
        Nl[l][index][index]=sigma_list[index]**2/(n_eff_list[index]*3600*57**2)
Nl_rot=np.swapaxes(Nl,0,1)
Nl_rot=np.swapaxes(Nl_rot,1,2)



Dl=[[] for i in range(number_of_spectrum)]
R3list=[]
R2list=[]
Rlist = []
for l_index in range(len(l_list)):
    Cl_matrix_para=Cl[l_index]
    Nl_matrix_para=Nl[l_index]
    

    #eigenvector_Nl=LA.eig(Nl_matrix_para)[1].transpose()
    R2=np.zeros(shape=(number_of_spectrum,number_of_spectrum))
    for i in range(number_of_spectrum):
        R2[i][i]=1/np.sqrt(Nl_matrix_para[i][i])
    before3=np.matmul(np.matmul(R2,Cl_matrix_para),R2)
    
    eigenValues, eigenVectors = LA.eig(before3)
    idx = eigenValues.argsort()[::-1]   
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    
    R3=eigenVectors.transpose()
    #R3=LA.eig(before3)[1].transpose()
    for i in range(number_of_spectrum):
        if abs(np.max(R3[i]))<abs(np.min(R3[i])):
            R3[i]=-R3[i]
    R3list.append(R3)
    R2list.append(R2)
    Rlist.append(np.matmul(R3,R2))
    #R.append(R3.transpose())
    #R.append(np.matmul(R2,R3).transpose())
    Cy=np.matmul(np.matmul(R3,before3),np.linalg.inv(R3))
    #print Cy
    Diag_list_Cy=np.diag(Cy)  

    Mode_list=np.sort(Diag_list_Cy)
    for mode_rank in range(len(Diag_list_Cy)):
        #Dl[mode_rank].append(Diag_list_Cy[mode_rank])
        Dl[mode_rank].append(Mode_list[mode_rank])


low_index=161
high_index=280
R_part=Rlist[low_index:high_index]
ep_bar=[]
epl=np.swapaxes(R_part,0,1)
epl=np.swapaxes(epl,1,2)

ep_bar=np.average(epl,axis=2,weights=l_list[low_index:high_index])
print np.shape(ep_bar)
print l_list[low_index],l_list[high_index]


def generate_vector_transformation(matrix_transformation):
    dimension = len(matrix_transformation)
    ans = np.zeros(shape=(dimension**2))
    off_dia_list=[1,2,3,4,6,7,8,9,11,12,13,14]
    for i in range(len(ans)):
        ans[i] = matrix_transformation[i//dimension]*matrix_transformation[i%dimension]*(1+float(i in off_dia_list))

    indexlist=[0,4,7,9]
    for i in range(dimension):
        for j in range(i):
            index = indexlist[i]
            ans = np.delete(ans,index)
    ans = ans/np.linalg.norm(ans)
    print(ans)
    return ans
#print Cl[0]
vector_Cl0 = np.zeros(shape=(10,))
vector_Cl0[0:4] = Cl[0][0][0:]
vector_Cl0[4:7] = Cl[0][1][1:]
vector_Cl0[7:9] = Cl[0][2][2:]
vector_Cl0[9:10] = Cl[0][3][3]

new_weight=generate_vector_transformation(np.array(ep_bar[0]))
print new_weight

Cl0_SN = np.dot(new_weight,vector_Cl0)
norm_weight = new_weight/np.linalg.norm(new_weight)

transformation_matrix = np.zeros(shape=(40,400))
repeat_weight = np.repeat(norm_weight,20)
for i in range(transformation_matrix.shape[0]):
    angle_index = i%20
    if i<20:
        for j in range(transformation_matrix.shape[1]/2):
            bin_index= (j%200)/20
            if j%20==angle_index:
                transformation_matrix[i][j]=norm_weight[bin_index]
    else:
        for j in range(transformation_matrix.shape[1]/2,transformation_matrix.shape[1]):
            bin_index= (j%200)/20
            if j%20==angle_index:
                transformation_matrix[i][j]=norm_weight[bin_index]

#plt.show()

xip_angle_min = 5
xip_angle_max = 20
xim_angle_min = 16
xim_angle_max = 20

remove_row_list = np.array(range(0,xip_angle_min)+range(xip_angle_max,20)+range(20,20+xim_angle_min)+range(xim_angle_max+20,40))
print remove_row_list


transformation_matrix=np.delete(transformation_matrix,remove_row_list,axis=0)

print transformation_matrix.shape,"here" 


#save weight

np.savetxt('Parameters'+'/weight.txt',transformation_matrix)

#calculate data vecter

data_txt = open('single_cosmology/test_output_single/data_vector/2pt_data.txt')
data = data_txt.readlines()[1:]
data = np.array([float(d) for d in data])
print data.shape

plt.plot(np.array(range(0,len(data))),data)
plt.show()

compressed_data = transformation_matrix.dot(data)
print(compressed_data.shape)

np.savetxt('Parameters'+'/data_c.txt',compressed_data)


def read_dat_covariance(filename):
    cov_txt = open(filename)
    data = cov_txt.readlines()
    last_line = data[0].split("\t")[:-1]

    shape = (len(data),len(last_line))
    cov_array = np.zeros(shape)
    for row in range(len(data)):
        line_list = data[row].split("\t")[:-1]
        for column in range(len(line_list)):
            cov_array[row][column]=float(line_list[column])
    header = fits.Header()

    header["COVDATA"]=True
    header["EXTNAME"]="COVMAT"
    cov_hdu = fits.ImageHDU(data=cov_array,header=header)

    return cov_hdu

def read_txt_covariance(filename):
    cov_txt = open(filename)
    data = cov_txt.readlines()
    last_line = data[-1].split(" ")
    last_line[0] = int(last_line[0])
    last_line[1] = int(last_line[1])
    shape = (last_line[0]+1,last_line[1]+1)
    cov_array = np.zeros(shape)
    print cov_array[0][0]
    for line in data:
        line_list = line.split(" ")
        row = int(line_list[0])
        column = int(line_list[1])
        correlation = float(line_list[2])
        cov_array[row][column] = correlation
    header = fits.Header()

    header["COVDATA"]=True
    header["EXTNAME"]="COVMAT"
    header["STRT_0"]=0
    header["NAME_0"]="xip"
    header["STRT_1"]=200
    header["NAME_1"]="xim"
    header["STRT_2"]=400
    header["NAME_2"]="gammat"
    header["STRT_3"]=800
    header["NAME_3"]="wtheta"
    cov_hdu = fits.ImageHDU(data=cov_array,header=header)

    return cov_hdu


# txt_filename = "cov_y3_mcal_emu_final.txt"
# cov_hdu2 = read_txt_covariance(txt_filename)

fits_filename = "2pt_NG_mcal_1110.fits"
cov_hdu2 = fits.open(fits_filename)[1]
original_cosmolike = cov_hdu2.data[0:400,0:400]
compressed_cosmolike = np.matmul(transformation_matrix,np.matmul(original_cosmolike,transformation_matrix.T))

np.linalg.inv(compressed_cosmolike)

np.savetxt('Parameters'+'/cov_c.txt',compressed_cosmolike)





