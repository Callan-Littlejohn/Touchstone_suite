# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 13:57:51 2024

@author: FTICR_Kool_Kidz_PC
"""


import numpy as np
import array
import matplotlib.pyplot as plot

sampling_rate=2**20
duration=0.1
no_zerofills=1

class urqrd(): #  class structure makes life a bit easier for long term data managment
    def __init__(self,data,rank,iterations): 
        self.data=data # makes things easier for data management, original code makes 3 copies 
        order=int(self.data.size/2) # order is availlable to change but for our purposes this works
        N=len(data)-order+1 # only useful in cases where order doesnt equal hald data size
        for i in range(iterations): # genrerally their benchmarking uses lots of iterations
            random_numbers = np.random.normal(size=(N,rank)) # this makes an array of random numbers with rank number of rows of N size
            first_order_q, data_multiplied_by_random_q=self.urqrd_core(random_numbers) #gets the first order q, a component of data, and the data convolved with this component
            self.data=self.FastHankel_2dt(first_order_q,data_multiplied_by_random_q)# main getback function
        if self.data.dtype =="float": #Kludge in spike that has been left in, this might not be needed
            self.data=np.real(self.data) # takes the real component
            
    
    def urqrd_core(self,random_numbers): # urqrd base function
        data_multiplied_by_random_numbers =  self.FastHankel_prod_mat_mat(self.data,random_numbers) # convolves data with random numbers
        first_order_q,r=np.linalg.qr(data_multiplied_by_random_numbers) # standard numpy function
        del(r) # sems a bit pointless
        data_multiplied_by_random_q=self.FastHankel_prod_mat_mat(self.data.conj(),first_order_q).conj().T # convolve data with componernt of itself
        return first_order_q, data_multiplied_by_random_q # return component of data and data multiplied by component of itself


    def FastHankel_prod_mat_mat(self,data,random_numbers):
        N,rank=random_numbers.shape # Rakes rank from first call and N is lendata- order
        length_of_data=len(data) # in essemce
        N_star=length_of_data-N+1# seems unnecesary 
        d=np.zeros(shape=(N_star,rank),dtype=complex) # zeros 2d array. seems unnecessary
        for r in range(rank): #sweep through rows
            column_of_random_nums=random_numbers[:,r] # take columns out
            d[:,r]=self.FastHankel_prod_mat_vec(column_of_random_nums,data) #input columns into data
        return d
            
    def FastHankel_prod_mat_vec(self,column_of_random_numbers,data): # convolving matrix with vector
        column_of_random_numbers_zerofilled=np.concatenate((np.zeros((len(data)-len(column_of_random_numbers))),column_of_random_numbers[::-1])) # zerofill random numbers
        data_fft,random_fft=np.fft.fft(data),np.fft.fft(column_of_random_numbers_zerofilled) # bring both datasets into frequency space
        product_of_two_ffts=data_fft*random_fft#multioly both frequency space vectors together
        resultant_signal=np.fft.ifft(product_of_two_ffts) #return combined dataset into time domain
        return np.roll(resultant_signal,+1)[:(len(data)-len(column_of_random_numbers)+1)]# roll data then retunrn up to halfway point
            
    def FastHankel_2dt(self,first_order_q,data_multiplied_by_random_q):
        #print(type(first_order_q)) #debug statement on my part
        data_len_half,rank=first_order_q.shape # seems lazy
        rank,data_len_half_plus1=data_multiplied_by_random_q.shape # nearly pointless double variable definition
        sum_of_all_data=np.zeros((len(self.data)),dtype=complex) # yet another 'lets make a large dataset full of zeroes
        for r in range(rank): # essentially flattens the data
            row=data_multiplied_by_random_q[r,:] # pointless but whatever
            zerofilled_first_order_q=np.concatenate((np.zeros(data_len_half_plus1-1),first_order_q[:,r],np.zeros(data_len_half_plus1-1))) # zerofill on both sides
            zerofilled_first_order_q_multi_data=self.FastHankel_prod_mat_vec(row[::-1],zerofilled_first_order_q) # convolvedata with dataconvolved with fragment of itself
            sum_of_all_data+=zerofilled_first_order_q_multi_data  #flattens array into 1D of row lengtj
        return (sum_of_all_data*self.vec_mean(data_len_half,(len(self.data)))) # returns that multiplied by the diagonal
        
    def vec_mean(self,M,L): # seems to find the diagonal
        '''
        Vector for calculating the mean from the sum on the antidiagonal.
        data = vec_sum*vec_mean
        '''
        vec_prod_diag = [1/float((i+1)) for i in range(M)]
        vec_prod_middle = [1/float(M) for i in range(L-2*M)]
        vec_mean_prod_tot = vec_prod_diag + vec_prod_middle + vec_prod_diag[::-1]
        return np.array(vec_mean_prod_tot)

class MAGIC():
    def __init__(self,data,rank):
        data2=[]
        for i in range(rank):
            #d=data**2
            data2.append(np.fft.irfft(np.fft.rfft(np.roll(data,+i))*np.sqrt(np.fft.rfft(np.roll(data,-i)))))
        data3=[0*len(data)]
        for i in data2:
            self.data3=data3+i

class MAGICnoroll():
    def __init__(self,data,rank):
        data2=[]
        for i in range(rank):
            #d=data**2
            data2.append(np.fft.irfft(np.fft.rfft(data)**2))
        data3=[0*len(data)]
        for i in data2:
            self.data3=data3+i        

class MAGICnosquare():
    def __init__(self,data,rank):
        data2=[]
        for i in range(rank):
            #d=data**2
            data2.append(np.fft.irfft(np.fft.rfft(np.roll(data,+i))))
        data3=[0*len(data)]
        for i in data2:
            self.data3=data3+i         
    
         
def import_data(filename,n_scans,s_size): # taken directly from spike
    data=[]
    with open(filename,"rb") as f:
        for i in range(n_scans):
            indiv_scan=f.read(4*s_size)
            indiv_scan=array.array("l",indiv_scan)
            data.append(indiv_scan)
    return data
def gen_signal(freq,samplrate,duration, amplitude):
    x=np.linspace(0.0,duration,int(duration*samplrate),endpoint=False)
    freqs=x*freq
    y=amplitude*np.sin((2*np.pi)*freqs)
    return x,y
# a,b=570,630
# #data=import_data("H:/solarix 12T/20220107/20220107_HPmix_sod_000005.d/fid",1,2**23)[0]
# #data=np.array(data)

# a3,khz5=gen_signal(5000,sampling_rate,duration,1)
# khz5p1=khz5+1
# a2,khz4=gen_signal(4000,sampling_rate,duration,1)
# a1,khz3=gen_signal(3000,sampling_rate,duration,1)
# data=khz5+khz3+khz4
# for i in range(no_zerofills):
#     data=np.pad(data,(0,len(data)),"constant")
# #MADIFC=MADIC(data,1)
# #data=MADIFC.data3
# urqrd_time=urqrd(data,25,1)
# data=urqrd_time.data

# print(len(data))
# plot.plot(abs(np.fft.rfft(data)[a:b]),linewidth=0.5,color="k")
# plot.title("URQRD denoising")

