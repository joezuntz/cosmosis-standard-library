import numpy as np
from scipy.signal import butter, lfilter, filtfilt, lfilter_zi
import copy 
from P_extend import k_extend 

def compare(A,B):
	return np.sum(A-B)


def filter_highk(k,P_in,start,end):
	P_out=copy.deepcopy(P_in)
	# start = k where you want to start cross fad
	# end = k where you want to end the cross fad 

	# filtering specific 
	k_start=start; k_end=end
	id1=np.where( k > k_end)[0]
	id2=np.where( k <= k_end)[0]
	id3=np.where( k > k_start)[0]
	#id4=np.where( k <= k_start)[0]
	id4=np.where( (k > k_start) & ( k<= k_end))[0]
	
	order=6; wn=.1
	B,A=butter(order,wn, btype='low', analog=False)
	
	theta=np.linspace(1,0,id4.size)
	W_fad=theta - 1/2./np.pi*np.sin(2*np.pi*theta) 
	filt_pad=id3.size
	# end filtering specific 
	
	def zero_phase(sig):
		
		sig=np.pad(sig,(filt_pad,filt_pad), 'constant', constant_values=(0, 0))
		#zi=lfilter_zi(B,A)
		#x,_=lfilter(B,A,sig, zi=zi*sig[0])
		x=lfilter(B,A,sig)
		#y,_=lfilter(B,A,x,zi=zi*x[0])
		y=lfilter(B,A,x[::-1])
		y=y[::-1]
		#return y
		return y[filt_pad:id3.size+filt_pad]
	
		
	P_smoothed=zero_phase(P_out[id3])
	P_patch=P_out[id4]*W_fad
	P_out[id3]=P_smoothed
	P_out[id4]=P_patch+(1-W_fad)*P_out[id4]
	
	return P_out

def filter_lowk(k,P_in,start,end):

	P_out=copy.deepcopy(P_in)
	# start = k where you want to start cross fad
	# end = k where you want to end the cross fad 

	# filtering specific 
	k_start=start; k_end=end
	id1=np.where( k > k_end)[0]
	id2=np.where( k <= k_end)[0]
	id3=np.where( k  < end)[0]
	id4=np.where( (k > k_start) & ( k<= k_end))[0]
	
	order=6; wn=.1
	B,A=butter(order,wn, btype='low', analog=False)
	
	theta=np.linspace(1,0,id4.size)
	theta=theta[::-1]
	W_fad=theta - 1/2./np.pi*np.sin(2*np.pi*theta) 
	filt_pad=id3.size
	# end filtering specific 
		
	def zero_phase(sig):
		
		sig=np.pad(sig,(filt_pad,filt_pad), 'constant', constant_values=(0, 0))
		#zi=lfilter_zi(B,A)
		#x,_=lfilter(B,A,sig, zi=zi*sig[0])
		x=lfilter(B,A,sig)
		#y,_=lfilter(B,A,x,zi=zi*x[0])
		y=lfilter(B,A,x[::-1])
		y=y[::-1]
		#return y
		return y[filt_pad:id3.size+filt_pad]
	
		
	P_smoothed=zero_phase(P_out[id3])
	P_patch=P_out[id4]*W_fad
	P_out[id3]=P_smoothed
	P_out[id4]=P_patch+(1-W_fad)*P_out[id4]
	
	return P_out

def BW_filter(P_in,order=3,nf=.01):
		
	print 'at butter, freq=', nf	
		
	B, A = butter(order, nf, 'low')
	
	sig_ff = filtfilt(B, A, P_in, padlen=200)

	return sig_ff

if __name__=="__main__":

	d=np.loadtxt('Pk_Planck15.dat')
	k=d[:,0]; P0=d[:,1]

	import copy
	test=copy.deepcopy(P0)
	
	low_extrap=-4
	high_extrap=5
	EK=k_extend(k,low_extrap,high_extrap)
	k=EK.extrap_k()
	
	
	P0=EK.extrap_P_low(P0)	
	P0=EK.extrap_P_high(P0)
	
	
	P1=filter_highk(k,P0,1,5)
	P2=filter_lowk(k,P0,.01,.05)
	
	k,P1=EK.PK_orginal(P1)
	k,P2=EK.PK_orginal(P2)
	k,P0=EK.PK_orginal(P0)

	import matplotlib.pyplot as plt

	ax=plt.subplot(141)
	ax.set_xscale('log')
	#ax.set_ylim(.99,1.01)
	ax.set_yscale('log')
	#P3=BW_filter(P0)
	ax.plot(k[:-2],np.absolute(np.diff(P0,2)), label='orginal')
	ax.plot(k[:-2],np.absolute(np.diff(P1,2)), '--', label='high filtered')
	

	plt.grid()
	plt.legend()

	ax=plt.subplot(142)
	ax.set_xscale('log')
	ax.set_yscale('log')
	P1=filter_highk(k,P0,1,5)
	P2=filter_lowk(k,P0,.01,.05)
	
	ax.plot(k[:-2],np.absolute(np.diff(P0,2)), label='orginal')
	ax.plot(k[:-2],np.absolute(np.diff(P2,2)), '--', label='low filtered')
	

	plt.grid()
	plt.legend()

	ax=plt.subplot(143)
	ax.set_ylim(.99,1.01)
	ax.set_xscale('log')
	
	P1=filter_highk(k,P0,1,5)
	P2=filter_lowk(k,P0,.01,.05)
	
	ax.plot(k,P1/P0,  label='high filtered')
	ax.plot(k,P2/P0, '--', label='low filtered')
	

	plt.grid()
	plt.legend()
	
	ax=plt.subplot(144)
	ax.set_yscale('log')
	ax.set_xscale('log')
	
	P1=filter_highk(k,P0,1,5)
	P2=filter_lowk(k,P0,.01,.05)
	
	ax.plot(k,P0,  label='original')
	ax.plot(k,P1,  label='high filtered')
	ax.plot(k,P2, '--', label='low filtered')
	

	plt.grid()
	plt.legend()

	plt.show()