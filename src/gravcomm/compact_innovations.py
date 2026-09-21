"""Bounded-order causal Gaussian innovations, exact during finite startup."""
import numpy as np
from scipy.linalg import solve_triangular
class CompactInnovations:
 def __init__(self,receiver,amplitude,order=64):
  if order<1 or int(order)!=order:raise ValueError('Positive integer order required')
  self.order=order;self.receiver=receiver;self.amplitude=amplitude
  finite=receiver.finite_record(order+1,amplitude)
  self.startup=solve_triangular(finite.factor,np.eye(order+1),lower=True)
  self.taps=self.startup[-1,::-1].copy()
 def block_matrices(self,start,count):
  """Map the most recent min(start,order) samples and the current block."""
  previous=min(start,self.order);matrix=np.zeros((count,previous+count),complex)
  for j in range(count):
   index=start+j
   taps=self.startup[index,:index+1][::-1] if index<=self.order else self.taps
   for lag,value in enumerate(taps):
    column=previous+j-lag
    if column>=0:matrix[j,column]=value
  return matrix[:,:previous],matrix[:,previous:]
 def spectral_error(self,nfft=65536):
  f=np.fft.fftfreq(nfft,1/self.receiver.output_rate_hz)
  whitened=abs(np.fft.fft(self.taps,nfft))**2*self.receiver.output_psd(f,self.amplitude)
  return dict(max_psd_error=float(max(abs(whitened-1))),rms_psd_error=float(np.sqrt(np.mean((whitened-1)**2))))
