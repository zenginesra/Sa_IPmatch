function [filtered_signal]=bandpass_filter(signal,dt,fcutlow,fcuthigh,nOrder)

% This function applies a band-pass filter to a signal,which contains a
% low-pass filter and a high-pass filter in the frequency domain.

% Input Variables
% signal    : input time series
% dt        : time step of the signal
% fcutlow   : lower cutoff frequency (Hz)
% fcuthigh  : upper cutoff frequency (Hz)
% nOrder    : order of the filter

% Output Variables
% filtered_signal : bandpass-filtered time series

signal=reshape(signal,1,[]);
signal=signal(all(~isnan(signal),2),:);

for m=1:length(fcuthigh)
N = 2^nextpow2(length(signal));
Y = fft(signal,N);
df=1/(N*dt);

% lowpass filter
freq=0:df:N/2*df-df;
ff=freq./fcuthigh(m);

for i=2:length(ff)
    % butterworth 
    H(i)= 1./sqrt(((1+ff(i).^(2*nOrder))));
    Y(i)= Y(i).*H(i);
    Y(length(Y)-i+2)=Y(length(Y)-i+2).*H(i);
end

% highpass filter
freq=0:df:N/2*df-df;
ff=freq./fcutlow(m);

H=zeros(1,length(N));
for i=2:length(ff)
    
    H(i)= sqrt((ff(i).^(2*nOrder)))./sqrt(((1+ff(i).^(2*nOrder))));
    Y(i)= Y(i).*H(i);
    Y(length(Y)-i+2)=Y(length(Y)-i+2).*H(i);
end

Y1= (real(ifft(Y,N)))';
Y1=Y1(1:length(signal));
filtered_signal(m,:)=Y1;
end
