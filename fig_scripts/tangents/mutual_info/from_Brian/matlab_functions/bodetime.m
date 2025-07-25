function [H f] = bodetime(h, flag);
%function [H f] = bodetime(h, flag);
%Bode plots for time domain h
%plot if flag=1

tmp = fftshift(fft(h));
H = tmp(round(length(tmp)/2+1):end);
%f1 = 1/dt/2*linspace(0,1,length(H));
f = logspace(-4,4,length(H));

Filt = H;

if flag==1
    figure(1)
    subplot(211)
    loglog(f,abs(Filt).^2),title('Magnitude')
    subplot(212)
    semilogx(f,angle(Filt)*180/pi),title('Phase (deg)')
    xlabel('Frequency (Hz)')
end