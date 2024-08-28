function [Y_out,f_out]=fft_plot(data,Fs,x_axis_in,format_in,win)

% function []=fft_izris(y,Fs)
% data = vhodni podatek;
% Fs= vzor?evalna frekvenca
% x_axis is log or lin
% format is abs or dB or log
% def window is kaiser (N,33)
% WINDOW Window function gateway.
%   WINDOW(@WNAME,N) returns an N-point window of type specified
%   by the function handle @WNAME in a column vector.  @WNAME can
%   be any valid window function name, for example:
%
%   @bartlett       - Bartlett window.
%   @barthannwin    - Modified Bartlett-Hanning window. 
%   @blackman       - Blackman window.
%   @blackmanharris - Minimum 4-term Blackman-Harris window.
%   @bohmanwin      - Bohman window.
%   @chebwin        - Chebyshev window.
%   @flattopwin     - Flat Top window.
%   @gausswin       - Gaussian window.
%   @hamming        - Hamming window.
%   @hann           - Hann window.
%   @kaiser         - Kaiser window.
%   @nuttallwin     - Nuttall defined minimum 4-term Blackman-Harris window.
%   @parzenwin      - Parzen (de la Valle-Poussin) window.
%   @rectwin        - Rectangular window.
%   @taylorwin      - Taylor window.
%   @tukeywin       - Tukey window.
%   @triang         - Triangular window.

switch nargin
    case 2
        x_axis='lin';
        format='abs';
        L=length(data);
        w=kaiser(L,33);
    case 3
        format='abs';
        x_axis=x_axis_in;
        L=length(data);
        w=kaiser(L,33);
    case 4
        x_axis=x_axis_in;
        format=format_in;
        L=length(data);
        w=kaiser(L,33);        
    case 5
        x_axis=x_axis_in;
        format=format_in;
        L=length(data);
        w=window(win,L); 
        
        
    otherwise
        disp('Wrong number of input arrguments')
%         disp(nargin)
        return;
end


if iscolumn(data)
    y=data;
else
    y=data';
end


%%% win generation def is kaiser 
% L=length(data);
% w=kaiser(L,33);
% w=hanning(L);
y=y.*w/(sum(w)/(L+1)); %%%rescaling to proper out value
% y=y.*w;
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% NFFT=L;
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
Y_abs=2*abs(Y(1:NFFT/2+1));
norm=sum(Y_abs)/sum(data(1:end).^2)*length(data);	


% Y_abs=  Y_abs/ norm;     % norm power spectrum
Y_out=Y_abs;
f_out=f;
%         keyboard
% Plot single-sided amplitude spectrum.
figure()
if strcmp(x_axis,'lin')
    if strcmp(format,'dB')
        plot(f,db(Y_abs))
    elseif strcmp(format,'log')
        semiloxy(f,(Y_abs));
    else
        plot(f,(Y_abs))
    end
elseif strcmp(x_axis,'log')
    if strcmp(format,'dB')
        semilogx(f,db(Y_abs))
    elseif strcmp(format,'log')
        loglog(f,(Y_abs));
    else
        semilogx(f,(Y_abs))
    end  
else
    semilogx(f,db(Y_abs))
end
%ax=xlim();
%xlim([1 10e6])
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
grid;

