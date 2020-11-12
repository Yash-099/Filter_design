%%%%%%%%%%%%%%%%%%%%%%%%
% FIR Band pass filter %
%%%%%%%%%%%%%%%%%%%%%%%%

close all

%% Determining the sequence hFIR[n]
n = -25:1:25;
w1 = 63.4*2/330;
w2 = 91.4*2/330;
h = (sinc(n.*w2).*w2)-(w1.*sinc(n.*w1));
figure(1);
stem(h);
title('FIR sequence of bandpass filter');
grid;

%% Z-transform of the FIR 
z = tf('z');
Hzs = filt(h,1);
Hz = Hzs*z^(25);
[num,den] = tfdata(Hzs,'v');

%% Plotting the frequency response
[H,w] = freqz(num,den,1000);
figure(2);
plot(w*(330/(2*pi)),abs(H));
grid;
title('magnitude response of FIR bandpass filter')
fvtool(num,den);

% hold on;
% ws1 = 56.6;
% wp1 = 60.6;
% wp2 = 80.6;
% ws2 = 84.6;
% plot(w*(330/(2*pi)),xline(ws1));
% plot(w*(330/(2*pi)),xline(wp1));
% plot(w*(330/(2*pi)),xline(wp2));
% plot(w*(330/(2*pi)),xline(ws2));
% y1 = 0.15;
% y2 = 0.85;
% y3 = 1.15;
% plot(w*(330/(2*pi)),yline(y1));
% plot(w*(330/(2*pi)),yline(y2));
% plot(w*(330/(2*pi)),yline(y3));
% ylim([0,1.5]);
