%Butterworth Analog LPF parameters
Wc = 1.05;              %cut-off frequency
N = 10;                  %order 

%poles of Butterworth polynomial of degree 8 in the open CLHP 
%{
p1 = Wc*cos(pi/2 + pi/20) + i*Wc*sin(pi/2 + pi/20);
p2 = Wc*cos(pi/2 + pi/20) - i*Wc*sin(pi/2 + pi/20);
p3 = Wc*cos(pi/2 + pi/20+ pi/10) + i*Wc*sin(pi/2 + pi/20+pi/10);
p4 = Wc*cos(pi/2 + pi/20+ pi/10) - i*Wc*sin(pi/2 + pi/20+pi/10);
p5 = Wc*cos(pi/2 + pi/20+ 2*pi/10) + i*Wc*sin(pi/2 + pi/20+2*pi/10);
p6 = Wc*cos(pi/2 + pi/20+ 2*pi/10) - i*Wc*sin(pi/2 + pi/20+2*pi/10);
p7 = Wc*cos(pi/2 + pi/20+ 3*pi/10) + i*Wc*sin(pi/2 + pi/20+3*pi/10);
p8 = Wc*cos(pi/2 + pi/20+ 3*pi/10) - i*Wc*sin(pi/2 + pi/20+3*pi/10);
p9 = Wc*cos(pi/2 + pi/20+ 4*pi/10) + i*Wc*sin(pi/2 + pi/20+4*pi/10);
p10 = Wc*cos(pi/2 + pi/20+ 4*pi/10) - i*Wc*sin(pi/2 + pi/20+4*pi/10);
%}
p1 = -0.16426 + i*1.0371;
p2 = -0.16426 - i*1.0371;
p3 = -0.47669+ i*0.93556;
p4 = -0.47669- i*0.93556;
p5 = -0.74246+ i*0.74246;
p6 = -0.74246- i*0.74246;
p7 = -0.93556+ i*0.47669;
p8 = -0.93556- i*0.47669;
p9 = -1.0371+ i*0.16426;
p10 = -1.0371- i*0.16426;

%Band Edge specifications
fp1 = 67.4;
fs1 = 63.4;
fs2 = 91.4;
fp2 = 87.4;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 330;         
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi);
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

%Parameters for Bandstop Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10],Wc^N);   %TF with poles p1-p8 and numerator Wc^N and no zeroes
                                                        %numerator chosen to make the DC Gain = 1

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bsf(s) = analog_lpf((s*s + W0*W0)/(B*s));        %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog bsf
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete bsf
[nz, dz] = numden(discrete_bsf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                           %frequency response

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 100e3);
plot((f*2+2.5*10^4),abs(H))
grid