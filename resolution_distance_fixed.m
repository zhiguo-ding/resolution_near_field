%courtesy to Zhaolin Wang @ QMUL
%clc
clear all
%close all
Nvec = [11 101 1001   ];
theta_vec = [0 : 10 : 90];
%N =123; % number of antennas
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % wavelength
d = lambda / 2; % antenna spacing

ratio_dis1 = 0.5;
ratio_dis2 =  0.1;

for in = 1 : length(theta_vec)
%for in = 1 : length(Nvec)
    N = 513;%Nvec(in);
    D = (N-1)*d; %apture or size
    dis_Ray = 2*D^2/lambda;
    dis_rea = 0.62*sqrt(D^3/lambda);

    r0 =   10;%  dis_Ray*ratio_dis1; % distance
    theta = theta_vec(in) / 180 * pi;%45 / 180 * pi; % direction
    %theta =  80 / 180 * pi; % direction
    wx = beamfocusing(r0, theta, N, d, f); % beamforming vector
    w = wx/sqrt(N);

    rm = r0+20;%dis_Ray*ratio_dis2;
    h = beamfocusing(rm, theta+0/180*pi, N, d, f);
    actual(in)=abs(w'*h)^2/N

  

    %existing approximation
    beta=sqrt(N^2*d^2*(1-sin(theta)^2)/2/lambda*abs(1/r0-1/rm));
    cbeta = 0;
    stepb = beta/100;
    beta_vec = [0:stepb:beta];
    for i=1: length(beta_vec)
        t = beta_vec(i);
        cbeta = cbeta + cos(pi/2*t^2)*stepb;
    end
    sbeta = 0; 
    for i=1: length(beta_vec)
        t = beta_vec(i);
        sbeta = sbeta + sin(pi/2*t^2)*stepb;
    end    
    existing(in)=abs((cbeta+1i*sbeta)/beta)^2;


    %new approximation
    tau = (1-sin(theta)^2)*abs(1/ratio_dis1  -1/(ratio_dis1+ratio_dis2));
    approi(in)=1+pi^2*tau^2/24^2*(N+1)^2/(N-1)^2 -pi^2*tau^2/480*(N+1)*(3*N^2-4)/(N-1)^3;
    approi2(in)=1 - 13*pi^2*tau^2/2880;
    zz(in) = abs(approi2(in)-approi2(in));
end
%plot(Nvec, actual, Nvec, existing, Nvec, approi)
%semilogx(Nvec, actual, Nvec, existing, Nvec, approi)
plot(theta_vec, actual)%, theta_vec, existing)%, theta_vec, approi)%,  theta_vec, approi2)