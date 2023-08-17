%courtesy to Zhaolin Wang @ QMUL
%clc
clear all
%close all

N =513; % number of antennas
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % wavelength
d = lambda / 2; % antenna spacing
D = (N-1)*d; %apture or size
dis_Ray = 2*D^2/lambda;
dis_rea = 0.62*sqrt(D^3/lambda);

dmax = 800;

r = dis_rea*2% dis_Ray*0.2; % distance
theta = 45 / 180 * pi; % direction
w = beamfocusing(r, theta, N, d, f)/sqrt(N); % beamforming vector

 

% calculate beampattern
m = 1000;
X = linspace(0, dmax, m);
Y = linspace(0, dmax, m);
[X, Y] = meshgrid(X, Y);
[theta_all, r_all] = cart2pol(X, Y);

P = zeros(length(r_all), length(theta_all));
parfor i = 1:m
    for j = 1:m
        a = beamfocusing(r_all(i,j), theta_all(i,j), N, d, f);
        P(i,j) = abs(a' * w)^2;
        if P(i,j)>0.99999
            dfd=0;
        end
    end
end
P = P ./ max(max(P));

% plot beampattern
[X, Y] = pol2cart(theta_all, r_all);
figure; hold on; colormap jet; colorbar;
mesh(X,Y,P); 
view([90,-90]);
xlim([0,dmax]); ylim([0,dmax]);
xlabel('x (m)'); ylabel('y (m)');
