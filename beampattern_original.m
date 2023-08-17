%courtesy to Zhaolin Wang @ QMUL
clear all
%close all

N =257; % number of antennas
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % wavelength
d = lambda / 2; % antenna spacing

r = 5.5; % distance
theta = 45 / 180 * pi; % direction
w = beamfocusing(r, theta, N, d, f); % beamforming vector

% calculate beampattern
m = 1000;
X = linspace(0, 8, m);
Y = linspace(0, 8, m);
[X, Y] = meshgrid(X, Y);
[theta_all, r_all] = cart2pol(X, Y);

P = zeros(length(r_all), length(theta_all));
parfor i = 1:m
    for j = 1:m
        a = beamfocusing(r_all(i,j), theta_all(i,j), N, d, f);
        P(i,j) = abs(a' * w)^2;
    end
end
P = P ./ max(max(P));

% plot beampattern
[X, Y] = pol2cart(theta_all, r_all);
figure; hold on; colormap jet; colorbar;
mesh(X,Y,P); 
view([90,-90]);
xlim([0,8]); ylim([0,8]);
xlabel('x (m)'); ylabel('y (m)');
