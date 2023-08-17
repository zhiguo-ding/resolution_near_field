clear all
figure
lambda=0.01;%the users' density
N=129;%number of antennas
fc=28*10^9;%carrier frequency
lamb = 3*10^8/fc; %wavelength
d_element = lamb/2; %antenna spacing
d_ray = 2*((N-1)*d_element)^2/lamb; %Rayleigh distance 
R=0.5; %target data rate
noise = 10^(-110/10); %noise power
%PS = 1; %signal power
Psvec = [0:5:30];
Dmax = 1000;%
ct=5000000;
kU = 1; % k th neighbour
alphaN = 4/5; alphaL = 1 - alphaN; % NOMA power allcoation

theta1 = 45 / 180 * pi; % the legacy user's angel NOTE: if N=513, it should be 80
rL1 = 50;% 0.1* d_ray;% the distance of the legacy user
wx = beamfocusing(rL1, theta1, N, d_element, fc)/sqrt(N); % beamforming vector

for ip = 1 : length(Psvec)
    PS = 10^((Psvec(ip)-30)/10);%signal power
 
    %the parameters related to the system models
    eps1 = 2^R-1;
    eta1 = pi^2*(N+1)*(26*N^2-38)/5760/(N-1)^3*(1-sin(theta1)^2)^2*d_ray^2;
    eps2 = noise*eps1/(alphaN-alphaL*eps1);
    eta2 = PS*N*lamb^2/16/pi^2;
    eta3 = 1-eta1/rL1;
    eta4 = sqrt(eps2/eta2);
    eta5 = 2*eta1/(-eta3+sqrt(eta3^2+4*eta1*eta4))-rL1;
    
    %find the second largest root 
    p = [-eta1 2*1/rL1*eta1 (1-1/rL1^2*eta1) 0 -eps2/eta2];
    rzeros = roots(p);
    rzeros_order = sort(rzeros, 'descend');
    
    sum1 = 0; sum3=0;
    for i = 1 : ct
        
        mu0  = (Dmax-rL1)*lambda;  %Dmax-rL1 is the length of the segement
        K = poissrnd(mu0); %the number of NOMA users on the segement
        loc = (Dmax-rL1)*rand(K,1)+rL1; %the locations of the K NOMA users
        rsort = sort(loc-rL1,'ascend'); %ordered distance to the legacy user
        if length(rsort)<kU %not enought users
            sum1 = sum1 +1;
            sum3 = sum3 +1;
        else
            r = rsort(kU); % the distance of the kth neighbour to the legacy user
            rN = r + rL1; % the distance of the kth neighbour to the base station
            bN = beamfocusing(rN, theta1, N, d_element, fc)/sqrt(N); % NOMA directional vector

            rnoma = log2(1+PS*N*lamb^2/16/pi^2/rN^2*alphaN*abs(bN'*wx)^2/...
                (PS*N*lamb^2/16/pi^2/rN^2*alphaL*abs(bN'*wx)^2+noise)); %NOMA user's data rate
            if rnoma<R
                sum1 = sum1 + 1;
            end
            %test
            xtest = eta2/rN^2*(1-eta1*(1/rL1-1/rN)^2);
            ytext = 1/rN*(1-eta1*(1/rL1-1/rN)^2);
            if r>1/rzeros_order(2)- 1/rL1 %xtest<eps2
                sum3 = sum3 + 1;
            end
        end
    end
    psim(ip) = sum1/ct;
    ptest(ip) = sum3/ct;

    %analysis

    
    pnot = 0;
    for i = 0 : kU-1
        pnot = pnot + exp(-lambda*(Dmax-rL1))/factorial(i)*lambda^i*(Dmax-rL1)^i;
    end
    if Dmax-rL1<1/rzeros_order(2)- rL1
        pana(ip) = pnot ;
    else        
        temp = min(max(0,1/rzeros_order(1)- rL1), Dmax-rL1)
        pana(ip) = pnot + (1-pnot)*(gammainc(lambda*(Dmax-rL1),kU)...
            -gammainc(lambda*(max(0,1/rzeros_order(2)- rL1)),kU))%...
            %+gammainc(lambda*(temp),kU));
    end
end

semilogy(Psvec,psim,Psvec,pana)
    

