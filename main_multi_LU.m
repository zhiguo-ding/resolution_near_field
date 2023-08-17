clear all

lambda=0.01;%the users' density
N=513;%number of antennas
M = 50;%number of legacy users
fc=28*10^9;%carrier frequency
lamb = 3*10^8/fc; %wavelength
d_element = lamb/2; %antenna spacing
d_ray = 2*((N-1)*d_element)^2/lamb; %Rayleigh distance 
R=0.5; %target data rate
noise = 10^(-110/10); %noise power
%PS = 1; %signal power
Psvec = [0:5:30];
Dmax = 100;%
ct=500;
kU = 1; % k th neighbour
alphaN = 4/5; alphaL = 1 - alphaN; % NOMA power allcoation

wx_LU = [];
theta_m = linspace(-90,90,M);% 
rL_m = 50*ones(M,1);
for m = 1: length(theta_m)
    theta1 = theta_m(m) / 180 * pi; % the legacy user's angel
    rL = rL_m(m);% 0.1* d_ray;% the distance of the legacy user
    wx = beamfocusing(rL, theta1, N, d_element, fc)/sqrt(N); % beamforming vector
    wx_LU = [wx_LU wx];
end

for ip =  1: length(Psvec)
    PS = 10^((Psvec(ip)-30)/10);%signal power
 
    %the parameters related to the system models
    eps1 = 2^R-1;
    eta1 = pi^2*(N+1)*(26*N^2-38)/5760/(N-1)^3*(1-sin(theta1)^2)^2*d_ray^2;
    eps2 = noise*eps1/(alphaN-alphaL*eps1);
    eta2 = PS*N*lamb^2/16/pi^2;
    eta4 = sqrt(eps2/eta2);
     
    sum1 = 0; sum2=0;
    for i = 1 : ct        
        K = 100;% the number of NOMA users
        NF_loc=[]; distance = []; angle=[];  theta_noma = [];
        while size(NF_loc,1)<K
            x_loc = [Dmax*rand(1,1) sign(randn)*Dmax*rand(1,1)];
            if sqrt(x_loc*x_loc')<Dmax  
                NF_loc = [NF_loc; x_loc];
                distance = [distance; max(1,sqrt(x_loc*x_loc'))];
                %step 1 - group users
                theta_temp = atan(NF_loc(end,2)/NF_loc(end,1))*180/pi;
                theta_noma = [theta_noma theta_temp];
                [xt,indx] = min(abs(theta_temp-theta_m));
                angle = [angle;indx]; %store the index, not actual angel
            end
        end

        %step 2 - schedule a single NOMA user from each group 
        index_schedule = [];
        for m= 1 : M
            index_mb = find(angle==m);%the index of users on this beam
            nub_mb = length(index_mb); %the number of NOMA users on this beam
            if nub_mb==0
                index_schedule = [index_schedule 0]; 
            else
                xt = NF_loc(index_mb,:)-[rL_m(m)*cos(theta_m(m)) rL_m(m)*sin(theta_m(m))];
                [yt, indexu] = min(sqrt(sum(xt.^2,2)));
                index_schedule = [index_schedule index_mb(indexu)]; 
            end
        end
           
         
        if index_schedule(1)==0 %no NOMA user on beam 1
            sum1 = sum1 + 1;
            sum2 = sum2 + 1;
        else
            % the channel vector of NOMA user on beam 1 - index_schedule(1) 
            dN = distance(index_schedule(1));
            bN = beamfocusing(dN, ...
                theta_noma(index_schedule(1)), N, d_element, fc)/sqrt(N); % beamforming vector;%
            %find the inter-beam interference
            interf = 0;
            for im = 2: M
                interf = interf + PS*N*lamb^2/16/pi^2/dN^2*abs(bN'*wx_LU(:,im))^2;
            end
        
            %the NOMA user's data rate
            rnoma = log2(1+PS*N*lamb^2/16/pi^2/dN^2*alphaN*abs(bN'*wx_LU(:,1))^2/...
                (PS*N*lamb^2/16/pi^2/dN^2*alphaL*abs(bN'*wx_LU(:,1))^2+interf+noise)); %NOMA user's data rate
            if rnoma<R
                sum1 = sum1 + 1;
            end
        end
         
    end
    psim(ip) = sum1/ct;
    pempty(ip) = sum2/ct;
 
 
end

semilogy(Psvec,psim, Psvec,pempty )
    

