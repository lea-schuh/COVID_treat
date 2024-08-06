function logL = loglikelihood_single_Greek(data1, time, par, treat, t_hosp, data2, data3)

%specify the free parameters
pB = 10^par(1); %immune response on infection rate - negative
pV = 10^par(2);
% dB = 10^par(3);
if treat == 1
    t_shift = floor(10^par(3));
    alpha = 10^par(4);
    t_treat = t_hosp+t_shift;
    sigma = 10^par(5);
    % t_shift = t_shift_t;
else
    t_shift = floor(10^par(3));
    alpha = 1;
    t_treat = 0;
    sigma = 10^par(4);
    % t_shift = t_shift_t;
end

%determine fixed parameter values for i_data = 1 (nasal) and i_data = 2
%(saliva)
S0 = 8*10^7; %total number of epithelial cells in nose at t=0, Ke et al., 2022
dN = 1/11; %death rate of all target cells, Tomasetti et al., 2017
pN = S0*dN; %production of new epithelial cells
b0 = 4.92*10^(-9); %infectivity rate, Ke et al., 2022
dI = 2.45; %death of infected cells, Ke et al., 2022
dV = 10; %deactivation virus, Ke et al., 2022
dB = 0;

%determine the initial values and B_thres for simulation
y0 = [S0, 1, 0, 0];
B_thres = 1-dI*dV/(b0*S0*(pV-dI));
tspan = [0:100]; %time span of solving ODE
options = odeset('NonNegative',[1,2,3,4]); %specify non-negative values
[t,yT] = ode45(@(t,y) odefcn_single_infection_S_R0(t,y,b0,dI,pV,dV,pN,dN,pB,dB,B_thres,alpha,t_treat,S0), tspan, y0,options);

% %evaluate the ODE at daily time points 0 to 20
% y = deval(sol,0:1:100);
% yT = y';

y_short = [];
ind = [];

%map the simulated time points to the observed time points of data
for i_time = 1:length(time)
    ind(i_time) = find(tspan-t_shift == time(i_time));
    y_short(i_time) = yT(ind(i_time),3);
end

%if values too small, fix at 1 (numerical problems)
y_short((y_short<1))=1;

%y_new = -(log(y_short/(1.441*10^14)))/(-0.685); %Ke 2021
y_new = -(log10(y_short)-11.35)/(-0.25); %Ke 2022

n = length(data1)+length(data2)+length(data3); %number of data points

%determine value of negative logL
logL = 0.5*n*log(2*pi*sigma^2)+0.5*nansum((data1-y_new).^2./sigma^2)+0.5*nansum((data2-y_new).^2./sigma^2)+0.5*nansum((data3-y_new).^2./sigma^2);

end