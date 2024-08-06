clearvars;
clc;

%load data set
data = readtable('Data_Greek_min_text.xlsx');

for i_patient = 1:size(data,1)
    for i_date = 1:10
        if ~(data.(sprintf('Date%d',i_date))(i_patient) == "")
            Dates{i_patient}(i_date) = data.(sprintf('Date%d',i_date))(i_patient);
            Vals{i_patient}(i_date) = data.(sprintf('Val%d',i_date))(i_patient);
        end
    end
    if ~isnan(data.('Treatment')(i_patient))
        Treat(i_patient) =  data.('Treatment')(i_patient);
    else
        Treat(i_patient) = 0;
    end
    Hosp(i_patient) = data.('Hospitalization')(i_patient);
    Hosp_dates{i_patient} = data.('DateOfHospitalization')(i_patient);
    Death(i_patient) = data.('Death')(i_patient);
    OnsetSympt(i_patient) = data.('OnsetSymptoms_daysBeforeDiagnosis_')(i_patient);
    ID(i_patient) = data.('ID')(i_patient);
end

icount = 1;
for i_patient = 1:size(data,1)
    if length(Dates{i_patient}) >= 4
        Dates_set{icount} = Dates{i_patient};
        Val_set{icount} = Vals{i_patient};
        Treat_set(icount) = Treat(i_patient);
        Hosp_set(icount) = Hosp(i_patient);
        Hosp_dates_set{icount} = Hosp_dates{i_patient};
        Death_set(icount) = Death(i_patient);
        OnsetSympt_set(icount) = OnsetSympt(i_patient);
        ID_set(icount) = ID(i_patient);
        icount = icount+1;
    end
end

days_month = [31,28,31,30,31,30,31,31,30,31,30,31];

for i_set = 1:size(Dates_set,2)
    clear month day time_new time val month_hosp day_hosp month_data1 day_data1 val1 val3
    for i_dates = 1:length(Dates_set{i_set})
        month(i_dates) = str2double(extractBefore(Dates_set{i_set}(i_dates),'/'));
        day(i_dates) = str2double(extractBetween(Dates_set{i_set}(i_dates),'/','/'));
        if i_dates > 1
            if month(i_dates) == month(i_dates-1)
                time_new = day(i_dates)-day(i_dates-1);
            elseif month(i_dates) == month(i_dates-1)+1
                time_new = days_month(month(i_dates-1))-day(i_dates-1) + day(i_dates);
            elseif month(i_dates) == 1 && month(i_dates-1) == 12
                time_new = days_month(month(i_dates-1))-day(i_dates-1) + day(i_dates);
            end
            time(i_dates) = time(i_dates-1)+time_new;
        else
            time(i_dates) = 0;
        end

        if Val_set{i_set}(i_dates) == ""
            val(i_dates) = NaN;
            val1(i_dates) = NaN;
            val3(i_dates) = NaN;
        else
            if strcmp(Val_set{i_set}(i_dates),'NEG (-)')
                % val(i_dates) = 42;
                val(i_dates) = 35; %according to Ioanna
                val1(i_dates) = 35; %according to Ioanna
                val3(i_dates) = 35; %according to Ioanna
            else
                val(i_dates) = str2double(extractBetween(Val_set{i_set}(i_dates),'/','/'));
                val1(i_dates) = str2double(extractBefore(Val_set{i_set}(i_dates),'/'));
                val3(i_dates) = str2double(extractAfter((extractAfter(Val_set{i_set}(i_dates),'/')),'/'));
            end
        end

    end

    T{i_set} = time;
    CT{i_set} = val;
    CT1{i_set} = val1;
    CT3{i_set} = val3;

    if Hosp_set(i_set) == 1
        if (Hosp_dates_set{i_set} == "") == 0
            month_hosp = str2double(extractBefore(Hosp_dates_set{i_set},'/'));
            day_hosp = str2double(extractBetween(Hosp_dates_set{i_set},'/','/'));
            month_data1 = str2double(extractBefore(Dates_set{i_set}(1),'/'));
            day_data1 = str2double(extractBetween(Dates_set{i_set}(1),'/','/'));
            if month_hosp == month_data1 %hospitalized in same month as first measurement
                T_hosp(i_set) = day_hosp-day_data1;
            elseif month_hosp == month_data1+1 %if hospitalization in month after first measurement
                T_hosp(i_set)  = days_month(month_data1)-day_data1 + day_hosp;
            elseif month_hosp == 1 && month_data1 == 12 %over year
                T_hosp(i_set)  = days_month(month_data1)-day_data1 + day_hosp;
            elseif month_hosp == month_data1-1
                T_hosp(i_set)  = -(day_data1 + days_month(month_hosp)-day_hosp);
            end
        else
            T_hosp(i_set) = NaN;
        end
    else
        T_hosp(i_set) = NaN;
    end

end

icount = 1;
for i_t = 1:length(T)
    if issorted(T{i_t})
        if sum(T{i_t} >= 0) == length(T{i_t})
            % if ~isnan(CT{i_t})
            T_check{icount} = T{i_t};
            CT_check{icount} = CT{i_t};
            CT1_check{icount} = CT1{i_t};
            CT3_check{icount} = CT3{i_t};
            Treat_check(icount) = Treat_set(i_t);
            T_hosp_check(icount) = T_hosp(i_t);
            Death_check(icount) = Death_set(i_t);
            OnsetSympt_check(icount) = OnsetSympt_set(i_t);
            ID_check(icount) = ID_set(i_t);
            icount = icount + 1;
            % end
        end
    end
end

%
inds_treat = find(Treat_check == 1);
inds_wo_treat = find(Treat_check == 0);

inds_hosp = find(~isnan(T_hosp_check));
inds_wo_hosp = find(isnan(T_hosp_check));

inds_dead = find(Death_check == 1);
inds_wo_dead = find(Death_check == 0);

inds_treat_hosp = intersect(inds_treat,inds_hosp);

Ind = [{inds_treat_hosp},0,{inds_wo_treat},0,{inds_hosp},[1:37,39:length(ID_check)],[1:37,39:length(ID_check)],[1:37,39:length(ID_check)]];

i_treat = 1; %with treatment
%i_treat = 3; %without treatment
%i_treat = 5; %both
%i_treat = 6; %untreated vs treated
%i_treat = 7; %untreated unhospitalized vs treated hospitalized (non-severe vs severe)
%i_treat = 8; %unhosp vs all hosp

if i_treat == 1
    load('sol_treat_3_t0')
    sol_notreat = load('sol_treat_no_treat_3_t0');
    t = tiledlayout(3,7,'TileSpacing','Compact','Padding','Compact');
elseif i_treat == 3
    load('sol_wo_treat_3_t0')
    t = tiledlayout(10,7,'TileSpacing','Compact','Padding','Compact');
else
    sol1 = load('sol_treat_3_t0');
    sol2 = load('sol_wo_treat_3_t0');
    for i = 1:length(sol2.sol)
        if isempty(sol2.sol{i})
            sol{i} = sol1.sol{i};
        else
            sol{i} = sol2.sol{i};
        end
    end
end


%fit over every patient seperately
for ID_opt = Ind{i_treat}

    clearvars -except C Ind i_treat ID_check i_ind data ID_opt S0 dN pN...
        b0 dI dV dB T_check CT_check CT1_check CT3_check...
        Treat_check T_hosp_check Death_check sol inds_treat_hosp inds_wo_treat inds_treat...
        inds_hosp inds_wo_hosp inds_dead inds_wo_dead OnsetSympt_check MSE Ind max_y max_V y_new_50 i_treat...
        ME sigma MSE_data alpha t_undetectV B_thres B_thres_dyn BgeqBthres BgeqBthres_dyn minS CN_sympt...
        pB pV t_shift areaVL areaVLtopeak areaVLafterpeak ind_max_y BgeqBthres_dyn_afterpeak areaBgeqBthres_dyn_afterpeak...
        R areaR max_t meanVLafterpeak Y_all sol_notreat BIC_wtreat...
        B_thres_notreat pV_notreat pB_notreat sigma_notreat t_shift_notreat BIC_notreat ME_notreat...
        t_detect meanVL meanVLtopeak K_all M_all CN_hosp

    ID_opt

    T_check{ID_opt} = T_check{ID_opt} + OnsetSympt_check(ID_opt);
    T_hosp_check(ID_opt) = T_hosp_check(ID_opt) + OnsetSympt_check(ID_opt);
    OnsetSympt_check(ID_opt) = 0;

    %set free parameters to best estimates from optimzation
    if i_treat == 1
        pB(ID_opt) = sol{ID_opt}.P_t(1); %immune response on infection rate
        pV(ID_opt) = sol{ID_opt}.P_t(2);
        % dB = 10^par_opt(ind(1),3);
        t_shift(ID_opt) = floor(sol{ID_opt}.P_t(3));
        alpha(ID_opt) = sol{ID_opt}.P_t(4);
        sigma(ID_opt) = sol{ID_opt}.P_t(5);

        pB_notreat(ID_opt) = sol_notreat.sol{ID_opt}.P(1); %immune response on infection rate
        pV_notreat(ID_opt) = sol_notreat.sol{ID_opt}.P(2);
        % dB = 10^par_opt(ind(1),3);
        t_shift_notreat(ID_opt) = floor(sol_notreat.sol{ID_opt}.P(3));
        sigma_notreat(ID_opt) = sol_notreat.sol{ID_opt}.P(4);
    elseif i_treat == 3
        pB(ID_opt) = sol{ID_opt}.P(1); %immune response on infection rate
        pV(ID_opt) = sol{ID_opt}.P(2);
        % dB = 10^par_opt(ind(1),3);
        t_shift(ID_opt) = floor(sol{ID_opt}.P(3));
        sigma(ID_opt) = sol{ID_opt}.P(4);
    else
        if isfield(sol{ID_opt},'P')
            pB(ID_opt) = sol{ID_opt}.P(1); %immune response on infection rate
            pV(ID_opt) = sol{ID_opt}.P(2);
            % dB = 10^par_opt(ind(1),3);
            t_shift(ID_opt) = floor(sol{ID_opt}.P(3));
            sigma(ID_opt) = sol{ID_opt}.P(4);
        else
            pB(ID_opt) = sol{ID_opt}.P_t(1); %immune response on infection rate
            pV(ID_opt) = sol{ID_opt}.P_t(2);
            % dB = 10^par_opt(ind(1),3);
            t_shift(ID_opt) = floor(sol{ID_opt}.P_t(3));
            alpha(ID_opt) = sol{ID_opt}.P_t(4);
            sigma(ID_opt) = sol{ID_opt}.P_t(5);
        end
    end

    S0 = 8*10^7; %total number of epithelial cells in nose at t=0, Ke et al., 2022
    dN = 1/11; %death rate of all target cells, Tomasetti et al., 2017
    pN = S0*dN; %production of new epithelial cells
    b0 = 4.92*10^(-9); %infectivity rate, Ke et al., 2022
    dI = 2.45; %death of infected cells, Ke et al., 2022
    dV = 10; %deactivation virus, Ke et al., 2022
    dB = 0;

    %plot the results for short- and long-term dynamics
    tspan = 0:0.1:100; %long term

    %determine the initial values and B_thres for simulation
    y0 = [S0, 1, 0, 0];
    B_thres(ID_opt) = 1-dI*dV/(b0*S0*(pV(ID_opt)-dI));
    options = odeset('NonNegative',[1,2,3,4]); %specify non-negative values
    if i_treat == 1
        [t,y] = ode45(@(t,y) odefcn_single_infection_S_R0(t,y,b0,dI,pV(ID_opt),dV,pN,dN,pB(ID_opt),dB,B_thres(ID_opt),alpha(ID_opt),T_hosp_check(ID_opt)+t_shift(ID_opt),S0), tspan, y0, options);
        B_thres_notreat(ID_opt) = 1-dI*dV/(b0*S0*(pV_notreat(ID_opt)-dI));
        [t_notreat,y_notreat] = ode45(@(t,y) odefcn_single_infection_S_R0(t,y,b0,dI,pV_notreat(ID_opt),dV,pN,dN,pB_notreat(ID_opt),dB,B_thres_notreat(ID_opt),1,0,S0), tspan, y0,options);
    elseif i_treat == 3
        [t,y] = ode45(@(t,y) odefcn_single_infection_S_R0(t,y,b0,dI,pV(ID_opt),dV,pN,dN,pB(ID_opt),dB,B_thres(ID_opt),1,0,S0), tspan, y0,options);
    else
        if isfield(sol{ID_opt},'P')
           [t,y] = ode45(@(t,y) odefcn_single_infection_S_R0(t,y,b0,dI,pV(ID_opt),dV,pN,dN,pB(ID_opt),dB,B_thres(ID_opt),1,0,S0), tspan, y0,options);
        else
           [t,y] = ode45(@(t,y) odefcn_single_infection_S_R0(t,y,b0,dI,pV(ID_opt),dV,pN,dN,pB(ID_opt),dB,B_thres(ID_opt),alpha(ID_opt),T_hosp_check(ID_opt)+t_shift(ID_opt),S0), tspan, y0, options);
        end
    end

    if isempty(find(y(1:501,4)>B_thres(ID_opt),1,'first'))
        BgeqBthres(ID_opt) = 50;
    else
        BgeqBthres(ID_opt) = (find(y(1:501,4)>B_thres(ID_opt),1,'first')-1)*0.1;
    end

    % tmax = floor(tspan(find(max(y(:,3)) == y(:,3))));
    B_thres_dyn(:,ID_opt) = 1-dI.*dV./(b0.*y(:,1).*(pV(ID_opt)-dI));
    % BgeqBthres_dyn(ID_opt) = length(find(y(1:501,4)>B_thres_dyn(1:501,ID_opt)));
    %BgeqBthres_dyn(ID_opt) = sum(y(tmax:501,4)-B_thres_dyn(tmax:501,ID_opt));

    % if isempty(find(y(1:501,4)<B_thres_dyn(1:501,ID_opt),1,'last'))
    %     % BgeqBthres_dyn(ID_opt) = 50;
    % else
    % BgeqBthres_dyn(ID_opt) = length(find(y(1:501,4)>B_thres_dyn(1:501,ID_opt)));
    BgeqBthres_dyn(ID_opt)  = (find(y(1:501,4)<B_thres_dyn(1:501,ID_opt),1,'last')-1)*0.1;
    % end

    R = (pV(ID_opt).*b0.*(1-y(1:501,4)).*y(1:501,1))./(dI.*(dV+b0.*(1-y(1:501,4)).*y(1:501,1)));

    if i_treat == 1
         nexttile
         if ismember(ID_opt,intersect(inds_treat_hosp,inds_wo_dead))
            % plot(tspan(1:501),y(1:501,4),'-','Color','k')
            plot(tspan(1:501),(pV(ID_opt).*b0.*(1-y(1:501,4)).*y(1:501,1))./(dI.*(dV+b0.*(1-y(1:501,4)).*y(1:501,1))),'-','Color','k');
            hold on
            yline(1)
            % plot(tspan(1:501),B_thres_dyn(1:501,ID_opt))
        elseif ismember(ID_opt,intersect(inds_treat_hosp,inds_dead))
            % plot(tspan(1:501),y(1:501,4),'-','Color','b')
            % hold on
            % plot(tspan(1:501),B_thres_dyn(1:501,ID_opt))
            plot(tspan(1:501),(pV(ID_opt).*b0.*(1-y(1:501,4)).*y(1:501,1))./(dI.*(dV+b0.*(1-y(1:501,4)).*y(1:501,1))),'-','Color','b');
            hold on
            yline(1)
        end
    elseif i_treat == 3
        nexttile
        if ismember(ID_opt,intersect(inds_wo_treat,inds_wo_hosp))
            % plot(tspan(1:501),y(1:501,4),'-','Color','k')
            % hold on
            % plot(tspan(1:501),B_thres_dyn(1:501,ID_opt))

            plot(tspan(1:501),(pV(ID_opt).*b0.*(1-y(1:501,4)).*y(1:501,1))./(dI.*(dV+b0.*(1-y(1:501,4)).*y(1:501,1))),'-','Color','k');
            % hold on
            yline(1)
        elseif ismember(ID_opt,intersect(inds_wo_treat,inds_hosp))
            % plot(tspan(1:501),y(1:501,4),'-','Color','b')
            % hold on
            % plot(tspan(1:501),B_thres_dyn(1:501,ID_opt))

            plot(tspan(1:501),(pV(ID_opt).*b0.*(1-y(1:501,4)).*y(1:501,1))./(dI.*(dV+b0.*(1-y(1:501,4)).*y(1:501,1))),'-','Color','b');
            % hold on
            yline(1)
        end
    end
    % ylim([0,1])
     ylim([0,5])

    minS(ID_opt) = min(y(:,1)./S0);

    time = T_check{ID_opt};

    max_t(ID_opt) = t(find(y(:,3) == max(y(:,3))));
    %max_V(ID_opt) = (log10(max(y(:,3)))-11.35)/(-0.25);
    max_V(ID_opt) = y(find(y(:,3) == max(y(:,3))),4);
    % max_V(ID_opt) = R(find(y(:,3) == max(y(:,3))));

    ind_max_y(ID_opt) = find(y(:,3) == max(y(:,3)));

    % R = (pV(ID_opt).*b0.*(1-y(1:501,4)).*y(1:501,1))./(dI.*(dV+b0.*(1-y(1:501,4)).*y(1:501,1)));
    R_afterpeak = R(ind_max_y(ID_opt):end);
    areaR(ID_opt) = sum(R_afterpeak(R_afterpeak>1)*0.1);

    %map the simulated time points to the observed time points of data
    for i_time = 1:length(time)
        ind = find(tspan-t_shift(ID_opt) == time(i_time));
        y_short(i_time) = y(ind,3);
    end

    %if values too small, fix at 1 (numerical problems)
    y_short((y_short<1))=1;

    %y_new = -(log(y_short/(1.441*10^14)))/(-0.685); %Ke 2021
    y_new = -(log10(y_short)-11.35)/(-0.25); %Ke 2022

    if i_treat == 1
        %map the simulated time points to the observed time points of data
        for i_time = 1:length(time)
            ind = find(tspan-t_shift_notreat(ID_opt) == time(i_time));
            y_short_notreat(i_time) = y_notreat(ind,3);
        end

        %if values too small, fix at 1 (numerical problems)
        y_short_notreat((y_short_notreat<1))=1;

        %y_new = -(log(y_short/(1.441*10^14)))/(-0.685); %Ke 2021
        y_new_notreat = -(log10(y_short_notreat)-11.35)/(-0.25);
    end

    Y_all(:,ID_opt) = -(log10(y(1:501,3))-11.35)/(-0.25);

    y_new_50(ID_opt) = -(log10(y(501,3))-11.35)/(-0.25);
    % if isempty(find(-(log10(y(151:501,3))-11.35)/(-0.25)<-35,1,'first'))
    %     t_undetectV(ID_opt) = 50;
    % else
    %     t_undetectV(ID_opt) = tspan(151+find(-(log10(y(151:501,3))-11.35)/(-0.25)<-35,1,'first'));
    % end

    VL_all = y(2:501,4);
    VL_topeak = y(2:ind_max_y(ID_opt),4);
    VL_afterpeak = y(ind_max_y(ID_opt):501,4);
    areaVL(ID_opt) = sum(VL_all)*0.1;
    areaVLtopeak(ID_opt) = sum(VL_topeak)*0.1;
    areaVLafterpeak(ID_opt) = sum(VL_afterpeak)*0.1;
    meanVL(ID_opt) = areaVL(ID_opt)*10/length(VL_all);
    meanVLtopeak(ID_opt) = areaVLtopeak(ID_opt)*10/length(VL_topeak);
    meanVLafterpeak(ID_opt) = areaVLafterpeak(ID_opt)*10/length(VL_afterpeak);

    % VL_all = R(2:501);
    % VL_topeak = R(2:ind_max_y(ID_opt));
    % VL_afterpeak = R(ind_max_y(ID_opt):501);
    % areaVL(ID_opt) = sum(VL_all-1)*0.1;
    % areaVLtopeak(ID_opt) = sum(VL_topeak-1)*0.1;
    % areaVLafterpeak(ID_opt) = sum(VL_afterpeak-1)*0.1;
    % meanVL(ID_opt) = areaVL(ID_opt)*10/length(VL_all)+1;
    % meanVLtopeak(ID_opt) = areaVLtopeak(ID_opt)*10/length(VL_topeak)+1;
    % meanVLafterpeak(ID_opt) = areaVLafterpeak(ID_opt)*10/length(VL_afterpeak)+1;

    if isempty(find(y(ind_max_y(ID_opt):501,4)<B_thres_dyn(ind_max_y(ID_opt):501,ID_opt),1,'first'))
        BgeqBthres_dyn_afterpeak(ID_opt)  = 0;
    else
        BgeqBthres_dyn_afterpeak(ID_opt)  = (find(y(ind_max_y(ID_opt):501,4)<B_thres_dyn(ind_max_y(ID_opt):501,ID_opt),1,'first')-1)*0.1;
    end

    diff_B = B_thres_dyn(ind_max_y(ID_opt):501,ID_opt)-y(ind_max_y(ID_opt):501,4);
    % areaBgeqBthres_dyn_afterpeak(ID_opt) = sum(diff_B(diff_B>0)*0.1);
    areaBgeqBthres_dyn_afterpeak(ID_opt) = sum((diff_B>0));
    % areaBgeqBthres_dyn_afterpeak(ID_opt) = sum(diff_B*0.1);

    CN_sympt(ID_opt) = y(find(tspan==t_shift(ID_opt)),4);
     % CN_sympt(ID_opt) = R(find(tspan==t_shift(ID_opt)));
    if ismember(i_treat,[1,6])
        CN_hosp(ID_opt) = y(find(tspan==t_shift(ID_opt)+T_hosp_check(ID_opt)),4);
    end
    %CN_sympt(ID_opt) = y(find(tspan==t_shift),3);

    t_detect(ID_opt) = (find(Y_all(:,ID_opt)>-35,1,'first')-1)*0.1;

    n = sum(~isnan(CT_check{ID_opt}))+sum(~isnan(CT1_check{ID_opt}))+sum(~isnan(CT3_check{ID_opt})); %number of data points

    %determine value of negative logL
    % logL = 0.5*n*log(2*pi*sigma^2)+0.5*nansum((data1-y_new).^2./sigma^2);
    MSE(ID_opt) = 1/n*(nansum((-CT_check{ID_opt}-y_new).^2)+nansum((-CT1_check{ID_opt}-y_new).^2)+nansum((-CT3_check{ID_opt}-y_new).^2));
    ME(ID_opt) = 1/n*(nansum(abs(-CT_check{ID_opt}-y_new))+nansum(abs(-CT1_check{ID_opt}-y_new))+nansum(abs(-CT3_check{ID_opt}-y_new)));

    MSE_data(ID_opt) = 1/(sum(~isnan([(-CT_check{ID_opt}+CT1_check{ID_opt}),(-CT_check{ID_opt}+CT3_check{ID_opt})])))*(nansum(abs(-CT_check{ID_opt}+CT1_check{ID_opt}))+nansum(abs(-CT_check{ID_opt}+CT3_check{ID_opt})));
    
     if i_treat == 1

         logL_wtreat = 0.5*n*log(2*pi*sigma(ID_opt)^2)+0.5*nansum((-CT1_check{ID_opt}-y_new).^2./sigma(ID_opt)^2)+...
             0.5*nansum((-CT_check{ID_opt}-y_new).^2./sigma(ID_opt)^2)+0.5*nansum((-CT3_check{ID_opt}-y_new).^2./sigma(ID_opt)^2);

         BIC_wtreat(ID_opt) = 5*log(n)+2*logL_wtreat;

         logL_notreat = 0.5*n*log(2*pi*sigma_notreat(ID_opt)^2)+0.5*nansum((-CT1_check{ID_opt}-y_new_notreat).^2./sigma_notreat(ID_opt)^2)+...
             0.5*nansum((-CT_check{ID_opt}-y_new_notreat).^2./sigma_notreat(ID_opt)^2)+0.5*nansum((-CT3_check{ID_opt}-y_new_notreat).^2./sigma_notreat(ID_opt)^2);

         BIC_notreat(ID_opt) = 4*log(n)+2*logL_notreat;

         ME_notreat(ID_opt) = 1/n*(nansum(abs(-CT_check{ID_opt}-y_new_notreat))+nansum(abs(-CT1_check{ID_opt}-y_new_notreat))+nansum(abs(-CT3_check{ID_opt}-y_new_notreat)));
       
     end
end

if i_treat == 1
    inds1 = intersect(inds_treat_hosp,inds_wo_dead);
    %inds1 = intersect(inds1,find(alpha>0 & alpha<0.9));
    inds2 = intersect(inds_treat_hosp,inds_dead);
    %inds2 = intersect(inds2,find(alpha>0 & alpha<0.9));
elseif i_treat == 3
    inds1 = intersect(inds_wo_treat,inds_wo_hosp);
    inds2 = intersect(inds_wo_treat,inds_hosp);
elseif i_treat == 5
    inds1 = intersect(intersect(inds_wo_treat,inds_hosp),inds_wo_dead);
    inds2 = intersect(inds_treat_hosp,inds_wo_dead);
    %inds2 = intersect(inds2,find(alpha>0 & alpha<0.9));
    %inds2 = intersect(inds2,find(alpha>0.9));
elseif i_treat == 6
    inds1 = inds_wo_treat;
    inds2 = inds_treat_hosp;
elseif i_treat == 7
    inds1 = intersect(inds_wo_treat,inds_wo_hosp);
    inds2 = inds_treat_hosp;
else
    inds1 = intersect(inds_wo_treat,inds_wo_hosp);
    inds2 = union(inds_treat_hosp,intersect(inds_wo_treat,inds_hosp));
end

MSE1 = MSE(inds1);
MSE2 = MSE(inds2);
% mediantest(MSE_wo_nohosp,MSE_wo_hosp)

t_detect1 = t_detect(inds1);
t_shift1 = t_shift(inds1);
pV1 = pV(inds1);
pB1 = pB(inds1);
if i_treat == 1
    alpha1 = alpha(inds1);
    T_hosp_check1 = T_hosp_check(inds1);
end
if ismember(i_treat, [5,6])
    T_hosp_check1 = T_hosp_check(inds1);
end
peak_VL1 = max_V(inds1);
tpeak_VL1 = max_t(inds1);
V50_1 = y_new_50(inds1);
V50_1_bin = sum(-V50_1<35);
areaVL1 = areaVL(inds1);
areaVLtopeak1 = areaVLtopeak(inds1);
areaVLafterpeak1 = areaVLafterpeak(inds1);
meanVL1 = meanVL(inds1);
meanVLafterpeak1 = meanVLafterpeak(inds1);
meanVLtopeak1 = meanVLtopeak(inds1);
areaR1 = areaR(inds1);
BgeqBthres_dyn1 = BgeqBthres_dyn(inds1);
BgeqBthres_dyn_afterpeak1 = BgeqBthres_dyn_afterpeak(inds1);
areaBgeqBthres_dyn_afterpeak1 = areaBgeqBthres_dyn_afterpeak(inds1);
CN_sympt1 = CN_sympt(inds1);
if ismember(i_treat,[1,6])
    CN_hosp1 = CN_hosp(inds1);
end
sigma1 = sigma(inds1);
MSE_data1 = MSE_data(inds1);

t_detect2 = t_detect(inds2);
t_shift2 = t_shift(inds2);
pV2 = pV(inds2);
pB2 = pB(inds2);
if i_treat == 1
    alpha2 = alpha(inds2);
    T_hosp_check2 = T_hosp_check(inds2);
end
if ismember(i_treat, [5,6])
    T_hosp_check2 = T_hosp_check(inds2);
end
peak_VL2 = max_V(inds2);
tpeak_VL2 = max_t(inds2);
areaVL2 = areaVL(inds2);
areaVLtopeak2 = areaVLtopeak(inds2);
areaVLafterpeak2 = areaVLafterpeak(inds2);
meanVL2 = meanVL(inds2);
meanVLafterpeak2 = meanVLafterpeak(inds2);
meanVLtopeak2 = meanVLtopeak(inds2);
areaR2 = areaR(inds2);
V50_2 = y_new_50(inds2);
V50_2_bin = sum(-V50_2<35);
BgeqBthres_dyn2 = BgeqBthres_dyn(inds2);
BgeqBthres_dyn_afterpeak2 = BgeqBthres_dyn_afterpeak(inds2);
areaBgeqBthres_dyn_afterpeak2 = areaBgeqBthres_dyn_afterpeak(inds2);
CN_sympt2 = CN_sympt(inds2);
if ismember(i_treat,[1,6])
    CN_hosp2 = CN_hosp(inds2);
end
sigma2 = sigma(inds2);
MSE_data2 = MSE_data(inds2);

sprintf('median MSE = %d', median(ME(ME>0)))
[h,p,ks2stat] = kstest2(MSE1,MSE2);
sprintf('p-value MSE differences = %d', p)

sprintf('median value of sigma for pop1 = %d', median(sigma1))
sprintf('median value of sigma for pop2= %d', median(sigma2))
[h,p,ks2stat] = kstest2(sigma1,sigma2);
sprintf('p-value differences in sigma = %d', p)

[h,p,ks2stat] = kstest2(MSE_data1,MSE_data2);
sprintf('median value of absolute differences across replicates for pop1 = %d', median(MSE_data1))
sprintf('median value of absolute differences across replicates for pop2 = %d',median(MSE_data2))
sprintf('p-value of absolute differences across replicates = %d', p)

%[h,p,ks2stat] = kstest2(peak_VL1,peak_VL2);
%sprintf('p-value differences in peak VLs= %d', p)
%[h,p,ks2stat] = kstest2(tpeak_VL1,tpeak_VL2);
%sprintf('p-value differences in timing of peak VLs= %d', p)
%[h,p,ks2stat] = kstest2(areaVLtopeak1,areaVLtopeak2);
%sprintf('p-value differences in area of CN values to peak value above detection threshold = %d', p)
%[h,p,ks2stat] = kstest2(areaVLafterpeak1,areaVLafterpeak2);
%sprintf('p-value differences in area of CN values after peak value above detection threshold = %d', p)
%[h,p,ks2stat] = kstest2(areaR1,areaR2);
%sprintf('p-value differences in area of R after peak value = %d', p)
%[h,p,ks2stat] = kstest2(BgeqBthres_dyn1,BgeqBthres_dyn2);
%sprintf('p-value differences in last time of B below the dynamic detection threshold = %d', p)
%[h,p,ks2stat] = kstest2(BgeqBthres_dyn_afterpeak1,BgeqBthres_dyn_afterpeak2);
%sprintf('p-value differences first time of B below the dynamic detection threshold after peak VL = %d', p)
%[h,p,ks2stat] = kstest2(areaBgeqBthres_dyn_afterpeak1,areaBgeqBthres_dyn_afterpeak2);
%sprintf('p-value differences area of between B and the dynamic detection threshold after peak VL = %d', p)
% [h,p,ks2stat] = kstest2(t_shift1,t_shift2);
% sprintf('p-value differences in incubation times = %d', p)
% [h,p,ks2stat] = kstest2(CN_sympt1,CN_sympt2);
% sprintf('p-value differences in CN value at symptom onset = %d', p)
% [h,p,ks2stat] = kstest2(T_hosp_check1,T_hosp_check2);
% sprintf('p-value differences in tsympt to thosp = %d', p)
% [h,p,ks2stat] = kstest2(pV1,pV2);
% sprintf('p-value differences in pV = %d', p)
% [h,p,ks2stat] = kstest2(log10(pB1),log10(pB2));
% sprintf('p-value differences in log10 pB = %d', p)
% if i_treat == 1
%     [h,p,ks2stat] = kstest2(alpha1,alpha2);
%     sprintf('p-value differences in alpha = %d', p)
% end
sprintf('percentage of pop1 with undetectable VL at 50 dpi = %d', 1-V50_1_bin/length(V50_1))
sprintf('percentage of pop2 with undetectable VL at 50 dpi = %d', 1-V50_2_bin/length(V50_2))
% [h,p,ks2stat] = kstest2(areaVL1,areaVL2);
% sprintf('p-value differences in area of CN values above detection threshold = %d', p)

if i_treat == 6
    pop11 = intersect(inds_wo_treat,inds_wo_hosp);
    pop12 = intersect(inds_wo_treat,inds_hosp);
    for h = 1:length(pop11)
        inds1a(h) = find(pop11(h) ==inds1);
    end
    for l = 1:length(pop12)
        inds1b(l) = find(pop12(l) ==inds1);
    end

    pop21 = intersect(inds_treat_hosp,inds_wo_dead);
    pop22 = intersect(inds_treat_hosp,inds_dead);
    for h = 1:length(pop21)
        inds2a(h) = find(pop21(h) ==inds2);
    end
    for l = 1:length(pop22)
        inds2b(l) = find(pop22(l) ==inds2);
    end
elseif i_treat == 7
    pop21 = intersect(inds_treat_hosp,inds_wo_dead);
    pop22 = intersect(inds_treat_hosp,inds_dead);
    for h = 1:length(pop21)
        inds2a(h) = find(pop21(h) ==inds2);
    end
    for l = 1:length(pop22)
        inds2b(l) = find(pop22(l) ==inds2);
    end
elseif i_treat == 8
    pop21 = intersect(inds_wo_treat,inds_hosp);
    pop22 = intersect(inds_treat_hosp,inds_wo_dead);
    pop23 = intersect(inds_treat_hosp,inds_dead);
    for k = 1:length(pop21)
        inds2a(k) = find(pop21(k) ==inds2);
    end
    for h = 1:length(pop22)
        inds2b(h) = find(pop22(h) ==inds2);
    end
    for l = 1:length(pop23)
        inds2c(l) = find(pop23(l) ==inds2);
    end
end

i1 = 1;
P1(i1,:) = t_detect1;
i1 = i1 + 1;
P1(i1,:) = t_shift1;
i1 = i1 + 1;
if ismember(i_treat, [1,6])
    P1(i1,:) = T_hosp_check1+t_shift1;
    i1 = i1 + 1;
end
P1(i1,:) = max_t(inds1);
i1 = i1 + 1;
if ismember(i_treat, [1,6])
    P1(i1,:) = T_hosp_check1;
    i1 = i1 + 1;
end
P1(i1,:) = max_t(inds1) - t_shift1;
i1 = i1 + 1;
if i_treat == 5
    P1(i1,:) = T_hosp_check1+t_shift1;
    i1 = i1 + 1;
    P1(i1,:) = T_hosp_check1;
    i1 = i1 + 1;
end

P1(i1,:) = CN_sympt1;
i1 = i1 + 1;
if ismember(i_treat, [1,6])
    P1(i1,:) = CN_hosp1;
    i1 = i1 + 1;
end
P1(i1,:) = max_V(inds1);
i1 = i1 + 1;

P1(i1,:) = areaVL1;
i1 = i1 + 1;
P1(i1,:) = areaVLtopeak1;
i1 = i1 + 1;
P1(i1,:) = areaVLafterpeak1;
i1 = i1 + 1;

P1(i1,:) = meanVL1;
i1 = i1 + 1;
P1(i1,:) = meanVLtopeak1;
i1 = i1 + 1;
P1(i1,:) = meanVLafterpeak1;
i1 = i1 + 1;

P1(i1,:) = pV1;
i1 = i1 + 1;
P1(i1,:) = log10(pB1);
i1 = i1 + 1;
if i_treat == 1
    P1(i1,:) = 1-alpha1;
    i1 = i1 + 1;
end

i2 = 1;
P2(i2,:) = t_detect2;
i2 = i2 + 1;
P2(i2,:) = t_shift2;
i2 = i2 + 1;
if ismember(i_treat, [1,6])
    P2(i2,:) = T_hosp_check2+t_shift2;
    i2 = i2 + 1;
end
P2(i2,:) = max_t(inds2);
i2 = i2 + 1;
if ismember(i_treat, [1,6])
    P2(i2,:) = T_hosp_check2;
    i2 = i2 + 1;
end
P2(i2,:) = max_t(inds2) - t_shift2;
i2 = i2 + 1;
if i_treat == 5
    P2(i2,:) = T_hosp_check2+t_shift2;
    i2 = i2 + 1;
    P2(i2,:) = T_hosp_check2;
    i2 = i2 + 1;
end

P2(i2,:) = CN_sympt2;
i2 = i2 + 1;
if ismember(i_treat, [1,6])
    P2(i2,:) = CN_hosp2;
    i2 = i2 + 1;
end
P2(i2,:) = max_V(inds2);
i2 = i2 + 1;

P2(i2,:) = areaVL2;
i2 = i2 + 1;
P2(i2,:) = areaVLtopeak2;
i2 = i2 + 1;
P2(i2,:) = areaVLafterpeak2;
i2 = i2 + 1;

P2(i2,:) = meanVL2;
i2 = i2 + 1;
P2(i2,:) = meanVLtopeak2;
i2 = i2 + 1;
P2(i2,:) = meanVLafterpeak2;
i2 = i2 + 1;

P2(i2,:) = pV2;
i2 = i2 + 1;
P2(i2,:) = log10(pB2);
i2 = i2 + 1;
if i_treat == 1
    P2(i2,:) = 1-alpha2;
    i2 = i2 + 1;
end

if i_treat == 1
    [c1,p1] = corr(T_hosp_check1',alpha1');
    sprintf('pval corr pop1: %d', p1)
    [c2,p2] = corr(T_hosp_check2',alpha2');
    sprintf('pval corr pop2: %d', p2)
end

if i_treat == 1
    inds_par = [{1:6},{7:9},{10:12},{13:15},16,17,18];
elseif i_treat == 3
    inds_par = [{1:4},{5:6},{7:9},{10:12},13,14];
elseif i_treat == 6
    inds_par = [{1:6},{7:9},{10:12},{13:15},16,17,18];
end

sD = 15;
C = [{[95,211,141]./255}, {[95,95,211]./255}, {[255,127,42]./255}, {[204,0,255]./255}, {[85,85,255]./255}, {[95,211,141]./255}, {[0,0,0]}, {[200,200,200]./255}];
M = {'o','^'};

FA = 0.5;
LW = 1; %linewidth

figure
t = tiledlayout(6,7,'TileSpacing','Compact','Padding','Compact');
n1 = normrnd(0,0.1,size(P1,2),1);
n2 = normrnd(0,0.1,size(P2,2),1);

for i = 1:length(inds_par)
    i_count = 1;  
        nexttile([6,1])
    for j = inds_par{i}
       
        percentile11 = prctile(P1(j,:),25);
        percentile12 = prctile(P1(j,:),75);

        plot([min(P1(j,:)),percentile11],[-((i_count-1)*2+1),-((i_count-1)*2+1)],'-', 'Color', C{7}, 'LineWidth', LW);
        hold on
        plot([percentile12,max(P1(j,:))],[-((i_count-1)*2+1),-((i_count-1)*2+1)],'-', 'Color', C{7}, 'LineWidth', LW);
        hold on
        rectangle('Position',[percentile11,-((i_count-1)*2+1)-0.2,percentile12-percentile11,0.4],'FaceColor',[C{i_treat}, 0.25],'EdgeColor','k')
        hold on
        plot([median(P1(j,:)),median(P1(j,:))],-[((i_count-1)*2+1)-0.2,((i_count-1)*2+1)+0.2],'-', 'Color', C{7}, 'LineWidth', LW)
        hold on


        % if ismember(i_treat,[1,6])
        %     alpha_bin = double(alpha(inds1)>0 & alpha(inds1)<0.9);
        %     alpha_bin(alpha_bin==0) = 2;
        %     k = alpha_bin;
        % 
        %     for j_data = 1:size(P1,2)
        %         scatter(P1(j,j_data),-((i_count-1)*2+1)+n1(j_data),M{k(j_data)},'MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{i_treat},'MarkerFaceAlpha',FA)
        %         hold on
        %     end
        % 
        %     plot([median(P1(j,find(k==1))),median(P1(j,find(k==1)))],-[((i_count-1)*2+1)-0.2,((i_count-1)*2+1)+0.2],'-', 'Color', C{8}, 'LineWidth', LW)
        %     hold on
        % else
             scatter(P1(j,:),-((i_count-1)*2+1)+n1(:),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{i_treat},'MarkerFaceAlpha',FA)
             hold on
        % end

        percentile21 = prctile(P2(j,:),25);
        percentile22 = prctile(P2(j,:),75);

        plot([min(P2(j,:)),percentile21],[-((i_count-1)*2+1.5),-((i_count-1)*2+1.5)],'-', 'Color', C{7}, 'LineWidth', LW);
        hold on
        plot([percentile22,max(P2(j,:))],[-((i_count-1)*2+1.5),-((i_count-1)*2+1.5)],'-', 'Color', C{7}, 'LineWidth', LW);
        hold on
        rectangle('Position',[percentile21,-((i_count-1)*2+1.5)-0.2,percentile22-percentile21,0.4],'FaceColor',[C{i_treat+1}, 0.25],'EdgeColor','k')
        hold on
        plot([median(P2(j,:)),median(P2(j,:))],-[((i_count-1)*2+1.5)-0.2,((i_count-1)*2+1.5)+0.2],'-', 'Color', C{7}, 'LineWidth', LW)
        hold on

        % if ismember(i_treat,[1,6])
        %     alpha_bin = double(alpha(inds2)>0 & alpha(inds2)<0.9);
        %     alpha_bin(alpha_bin==0) = 2;
        %     k = alpha_bin;
        % 
        %     for j_data = 1:size(P2,2)
        %         scatter(P2(j,j_data),-((i_count-1)*2+1.5)+n1(j_data),M{k(j_data)},'MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{i_treat+1},'MarkerFaceAlpha',FA)
        %         hold on
        %     end
        % 
        %     plot([median(P2(j,find(k==1))),median(P2(j,find(k==1)))],-[((i_count-1)*2+1.5)-0.2,((i_count-1)*2+1.5)+0.2],'-', 'Color', C{8}, 'LineWidth', LW)
        %     %scatter(1.1,median(P2(i,find(k==1))),'o','MarkerEdgeColor','none','SizeData',sD+5,'MarkerFaceColor',C{8},'MarkerFaceAlpha',1)
        %     hold on
        % 
        % else
             scatter(P2(j,:),-((i_count-1)*2+1.5)+n2(:),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{i_treat+1},'MarkerFaceAlpha',FA)
        % end

        plot([median(P2(j,:)),median(P2(j,:))],-[((i_count-1)*2+1.5)-0.2,((i_count-1)*2+1.5)+0.2],'-', 'Color', C{7}, 'LineWidth', LW)
        i_count = i_count + 1;
        % 
        % if i_treat == 6
        %     if ~ismember(i,[7,8])
        %         [h,p,ks2stat] = kstest2(P1(i,:),P2(i,:));
        %         sprintf('kstest2 of par_%d: %d', i, p)
        %         sprintf('mediantest of par_%d: %d', i, mediantest(P1(i,:),P2(i,:)))
        %         M_all(i) = mediantest(P1(i,:),P2(i,:));
        %     else
        %         [h,p,ks2stat] = kstest2(P1(i,inds1b),P2(i,:));
        %         sprintf('kstest2 of par_%d: %d', i, p)
        %         sprintf('mediantest of par_%d: %d', i, mediantest(P1(i,inds1b),P2(i,:)))
        %         M_all(i) = mediantest(P1(i,inds1b),P2(i,:));
        %     end
        %     K_all(i) = p;
        %     sprintf('medians of par_%d: pop1 %d vs pop2 %d', i, median(P1(i,:)), median(P2(i,:)))
        % else
            [h,p,ks2stat] = kstest2(P1(j,:),P2(j,:));
            sprintf('kstest2 of par_%d: %d', j, p)
            sprintf('mediantest of par_%d: %d', j, mediantest(P1(j,:),P2(j,:)))
            sprintf('medians of par_%d: pop1 %d vs pop2 %d', j, median(P1(j,:)), median(P2(j,:)))

            K_all(j) = p;
            M_all(j) = mediantest(P1(j,:),P2(j,:));
        % end

    end

    ylim([-12,0])
    if ismember(i,[2,3,4,5])
        xlim([0,Inf])
    end
    if i == 7
        xlim([0,1])
    end

    yticklabels({})
    box off

end

% for i = 1:size(P1,1)
%     nexttile
%     if i_treat == 1
%         if i == 9
%             nexttile
%             % nexttile
%         end
%     elseif i_treat == 3
%         if i == 7
%             nexttile
%             nexttile
%         end
%     elseif i_treat >= 7
%         if i == 8
%             nexttile
%             % nexttile
%         end
%     end
%     if i_treat == 1    
%         alpha_bin = double(alpha(inds1)>0 & alpha(inds1)<0.9);
%         alpha_bin(alpha_bin==0) = 2;
%         k = alpha_bin;
%         for j = 1:size(P1,2)
%             scatter(1+n1(j),P1(i,j),M{k(j)},'MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{i_treat},'MarkerFaceAlpha',FA)
%             hold on
%         end
%         plot([1-0.02,1+0.02],[median(P1(i,find(k==1))),median(P1(i,find(k==1)))],'-', 'Color', C{8}, 'LineWidth', LW)
%         hold on
%         plot([1-0.02,1+0.02],[median(P1(i,:)),median(P1(i,:))],'-', 'Color', C{7}, 'LineWidth', LW)
%         hold on
% 
%         alpha_bin = double(alpha(inds2)>0 & alpha(inds2)<0.9);
%         alpha_bin(alpha_bin==0) = 2;
%         k = alpha_bin;
%         for j = 1:size(P2,2)
%             scatter(1.1+n1(j),P2(i,j),M{k(j)},'MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{i_treat+1},'MarkerFaceAlpha',FA)
%             hold on
%         end
%         plot([1.1-0.02,1.1+0.02],[median(P2(i,find(k==1))),median(P2(i,find(k==1)))],'-', 'Color', C{8}, 'LineWidth', LW)
%         %scatter(1.1,median(P2(i,find(k==1))),'o','MarkerEdgeColor','none','SizeData',sD+5,'MarkerFaceColor',C{8},'MarkerFaceAlpha',1)
%         hold on
%         plot([1.1-0.02,1.1+0.02],[median(P2(i,:)),median(P2(i,:))],'-', 'Color', C{7}, 'LineWidth', LW)
%     elseif i_treat == 3
%         scatter(1+n1,P1(i,:),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{i_treat},'MarkerFaceAlpha',FA)
%         hold on
%         plot([1-0.02,1+0.02],[median(P1(i,:)),median(P1(i,:))],'-', 'Color', C{7}, 'LineWidth', LW)
%         hold on
%         scatter(1.1+n2,P2(i,:),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{i_treat+1},'MarkerFaceAlpha',FA)
%         hold on
%         plot([1.1-0.02,1.1+0.02],[median(P2(i,:)),median(P2(i,:))],'-', 'Color', C{7}, 'LineWidth', LW)
%     elseif i_treat == 5
%         scatter(1+n1,P1(i,:),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{i_treat},'MarkerFaceAlpha',FA)
%         hold on
%         plot([1-0.02,1+0.02],[median(P1(i,:)),median(P1(i,:))],'-', 'Color', C{7}, 'LineWidth', LW)
%         hold on
%         alpha_bin = double(alpha(inds2)>0 & alpha(inds2)<0.9);
%         alpha_bin(alpha_bin==0) = 2;
%         k = alpha_bin;
%         for j = 1:size(P2,2)
%             scatter(1.1+n2(j),P2(i,j),M{k(j)},'MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{i_treat+1},'MarkerFaceAlpha',FA)
%             hold on
%         end
%         plot([1.1-0.02,1.1+0.02],[median(P2(i,find(k==2))),median(P2(i,find(k==2)))],'-','Color',C{8},'LineWidth',LW)
%         hold on
%         scatter(1.1,median(P2(i,find(k==1))),'o','MarkerEdgeColor','none','SizeData',sD+5,'MarkerFaceColor',C{8},'MarkerFaceAlpha',1)
%     elseif i_treat == 6
%         if ~ismember(i,[7,8])
%             scatter(1+n1(inds1a),P1(i,inds1a),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{3},'MarkerFaceAlpha',FA)
%             hold on
%         end
%         scatter(1+n1(inds1b),P1(i,inds1b),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{4},'MarkerFaceAlpha',FA)
%         hold on
%         if ~ismember(i,[7,8])
%             plot([1-0.02,1+0.02],[median(P1(i,:)),median(P1(i,:))],'-', 'Color', C{7}, 'LineWidth', LW)
%             hold on
%         else
%             plot([1-0.02,1+0.02],[median(P1(i,inds1b)),median(P1(i,inds1b))],'-', 'Color', C{7}, 'LineWidth', LW)
%             hold on
%         end
%         scatter(1.1+n2(inds2a),P2(i,inds2a),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{1},'MarkerFaceAlpha',FA)
%         hold on
%         scatter(1.1+n2(inds2b),P2(i,inds2b),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{2},'MarkerFaceAlpha',FA)
%         hold on
%         plot([1.1-0.02,1.1+0.02],[median(P2(i,:)),median(P2(i,:))],'-', 'Color', C{7}, 'LineWidth', LW)
%     elseif i_treat == 7
%         scatter(1+n1,P1(i,:),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{3},'MarkerFaceAlpha',FA)
%         hold on
%         plot([1-0.02,1+0.02],[median(P1(i,:)),median(P1(i,:))],'-', 'Color', C{7}, 'LineWidth', LW)
%         hold on
%         scatter(1.1+n2(inds2a),P2(i,inds2a),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{1},'MarkerFaceAlpha',FA)
%         hold on
%         scatter(1.1+n2(inds2b),P2(i,inds2b),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{2},'MarkerFaceAlpha',FA)
%         hold on
%         plot([1.1-0.02,1.1+0.02],[median(P2(i,:)),median(P2(i,:))],'-', 'Color', C{7}, 'LineWidth', LW)
%     else
%         scatter(1+n1,P1(i,:),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{3},'MarkerFaceAlpha',FA)
%         hold on
%         plot([1-0.02,1+0.02],[median(P1(i,:)),median(P1(i,:))],'-', 'Color', C{7}, 'LineWidth', LW)
%         hold on
%         scatter(1.1+n2(inds2a),P2(i,inds2a),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{4},'MarkerFaceAlpha',FA)
%         hold on
%         scatter(1.1+n2(inds2b),P2(i,inds2b),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{1},'MarkerFaceAlpha',FA)
%         hold on
%         scatter(1.1+n2(inds2c),P2(i,inds2c),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{2},'MarkerFaceAlpha',FA)
%         hold on
%         plot([1.1-0.02,1.1+0.02],[median(P2(i,:)),median(P2(i,:))],'-', 'Color', C{7}, 'LineWidth', LW)
%     end
%     if i_treat == 1
%         if ~ismember(i,[6,16])
%             ylim([0,Inf])
%         end
%     elseif i_treat == 3
%         if ~ismember(i,[6,14])
%             ylim([0,Inf])
%         end
%     elseif i_treat == 5
%         if ~ismember(i,[6,16])
%             ylim([0,Inf])
%         end
%     elseif i_treat == 6
%         if ~ismember(i,[6,16])
%             ylim([0,Inf])
%         end
%     else
%         if ~ismember(i,[6,15])
%             ylim([0,Inf])
%         end
%     end
%     hold on
% 
%     xticks([1,1.1])
%     xticklabels({})
%     xlim([1-max(abs(n1))-0.01,1.1+max(abs(n2))+0.01])
% 
%     box off

    % if i_treat == 6
    %     if ~ismember(i,[7,8])
    %         [h,p,ks2stat] = kstest2(P1(i,:),P2(i,:));
    %         sprintf('kstest2 of par_%d: %d', i, p)
    %         sprintf('mediantest of par_%d: %d', i, mediantest(P1(i,:),P2(i,:)))
    %         M_all(i) = mediantest(P1(i,:),P2(i,:));
    %     else
    %         [h,p,ks2stat] = kstest2(P1(i,inds1b),P2(i,:));
    %         sprintf('kstest2 of par_%d: %d', i, p)
    %         sprintf('mediantest of par_%d: %d', i, mediantest(P1(i,inds1b),P2(i,:)))
    %         M_all(i) = mediantest(P1(i,inds1b),P2(i,:));
    %     end
    %     K_all(i) = p;
    %     sprintf('medians of par_%d: pop1 %d vs pop2 %d', i, median(P1(i,:)), median(P2(i,:)))
    % else
    %     [h,p,ks2stat] = kstest2(P1(i,:),P2(i,:));
    %     sprintf('kstest2 of par_%d: %d', i, p)
    %     sprintf('mediantest of par_%d: %d', i, mediantest(P1(i,:),P2(i,:)))
    %     sprintf('medians of par_%d: pop1 %d vs pop2 %d', i, median(P1(i,:)), median(P2(i,:)))
    % 
    %     K_all(i) = p;
    %     M_all(i) = mediantest(P1(i,:),P2(i,:));
    % end
    % 
% end

clear figure_property;
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
figure_property.Width= '10'; % Figure width on canvas
figure_property.Height= '2.5'; % Figure height on canvas
figure_property.Units= 'inches';
figure_property.Color= 'rgb';
figure_property.Background= 'w';
figure_property.FixedfontSize= '8';
figure_property.ScaledfontSize= 'auto';
figure_property.FontMode= 'fixed';
figure_property.FontSizeMin= '8';
figure_property.FixedLineWidth= '1';
figure_property.ScaledLineWidth= 'auto';
figure_property.LineMode= 'none';
figure_property.LineWidthMin= '1';
figure_property.FontName= 'Arial';% Might want to change this to something that is available
figure_property.FontWeight= 'auto';
figure_property.FontAngle= 'auto';
figure_property.FontEncoding= 'latin1';
figure_property.PSLevel= '3';
figure_property.Renderer= 'painters';
figure_property.Resolution= '600';
figure_property.LineStyleMap= 'none';
figure_property.ApplyStyle= '0';
figure_property.Bounds= 'tight';
figure_property.LockAxes= 'off';
figure_property.LockAxesTicks= 'off';
figure_property.ShowUI= 'off';
figure_property.SeparateText= 'off';
chosen_figure=gcf;
set(chosen_figure,'PaperUnits','inches');
set(chosen_figure,'PaperPositionMode','auto');
set(chosen_figure,'PaperSize',[8.27,11.69]); % Canvas Size
set(chosen_figure,'Units','inches');
% if i_treat == 1
%     hgexport(gcf,'Par_analysis_treat.pdf',figure_property); %Set desired file name
% elseif i_treat == 3
%     hgexport(gcf,'Par_analysis_hosp.pdf',figure_property); %Set desired file name
% elseif i_treat == 5
%     hgexport(gcf,'Par_analysis_both.pdf',figure_property); %Set desired file name
% elseif i_treat == 6
%     hgexport(gcf,'Par_analysis_all.pdf',figure_property); %Set desired file name
% elseif i_treat == 7
%     hgexport(gcf,'Par_analysis_extremes.pdf',figure_property); %Set desired file name
% else
%      hgexport(gcf,'Par_analysis_severe.pdf',figure_property); %Set desired file name
% end

% S(inds1a) = 1;
% S(inds1b) = 2;
% T(inds2a) = 3;
% T(inds2b) = 4;
% species = [S';T'];
% 
% X = [P1([1:6,9:size(P1,1)],:)';P2([1:6,9:size(P2,1)],:)'];
% 
% alg = {'euclidean' 'seuclidean' 'fasteuclidean' 'fastseuclidean' 'cityblock'...
%     'chebychev' 'minkowski' 'cosine' 'correlation' 'spearman'...
%     'hamming' 'jaccard'};
% p = 3;
% p2 = 38;
% % 
% figure
% % gscatter(X(:,10),X(:,11),species',[C{3};C{4};C{1};C{2}],'.',10,'off')
% gscatter(X(1:size(P1,2),13),X(1:size(P1,2),14),species(1:size(P1,2)),[C{3};C{4};C{1};C{2}],'.',10,'off')
% hold on
% alpha2 = alpha(inds2);
% gscatter(X(size(P1,2)+find(1-alpha2 > 0.1 & 1-alpha2 < 1),13),X(size(P1,2)+find(1-alpha2 > 0.1 & 1-alpha2 < 1),14),species(size(P1,2)+find(1-alpha2 > 0.1 & 1-alpha2 < 1)),[C{1};C{2}],'.',10,'off')
% hold on
% gscatter(X(size(P1,2)+find(1-alpha2 < 0.1),13),X(size(P1,2)+find(1-alpha2 < 0.1),14),species(size(P1,2)+find(1-alpha2 < 0.1)),[C{1};C{2}],'^',5,'off')
% box off
% 
% % 
% % clear figure_property;
% % figure_property.units = 'inches';
% % figure_property.format = 'pdf';
% % figure_property.Preview= 'none';
% % figure_property.Width= '2'; % Figure width on canvas
% % figure_property.Height= '2'; % Figure height on canvas
% % figure_property.Units= 'inches';
% % figure_property.Color= 'rgb';
% % figure_property.Background= 'w';
% % figure_property.FixedfontSize= '8';
% % figure_property.ScaledfontSize= 'auto';
% % figure_property.FontMode= 'fixed';
% % figure_property.FontSizeMin= '8';
% % figure_property.FixedLineWidth= '1';
% % figure_property.ScaledLineWidth= 'auto';
% % figure_property.LineMode= 'none';
% % figure_property.LineWidthMin= '1';
% % figure_property.FontName= 'Arial';% Might want to change this to something that is available
% % figure_property.FontWeight= 'auto';
% % figure_property.FontAngle= 'auto';
% % figure_property.FontEncoding= 'latin1';
% % figure_property.PSLevel= '3';
% % figure_property.Renderer= 'painters';
% % figure_property.Resolution= '600';
% % figure_property.LineStyleMap= 'none';
% % figure_property.ApplyStyle= '0';
% % figure_property.Bounds= 'tight';
% % figure_property.LockAxes= 'off';
% % figure_property.LockAxesTicks= 'off';
% % figure_property.ShowUI= 'off';
% % figure_property.SeparateText= 'off';
% % chosen_figure=gcf;
% % set(chosen_figure,'PaperUnits','inches');
% % set(chosen_figure,'PaperPositionMode','auto');
% % set(chosen_figure,'PaperSize',[8.27,11.69]); % Canvas Size
% % set(chosen_figure,'Units','inches');
% % % hgexport(gcf,'Par_analysis_par_comp.pdf',figure_property); %Set desired file name
% % 
% 
% %Y = tsne(X,Perplexity = p);
% for p = 3
%     figure
%     itile = 1;
%     for p2 = 41
%         rng default % for reproducibility
%         % subplot(3,4,itile)
%         Y = tsne(X,'Algorithm','exact','Distance',alg{p},Perplexity=p2);
%         % subplot(5,10,p)
%         gscatter(Y(1:size(P1,2),1),Y(1:size(P1,2),2),species(1:size(P1,2)),[C{3};C{4};C{1};C{2}],'.',10,'off')
%         hold on
%         alpha2 = alpha(inds2);
%         gscatter(Y(size(P1,2)+find(1-alpha2 > 0.1 & 1-alpha2 < 1),1),Y(size(P1,2)+find(1-alpha2 > 0.1 & 1-alpha2 < 1),2),species(size(P1,2)+find(1-alpha2 > 0.1 & 1-alpha2 < 1)),[C{1};C{2}],'.',10,'off')
%         hold on
%         gscatter(Y(size(P1,2)+find(1-alpha2 < 0.1),1),Y(size(P1,2)+find(1-alpha2 < 0.1),2),species(size(P1,2)+find(1-alpha2 < 0.1)),[C{1};C{2}],'^',5,'off')
%         box off
%         itile = itile + 1;
%     end
% end
% 
% % clear figure_property;
% % figure_property.units = 'inches';
% % figure_property.format = 'pdf';
% % figure_property.Preview= 'none';
% % figure_property.Width= '2'; % Figure width on canvas
% % figure_property.Height= '2'; % Figure height on canvas
% % figure_property.Units= 'inches';
% % figure_property.Color= 'rgb';
% % figure_property.Background= 'w';
% % figure_property.FixedfontSize= '8';
% % figure_property.ScaledfontSize= 'auto';
% % figure_property.FontMode= 'fixed';
% % figure_property.FontSizeMin= '8';
% % figure_property.FixedLineWidth= '1';
% % figure_property.ScaledLineWidth= 'auto';
% % figure_property.LineMode= 'none';
% % figure_property.LineWidthMin= '1';
% % figure_property.FontName= 'Arial';% Might want to change this to something that is available
% % figure_property.FontWeight= 'auto';
% % figure_property.FontAngle= 'auto';
% % figure_property.FontEncoding= 'latin1';
% % figure_property.PSLevel= '3';
% % figure_property.Renderer= 'painters';
% % figure_property.Resolution= '600';
% % figure_property.LineStyleMap= 'none';
% % figure_property.ApplyStyle= '0';
% % figure_property.Bounds= 'tight';
% % figure_property.LockAxes= 'off';
% % figure_property.LockAxesTicks= 'off';
% % figure_property.ShowUI= 'off';
% % figure_property.SeparateText= 'off';
% % chosen_figure=gcf;
% % set(chosen_figure,'PaperUnits','inches');
% % set(chosen_figure,'PaperPositionMode','auto');
% % set(chosen_figure,'PaperSize',[8.27,11.69]); % Canvas Size
% % set(chosen_figure,'Units','inches');
% % hgexport(gcf,'Par_analysis_tsne.pdf',figure_property); %Set desired file name
% % 
