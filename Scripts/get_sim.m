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
    Sex(i_patient) = data.('Sex')(i_patient);
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
        Sex_set(icount) = Sex(i_patient);
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
            Sex_check(icount) = Sex_set(i_t);
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

load('sol_treat_3_t0')

i_count = 1;
for ID_opt = inds_treat_hosp(1:5)
    j_count = 1;
    figure
    for i_alpha = 0:0.1:1

        % for i_alpha = 0
        %     for ID_opt = inds_treat_hosp(1)

        pB(ID_opt) = sol{ID_opt}.P_t(1); %immune response on infection rate
        pV(ID_opt) = sol{ID_opt}.P_t(2);
        % dB = 10^par_opt(ind(1),3);
        t_shift(ID_opt) = floor(sol{ID_opt}.P_t(3));
        % alpha(ID_opt) = sol{ID_opt}.P_t(4);
        sigma(ID_opt) = sol{ID_opt}.P_t(5);

        T_check_new = T_check{ID_opt} + OnsetSympt_check(ID_opt)+t_shift(ID_opt);
        T_check_new = T_check_new(find(T_check_new<=50));

        T_hosp_check_new = T_hosp_check(ID_opt) + OnsetSympt_check(ID_opt)+t_shift(ID_opt);

        S0 = 8*10^7; %total number of epithelial cells in nose at t=0, Ke et al., 2022
        dN = 1/11; %death rate of all target cells, Tomasetti et al., 2017
        pN = S0*dN; %production of new epithelial cells
        b0 = 4.92*10^(-9); %infectivity rate, Ke et al., 2022
        dI = 2.45; %death of infected cells, Ke et al., 2022
        dV = 10; %deactivation virus, Ke et al., 2022
        dB = 0;

        %plot the results for short- and long-term dynamics
        tspan = 0:0.1:50; %long term

        %determine the initial values and B_thres for simulation
        y0 = [S0, 1, 0, 0];
        B_thres(ID_opt) = 1-dI*dV/(b0*S0*(pV(ID_opt)-dI));
        options = odeset('NonNegative',[1,2,3,4]); %specify non-negative values
        [t,y] = ode45(@(t,y) odefcn_single_infection_S_R0(t,y,b0,dI,pV(ID_opt),dV,pN,dN,pB(ID_opt),dB,B_thres(ID_opt),i_alpha,T_hosp_check(ID_opt)+t_shift(ID_opt),S0), tspan, y0, options);

        y_short = y(:,3);
        
        y_short_after_treat = y_short(find(t == T_hosp_check_new):end);
        ind_last_after_treat = find(y_short_after_treat<1,1,'last');

        % y_short(find(t == T_hosp_check_new)+ind_last_after_treat:501) = 0;
        %if values too small, fix at 1 (numerical problems)
        y_short(y_short<1)=1;

        y_short = -(log10(y_short)-11.35)/(-0.25);
        % 
        % plot(t,y_short,'-')
        % hold on
        plot(t,y(:,4),'-')
        hold on
        yline( B_thres(ID_opt))

        for i_tpt = 1:length(T_check_new)
            r = normrnd(0,sigma(ID_opt),1,3);
            D_sim{i_count}{j_count}(:,i_tpt) = y_short(find(t==T_check_new(i_tpt)))+r;
        end
        j_count = j_count +1;
    end
    i_count = i_count + 1;
end

% save('D_sim','D_sim')



