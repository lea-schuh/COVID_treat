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
    clear month day time_new time val month_hosp day_hosp month_data1 day_data1
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
        else
            if strcmp(Val_set{i_set}(i_dates),'NEG (-)')
                % val(i_dates) = 42;
                val(i_dates) = 35;
            else
                val(i_dates) = str2double(extractBetween(Val_set{i_set}(i_dates),'/','/'));
            end
        end

    end

    T{i_set} = time;
    CT{i_set} = val;

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

% for i_checked = 1:length(T_check)
% for i_checked = 1
%     figure
%     if Treat_check(i_checked) == 0
%         plot(T_check{i_checked},-CT_check{i_checked},'o', 'Color', 'k','Markersize',2,'MarkerFaceColor','k')
%     else
%         plot(T_check{i_checked},-CT_check{i_checked},'o', 'Color', 'r','Markersize',2,'MarkerFaceColor','r')
%     end
%     set(gcf, 'Position',  [100, 100, 250, 250])
% end

%determine fixed parameter values for i_data = 1 (nasal) and i_data = 2
%(saliva) 
%nasal
S0 = 8*10^7; %total number of epithelial cells in nose at t=0, Ke et al., 2022
dN = 1/11; %death rate of all target cells, Tomasetti et al., 2017
pN = S0*dN; %production of new epithelial cells
b0 = 4.92*10^(-9); %infectivity rate, Ke et al., 2022
dI = 2.45; %death of infected cells, Ke et al., 2022
dV = 10; %deactivation virus, Ke et al., 2022
dB = 0;

inds_treat = find(Treat_check == 1);
inds_wo_treat = find(Treat_check == 0);

i_ind = 1;

% load('sol_wo_treat')
load('sol_treat')
treat = 1;

% 
% w = 10;
% h = 7;

w = 5;
h = 4;


figure
t = tiledlayout(h,w,'TileSpacing','Compact','Padding','Compact');

%fit over every patient seperately
% for ID_opt = 1:length(T_check)
% for ID_opt = inds_treat([2,7,8,9,12,13,16,17,18])
for ID_opt = inds_treat
% for ID_opt = inds_wo_treat

    nexttile

    clearvars -except ID_check i_ind data ID_opt S0 dN pN b0 dI dV dB T_check...
        CT_check Treat_check T_hosp_check Death_check sol inds_treat inds_wo_treat...
        OnsetSympt_check h w treat

    ID_opt 


    if treat == 0
    %set free parameters to best estimates from optimzation
    pB = sol{ID_opt}.P(1); %immune response on infection rate
    pV = sol{ID_opt}.P(2);
    t_shift = floor(sol{ID_opt}.P(3));

    %plot the results for short- and long-term dynamics
    tspan = 0:0.1:max(cell2mat(T_check))+t_shift; %long term

    %determine the initial values and B_thres for simulation
    y0 = [S0, 1, 0, 0];
    B_thres = 1-dI*dV/(b0*S0*(pV-dI));
    options = odeset('NonNegative',[1,2,3,4]); %specify non-negative values
    [t,y] = ode45(@(t,y) odefcn_single_infection_S_R0(t,y,b0,dI,pV,dV,pN,dN,pB,dB,B_thres,1,0), tspan, y0,options);

    y_short = y(:,3);
    %if values too small, fix at 1 (numerical problems)
    y_short(y_short<1)=1;
    %plot(t,-(log(y(:,3)/(1.441*10^14)))/(-0.685),'Color','k','Linewidth',2) %Ke 2021
    plot(t,-(log10(y_short)-11.35)/(-0.25),'Color','k','Linewidth',1) %Ke 2022 nasal
    hold on
 
    else

    %set free parameters to best estimates from optimzation
    pB_t = sol{ID_opt}.P_t(1); %immune response on infection rate
    pV_t = sol{ID_opt}.P_t(2);
    t_shift = floor(sol{ID_opt}.P(3));
    alpha_t = sol{ID_opt}.P_t(3);
    treat_t = sol{ID_opt}.P_t(4);

    %plot the results for short- and long-term dynamics
    tspan = 0:0.1:max(cell2mat(T_check))+t_shift; %long term

    %determine the initial values and B_thres for simulation
    y0_t = [S0, 1, 0, 0];
    B_thres_t = 1-dI*dV/(b0*S0*(pV_t-dI));
    options = odeset('NonNegative',[1,2,3,4]); %specify non-negative values
    [t,y] = ode45(@(t,y) odefcn_single_infection_S_R0(t,y,b0,dI,pV_t,dV,pN,dN,pB_t,dB,B_thres_t,alpha_t,treat_t), tspan, y0_t, options);

    y_short_t = y(:,3);
    %if values too small, fix at 1 (numerical problems)
    y_short_t(y_short_t<1)=1;
    %plot(t,-(log(y(:,3)/(1.441*10^14)))/(-0.685),'Color','k','Linewidth',2) %Ke 2021
    plot(t,-(log10(y_short_t)-11.35)/(-0.25),'Color','r','Linewidth',1) %Ke 2022 nasal
    hold on

    end

    if Treat_check(ID_opt) == 0
        plot(T_check{ID_opt}+t_shift, -CT_check{ID_opt}, 'o', 'Color','k','Markersize',4,'MarkerFaceColor','k') %data
    else
        plot(T_check{ID_opt}+t_shift, -CT_check{ID_opt}, 'o', 'Color','r','Markersize',4,'MarkerFaceColor','r') %data
    end
    hold on

    plot([T_hosp_check(ID_opt)+t_shift,T_hosp_check(ID_opt)+t_shift],[-50,-10],':','Color','k','Linewidth',1)
    hold on
    %plot([0,max(t)],[-42,-42],'--','Color','k','Linewidth',1)
    plot([0,max(t)],[-35,-35],'--','Color','k','Linewidth',1)
    hold on
    plot(T_check{ID_opt}(1)+t_shift-OnsetSympt_check(ID_opt), -50, '*', 'Color', 'k','Markersize',5) %data
    if Death_check(ID_opt) == 0
        title(sprintf('#%s',ID_check{ID_opt}),'FontWeight','normal')
    else
        title(sprintf('#%s*',ID_check{ID_opt}),'FontWeight','normal')
    end

    ylim([-50,-10])
    xlim([0,50])

    % set(gcf, 'Position',  [100, 100, 250, 600])
    if rem(i_ind,w)==1
        yticks([-50,-30,-10]);
    else
        yticks([-50,-30,-10]);
        yticklabels({});
    end

    if i_ind>((h-1)*w)
        xticks([0,20,40]);
    else
        xticks([0,20,40]);
        xticklabels({});
    end

    box off
    
    i_ind = i_ind+1;
end

% hgexport(gcf,'Figure_Fits_wo_treat.pdf',figure_property); %Set desired file name
% hgexport(gcf,'Figure_Fits_treat.pdf',figure_property); %Set desired file name

