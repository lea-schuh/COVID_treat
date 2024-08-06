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
                val(i_dates) = 42;
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

inds_treat = find(Treat_check == 1);
inds_wo_treat = find(Treat_check == 0);

i_ind = 1;

t = tiledlayout(4,5,'TileSpacing','Compact','Padding','Compact');
%fit over every patient seperately
% for ID_opt = 1:length(T_check)
for ID_opt = inds_treat
% for ID_opt = inds_wo_treat(1:20)

    nexttile

    clearvars -except ID_check i_ind data ID_opt S0 dN pN b0 dI dV dB T_check CT_check Treat_check T_hosp_check Death_check sol inds_treat inds_wo_treat OnsetSympt_check

    ID_opt 

    %number of single optimzation runs
    n_start = 10;

    %determine lower (lb) and upper (ub) parameter boundaries
    lb = [-10,2,-4,log10(2+OnsetSympt_check(ID_opt)),-2]; %pB, pV, dB, t_shift, sigma
    ub = [-4,3,0,log10(20),1];

    %use latin hypercube sampling to sample n_start initial values

    par0_all = lhsdesign_modified(n_start,lb,ub); %latinhypercube sampled startign parameters

    %optimize for each n_start
    for i = 1:n_start
        [par_opt(i,:),fval(i)] = fmincon(@(par)loglikelihood_single_Greek(-CT_check{ID_opt},T_check{ID_opt},par),par0_all(i,:),[],[],[],[],lb,ub);
        % BIC(i) = size(par0_all,2)*log(size(data_ID{ID_opt},1))+2*fval(i);
    end

    %sort all n_starts according to best runs (highest fval)
    [fval_sorted,ind] = sort(fval);
    sol{ID_opt}.fval_sorted = fval_sorted;
    % sol{ID_opt}.BIC_sorted= BIC(ind);

    %set free parameters to best estimates from optimzation
    pB = 10^par_opt(ind(1),1); %immune response on infection rate
    pV = 10^par_opt(ind(1),2);
    dB = 10^par_opt(ind(1),3);
    t_shift = ceil(10^par_opt(ind(1),4));

    sol{ID_opt}.P = 10.^par_opt(ind(1),:);
    sol{ID_opt}.Pall = 10.^par_opt(ind,:);

    %save sol structure
    % save(sprintf('sol_%d_wo_dB',i_data),'sol');

    %plot the results for short- and long-term dynamics
    tspan = 0:0.1:max(cell2mat(T_check))+t_shift; %long term

    %determine the initial values and B_thres for simulation
    y0 = [S0, 1, 0, 0];
    B_thres = 1-dI*dV/(b0*S0*(pV-dI));
    options = odeset('NonNegative',[1,2,3,4]); %specify non-negative values
    [t,y] = ode45(@(t,y) odefcn_single_infection_S_R0(t,y,b0,dI,pV,dV,pN,dN,pB,dB,B_thres), tspan, y0,options);

    % figure
    % subplot(size(y,2)+1,1,1)
    % plot(t,y(:,1),'Color','k','Linewidth',2) %S
    % 
    % subplot(size(y,2)+1,1,2)
    % plot(t,y(:,2),'Color','k','Linewidth',2) %I
    % 
    % subplot(size(y,2)+1,1,3)
    % plot(t,y(:,3),'Color','k','Linewidth',2) %V
    % 
    % subplot(size(y,2)+1,1,4)
    if Treat_check(ID_opt) == 0
        plot(T_check{ID_opt}+t_shift, -CT_check{ID_opt}, 'o', 'Color','k','Markersize',2,'MarkerFaceColor','k') %data
    else
        plot(T_check{ID_opt}+t_shift, -CT_check{ID_opt}, 'o', 'Color','r','Markersize',2,'MarkerFaceColor','r') %data
    end
    hold on
    y_short = y(:,3);
    %if values too small, fix at 1 (numerical problems)
    y_short(y_short<1)=1;
    %plot(t,-(log(y(:,3)/(1.441*10^14)))/(-0.685),'Color','k','Linewidth',2) %Ke 2021
    plot(t,-(log10(y_short)-11.35)/(-0.25),'Color','k','Linewidth',1) %Ke 2022 nasal
    hold on
    plot([T_hosp_check(ID_opt)+t_shift,T_hosp_check(ID_opt)+t_shift],[-50,-10],':','Color','k','Linewidth',1)
    hold on
    plot([0,max(t)],[-42,-42],'--','Color','k','Linewidth',1)
    hold on
    plot(T_check{ID_opt}(1)+t_shift-OnsetSympt_check(ID_opt), -50, '*', 'Color', 'k','Markersize',5) %data
    if Death_check(ID_opt) == 0
        title(sprintf('#%s',ID_check{ID_opt}),'FontWeight','normal')
    else
        title(sprintf('#%s*',ID_check{ID_opt}),'FontWeight','normal')
    end
    ylim([-50,-10])
    xlim([0,50])

    % subplot(size(y,2)+1,1,5)
    % plot(t,y(:,4), '-', 'Color', 'k','Linewidth',2) %B
    % hold on
    % plot(t,1-(dI*dV/(b0*S0*(pV-dI)))*ones(length(t),1),'--','Color', 'k','Linewidth',2) %B_thres
    % ylim([0,1])
    % hold on

    % set(gcf, 'Position',  [100, 100, 250, 600])
    if ismember(i_ind,[1,6,11,16])
        yticks([-50,-30,-10]);
    else
        yticks([-50,-30,-10]);
        yticklabels({});
    end

    if i_ind == 11
       ylabel('                                     Ct values')
    end

    if i_ind == 18
       xlabel('Days since infection')
    end

    if ismember(i_ind,[16,17,18,19,20])
        xticks([0,20,40]);
    else
        xticks([0,20,40]);
        xticklabels({});
    end

    box off
    
    i_ind = i_ind+1;
end

clear figure_property;
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
figure_property.Width= '6'; % Figure width on canvas
figure_property.Height= '3'; % Figure height on canvas
figure_property.Units= 'inches';
figure_property.Color= 'rgb';
figure_property.Background= 'w';
figure_property.FixedfontSize= '8';
figure_property.ScaledfontSize= 'auto';
figure_property.FontMode= 'scaled';
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
set(chosen_figure,'PaperSize',[str2num(figure_property.Width) str2num(figure_property.Height)]); % Canvas Size
set(chosen_figure,'Units','inches');
% hgexport(gcf,'Figure_Fits_wo_treat.pdf',figure_property); %Set desired file name
hgexport(gcf,'Figure_Fits_treat.pdf',figure_property); %Set desired file name
