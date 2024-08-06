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
                val(i_dates) = 35; %according to Ioanna
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

load('sol_wo_treat_3_t0')
%
inds_treat = find(Treat_check == 1);
inds_wo_treat = find(Treat_check == 0);

inds_hosp = find(~isnan(T_hosp_check));
inds_wo_hosp = find(isnan(T_hosp_check));

inds_dead = find(Death_check == 1);
inds_wo_dead = find(Death_check == 0);

icount = 1;
jcount = 1;
kcount = 1;
lcount = 1;
mcount = 1;

for k = inds_wo_treat
    %hospitalization
    if ismember(k,inds_hosp)
        P_hosp(icount,:) = sol{k}.P;
        icount = icount +1;
    else
        P_wo_hosp(jcount,:) = sol{k}.P;
        jcount = jcount +1;
    end

    %death
    if ismember(k,inds_dead)
        P_dead(kcount,:) = sol{k}.P;
        kcount = kcount +1;
    else
        P_wo_dead(lcount,:) = sol{k}.P;
        lcount = lcount +1;
    end

    mcount = mcount +1;

end

P_hosp(:,3) = floor(P_hosp(:,3));
P_wo_hosp(:,3) = floor(P_wo_hosp(:,3));
P_dead(:,3) = floor(P_dead(:,3));
P_wo_dead(:,3) = floor(P_wo_dead(:,3));

sD = 20;
C = [{[95,211,141]./255}, {[255,0,102]./255}, {[255,127,42]./255}, {[204,0,255]./255}, {[200,200,200]./255}, {[0,0,0]}];

FA = 0.5;
LW = 1; %linewidth

figure
t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
n1 = normrnd(0,0.01,length(P_wo_hosp),1);
n2 = normrnd(0,0.01,length(P_hosp),1);
for i = 1:3
    subplot(1,3,i)
    if i == 1
        scatter(1+n1,log10(P_wo_hosp(:,i)),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{3},'MarkerFaceAlpha',FA)
        hold on
        plot([1-0.02,1+0.02],[median(log10(P_wo_hosp(:,i))),median(log10(P_wo_hosp(:,i)))],'-', 'Color', C{6}, 'LineWidth', LW)
        hold on
        scatter(1.1+n2,log10(P_hosp(:,i)),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{4},'MarkerFaceAlpha',FA)
        hold on
        plot([1.1-0.02,1.1+0.02],[median(log10(P_hosp(:,i))),median(log10(P_hosp(:,i)))],'-', 'Color', C{6}, 'LineWidth', LW)
    else
        scatter(1+n1,P_wo_hosp(:,i),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{3},'MarkerFaceAlpha',FA)
        hold on
        plot([1-0.02,1+0.02],[median(P_wo_hosp(:,i)),median(P_wo_hosp(:,i))],'-', 'Color', C{6}, 'LineWidth', LW)
        hold on
        scatter(1.1+n2,P_hosp(:,i),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{4},'MarkerFaceAlpha',FA)
        hold on
        plot([1.1-0.02,1.1+0.02],[median(P_hosp(:,i)),median(P_hosp(:,i))],'-', 'Color', C{6}, 'LineWidth', LW)
        ylim([0,Inf])
    end
    hold on

    xticks([1,1.1])
    xticklabels({})
    xlim([1-max(abs(n1))-0.01,1.1+max(abs(n2))+0.01])

    box off

    if i == 1
        sprintf('median test of par_%d: %d', i, mediantest(log10(P_wo_hosp(:,i)),log10(P_hosp(:,i))))
        sprintf('kstest2 of par_%d: %d', i, kstest2(log10(P_wo_hosp(:,i)),log10(P_hosp(:,i))))
        sprintf('medians of par_%d: hosp %d vs non-hosp %d', i, median(log10(P_wo_hosp(:,i))), median(log10(P_hosp(:,i))))
    else
        sprintf('median test of par_%d: %d', i, mediantest(P_wo_hosp(:,i),P_hosp(:,i)))
        sprintf('kstest2 of par_%d: %d', i, kstest2(P_wo_hosp(:,i),P_hosp(:,i)))
        sprintf('medians of par_%d: hosp %d vs non-hosp %d', i, median(P_wo_hosp(:,i)), median(P_hosp(:,i)))
    end
end


clear figure_property;
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
figure_property.Width= '3'; % Figure width on canvas
figure_property.Height= '1.5'; % Figure height on canvas
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
hgexport(gcf,'Par_analysis_hosp.pdf',figure_property); %Set desired file name


% subplot(1,5,4)
% scatter(ones(size(Tsymp_hosp,1),1),P_hosp(:,3)-Tsymp_hosp')
% hold on
% scatter(2*ones(size(Tsymp_wo_hosp,1),1),P_wo_hosp(:,3)-Tsymp_wo_hosp')
%
% mediantest(P_hosp(:,3)-Tsymp_hosp',P_wo_hosp(:,3)-Tsymp_wo_hosp')
%
% subplot(1,5,5)
% scatter(ones(size(Tsymp_hosp,1),1),Tsymp_hosp')
% hold on
% scatter(2*ones(size(Tsymp_wo_hosp,1),1),Tsymp_wo_hosp')

% mediantest(Tsymp_hosp',Tsymp_wo_hosp')


inds_treat_hosp = intersect(inds_treat,inds_hosp);

%with treatment
load('sol_treat_3_t0')

icount = 1;
jcount = 1;
kcount = 1;
lcount = 1;
mcount = 1;

for k = inds_treat_hosp
    %hospitalization
    % if ismember(k,inds_hosp)
    %     P_hosp_treat(icount,:) = [sol{k}.P_t(1:4)];
    %     Tsymp_hosp_treat(icount) = OnsetSympt_check(k);
    %     Thosp_hosp_treat(icount) = -T_hosp_check(k);
    %     icount = icount +1;
    % else
    %     % P_wo_hosp_treat(jcount,:) = [sol{k}.P_t(1:4)];
    %     % Tsymp_wo_hosp_treat(jcount) = OnsetSympt_check(k);
    %     % Thosp_wo_hosp_treat(jcount) = -T_hosp_check(k);
    %     % jcount = jcount +1;
    % end

    %death
    if ismember(k,inds_dead)
        P_dead_treat(kcount,:) = [sol{k}.P_t(1:4)];
        %P_dead_treat(kcount,:) = [sol{k}.P(1:4)];
        Tsymp_dead_treat(kcount) = OnsetSympt_check(k);
        Thosp_dead_treat(kcount) = -T_hosp_check(k);
        kcount = kcount +1;
    else
        P_wo_dead_treat(lcount,:) = [sol{k}.P_t(1:4)];
        %P_wo_dead_treat(lcount,:) = [sol{k}.P(1:4)];
        Tsymp_wo_dead_treat(lcount) = OnsetSympt_check(k);
        Thosp_wo_dead_treat(lcount) = -T_hosp_check(k);
        lcount = lcount +1;
    end

    % MSE_treat(mcount) = sol{k}.MSE(1);
    % MSE_t_treat(mcount) = sol{k}.MSE_t(1);

    mcount = mcount+1;

end

% P_hosp_treat(:,1) = floor(P_hosp_treat(:,1));
% P_wo_hosp_treat(:,1) = floor(P_wo_hosp_treat(:,1));
P_dead_treat(:,3) = floor(P_dead_treat(:,3));
P_wo_dead_treat(:,3) = floor(P_wo_dead_treat(:,3));

figure
t = tiledlayout(1,5,'TileSpacing','Compact','Padding','Compact');
n1 = normrnd(0,0.01,length(P_wo_dead_treat),1);
n2 = normrnd(0,0.01,length(P_dead_treat),1);
for i = 1:4

    subplot(1,5,i)
    if i == 1
        scatter(1+n1,log10(P_wo_dead_treat(:,i)),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{1},'MarkerFaceAlpha',FA)
        hold on
        plot([1-0.02,1+0.02],[median(log10(P_wo_dead_treat(:,i))),median(log10(P_wo_dead_treat(:,i)))],'-', 'Color', C{6}, 'LineWidth', LW)
        hold on
        scatter(1.1+n2,log10(P_dead_treat(:,i)),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{2},'MarkerFaceAlpha',FA)
        plot([1.1-0.02,1.1+0.02],[median(log10(P_dead_treat(:,i))),median(log10(P_dead_treat(:,i)))],'-', 'Color', C{6}, 'LineWidth', LW)

    else
        %for i == 1: time from infection to first measurement
        scatter(1+n1,P_wo_dead_treat(:,i),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{1},'MarkerFaceAlpha',FA)
        hold on
        plot([1-0.02,1+0.02],[median(P_wo_dead_treat(:,i)),median(P_wo_dead_treat(:,i))],'-', 'Color', C{6}, 'LineWidth', LW)
        hold on
        scatter(1.1+n2,P_dead_treat(:,i),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{2},'MarkerFaceAlpha',FA)
        plot([1.1-0.02,1.1+0.02],[median(P_dead_treat(:,i)),median(P_dead_treat(:,i))],'-', 'Color', C{6}, 'LineWidth', LW)
        ylim([0,Inf])
    end

    box off

    xticks([1,1.1])
    xticklabels({})
    xlim([1-max(abs(n1))-0.01,1.1+max(abs(n2))+0.01])

    if i == 1
        sprintf('median test of par_%d: %d', i, mediantest(log10(P_wo_dead_treat(:,i)),log10(P_dead_treat(:,i))))
        sprintf('kstest2 of par_%d: %d', i, kstest2(log10(P_wo_dead_treat(:,i)),log10(P_dead_treat(:,i))))
        sprintf('medians of par_%d: hosp %d vs non-hosp %d', i, median(log10(P_wo_dead_treat(:,i))), median(log10(P_dead_treat(:,i))))
    else
        sprintf('median test of par_%d: %d', i, mediantest(P_wo_dead_treat(:,i),P_dead_treat(:,i)))
        sprintf('kstest2 of par_%d: %d', i, kstest2(P_wo_dead_treat(:,i),P_dead_treat(:,i)))
        sprintf('medians of par_%d: hosp %d vs non-hosp %d', i, median(P_wo_dead_treat(:,i)), median(P_dead_treat(:,i)))
    end
end

%time from symptom onset to hospitalization
subplot(1,5,5)
scatter(1+n1(~isnan(Thosp_wo_dead_treat')),Tsymp_wo_dead_treat(~isnan(Thosp_wo_dead_treat'))-Thosp_wo_dead_treat,'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{1},'MarkerFaceAlpha',FA)
hold on
plot([1-0.02,1+0.02],[median(Tsymp_wo_dead_treat(~isnan(Thosp_wo_dead_treat'))-Thosp_wo_dead_treat),...
    median(Tsymp_wo_dead_treat(~isnan(Thosp_wo_dead_treat'))-Thosp_wo_dead_treat)],'-', 'Color', C{6}, 'LineWidth', LW)
hold on
scatter(1.1+n2(~isnan(Thosp_dead_treat')),Tsymp_dead_treat(~isnan(Thosp_dead_treat'))-Thosp_dead_treat(~isnan(Thosp_dead_treat')),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{2},'MarkerFaceAlpha',FA)
plot([1.1-0.02,1.1+0.02],[median(Tsymp_dead_treat(~isnan(Thosp_dead_treat'))-Thosp_dead_treat),...
    median(Tsymp_dead_treat(~isnan(Thosp_dead_treat'))-Thosp_dead_treat)],'-', 'Color', C{6}, 'LineWidth', LW)
hold on
ylim([0,Inf])

xticks([1,1.1])
xticklabels({})
xlim([1-max(abs(n1))-0.01,1.1+max(abs(n2))+0.01])

sprintf('median test of par_%d: %d', 5, mediantest(Tsymp_wo_dead_treat(~isnan(Thosp_wo_dead_treat'))-Thosp_wo_dead_treat(~isnan(Thosp_wo_dead_treat')),Tsymp_dead_treat(~isnan(Thosp_dead_treat'))-Thosp_dead_treat(~isnan(Thosp_dead_treat'))))
sprintf('kstest2 of par_%d: %d', 5, kstest2(Tsymp_wo_dead_treat(~isnan(Thosp_wo_dead_treat'))-Thosp_wo_dead_treat(~isnan(Thosp_wo_dead_treat')),Tsymp_dead_treat(~isnan(Thosp_dead_treat'))-Thosp_dead_treat(~isnan(Thosp_dead_treat'))))
sprintf('medians of par_%d: dead %d vs non-dead %d', 5, median(Tsymp_wo_dead_treat(~isnan(Thosp_wo_dead_treat'))-Thosp_wo_dead_treat(~isnan(Thosp_wo_dead_treat'))), median(Tsymp_dead_treat(~isnan(Thosp_dead_treat'))-Thosp_dead_treat(~isnan(Thosp_dead_treat'))))

clear figure_property;
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
figure_property.Width= '6'; % Figure width on canvas
figure_property.Height= '1.5'; % Figure height on canvas
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
hgexport(gcf,'Par_analysis_treat_dead.pdf',figure_property); %Set desired file name


% [val,ind] = sort(MSE_t./MSE)
%
% [val,ind] = sort(MSE_t_treat./MSE_treat)

% figure
% scatter(MSE_treat,MSE_t_treat)
% hold on
% plot(0:max(max(MSE_treat),max(MSE_t_treat))+10,0:max(max(MSE_treat),max(MSE_t_treat))+10,'-')
% hold on
% plot(0:max(max(MSE_treat),max(MSE_t_treat))+10,(0:max(max(MSE_treat),max(MSE_t_treat))+10)*0.9,'--')
% hold on
% plot(0:max(max(MSE_treat),max(MSE_t_treat))+10,(0:max(max(MSE_treat),max(MSE_t_treat))+10)*0.8,':')
% hold on
% plot(0:max(max(MSE_treat),max(MSE_t_treat))+10,(0:max(max(MSE_treat),max(MSE_t_treat))+10)*0.7,'--')
% hold on
% plot(0:max(max(MSE_treat),max(MSE_t_treat))+10,(0:max(max(MSE_treat),max(MSE_t_treat))+10)*0.6,':')
% hold on
% plot(0:max(max(MSE_treat),max(MSE_t_treat))+10,(0:max(max(MSE_treat),max(MSE_t_treat))+10)*0.5,'--')

inds_wo_treat_hosp = intersect(inds_wo_treat,inds_hosp);

%compare untreated vs treated
load('sol_wo_treat_3_t0')
icount = 1;
% for k = inds_wo_treat_hosp
%     P_wo_treat(icount,:) = sol{k}.P;
%     icount = icount + 1;
% end
for k = inds_wo_treat
    P_wo_treat(icount,:) = sol{k}.P;
    icount = icount + 1;
end

P_wo_treat(:,3) = floor(P_wo_treat(:,3));

load('sol_treat_3_t0')
jcount = 1;
for k = inds_treat_hosp
    P_treat(jcount,:) = sol{k}.P_t;
    jcount = jcount + 1;
end

P_treat(:,3) = floor(P_treat(:,3));

figure
t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
n1 = normrnd(0,0.01,length(P_wo_hosp),1);
n2 = normrnd(0,0.01,length(P_hosp),1);
n3 = normrnd(0,0.01,length(P_wo_dead_treat),1);
n4 = normrnd(0,0.01,length(P_dead_treat),1);
for i = 1:3
    subplot(1,3,i)
    if i == 1
        scatter(1+n1,log10(P_wo_hosp(:,i)),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{3},'MarkerFaceAlpha',FA)
        hold on
        scatter(1+n2,log10(P_hosp(:,i)),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{4},'MarkerFaceAlpha',FA)
        hold on
        plot([1-0.02,1+0.02],[median(log10(P_wo_treat(:,i))),median(log10(P_wo_treat(:,i)))],'-', 'Color', C{6}, 'LineWidth', LW)
        hold on
        scatter(1.1+n3,log10(P_wo_dead_treat(:,i)),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{1},'MarkerFaceAlpha',FA)
        hold on
        scatter(1.1+n4,log10(P_dead_treat(:,i)),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{2},'MarkerFaceAlpha',FA)
        hold on
        plot([1.1-0.02,1.1+0.02],[median(log10(P_treat(:,i))),median(log10(P_treat(:,i)))],'-', 'Color', C{6}, 'LineWidth', LW)
    else
        scatter(1+n1,P_wo_hosp(:,i),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{3},'MarkerFaceAlpha',FA)
        hold on
        scatter(1+n2,P_hosp(:,i),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{4},'MarkerFaceAlpha',FA)
        hold on
        plot([1-0.02,1+0.02],[median(P_wo_treat(:,i)),median(P_wo_treat(:,i))],'-', 'Color', C{6}, 'LineWidth', LW)
        hold on
        scatter(1.1+n3,P_wo_dead_treat(:,i),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{1},'MarkerFaceAlpha',FA)
        hold on
        scatter(1.1+n4,P_dead_treat(:,i),'o','MarkerEdgeColor','none','SizeData',sD,'MarkerFaceColor',C{2},'MarkerFaceAlpha',FA)
        hold on
        plot([1.1-0.02,1.1+0.02],[median(P_treat(:,i)),median(P_treat(:,i))],'-', 'Color', C{6}, 'LineWidth', LW)
        ylim([0,Inf])
    end
    hold on

    xticks([1,1.1])
    xticklabels({})
    xlim([1-max(abs(n1))-0.01,1.1+max(abs(n2))+0.01])

    box off

    if i == 1
        sprintf('median test of par_%d: %d', i, mediantest(log10(P_wo_treat(:,i)),log10(P_treat(:,i))))
        sprintf('kstest2 of par_%d: %d', i, kstest2(log10(P_wo_treat(:,i)),log10(P_treat(:,i))))
        sprintf('medians of par_%d: hosp %d vs non-hosp %d', i, median(log10(P_wo_treat(:,i))), median(log10(P_treat(:,i))))
    else
        sprintf('median test of par_%d: %d', i, mediantest(P_wo_treat(:,i),P_treat(:,i)))
        sprintf('kstest2 of par_%d: %d', i, kstest2(P_wo_treat(:,i),P_treat(:,i)))
        sprintf('medians of par_%d: hosp %d vs non-hosp %d', i, median(P_wo_treat(:,i)), median(P_treat(:,i)))
    end
end

clear figure_property;
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
figure_property.Width= '3'; % Figure width on canvas
figure_property.Height= '1.5'; % Figure height on canvas
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
hgexport(gcf,'Par_analysis_wotreat_treat.pdf',figure_property); %Set desired file name