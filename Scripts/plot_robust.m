clearvars;
clc;

%load data set
data = readtable('Data_Greek_min_text.xlsx');

%extract the CT values (Vals), the dates at which the CT values where taken (Dates),
%whether an individual was treated (Treat), hospitalized (Hosp), or died (Death).
%extract the days of symptom onset relative to the first CT value (OnsetSympt)
%and anonymized ID of individual (ID).
%if the individual was hospitalized, also extract the date of hospitalization (Hosp_dates).
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

%find the subset of indiviudals for whcih we have at least 4 CT
%measurements (model cohort)
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

%convert the dates to days relative to the day of symptom onset

%number of days in a month for a regular year
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
                val(i_dates) = 37; %detection threshold, according to Ioanna
                val1(i_dates) = 37; %detection threshold, according to Ioanna
                val3(i_dates) = 37; %detection threshold, according to Ioanna
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

%sanity check: are the dates of CT measurements in order? are there the
%same number of CT measurements as dates?
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

%extract the subset of indiviudals that were treated/untreated
inds_treat = find(Treat_check == 1);
inds_wo_treat = find(Treat_check == 0);

%exract the subset of indiviudals that were hospitalized/non-hospitalzied
inds_hosp = find(~isnan(T_hosp_check));
inds_wo_hosp = find(isnan(T_hosp_check));

%exract the subset of indiviudals that survived/died
inds_dead = find(Death_check == 1);
inds_wo_dead = find(Death_check == 0);

%get indices of all 4 classes
inds_wo_treat_hosp = intersect(inds_wo_treat,inds_hosp);
inds_wo_treat_wo_hosp = intersect(inds_wo_treat,inds_wo_hosp);
inds_treat_hosp_wo_dead = intersect(intersect(inds_treat,inds_hosp), inds_wo_dead);
inds_treat_hosp_dead = intersect(intersect(inds_treat,inds_hosp), inds_dead);

C = [{[95,211,141]./255}, {[95,95,211]./255}, {[209,94,94]./255}, {[187,94,209]./255}, {[200,200,200]./255}, {[0,0,0]}];

sol_treat = load('sol_treat_3_t0');
sol_wo_treat = load('sol_wo_treat_3_t0');

% T1 = [];
% T2 = [];
% T3 = [];
% T4 = [];
i_top = 3;

for i_par = 1:4
    figure;
    i_count = 1;

    for i_ind = 1:length(inds_wo_treat_wo_hosp)

        yline(-i_count,'Color', C{5});
        hold on
        if i_par < 4
            if i_par == 1
                plot(log10(sol_wo_treat.sol{inds_wo_treat_wo_hosp(i_ind)}.Pall(1:i_top,i_par)),-i_count,'.','Color', C{1});
                hold on
                plot([min(log10(sol_wo_treat.sol{inds_wo_treat_wo_hosp(i_ind)}.Pall(1:i_top,i_par))),max(log10(sol_wo_treat.sol{inds_wo_treat_wo_hosp(i_ind)}.Pall(1:i_top,i_par)))],[-i_count,-i_count],'-','Color', C{1})
                hold on
                plot(log10(sol_wo_treat.sol{inds_wo_treat_wo_hosp(i_ind)}.Pall(1,i_par)),-i_count,'.','Color', C{1},'Markersize',9);
                hold on
            else
                plot(sol_wo_treat.sol{inds_wo_treat_wo_hosp(i_ind)}.Pall(1:i_top,i_par),-i_count,'.','Color', C{1});
                hold on
                plot([min(sol_wo_treat.sol{inds_wo_treat_wo_hosp(i_ind)}.Pall(1:i_top,i_par)),max(sol_wo_treat.sol{inds_wo_treat_wo_hosp(i_ind)}.Pall(1:i_top,i_par))],[-i_count,-i_count],'-','Color', C{1})
                hold on
                plot((sol_wo_treat.sol{inds_wo_treat_wo_hosp(i_ind)}.Pall(1,i_par)),-i_count,'.','Color', C{1},'Markersize',9);
                hold on
            end
        end
        i_count = i_count+1;

        % T1 = [T1,T_check{inds_wo_treat_wo_hosp(i_ind)} + OnsetSympt_check(inds_wo_treat_wo_hosp(i_ind))];
    end

    for i_ind = 1:length(inds_wo_treat_hosp)

        yline(-i_count,'Color', C{5});
        hold on
        if i_par < 4
            if i_par == 1
                plot(log10(sol_wo_treat.sol{inds_wo_treat_hosp(i_ind)}.Pall(1:i_top,i_par)),-i_count,'.','Color', C{2});
                hold on
                plot([min(log10(sol_wo_treat.sol{inds_wo_treat_hosp(i_ind)}.Pall(1:i_top,i_par))),max(log10(sol_wo_treat.sol{inds_wo_treat_hosp(i_ind)}.Pall(1:i_top,i_par)))],[-i_count,-i_count],'-','Color', C{2})
                hold on
                plot(log10(sol_wo_treat.sol{inds_wo_treat_hosp(i_ind)}.Pall(1,i_par)),-i_count,'.','Color', C{2},'Markersize',9);
                hold on
            else
                plot(sol_wo_treat.sol{inds_wo_treat_hosp(i_ind)}.Pall(1:i_top,i_par),-i_count,'.','Color', C{2});
                hold on
                plot([min(sol_wo_treat.sol{inds_wo_treat_hosp(i_ind)}.Pall(1:i_top,i_par)),max(sol_wo_treat.sol{inds_wo_treat_hosp(i_ind)}.Pall(1:i_top,i_par))],[-i_count,-i_count],'-','Color', C{2})
                hold on
                plot((sol_wo_treat.sol{inds_wo_treat_hosp(i_ind)}.Pall(1,i_par)),-i_count,'.','Color', C{2},'Markersize',9);
                hold on
            end
        end
        i_count = i_count+1;

        % T2 = [T2,T_check{inds_wo_treat_hosp(i_ind)} + OnsetSympt_check(inds_wo_treat_hosp(i_ind))];
    end

    for i_ind = 1:length(inds_treat_hosp_wo_dead)

        yline(-i_count,'Color', C{5});
        hold on
        if i_par == 1
            plot(log10(sol_treat.sol{inds_treat_hosp_wo_dead(i_ind)}.Pall_t(1:i_top,i_par)),-i_count,'.','Color', C{4});
            hold on
            plot([min(log10(sol_treat.sol{inds_treat_hosp_wo_dead(i_ind)}.Pall_t(1:i_top,i_par))),max(log10(sol_treat.sol{inds_treat_hosp_wo_dead(i_ind)}.Pall_t(1:i_top,i_par)))],[-i_count,-i_count],'-','Color', C{4})
            hold on
            plot(log10(sol_treat.sol{inds_treat_hosp_wo_dead(i_ind)}.Pall_t(1,i_par)),-i_count,'.','Color', C{4},'Markersize',9);
            hold on
        else
            plot(sol_treat.sol{inds_treat_hosp_wo_dead(i_ind)}.Pall_t(1:i_top,i_par),-i_count,'.','Color', C{4});
            hold on
            plot([min(sol_treat.sol{inds_treat_hosp_wo_dead(i_ind)}.Pall_t(1:i_top,i_par)),max(sol_treat.sol{inds_treat_hosp_wo_dead(i_ind)}.Pall_t(1:i_top,i_par))],[-i_count,-i_count],'-','Color', C{4})
            hold on
            plot((sol_treat.sol{inds_treat_hosp_wo_dead(i_ind)}.Pall_t(1,i_par)),-i_count,'.','Color', C{4},'Markersize',9);
            hold on
        end
        i_count = i_count+1;

        % T3 = [T3,T_check{inds_treat_hosp_wo_dead(i_ind)} + OnsetSympt_check(inds_treat_hosp_wo_dead(i_ind))];
    end

    for i_ind = 1:length(inds_treat_hosp_dead)

        yline(-i_count,'Color', C{5});
        hold on
        if i_par == 1
            plot(log10(sol_treat.sol{inds_treat_hosp_dead(i_ind)}.Pall_t(1:i_top,i_par)),-i_count,'.','Color', C{3});
            hold on
            plot([min(log10(sol_treat.sol{inds_treat_hosp_dead(i_ind)}.Pall_t(1:i_top,i_par))),max(log10(sol_treat.sol{inds_treat_hosp_dead(i_ind)}.Pall_t(1:i_top,i_par)))],[-i_count,-i_count],'-','Color', C{3})
            hold on
            plot(log10(sol_treat.sol{inds_treat_hosp_dead(i_ind)}.Pall_t(1,i_par)),-i_count,'.','Color', C{3},'Markersize',9);
            hold on
        else
            plot(sol_treat.sol{inds_treat_hosp_dead(i_ind)}.Pall_t(1:i_top,i_par),-i_count,'.','Color', C{3});
            hold on
            plot([min(sol_treat.sol{inds_treat_hosp_dead(i_ind)}.Pall_t(1:i_top,i_par)),max(sol_treat.sol{inds_treat_hosp_dead(i_ind)}.Pall_t(1:i_top,i_par))],[-i_count,-i_count],'-','Color', C{3})
            hold on
            plot((sol_treat.sol{inds_treat_hosp_dead(i_ind)}.Pall_t(1,i_par)),-i_count,'.','Color', C{3},'Markersize',9);
            hold on
        end
        i_count = i_count+1;

        % T4 = [T4,T_check{inds_treat_hosp_dead(i_ind)} + OnsetSympt_check(inds_treat_hosp_dead(i_ind))];
    end
    set(gcf, 'Position',  [100, 100, 200, 450])
    saveas(gcf,sprintf('Figure_Robust%d.pdf',i_par)); %Set desired file name
end
% saveas(gcf,'Figure_Timing.pdf'); %Set desired file name


figure;
    i_count = 1;

    for i_ind = 1:length(inds_wo_treat_wo_hosp)

        yline(-i_count,'Color', C{5});
        hold on
            plot(sol_wo_treat.sol{inds_wo_treat_wo_hosp(i_ind)}.fval_sorted(1:i_top),-i_count,'.','Color', C{1});
            hold on
            plot([min(sol_wo_treat.sol{inds_wo_treat_wo_hosp(i_ind)}.fval_sorted(1:i_top)),max(sol_wo_treat.sol{inds_wo_treat_wo_hosp(i_ind)}.fval_sorted(1:i_top))],[-i_count,-i_count],'-','Color', C{1})
            hold on
        i_count = i_count+1;

        % T1 = [T1,T_check{inds_wo_treat_wo_hosp(i_ind)} + OnsetSympt_check(inds_wo_treat_wo_hosp(i_ind))];
    end

    for i_ind = 1:length(inds_wo_treat_hosp)

        yline(-i_count,'Color', C{5});
        hold on
            plot(sol_wo_treat.sol{inds_wo_treat_hosp(i_ind)}.fval_sorted(1:i_top),-i_count,'.','Color', C{2});
            hold on
            plot([min(sol_wo_treat.sol{inds_wo_treat_hosp(i_ind)}.fval_sorted(1:i_top)),max(sol_wo_treat.sol{inds_wo_treat_hosp(i_ind)}.fval_sorted(1:i_top))],[-i_count,-i_count],'-','Color', C{2})
            hold on
        i_count = i_count+1;

        % T2 = [T2,T_check{inds_wo_treat_hosp(i_ind)} + OnsetSympt_check(inds_wo_treat_hosp(i_ind))];
    end

    for i_ind = 1:length(inds_treat_hosp_wo_dead)

        yline(-i_count,'Color', C{5});
        hold on
        plot(sol_treat.sol{inds_treat_hosp_wo_dead(i_ind)}.fval_sorted_t(1:i_top),-i_count,'.','Color', C{4});
        hold on
        plot([min(sol_treat.sol{inds_treat_hosp_wo_dead(i_ind)}.fval_sorted_t(1:i_top)),max(sol_treat.sol{inds_treat_hosp_wo_dead(i_ind)}.fval_sorted_t(1:i_top))],[-i_count,-i_count],'-','Color', C{4})
        hold on
        i_count = i_count+1;

        % T3 = [T3,T_check{inds_treat_hosp_wo_dead(i_ind)} + OnsetSympt_check(inds_treat_hosp_wo_dead(i_ind))];
    end

    for i_ind = 1:length(inds_treat_hosp_dead)

        yline(-i_count,'Color', C{5});
        hold on
        plot(sol_treat.sol{inds_treat_hosp_dead(i_ind)}.fval_sorted_t(1:i_top),-i_count,'.','Color', C{3});
        hold on
        plot([min(sol_treat.sol{inds_treat_hosp_dead(i_ind)}.fval_sorted_t(1:i_top)),max(sol_treat.sol{inds_treat_hosp_dead(i_ind)}.fval_sorted_t(1:i_top))],[-i_count,-i_count],'-','Color', C{3})
        hold on
        i_count = i_count+1;

        % T4 = [T4,T_check{inds_treat_hosp_dead(i_ind)} + OnsetSympt_check(inds_treat_hosp_dead(i_ind))];
    end
    set(gcf, 'Position',  [100, 100, 200, 450])
saveas(gcf,'Figure_RobustLL.pdf');
% h2 = figure;
% [f1,xi1] = ksdensity(T1,'Bandwidth',2);
% area(xi1,f1,'EdgeColor','none','FaceColor',C{1},'FaceAlpha',0.75)
% hold on
% [f2,xi2] = ksdensity(T2,'Bandwidth',2);
% area(xi2,f2,'EdgeColor','none','FaceColor',C{2},'FaceAlpha',0.75)
% hold on
% [f3,xi3] = ksdensity(T3,'Bandwidth',2);
% area(xi3,f3,'EdgeColor','none','FaceColor',C{4},'FaceAlpha',0.75)
% hold on
% [f4,xi4] = ksdensity(T4,'Bandwidth',2);
% area(xi4,f4,'EdgeColor','none','FaceColor',C{3},'FaceAlpha',0.75)
%
% set(h2, 'Position',  [100, 100, 200, 150])
% saveas(gcf,'Figure_Timing2.pdf'); %Set desired file name
%
