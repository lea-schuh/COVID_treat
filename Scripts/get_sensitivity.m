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
                val(i_dates) = 35; %detection threshold, according to Ioanna
                val1(i_dates) = 35; %detection threshold, according to Ioanna
                val3(i_dates) = 35; %detection threshold, according to Ioanna
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

for i_cohort = 1:4
    figure
    if i_cohort < 3
        t = tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact');
    else
        t = tiledlayout(3,4,'TileSpacing','Compact','Padding','Compact');
    end
    for i_ind = 1:3

        if i_cohort == 1
            ID_opt = inds_wo_treat_wo_hosp(i_ind);
        elseif i_cohort == 2
            ID_opt = inds_wo_treat_hosp(i_ind);
        elseif i_cohort == 3
            ID_opt = inds_treat_hosp_wo_dead(i_ind);
        else
            ID_opt = inds_treat_hosp_dead(i_ind);
        end
        %make time points relative to symptom onset
        T_check{ID_opt} = T_check{ID_opt} + OnsetSympt_check(ID_opt);
        T_hosp_check(ID_opt) = T_hosp_check(ID_opt) + OnsetSympt_check(ID_opt);
        OnsetSympt_check(ID_opt) = 0;

        if i_cohort < 3
            par = log10(sol_wo_treat.sol{ID_opt}.P);
            max_var = 3;
        else
            par = log10(sol_treat.sol{ID_opt}.P_t);
            max_var = 4;
        end
        LL_norm= loglikelihood_single_Greek(-CT_check{ID_opt},T_check{ID_opt},par,0,[],-CT1_check{ID_opt},-CT3_check{ID_opt});

        for i_var = max_var:-1:1
            nexttile
            if i_cohort < 3
                par = log10(sol_wo_treat.sol{ID_opt}.P);
            else
                par = log10(sol_treat.sol{ID_opt}.P_t);
            end
            i_count = 1;
            for i_dev = [linspace(0.9,1.1,99)]
                if i_cohort < 3
                    par(i_var) = log10(i_dev*sol_wo_treat.sol{ID_opt}.P(i_var));
                    LL{i_ind}(i_var,i_count) = loglikelihood_single_Greek(-CT_check{ID_opt},T_check{ID_opt},par,0,[],-CT1_check{ID_opt},-CT3_check{ID_opt});
                else
                    par(i_var) = log10(i_dev*sol_treat.sol{ID_opt}.P_t(i_var));
                    if i_var == 4
                        par(i_var) = min(log10(i_dev*sol_treat.sol{ID_opt}.P_t(i_var)),log10(1));
                    end
                    LL{i_ind}(i_var,i_count) = loglikelihood_single_Greek(-CT_check{ID_opt},T_check{ID_opt},par,1,T_hosp_check(ID_opt),-CT1_check{ID_opt},-CT3_check{ID_opt});
                end 
                i_count = i_count+1;
            end
            if i_cohort < 3
                par = log10(sol_wo_treat.sol{ID_opt}.P);
                if i_var == 1
                    plot(log10(linspace(0.9,1.1,99)*sol_wo_treat.sol{ID_opt}.P(i_var)),LL{i_ind}(i_var,:),'-','Color',C{i_cohort},'LineWidth',1.5)
                    hold on
                    plot(log10(sol_wo_treat.sol{ID_opt}.P(i_var)),LL{i_ind}(i_var,50),'.','Color','k','Markersize',7)
                else
                    plot(linspace(0.9,1.1,99)*sol_wo_treat.sol{ID_opt}.P(i_var),LL{i_ind}(i_var,:),'-','Color',C{i_cohort},'LineWidth',1.5)
                    hold on
                    plot(sol_wo_treat.sol{ID_opt}.P(i_var),LL{i_ind}(i_var,50),'.','Color','k','Markersize',7)
                end
            else
                par = log10(sol_treat.sol{ID_opt}.P_t);
                if i_var == 1
                    plot(log10(linspace(0.9,1.1,99)*sol_treat.sol{ID_opt}.P_t(i_var)),LL{i_ind}(i_var,:),'-','Color',C{i_cohort},'LineWidth',1.5)
                    hold on
                    plot(log10(sol_treat.sol{ID_opt}.P_t(i_var)),LL{i_ind}(i_var,50),'.','Color','k','Markersize',7)
                else
                    plot(linspace(0.9,1.1,99)*sol_treat.sol{ID_opt}.P_t(i_var),LL{i_ind}(i_var,:),'-','Color',C{i_cohort},'LineWidth',1.5)
                    hold on
                    plot(sol_treat.sol{ID_opt}.P_t(i_var),LL{i_ind}(i_var,50),'.','Color','k','Markersize',7)
                end
            end

            box off
            
        end
    end
    set(gcf, 'Position',  [100, 100, 400, 300])
    saveas(gcf,sprintf('Figure_Sens%d.pdf',i_cohort));
end