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
    t_Hosp(i_patient) = data.('HospitalizationDuration')(i_patient);
    Death(i_patient) = data.('Death')(i_patient);
    ID(i_patient) = data.('ID')(i_patient);
    Age(i_patient) = data.('Age')(i_patient);
    Sex(i_patient) = data.('Sex')(i_patient);
    Variant(i_patient) = data.('Variant')(i_patient);
    Symptoms(i_patient) = data.('Symptoms')(i_patient);
    t_Symptoms(i_patient) = data.('SymptomDuration')(i_patient);
    ICU(i_patient) = data.('ICU')(i_patient);
    t_ICU(i_patient) = data.('ICUDuration')(i_patient);
    Intubation(i_patient) = data.('Intubation')(i_patient);
    t_Intubation(i_patient) = data.('IntubationDuration')(i_patient);
    Commorbity(i_patient) = data.('Commorbity')(i_patient);
    Vaccination(i_patient) = data.('Vaccination')(i_patient);
    prevInfection(i_patient) = data.('PreviousInfection')(i_patient);

    Cough(i_patient) = data.('Cough')(i_patient);
    t_Cough(i_patient) = data.('CoughDuration')(i_patient);
    Fever(i_patient) = data.('Fever')(i_patient);
    t_Fever(i_patient) = data.('FeverDuration')(i_patient);
    SoreThroat(i_patient) = data.('SoreThroat')(i_patient);
    NasalCongestion(i_patient) = data.('NasalCongestion')(i_patient);
    t_NasalCongestion(i_patient) = data.('NasalCongestionDuration')(i_patient);
    Headache(i_patient) = data.('Headache')(i_patient);
    Fatigue(i_patient) = data.('Fatigue')(i_patient);
    TasteDisorders(i_patient) = data.('TasteDisorders')(i_patient);
    t_TasteDisorders(i_patient) = data.('TasteDisordersDuration')(i_patient);
    SmellDisorders(i_patient) = data.('SmellDisorders')(i_patient);
    t_SmellDisorders(i_patient) = data.('SmellDisordersDuration')(i_patient);
    GISymptoms(i_patient) = data.('GISymptoms')(i_patient);
    AppetiteDisorders(i_patient) = data.('AppetiteDisorders')(i_patient);
    Dyspnoea(i_patient) = data.('Dyspnoea')(i_patient);
    Pneumonia(i_patient) = data.('Pneumonia')(i_patient);

end

icount = 1;
for i_patient = 1:size(data,1)
    if length(Dates{i_patient}) >= 4
        Dates_set{icount} = Dates{i_patient};
        Val_set{icount} = Vals{i_patient};
        ID_set(icount) = ID(i_patient);
        inds_set(icount) = i_patient;
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

    end

    T{i_set} = time;

end

icount = 1;
for i_t = 1:length(T)
    if issorted(T{i_t})
        if sum(T{i_t} >= 0) == length(T{i_t})
            % if ~isnan(CT{i_t})
            inds_check(icount) = inds_set(i_t);
            icount = icount + 1;
            % end
        end
    end
end

inds_treat = intersect(inds_check,find(Treat == 1));
inds_wo_treat = intersect(inds_check,find(Treat == 0));

inds_hosp = intersect(inds_check,find(Hosp == 1));
inds_wo_hosp = intersect(inds_check,find(Hosp == 0));

inds_dead = intersect(inds_check,find(Death == 1));
inds_wo_dead = intersect(inds_check,find(Death == 0));

inds_treat_hosp = intersect(inds_treat,inds_hosp);


inds_wotreat_wohosp = intersect(inds_wo_treat,inds_wo_hosp);
inds_wotreat_hosp = intersect(inds_wo_treat,inds_hosp);
inds_treat_hosp_survived = intersect(inds_treat_hosp,inds_wo_dead);
inds_treat_hosp_dead = intersect(inds_treat_hosp,inds_dead);

Inds{1} = 1:length(Age);
Inds{2} = inds_check;
Inds{3} = inds_wotreat_wohosp;
Inds{4} = inds_wotreat_hosp;
Inds{5} = inds_treat_hosp_survived;
Inds{6} = inds_treat_hosp_dead;

uVariant = unique(Variant);

for i_var = 1:length(uVariant)
    inds = find(strcmpi(Variant, uVariant(i_var))==1);
    if i_var == 1
        bin_Variant(inds) = 1;
    elseif ismember(i_var,[2,4,5])
        bin_Variant(inds) = 3;
    elseif i_var == 3
        bin_Variant(inds) = 2;
    else
        bin_Variant(inds) = 0;
    end
end

icount = 1;
for i_patient = 1:length(Age)
    if Age(i_patient) < 20
        bin_Age(i_patient) = 0;
    elseif Age(i_patient) >= 20 && Age(i_patient) < 40
        bin_Age(i_patient) = 1;
    elseif Age(i_patient) >= 40 && Age(i_patient) < 60
        bin_Age(i_patient) = 2;
    elseif Age(i_patient) >= 60 && Age(i_patient) < 80
        bin_Age(i_patient) = 3;
    else
        bin_Age(i_patient) = 4;
    end

end

labels2 = repmat({''},1,2);
labels4 = repmat({''},1,4);
labels5 = repmat({''},1,5);

Stat(1,:) = Sex;
Stat(2,:) = bin_Age;
Stat(3,:) = bin_Variant;
Stat(4,:) = Vaccination;
Stat(5,:) = prevInfection;
Stat(6,:) = Symptoms;
Stat(7,:) = Commorbity;
Stat(8,:) = Hosp;
Stat(9,:) = Treat;
Stat(10,:) = ICU;
Stat(11,:) = Intubation;
Stat(12,:) = Death;

Stat_sympt(1,:) = Cough;
Stat_sympt(2,:) = Fever;
Stat_sympt(3,:) = SoreThroat;
Stat_sympt(4,:) = NasalCongestion;
Stat_sympt(5,:) = Headache;
Stat_sympt(6,:) = Fatigue;
Stat_sympt(7,:) = TasteDisorders;
Stat_sympt(8,:) = SmellDisorders;
Stat_sympt(9,:) = GISymptoms;
%Stat_sympt(10,:) = AppetiteDisorders;
Stat_sympt(10,:) = Dyspnoea;
Stat_sympt(11,:) = Pneumonia;

newColors = [...
    95, 141,211;   %blue
    255, 85, 85;   %red
    170, 255, 204;  %ligh light green
    85, 255, 153; % light green
    0, 170, 68; %green
    0, 85, 34; %dark green
    0, 43, 17; %dark dark green
    255, 170, 238; %light pink
    255, 85, 221; %pink
    170, 0, 136; %pink pink
    85, 0, 68; %dark pink
    51, 0, 128; %purple
    255, 204, 0]./255; %yellow

figure
t = tiledlayout(length(Inds),size(Stat,1),'TileSpacing','compact');

% for i_stat = 1:size(Stat,2)
for i_set = 1:length(Inds)
    for i_stat = 1:size(Stat,1)
        S0(i_set,i_stat) = sum(Stat(i_stat,Inds{i_set})==0)./length(Stat(i_stat,Inds{i_set}));
        S1(i_set,i_stat) = sum(Stat(i_stat,Inds{i_set})==1)./length(Stat(i_stat,Inds{i_set}));

        nexttile
        if i_stat == 2
            h = pie([sum(Stat(i_stat,Inds{i_set})==0),sum(Stat(i_stat,Inds{i_set})==1),...
                sum(Stat(i_stat,Inds{i_set})==2),sum(Stat(i_stat,Inds{i_set})==3),sum(Stat(i_stat,Inds{i_set})==4)],labels5);
            patchHand = findobj(h, 'Type', 'Patch');
            patchHand(1).FaceColor = newColors(3,:);
            patchHand(2).FaceColor = newColors(4,:);
            patchHand(3).FaceColor = newColors(5,:);
            patchHand(4).FaceColor = newColors(6,:);
            patchHand(5).FaceColor = newColors(7,:);
            S2(i_set,i_stat) = sum(Stat(i_stat,Inds{i_set})==2)./length(Stat(i_stat,Inds{i_set}));
            S3(i_set,i_stat) = sum(Stat(i_stat,Inds{i_set})==3)./length(Stat(i_stat,Inds{i_set}));
            S4(i_set,i_stat) = sum(Stat(i_stat,Inds{i_set})==4)./length(Stat(i_stat,Inds{i_set}));
        elseif i_stat == 3
            h = pie([sum(Stat(i_stat,Inds{i_set})==0),sum(Stat(i_stat,Inds{i_set})==1),sum(Stat(i_stat,Inds{i_set})==2),sum(Stat(i_stat,Inds{i_set})==3)],labels4);
            patchHand = findobj(h, 'Type', 'Patch');
            patchHand(1).FaceColor = newColors(8,:);
            patchHand(2).FaceColor = newColors(9,:);
            patchHand(3).FaceColor = newColors(10,:);
            patchHand(4).FaceColor = newColors(11,:);
            S2(i_set,i_stat) = sum(Stat(i_stat,Inds{i_set})==2)./length(Stat(i_stat,Inds{i_set}));
            S3(i_set,i_stat) = sum(Stat(i_stat,Inds{i_set})==3)./length(Stat(i_stat,Inds{i_set}));
        else
            h = pie([sum(Stat(i_stat,Inds{i_set})==0),sum(Stat(i_stat,Inds{i_set})==1)],labels2);
            patchHand = findobj(h, 'Type', 'Patch');
            if i_stat == 1
                patchHand(1).FaceColor = newColors(1,:);
                patchHand(2).FaceColor = newColors(2,:);
            else
                patchHand(1).FaceColor = newColors(end-1,:);
                patchHand(2).FaceColor = newColors(end,:);
            end
        end
    end
end

median_t_Symptoms = median(t_Symptoms(find(Symptoms==1)));
prctile_t_Symptoms_25 = prctile(t_Symptoms(find(Symptoms==1)),25);
prctile_t_Symptoms_75 = prctile(t_Symptoms(find(Symptoms==1)),75);

median_t_Hosp = median(t_Hosp(find(Hosp==1)));
prctile_t_Hosp_25 = prctile(t_Hosp(find(Hosp==1)),25);
prctile_t_Hosp_75 = prctile(t_Hosp(find(Hosp==1)),75);

median_t_ICU = median(t_ICU(find(ICU==1)));
prctile_t_ICU_25 = prctile(t_ICU(find(ICU==1)),25);
prctile_t_ICU_75 = prctile(t_ICU(find(ICU==1)),75);

median_t_Intubation = median(t_Intubation(find(Intubation==1)));
prctile_t_Intubation_25 = prctile(t_Intubation(find(Intubation==1)),25);
prctile_t_Intubation_75 = prctile(t_Intubation(find(Intubation==1)),75);

clear figure_property;
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
%figure_property.Width= '8'; % Figure width on canvas
%figure_property.Height= '4'; % Figure height on canvas
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
% hgexport(gcf,'Data.pdf',figure_property); %Set desired file name

figure
t = tiledlayout(length(Inds),size(Stat_sympt,1),'TileSpacing','compact');
% for i_stat = 1:size(Stat,2)
for i_set = 1:length(Inds)
    for i_stat = 1:size(Stat_sympt,1)
        nexttile
        S0_sympt(i_set,i_stat) = sum(Stat_sympt(i_stat,Inds{i_set})==0)./length(Stat_sympt(i_stat,Inds{i_set}));
        S1_sympt(i_set,i_stat) = sum(Stat_sympt(i_stat,Inds{i_set})==1)./length(Stat_sympt(i_stat,Inds{i_set}));
        h = pie([sum(Stat_sympt(i_stat,Inds{i_set})==0),sum(Stat_sympt(i_stat,Inds{i_set})==1)],labels2);
        patchHand = findobj(h, 'Type', 'Patch');
        patchHand(1).FaceColor = newColors(end-1,:);
        patchHand(2).FaceColor = newColors(end,:);
    end
end

clear figure_property;
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
%figure_property.Width= '8'; % Figure width on canvas
%figure_property.Height= '4'; % Figure height on canvas
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
% hgexport(gcf,'Sympt.pdf',figure_property); %Set desired file name

figure
t = tiledlayout(3,1,'TileSpacing','compact');

nexttile
h = pie([52,17,10,9,1],labels5);
patchHand = findobj(h, 'Type', 'Patch');
patchHand(1).FaceColor = newColors(1,:);
patchHand(2).FaceColor = newColors(2,:);
patchHand(3).FaceColor = newColors(3,:);
patchHand(4).FaceColor = newColors(4,:);
patchHand(5).FaceColor = newColors(5,:);

nexttile
h = pie([52,17],labels2);
patchHand = findobj(h, 'Type', 'Patch');
patchHand(1).FaceColor = newColors(1,:);
patchHand(2).FaceColor = newColors(2,:);

nexttile
h = pie([10,9],labels2);
patchHand = findobj(h, 'Type', 'Patch');
patchHand(1).FaceColor = newColors(1,:);
patchHand(2).FaceColor = newColors(2,:);

clear figure_property;
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
%figure_property.Width= '8'; % Figure width on canvas
%figure_property.Height= '4'; % Figure height on canvas
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
% hgexport(gcf,'Stat_total.pdf',figure_property); %Set desired file name

%check whether subset representative 
x = norminv([0.025 0.975]);

for i_par = 1:size(Stat,1)
    B = binornd(89,sum(Stat(i_par,Inds{1})==0)./length(Stat(i_par,Inds{1})),10000,1)./89;
    UB(i_par) = mean(B)+x(2)*std(B);
    LB(i_par) = mean(B)+x(1)*std(B);
    if LB(i_par) <= sum(Stat(i_par,Inds{2})==0)/length(Stat(i_par,Inds{2})) && sum(Stat(i_par,Inds{2})==0)/length(Stat(i_par,Inds{2})) <= UB(i_par)
        dec(i_par) = 1;
    else
        dec(i_par) = 0;
    end
end

for i_par = 2:3
    for i_check = 0:4
        B = binornd(89,sum(Stat(i_par,Inds{1})==i_check)./length(Stat(i_par,Inds{1})),10000,1)./89;
        UB= mean(B)+x(2)*std(B);
        LB= mean(B)+x(1)*std(B);
        if LB <= sum(Stat(i_par,Inds{2})==i_check)/length(Stat(i_par,Inds{2})) && sum(Stat(i_par,Inds{2})==i_check)/length(Stat(i_par,Inds{2})) <= UB
            dec2(i_par,i_check+1) = 1;
        else
            dec2(i_par,i_check+1) = 0;
        end
    end
end


for i_subset = 3:6
    for i_par = 1:size(Stat,1)
        B = binornd(length(Inds{i_subset}),sum(Stat(i_par,Inds{2})==0)./length(Stat(i_par,Inds{2})),10000,1)./length(Inds{i_subset});
        UB(i_par) = mean(B)+x(2)*std(B);
        LB(i_par) = mean(B)+x(1)*std(B);
        if LB(i_par) <= sum(Stat(i_par,Inds{i_subset})==0)/length(Stat(i_par,Inds{i_subset})) && sum(Stat(i_par,Inds{i_subset})==0)/length(Stat(i_par,Inds{i_subset})) <= UB(i_par)
            dec3(i_subset-2,i_par) = 1;
        else
            dec3(i_subset-2,i_par) = 0;
        end
    end
end

for i_subset = 3:6
    clear dec4
    for i_par = 2:3
        for i_check = 0:4
            B = binornd(length(Inds{i_subset}),sum(Stat(i_par,Inds{2})==i_check)./length(Stat(i_par,Inds{2})),10000,1)./length(Inds{i_subset});
            UB= mean(B)+x(2)*std(B);
            LB= mean(B)+x(1)*std(B);
            if LB <= sum(Stat(i_par,Inds{i_subset})==i_check)/length(Stat(i_par,Inds{i_subset})) && sum(Stat(i_par,Inds{i_subset})==i_check)/length(Stat(i_par,Inds{i_subset})) <= UB
                dec4(i_par,i_check+1) = 1;
            else
                dec4(i_par,i_check+1) = 0;
            end
        end
    end
    dec4
end




figure
t = tiledlayout(length(Inds),size(Stat_sympt,1),'TileSpacing','compact');
% for i_stat = 1:size(Stat,2)
for i_set = 1:length(Inds)
    for i_stat = 1:size(Stat_sympt,1)
        nexttile
        S0_sympt(i_set,i_stat) = sum(Stat_sympt(i_stat,Inds{i_set}(Sex(Inds{i_set})==0))==0)./length(Stat_sympt(i_stat,Inds{i_set}(Sex(Inds{i_set})==0)));
        S1_sympt(i_set,i_stat) = sum(Stat_sympt(i_stat,Inds{i_set}(Sex(Inds{i_set})==0))==1)./length(Stat_sympt(i_stat,Inds{i_set}(Sex(Inds{i_set})==0)));
        h = pie([sum(Stat_sympt(i_stat,Inds{i_set}(Sex(Inds{i_set})==0))==0),sum(Stat_sympt(i_stat,Inds{i_set}(Sex(Inds{i_set})==0))==1)],labels2);
        patchHand = findobj(h, 'Type', 'Patch');
        patchHand(1).FaceColor = newColors(end-1,:);
        patchHand(2).FaceColor = newColors(end,:);
    end
end

figure
t = tiledlayout(length(Inds),size(Stat_sympt,1),'TileSpacing','compact');
% for i_stat = 1:size(Stat,2)
for i_set = 1:length(Inds)
    for i_stat = 1:size(Stat_sympt,1)
        nexttile
        S0_sympt(i_set,i_stat) = sum(Stat_sympt(i_stat,Inds{i_set}(Sex(Inds{i_set})==1))==0)./length(Stat_sympt(i_stat,Inds{i_set}(Sex(Inds{i_set})==1)));
        S1_sympt(i_set,i_stat) = sum(Stat_sympt(i_stat,Inds{i_set}(Sex(Inds{i_set})==1))==1)./length(Stat_sympt(i_stat,Inds{i_set}(Sex(Inds{i_set})==1)));
        h = pie([sum(Stat_sympt(i_stat,Inds{i_set}(Sex(Inds{i_set})==1))==0),sum(Stat_sympt(i_stat,Inds{i_set}(Sex(Inds{i_set})==1))==1)],labels2);
        patchHand = findobj(h, 'Type', 'Patch');
        patchHand(1).FaceColor = newColors(end-1,:);
        patchHand(2).FaceColor = newColors(end,:);
    end
end


