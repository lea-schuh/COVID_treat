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
                val(i_dates) = 37; %according to Ioanna
                val1(i_dates) = 37; %according to Ioanna
                val3(i_dates) = 37; %according to Ioanna
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

inds_hosp = find(~isnan(T_hosp_check));
inds_wo_hosp = find(isnan(T_hosp_check));

inds_dead = find(Death_check == 1);
inds_wo_dead = find(Death_check == 0);

inds_treat_hosp = intersect(inds_treat,inds_hosp);

inds_wo_treat_sort = [intersect(inds_wo_treat,inds_wo_hosp), intersect(inds_wo_treat,inds_hosp)];
inds_treat_hosp_sort = [intersect(inds_treat_hosp,inds_wo_dead), intersect(inds_treat_hosp,inds_dead)];

C = [{[95,211,141]./255}, {[95,95,211]./255}, {[255,127,42]./255}, {[204,0,255]./255}, {[200,200,200]./255}, {[0,0,0]}];

Ind = [{inds_treat_hosp_sort},{inds_wo_treat_sort}];

%for i_var = 1:2
for i_var = 1

    for i_treat = 1:2

        i_ind = 1;

        figure
        % % t = tiledlayout(7,10,'TileSpacing','Compact','Padding','Compact');
        if i_treat == 1
            t = tiledlayout(4,7,'TileSpacing','Compact','Padding','Compact');
            load('sol_treat_3_t0')
            sol_notreat = load('sol_treat_no_treat_3_t0');
        else
            t = tiledlayout(11,7,'TileSpacing','Compact','Padding','Compact');
            load('sol_wo_treat_3_t0')
        end


        %fit over every patient seperately
        % for ID_opt = 1:length(T_check)
        % for ID_opt = inds_treat([2,7,8,9,12,13,16,17,18])
        % for ID_opt = inds_treat_hosp
        for ID_opt = Ind{i_treat}
            % for ID_opt = inds_wo_treat(3)
            % for ID_opt = inds_wo_treat

            nexttile
            if i_treat == 1
                if find(ID_opt==Ind{i_treat}) == size(intersect(inds_treat_hosp,inds_wo_dead),2)+1
                    yticklabels({});
                    xticklabels({});
                    nexttile
                    yticklabels({});
                    xticklabels({});
                    nexttile
                    yticklabels({});
                    xticklabels({});
                    nexttile
                    yticklabels({});
                    xticklabels({});
                    nexttile
                    yticklabels({});
                    xticklabels({});
                end
            else
                if find(ID_opt==Ind{i_treat}) == size(intersect(inds_wo_treat,inds_wo_hosp),2)+1
                    yticklabels({});
                    xticklabels({});
                    nexttile
                    yticklabels({});
                    xticklabels({});
                    nexttile
                    yticklabels({});
                    xticklabels({});
                    nexttile
                    yticklabels({});
                    xticklabels({});
                    nexttile
                    yticklabels({});
                    xticklabels({});
                end
            end
            % figure

            clearvars -except C Ind i_treat ID_check i_ind data ID_opt S0 dN pN...
                b0 dI dV dB T_check CT_check CT1_check CT3_check...
                Treat_check T_hosp_check Death_check sol inds_treat_hosp inds_wo_treat...
                sol_notreat OnsetSympt_check inds_wo_hosp inds_treat inds_wo_treat inds_treat_hosp inds_wo_dead...
                i_var

            ID_opt

            T_check{ID_opt} = T_check{ID_opt} + OnsetSympt_check(ID_opt);
            T_hosp_check(ID_opt) = T_hosp_check(ID_opt) + OnsetSympt_check(ID_opt);
            OnsetSympt_check(ID_opt) = 0;

            %plot the results for short- and long-term dynamics
            tspan = 0:0.1:50; %long term

            %set free parameters to best estimates from optimzation
            if i_treat == 1
                pB = sol{ID_opt}.P_t(1); %immune response on infection rate
                pV = sol{ID_opt}.P_t(2);
                % dB = 10^par_opt(ind(1),3);
                t_shift = floor(sol{ID_opt}.P_t(3));
                alpha = sol{ID_opt}.P_t(4);
            else
                pB = sol{ID_opt}.P(1); %immune response on infection rate
                pV = sol{ID_opt}.P(2);
                % dB = 10^par_opt(ind(1),3);
                t_shift = floor(sol{ID_opt}.P(3));
            end

            if i_var == 1
                if i_treat == 1
                    rectangle('Position',[T_hosp_check(ID_opt)+t_shift,0,5,10],'FaceColor',[200./255,200./255,200./255, 0.5],'EdgeColor','none')
                    hold on
                end

                if isnan(T_hosp_check(ID_opt))
                else
                    xline(T_hosp_check(ID_opt)+t_shift,':','Color','k','Linewidth',1)
                end
                hold on


                % if Treat_check(ID_opt) == 0
                %     plot(T_check{ID_opt}+t_shift_t, -CT_check{ID_opt}, 'o', 'Color','k','Markersize',4,'MarkerFaceColor','k') %data
                % else
                %     plot(T_check{ID_opt}+t_shift_t, -CT_check{ID_opt}, 'o', 'Color','r','Markersize',4,'MarkerFaceColor','r') %data
                % end
                if i_treat == 1
                    if Death_check(ID_opt) == 0
                        scatter(T_check{ID_opt}+t_shift, -0.25*(CT_check{ID_opt})+11.35, 'o','MarkerEdgeColor', 'none','SizeData',12,'MarkerFaceColor',C{1},'MarkerFaceAlpha',0.5) %data
                        hold on
                        scatter(T_check{ID_opt}+t_shift, -0.25*(CT1_check{ID_opt})+11.35, 'o','MarkerEdgeColor', 'none','SizeData',12,'MarkerFaceColor',C{1},'MarkerFaceAlpha',0.5) %data
                        hold on
                        scatter(T_check{ID_opt}+t_shift, -0.25*(CT3_check{ID_opt})+11.35,'o','MarkerEdgeColor', 'none', 'SizeData',12,'MarkerFaceColor',C{1},'MarkerFaceAlpha',0.5) %data
                    else
                        scatter(T_check{ID_opt}+t_shift, -0.25*(CT_check{ID_opt})+11.35, 'o','MarkerEdgeColor', 'none','SizeData',12,'MarkerFaceColor',C{2},'MarkerFaceAlpha',0.5) %data
                        hold on
                        scatter(T_check{ID_opt}+t_shift, -0.25*(CT1_check{ID_opt})+11.35, 'o','MarkerEdgeColor', 'none','SizeData',12,'MarkerFaceColor',C{2},'MarkerFaceAlpha',0.5) %data
                        hold on
                        scatter(T_check{ID_opt}+t_shift, -0.25*(CT3_check{ID_opt})+11.35, 'o','MarkerEdgeColor', 'none','SizeData',12,'MarkerFaceColor',C{2},'MarkerFaceAlpha',0.5) %data
                    end
                else
                    if isnan(T_hosp_check(ID_opt))
                        scatter(T_check{ID_opt}+t_shift, -0.25*(CT_check{ID_opt})+11.35, 'o','MarkerEdgeColor', 'none','SizeData',12,'MarkerFaceColor',C{3},'MarkerFaceAlpha',0.5) %data
                        hold on
                        scatter(T_check{ID_opt}+t_shift, -0.25*(CT1_check{ID_opt})+11.35, 'o','MarkerEdgeColor', 'none','SizeData',12','MarkerFaceColor',C{3},'MarkerFaceAlpha',0.5) %data
                        hold on
                        scatter(T_check{ID_opt}+t_shift, -0.25*(CT3_check{ID_opt})+11.35, 'o','MarkerEdgeColor', 'none','SizeData',12,'MarkerFaceColor',C{3},'MarkerFaceAlpha',0.5) %data
                    else
                        scatter(T_check{ID_opt}+t_shift, -0.25*(CT_check{ID_opt})+11.35, 'o','MarkerEdgeColor', 'none','SizeData',12,'MarkerFaceColor',C{4},'MarkerFaceAlpha',0.5) %data
                        hold on
                        scatter(T_check{ID_opt}+t_shift, -0.25*(CT1_check{ID_opt})+11.35, 'o','MarkerEdgeColor', 'none','SizeData',12,'MarkerFaceColor',C{4},'MarkerFaceAlpha',0.5) %data
                        hold on
                        scatter(T_check{ID_opt}+t_shift, -0.25*(CT3_check{ID_opt})+11.35, 'o','MarkerEdgeColor', 'none','SizeData',12,'MarkerFaceColor',C{4},'MarkerFaceAlpha',0.5) %data
                    end
                end
                hold on

            % plot([T_hosp_check(ID_opt)+t_shift_t,T_hosp_check(ID_opt)+t_shift_t],[-40,-10],':','Color','k','Linewidth',1)
            % hold on
            % plot([0,max(t)],[-42,-42],'--','Color','k','Linewidth',1)
            plot([0,max(tspan)],[-0.25*(37)+11.35,-0.25*(37)+11.35],'--','Color','k','Linewidth',1)
            hold on
            plot(t_shift, 0, '*', 'Color', 'k','Markersize',5) %data
            else
                if isnan(T_hosp_check(ID_opt))
                else
                    xline(T_hosp_check(ID_opt)+t_shift,':','Color','k','Linewidth',1)
                end
                hold on
                plot(t_shift, 0, '*', 'Color', 'k','Markersize',5) %data
                hold on
            end

            %determine the initial values and B_thres for simulation
            y0 = [S0, 1, 0, 0];
            B_thres = 1-dI*dV/(b0*S0*(pV-dI));
            options = odeset('NonNegative',[1,2,3,4]); %specify non-negative values
            if i_treat == 1
                [t,y] = ode45(@(t,y) odefcn_single_infection_S_R0(t,y,b0,dI,pV,dV,pN,dN,pB,dB,B_thres,alpha,T_hosp_check(ID_opt)+t_shift,S0), tspan, y0, options);
            else
                [t,y] = ode45(@(t,y) odefcn_single_infection_S_R0(t,y,b0,dI,pV,dV,pN,dN,pB,dB,B_thres,1,0,S0), tspan, y0,options);
            end
            y_short = y(:,3);
            %if values too small, fix at 1 (numerical problems)
            y_short(y_short<1)=1;
            %plot(t,-(log(y(:,3)/(1.441*10^14)))/(-0.685),'Color','k','Linewidth',2) %Ke 2021
            % plot(t,-(log10(y_short)-11.35)/(-0.25),'Color','k','Linewidth',1) %Ke 2022 nasal
            if i_var == 1
                plot(t,log10(y_short),'Color',C{6},'Linewidth',1) %Ke 2022 nasal
            else
                plot(t,y(:,4),'Color',C{6},'Linewidth',1) %Ke 2022 nasal
                hold on
                xline(t(find(max(y(:,3)) == y(:,3))),'--','Color','k','Linewidth',1)
            end
            hold on

            if i_treat == 1
                pB = sol_notreat.sol{ID_opt}.P(1); %immune response on infection rate
                pV = sol_notreat.sol{ID_opt}.P(2);
                % dB = 10^par_opt(ind(1),3);
                t_shift = floor(sol_notreat.sol{ID_opt}.P(3));

                tspan = 0:0.1:100;

                y0 = [S0, 1, 0, 0];
                B_thres = 1-dI*dV/(b0*S0*(pV-dI));
                options = odeset('NonNegative',[1,2,3,4]); %specify non-negative values
                [t,y] = ode45(@(t,y) odefcn_single_infection_S_R0(t,y,b0,dI,pV,dV,pN,dN,pB,dB,B_thres,1,0,S0), tspan, y0,options);

                y_short = y(:,3);
                %if values too small, fix at 1 (numerical problems)
                y_short(y_short<1)=1;
                %plot(t,-(log(y(:,3)/(1.441*10^14)))/(-0.685),'Color','k','Linewidth',2) %Ke 2021
                % plot(t,-(log10(y_short)-11.35)/(-0.25),'Color','k','Linewidth',1) %Ke 2022 nasal
                plot(t+(floor(sol{ID_opt}.P_t(3))-t_shift),log10(y_short),'Color',C{5},'Linewidth',1) %Ke 2022 nasal
                hold on
            end

            % if Death_check(ID_opt) == 0
            %     title(sprintf('#%s',ID_check{ID_opt}),'FontWeight','normal')
            % else
            %     title(sprintf('#%sx',ID_check{ID_opt}),'FontWeight','normal')
            % end

            if i_var == 1
                ylim([0,10])
                xlim([0,50])
                xticks([0,20,40]);
                xticklabels({});
                yticks([0,5,10]);
                yticklabels({});
            else
                ylim([0,1])
                xlim([0,50])
                xticks([0,20,40]);
                xticklabels({});
                yticks([0,0.5,1]);
                yticklabels({});
            end

            % subplot(size(y,2)+1,1,5)
            % plot(t,y(:,4), '-', 'Color', 'k','Linewidth',2) %B
            % hold on
            % plot(t,1-(dI*dV/(b0*S0*(pV-dI)))*ones(length(t),1),'--','Color', 'k','Linewidth',2) %B_thres
            % ylim([0,1])
            % hold on

            % set(gcf, 'Position',  [100, 100, 250, 600])
            % if i_treat == 1
            % if ismember(i_ind,[1:7:28])
            %     yticks([-40,-30,-20,-10]);
            % else
            % yticks([-40,-30,-20,-10]);
            % yticklabels({});
            % end

            % if i_ind == 8
            %     ylabel('CT values')
            % end
            %
            % if i_ind == 18
            %     xlabel('Days since infection')
            % end

            % if ismember(i_ind,[13:19])
            %     xticks([0,20,40]);
            % else
            % xticks([0,20,40]);
            % xticklabels({});
            % end

            % else
            %
            %     if ismember(i_ind,[1:7:69])
            %         yticks([-40,-30,-20,-10]);
            %     else
            %         yticks([-40,-30,-20,-10]);
            %         yticklabels({});
            %     end
            %
            %     if i_ind == 29
            %         ylabel('CT values')
            %     end
            %
            %     if i_ind == 68
            %         xlabel('Days since infection')
            %     end
            %
            %     if ismember(i_ind,[63:69])
            %         xticks([0,20,40]);
            %     else
            %         xticks([0,20,40]);
            %         xticklabels({});
            %     end
            % end
            box off

            i_ind = i_ind+1;
        end

        clear figure_property;
        figure_property.units = 'inches';
        figure_property.format = 'pdf';
        figure_property.Preview= 'none';
        figure_property.Width= '8'; % Figure width on canvas
        if i_treat == 1
            figure_property.Height= '3'; % Figure height on canvas
        else
            figure_property.Height= '8.25'; % Figure height on canvas
        end
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
        if i_var == 1
            if i_treat == 1
                hgexport(gcf,'Figure_treat_3_t0.pdf',figure_property); %Set desired file name
            else
                hgexport(gcf,'Figure_wo_treat_3_t0.pdf',figure_property); %Set desired file name
            end
        else
            if i_treat == 1
                hgexport(gcf,'Figure_treat_3_t0_B.pdf',figure_property); %Set desired file name
            else
                hgexport(gcf,'Figure_wo_treat_3_t0_B.pdf',figure_property); %Set desired file name
            end
        end
    end
end