clear 
close all
% Matlab code to process mechanical stick slip experimental data (stress 
% and strain for comparison to acoustic data
% Louisa Brotherson - 11/06/21

%% read in files

addpath 'C:\\Users\\lbrot\\Documents\\Experiments'
addpath 'C:\\Users\\lbrot\\Documents\\Experiments\\Stick-slip\\Images';
addpath 'C:\\Users\\lbrot\\Documents\\Codes';
folder = 'C:\\Users\\lbrot\\Documents\\Experiments\\Stick-slip\\Images\\';
loading_velocity = 1.0; % in um/s (micro meters/s)
norm = 1000; % the lower the norm, the less points you'll have (will start at higher index)
stickslip = 1;
plots = 1;

file_root = 'C:\\Users\\lbrot\\Documents\\Experiments\\Stick-slip\\Processed\\';
ext  = 'Perspex_20_1.0_27_10_20';
file = [file_root ext '_processed.txt'];
Pc = 20;
grit = 0;
cmap = lines(8);

if grit == 120
    
    plotc = cmap(1,:);
   
elseif grit == 400

    plotc = cmap(2,:);
    
elseif grit == 800
    
    plotc = cmap(3,:);    
    
elseif grit == 1200

    plotc = cmap(4,:);
    
else
    
    plotc = cmap(1,:);
    
end

% process files and find stick slips 
[d] = friction_process(file,1);

if stickslip == 1

    [M1,M2,M3,stickslip_ind_start,stickslip_ind_end,breaking,rejected_stickslip] = find_stickslip_v4(d,5,0.75,0.97,0,5.65,0.01);
%     [stickslip_ind_start,stickslip_ind_end] = find_stickslip_v3(d,0.99,0,5.5,0.5);
    % calculate magntiude of stick slips (stress drop) in MPa and strain

    dtau = d(stickslip_ind_start,2)-d(stickslip_ind_end,2); % stress drop in MPa
    epsilon = abs(d(stickslip_ind_start,4))-abs(d(stickslip_ind_end,4)); % stress drop in strain units
    dF = d(stickslip_ind_start,6)-d(stickslip_ind_end,6); % force drop in kN
    ratio_tau = dtau./epsilon;
    ratio_F = dF./epsilon;

    %% plots

    close all

    figure
    plot(d(stickslip_ind_start(1)-norm:stickslip_ind_end(end),1),d(stickslip_ind_start(1)-norm:stickslip_ind_end(end),2),'Color',plotc)
    hold on 
    plot(d(stickslip_ind_start(:),1),d(stickslip_ind_start(:),2),'xr','MarkerSize',5,'LineWidth',1.2)
    hold on 
    plot(d(stickslip_ind_end(:),1),d(stickslip_ind_end(:),2),'xk','MarkerSize',5,'LineWidth',1.2)    
    xlabel('Load point displacement (mm)')
    ylabel('Shear stress (MPa)')
    set(gca,'FontSize',12) 
    saveas(gcf, [folder ext '\\shearstress_stickslips'], 'png')
    
    figure
    plot(d(stickslip_ind_start(1)-norm:stickslip_ind_end(end),1),d(stickslip_ind_start(1)-norm:stickslip_ind_end(end),2),'Color',plotc)
    hold on 
    xlabel('Load point displacement (mm)')
    ylabel('Shear stress (MPa)')
    set(gca,'FontSize',12) 
    saveas(gcf, [folder ext '\\shearstress'], 'png')

    figure
    plot(d(stickslip_ind_start(1)-norm:stickslip_ind_end(end),7),d(stickslip_ind_start(1)-norm:stickslip_ind_end(end),2),'Color',plotc)
    xlabel('Time (s)')
    ylabel('Shear stress (MPa)')
    set(gca,'FontSize',12) 
    saveas(gcf, [folder ext '\\time_shearstress'], 'png')

    % strain calibration

    [p1,S1] = polyfit(d(stickslip_ind_start(1)-norm:stickslip_ind_start(1),4),d(stickslip_ind_start(1)-norm:stickslip_ind_start(1),2),1);
    [y1,delta1]= polyval(p1,d(stickslip_ind_start(1)-norm:stickslip_ind_start(1),4),S1);
    [p2,S2] = polyfit(d(stickslip_ind_start(1)+1:stickslip_ind_end(end),4),d(stickslip_ind_start(1)+1:stickslip_ind_end(end),2),1);
    [y2,delta2] = polyval(p2,d(stickslip_ind_start(1)+1:stickslip_ind_end(end),4),S2);
    r1 = rsquare(d(stickslip_ind_start(1)-norm:stickslip_ind_start(1),2),y1);
    r2 = rsquare(d(stickslip_ind_start(1)+1:stickslip_ind_end(end),4),y2);

    figure
    h(1) = plot(d(stickslip_ind_start(1)-norm:stickslip_ind_start(1),4),d(stickslip_ind_start(1)-norm:stickslip_ind_start(1),2),'xb','MarkerSize',5);
    hold on  
    h(2) = plot(d(stickslip_ind_start(1)+1:stickslip_ind_end(end),4),d(stickslip_ind_start(1)+1:stickslip_ind_end(end),2),'xb','MarkerSize',5);
    hold on
    h(3) = plot(d(stickslip_ind_start(1)-norm:stickslip_ind_start(1),4),y1,'-k','LineWidth',1.2);
    hold on
    h(4:5) = plot(d(stickslip_ind_start(1)-norm:stickslip_ind_start(1),4),y1+2*delta1,'--k',d(stickslip_ind_start(1)-norm:stickslip_ind_start(1),4),y1-2*delta1,'--k','LineWidth',0.6);
    hold on
    h(6) = plot(d(stickslip_ind_start(1)+1:stickslip_ind_end(end),4),y2,'-r','LineWidth',1.2);
    hold on
    h(7:8) = plot(d(stickslip_ind_start(1)+1:stickslip_ind_end(end),4),y2+2*delta2,'--r',d(stickslip_ind_start(1)+1:stickslip_ind_end(end),4),y2-2*delta2,'--r','LineWidth',0.6);
    hold on
    xlabel(['Raw strain gauge signal, (', char(949) ')'])
    ylabel('Shear stress (MPa)')
    legend([h(1),h(3),h(6),h(4),h(7)],'Data','Linear fit 1 (loading)', 'Linear fit 2 (stick-slip)','95% confidence interval','95% confidence interval','Location','northwest')
    saveas(gcf, [folder ext '\\straincalibration'], 'png')

    s_file1 = [dtau dF epsilon ratio_tau ratio_F]; 
    s_file2 = [p1(1) p1(2) p2(1) p2(2) delta1(1) delta2(2) r1 r2];

    fid = fopen([file_root ext '_straincalibration.txt'],'w');
    fprintf(fid,'%5.5d %5.5d %5.5d %5.5d %5.5d %5.5d %5.5d \n',s_file1);
    fprintf(fid,'%5.5d %5.5d %5.5d %5.5d %5.5d %5.5d %5.5d \n',s_file2);
    fclose(fid);

    figure
    plot(d(stickslip_ind_start(1)-norm:stickslip_ind_end(end),1),d(stickslip_ind_start(1)-norm:stickslip_ind_end(end),4),'Color',plotc)
    hold on 
    plot(d(stickslip_ind_start(:),1),d(stickslip_ind_start(:),4),'xr','MarkerSize',5,'LineWidth',1.2)
    hold on 
    plot(d(stickslip_ind_end(:),1),d(stickslip_ind_end(:),4),'xk','MarkerSize',5,'LineWidth',1.2)
    xlabel('Load point displacement (mm)')
    ylabel(['Raw strain gauge signal, (', char(949) ')'])
    saveas(gcf, [folder ext '\\strainstress'], 'png')

    %% displacement calculation

    clear D F_hold
    triaxial_k = load('C:\Users\lbrot\Documents\Experiments\Stick-slip\Triaxial calibration\stiffness50.txt')';

    for i = 1:size(triaxial_k,1)

        for j = Pc/10

            D{i,j}= dF/triaxial_k(i,j); % in mm

        end

    end

    D = D(~cellfun(@isempty,D));
    D = cell2mat(D');
    loading = [D(:,1) D(:,3) D(:,5)];
    unloading = [D(:,2) D(:,4) D(:,6)];

    D_load = mean(loading,2);
    D_unload = mean(unloading,2);
    D_total = mean(D,2);

    for i = 1:length(D)

        F_hold(i,:) = [d(stickslip_ind_start(i),6) d(stickslip_ind_end(i),6)];

    end

    F_mean = mean(F_hold,2); %in kN
    M_0 = (F_mean*1000).*(D/1000); %in Nm
    M_0_total = (F_mean*1000).*(D_total/1000); %in Nm

    m_file = [F_mean'; D_total'; M_0_total'; dtau'];

    fid2 = fopen([file_root ext '_mechprop_calc.txt'],'w');
    fprintf(fid2,'%1.5d %1.5d %1.5d %1.5d \n',m_file);
    fclose(fid2);
    save([file_root ext  '_mech_process.mat'])

    d = d';

    fid3 = fopen([file_root ext '_data.txt'],'w');
    fprintf(fid3,'%5.5d %5.5d %5.5d %5.5d %5.5d %5.5d %5.5d \n',d);
    fclose(fid3);
    
 %% if plots == 1
 
 close all
    
if plots == 1
% plots with all   
    
    figure
    yyaxis left %plot left
    Lyax = gca;                                 % Left Axis Handle
    plot(d(1,:),d(2,:),'LineWidth',1.2,'Color',plotc)
    ylabel('Shear stress (MPa)')    
    yyaxis right
    plot(d(1,:),d(2,:),'-','LineWidth',1.2,'Color',plotc)    
    Ryax = gca;                                 % Right axis handle
    Lyaxt = get(Lyax, 'YTick');                 % Get Left axis ticks
    Lyaxfricco = round(10*(Lyaxt/mean(d(3,:))),1)/10;     % Convert to °C (or whatever you want)
    set(Ryax, 'YTick',Lyaxt,'YTickLabel',num2str(Lyaxfricco','%.2f'))    % New Y-tick values    
    xlim([0 5.5])   
    ylim([0 6])    
    xlabel('Load point displacement (mm)')
    ylabel('Friction coefficient')
    set(gca,'FontSize',12) 
    saveas(gcf, [folder ext '\\shearstress_all'], 'png')

    figure
    yyaxis left %plot left    
    Lyax = gca;                                 % Left Axis Handle
    plot(d(7,:),d(2,:),'LineWidth',1.2,'Color',plotc)
    ylabel('Shear stress (MPa)')    
    yyaxis right
    plot(d(7,:),d(2,:),'-','LineWidth',1.2,'Color',plotc)    
    Ryax = gca;                                 % Right axis handle
    Lyaxt = get(Lyax, 'YTick');                 % Get Left axis ticks
    Lyaxfricco = round(10*(Lyaxt/mean(d(3,:))),1)/10;     % Convert to °C (or whatever you want)
    set(Ryax, 'YTick',Lyaxt,'YTickLabel',num2str(Lyaxfricco','%.2f'))    % New Y-tick values
    ylim([0 6])    
    xlabel('Time (s)')
    ylabel('Friction coefficient')
    set(gca,'FontSize',12) 
    saveas(gcf, [folder ext '\\time_shearstress_all'], 'png') 
    
% plots with zooms - edit as you wish

    figure
    yyaxis left %plot left
    plot(d(1,:),d(2,:),'LineWidth',1.2,'Color',plotc)
    Lyax = gca;                                 % Left Axis Handle
    ylabel('Shear stress (MPa)')    
    yyaxis right
    plot(d(1,:),d(2,:),'-','LineWidth',1.2,'Color',plotc)    
    Ryax = gca;                                 % Right axis handle
    Lyaxt = get(Lyax, 'YTick');                 % Get Left axis ticks
    Lyaxfricco = round(10*(Lyaxt/mean(d(3,:))),1)/10;     % Convert to °C (or whatever you want)
    set(Ryax, 'YTick',Lyaxt,'YTickLabel',num2str(Lyaxfricco','%.2f'))    % New Y-tick values
    ylabel('Friction coefficient')      
    xlim([0.6 1.4])
    ylim([2.5 6])
    xlabel('Load point displacement (mm)')

    set(gca,'FontSize',12) 
    saveas(gcf, [folder ext '\\shearstress_zoom1'], 'png')

    figure
    yyaxis left %plot left
    plot(d(1,:),d(2,:),'LineWidth',1.2,'Color',plotc)
    Lyax = gca;                                 % Left Axis Handle
    ylabel('Shear stress (MPa)')    
    yyaxis right
    plot(d(1,:),d(2,:),'-','LineWidth',1.2,'Color',plotc)    
    Ryax = gca;                                 % Right axis handle
    Lyaxt = get(Lyax, 'YTick');                 % Get Left axis ticks
    Lyaxfricco = round(10*(Lyaxt/mean(d(3,:))),1)/10;     % Convert to °C (or whatever you want)
    set(Ryax, 'YTick',Lyaxt,'YTickLabel',num2str(Lyaxfricco','%.2f'))    % New Y-tick values
    ylabel('Friction coefficient')      
    xlim([0.8 1.6])
    ylim([3 6])
    xlabel('Load point displacement (mm)')
    ylabel('Friction coefficent')
    set(gca,'FontSize',12) 
    saveas(gcf, [folder ext '\\shearstress_zoom2'], 'png')
    
    figure
    yyaxis left %plot left
    plot(d(1,:),d(2,:),'LineWidth',1.2,'Color',plotc)
    Lyax = gca;                                 % Left Axis Handle
    ylabel('Shear stress (MPa)')    
    yyaxis right
    plot(d(1,:),d(2,:),'-','LineWidth',1.2,'Color',plotc)    
    Ryax = gca;                                 % Right axis handle
    Lyaxt = get(Lyax, 'YTick');                 % Get Left axis ticks
    Lyaxfricco = round(10*(Lyaxt/mean(d(3,:))),1)/10;     % Convert to °C (or whatever you want)
    set(Ryax, 'YTick',Lyaxt,'YTickLabel',num2str(Lyaxfricco','%.2f'))    % New Y-tick values
    ylabel('Friction coefficient')      
    xlim([2 2.7])
    ylim([2.5 5])
    xlabel('Load point displacement (mm)')
    ylabel('Shear stress (MPa)')
    set(gca,'FontSize',12) 
    saveas(gcf, [folder ext '\\shearstress_zoom3'], 'png')    
    
    figure
    yyaxis left %plot left
    plot(d(1,:),d(2,:),'LineWidth',1.2,'Color',plotc)
    Lyax = gca;                                 % Left Axis Handle
    ylabel('Shear stress (MPa)')    
    yyaxis right
    plot(d(1,:),d(2,:),'-','LineWidth',1.2,'Color',plotc)    
    Ryax = gca;                                 % Right axis handle
    Lyaxt = get(Lyax, 'YTick');                 % Get Left axis ticks
    Lyaxfricco = round(10*(Lyaxt/mean(d(3,:))),1)/10;     % Convert to °C (or whatever you want)
    set(Ryax, 'YTick',Lyaxt,'YTickLabel',num2str(Lyaxfricco','%.2f'))    % New Y-tick values
    ylabel('Friction coefficient')      
    xlim([3.16 3.64])
    ylim([2.5 4.5])
    xlabel('Load point displacement (mm)')
    set(gca,'FontSize',12) 
    saveas(gcf, [folder ext '\\shearstress_zoom4'], 'png')        
    
    figure
    yyaxis left %plot left
    plot(d(1,:),d(2,:),'LineWidth',1.2,'Color',plotc)
    Lyax = gca;                                 % Left Axis Handle
    ylabel('Shear stress (MPa)')    
    yyaxis right
    plot(d(1,:),d(2,:),'-','LineWidth',1.2,'Color',plotc)    
    Ryax = gca;                                 % Right axis handle
    Lyaxt = get(Lyax, 'YTick');                 % Get Left axis ticks
    Lyaxfricco = round(10*(Lyaxt/mean(d(3,:))),1)/10;     % Convert to °C (or whatever you want)
    set(Ryax, 'YTick',Lyaxt,'YTickLabel',num2str(Lyaxfricco','%.2f'))    % New Y-tick values
    ylabel('Friction coefficient')      
    xlim([5.2 5.4])
    ylim([2.5 4.5])
    xlabel('Load point displacement (mm)')
    set(gca,'FontSize',12) 
    saveas(gcf, [folder ext '\\shearstress_zoom5'], 'png')      
    
end
        
else   
    
    load = 300;
    unload = load/8;
    
    close all

    figure
    plot(d(:,1),d(:,2))
    xlabel('Load point displacement (mm)')
    ylabel('Shear stress (MPa)')
    set(gca,'FontSize',12) 
    saveas(gcf, [folder ext '\\shearstress'], 'png')

    figure
    plot(d(:,7),d(:,2))
    xlabel('Time (s)')
    ylabel('Shear stress (MPa)')
    set(gca,'FontSize',12) 
    saveas(gcf, [folder ext '\\time_shearstress'], 'png')

    % strain calibration

    [p1,S1] = polyfit(d(1:load,4),d(1:load,2),1);
    [y1,delta1]= polyval(p1,d(1:load,4),S1);
    [p2,S2] = polyfit(d(end-unload:end,4),d(end-unload:end,2),1);
    [y2,delta2] = polyval(p2,d(end-unload:end,4),S2);
    r1 = rsquare(d(1:load,2),y1);
    r2 = rsquare(d(end-unload:end,2),y2);

    figure
    h(1) = plot(d(1:load,4),d(1:load,2),'xb','MarkerSize',8);
    hold on  
    h(2) = plot(d(end-unload:end,4),d(end-unload:end,2),'xb','MarkerSize',8);
    hold on
    h(3) = plot(d(1:load,4),y1,'-k','LineWidth',1.5);
    hold on
    h(4:5) = plot(d(1:load,4),y1+2*delta1,'--k',d(1:load,4),y1-2*delta1,'--k','LineWidth',0.6);
    hold on
    h(6) = plot(d(end-unload:end,4),y2,'-r','LineWidth',1.5);
    hold on
    h(7:8) = plot(d(end-unload:end,4),y2+2*delta2,'--r',d(end-unload:end,4),y2-2*delta2,'--r','LineWidth',0.6);
    hold on
    xlabel(['Raw strain gauge signal, (', char(949) ')'])
    ylabel('Shear stress (MPa)')
    legend([h(1),h(3),h(6),h(4),h(7)],'Data','Linear fit 1 (loading)', 'Linear fit 2 (unloading)','95% confidence interval','95% confidence interval','Location','northwest')
    saveas(gcf, [folder ext '\\straincalibration'], 'png')

    figure
    plot(d(:,1),d(:,4))
    hold on 
%     plot(d(:,4),'xr','MarkerSize',5,'LineWidth',1.2)
    hold on 
%     plot(d(:,4),'xk','MarkerSize',5,'LineWidth',1.2)
    xlabel('Displacement (mm)')
    ylabel(['Raw strain gauge signal, (', char(949) ')'])
    saveas(gcf, [folder ext '\\strainstress'], 'png')
    
    s_file2 = [p1(1) p1(2) p2(1) p2(2) delta1(1) delta2(2) r1 r2];

    fid = fopen([file_root ext '_straincalibration.txt'],'w');
    fprintf(fid,'%5.5d %5.5d %5.5d %5.5d %5.5d %5.5d %5.5d \n',s_file2);
    fclose(fid);    
    
    
end