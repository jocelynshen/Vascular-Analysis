% This script analyzes data from the Angicart output TSV files at
% different intensity thresholds. It calculates the scaling exponents a and
% b using five different methods.

clc; close all; clear

% % Find .tsv files in directory (data for different thresholds)
TSVfiles=dir('*.tsv');      
s=size(TSVfiles);
for i=1:s(1)
    tsv=tdfread(TSVfiles(i).name);          % Import .tsv file
    new_name=TSVfiles(i).name;              % Name new variables
    assignin('base',['Im' new_name(6) new_name(8:12)],tsv)  % Assign values to variable names
end

% List chosen thresholds from 'local/process':

chosen_th=[0.28 0.34 0.68 0.65 0.76 0.88 0.40 0.63 0.65...
    0.78 0.70 0.31 0.67 0.70 0.64 0.34 0.50 0.51]*100;

% Define threshold samples for each image:

Im11_th=25:5:60;
Im12_th=25:5:80;
Im13_th=35:5:95;
Im14_th=30:5:95;
Im15_th=40:5:95;
Im16_th=40:5:95;
Im21_th=30:5:85;
Im22_th=35:5:95;
Im23_th=35:5:95;
Im24_th=35:5:95;
Im25_th=35:5:90;
Im26_th=20:5:50;
Im27_th=35:5:90;
Im31_th=35:5:95;
Im32_th=35:5:90;
Im33_th=15:5:60;
Im34_th=20:5:55;
Im35_th=20:5:60;

% Loop through all images:
for im_num=[11:16 21:27 31:35]
    
    % Values to save between thresholds:
    a_average_ALL_th=[]; b_average_ALL_th=[];
    a_average_1rem_th=[]; b_average_1rem_th=[];
    
    a_q_avg_th=[]; b_s_avg_th=[]; 
    
    a_PL_ALL_th=[]; b_PL_ALL_th=[];
    a_REG_ALL_th=[]; b_REG_ALL_th=[];
    a_HA_ALL_th=[]; b_HA_ALL_th=[];
    
    a_PL_1rem_th=[]; b_PL_1rem_th=[];
    a_REG_1rem_th=[]; b_REG_1rem_th=[];
    a_HA_1rem_th=[]; b_HA_1rem_th=[];
    
    a_count_ALL=[]; b_count_ALL=[];
    a_count_1rem=[]; b_count_1rem=[];
    
    per_sym=[];

    % Loop through all thresholds:
    for th_ind=1:length(eval(['Im' num2str(im_num) '_th']))
       
        % Set variables to columns of .tsv file:
        th_name=eval(['Im' num2str(im_num) '_th']);
        th_val=th_name(th_ind);
        var_name=eval(['Im' num2str(im_num) '_t' num2str(th_val)]);
        beta=var_name.beta;
        gamma=var_name.gamma;
        n_child=2;          % Assume n=2
        len=var_name.len;
        rad=var_name.rad;
        tips=var_name.tips;
        tag=var_name.tag;
        name=var_name.name;
        parent=var_name.parent;
        q=var_name.q;
        s=var_name.s;
        
        rad_th=eval(['Im' num2str(im_num) '_th']);

        %% CONSERVATION-BASED SCALING EXPONENTS
        
        ind_NA_q=[];
        for i=1:length(q)
            if q(i) == 'N'
                ind_NA_q=[ind_NA_q i];
            end
        end
        q(ind_NA_q,:)=[];
        q=str2num(q);
        a_q_avg=mean(1./q);
        a_q_avg_th=[a_q_avg_th a_q_avg];
        
        ind_NA_s=[];
        for i=1:length(s)
            if s(i) == 'N'
                ind_NA_s=[ind_NA_s i];
            end
        end
        s(ind_NA_s,:)=[];
        s=str2num(s);
        b_s_avg=mean(1./s);
        b_s_avg_th=[b_s_avg_th b_s_avg];
        
        %% RATIO-BASED SCALING EXPONENTS

        % Find NA rows in beta and gamma values
        ind_NA=1;
        gamma(ind_NA,:)='9';
        beta(ind_NA,:)='9';
        beta=str2num(beta);     % Convert to double class
        gamma=str2num(gamma);

        % ----- FOR ALL BETA AND GAMMA -----
        
        % Remove NA beta and gamma values:
        beta_NArem=beta; beta_NArem(ind_NA)=[];
        gamma_NArem=gamma; gamma_NArem(ind_NA)=[];
        
        % Calculate a_i and b_i for each node:
        a_i_ALL=-log(beta_NArem)/log(2);
        b_i_ALL=-log(gamma_NArem)/log(2);

        % Calculate arithmetic averages for a and b:
        a_ari_ALL=mean(a_i_ALL);            % Arithmetic
        b_ari_ALL=mean(b_i_ALL);

        % Save the arithmetic average values for this threshold:
        a_average_ALL_th=[a_average_ALL_th a_ari_ALL];
        b_average_ALL_th=[b_average_ALL_th b_ari_ALL];
        
        % Also save how many points are used to calculate a and b:
        a_count_ALL=[a_count_ALL length(a_i_ALL)];
        b_count_ALL=[b_count_ALL length(b_i_ALL)];
        
        % ----- FOR BETA<1 AND GAMMA<1 -----
        
        % Find beta and gamma values >1:
        ind_bad_beta=find(beta>=1); 
        ind_bad_gamma=find(gamma>=1);

        % Remove beta and gamma values >1
        beta_1rem=beta; beta_1rem(ind_bad_beta)=[];
        gamma_1rem=gamma; gamma_1rem(ind_bad_gamma)=[];
        
        % Calculate a_i and b_i for each node:
        a_i_1rem=-log(beta_1rem)/log(2);
        b_i_1rem=-log(gamma_1rem)/log(2);

        % Calculate arithmetic averages for a and b:
        a_ari_1rem=mean(a_i_1rem);            % Arithmetic
        b_ari_1rem=mean(b_i_1rem);

        % Save the arithmetic average values for this threshold:
        a_average_1rem_th=[a_average_1rem_th a_ari_1rem];
        b_average_1rem_th=[b_average_1rem_th b_ari_1rem];
        
        % Also save how many points are used to calculate a and b:
        a_count_1rem=[a_count_1rem length(a_i_1rem)];
        b_count_1rem=[b_count_1rem length(b_i_1rem)];
        
        
        %% DISTRIBUTION, HIERARCHICAL, REGRESSION-BASED EXPONENTS
        
        % ----- FOR ALL BETA AND GAMMA -----
        
        [a_PL,b_PL] = PowerLawFit(rad,len,[],1);
        
        if abs(a_PL)>5
            a_PL=NaN;
        end
        if abs(b_PL)>5
            b_PL=NaN;
        end
        
        a_PL_ALL_th=[a_PL_ALL_th a_PL]; b_PL_ALL_th=[b_PL_ALL_th b_PL];
        
        [a_HA,b_HA] = HierarchicalAveraging(rad,beta,len,gamma,n_child,ind_NA,[],1);
        
        if abs(a_HA)>5
            a_HA=NaN;
        end
        if abs(b_HA)>5
            b_HA=NaN;
        end
        
        a_HA_ALL_th=[a_HA_ALL_th a_HA]; b_HA_ALL_th=[b_HA_ALL_th b_HA];
       
        [a_REG,b_REG] = RegressionExponents(tips,rad,len,[]);
        
        if abs(a_REG)>5
            a_REG=NaN;
        end
        if abs(b_REG)>5
            b_REG=NaN;
        end
        
        a_REG_ALL_th=[a_REG_ALL_th a_REG]; b_REG_ALL_th=[b_REG_ALL_th b_REG];
        
        % ----- FOR BETA<1 AND GAMMA<1 -----
        
        % For a:
        
        [a_PL,b_PL] = PowerLawFit(rad,len,ind_bad_beta',1);
        if abs(a_PL)>5
            a_PL=NaN;
        end
        a_PL_1rem_th=[a_PL_1rem_th a_PL];
        
        [a_HA,b_HA] = HierarchicalAveraging(rad,beta,len,gamma,n_child,ind_NA,ind_bad_beta',1);
        if abs(a_HA)>5
            a_HA=NaN;
        end
        a_HA_1rem_th=[a_HA_1rem_th a_HA];
        
        [a_REG,b_REG] = RegressionExponents(tips,rad,len,ind_bad_beta');
        if abs(a_REG)>5
            a_REG=NaN;
        end
        a_REG_1rem_th=[a_REG_1rem_th a_REG];
        
        % For b:
        
        [a_PL,b_PL] = PowerLawFit(rad,len,ind_bad_gamma',1);
        if abs(b_PL)>5
            b_PL=NaN;
        end
        b_PL_1rem_th=[b_PL_1rem_th b_PL];
        
        [a_HA,b_HA] = HierarchicalAveraging(rad,beta,len,gamma,n_child,ind_NA,ind_bad_gamma',1);
        if abs(b_HA)>5
            b_HA=NaN;
        end
        b_HA_1rem_th=[b_HA_1rem_th b_HA];
       
        [a_REG,b_REG] = RegressionExponents(tips,rad,len,ind_bad_gamma');
        if abs(b_REG)>5
            b_REG=NaN;
        end
        b_REG_1rem_th=[b_REG_1rem_th b_REG];
        
    end
    
    %% PLOTTING RESULTS FOR EACH METHOD ACROSS THRESHOLDS
    
    % Retrieve actual image/patient #:
    im_num_string=num2str(im_num);
    switch im_num_string(1)
        case '1'
            actual_image_num=str2num(im_num_string(2));
        case '2'
            actual_image_num=6+str2num(im_num_string(2));
        case '3'
            actual_image_num=13+str2num(im_num_string(2));
    end  
    
    % The scaling exponents are plotted across thresholds for each
    % method of calculation:
    
    %   (1) Ratio-Based
    %   (2) Distribution-Based
    %   (3) Regression-Based
    %   (4) Hierarchical Averaging
    %   (5) Conservation-Based
    
    method={'average','PL','REG','HA'};
    method_name={'Ratio-Based','Distribution-Based','Regression-Based',...
        'Hierarchical Averaging','Conservation-Based'};
    
    for method_num=1:4
        
        ALL_name_a=eval(['a_' method{method_num} '_ALL_th']);
        ALL_name_b=eval(['b_' method{method_num} '_ALL_th']);
        
        rem1_name_a=eval(['a_' method{method_num} '_1rem_th']);
        rem1_name_b=eval(['b_' method{method_num} '_1rem_th']);
        
        figure; plot(th_name,ALL_name_a,'Color','b'); hold on
        for i=1:length(ALL_name_a)
            text(th_name(i),ALL_name_a(i), num2str(a_count_ALL(i)))
        end

        plot(th_name,ALL_name_b,'Color','r')
        for i=1:length(ALL_name_b)
            text(th_name(i),ALL_name_b(i), num2str(b_count_ALL(i)))
        end

        plot(th_name,rem1_name_a,'Color','b','LineStyle',':')
        for i=1:length(rem1_name_a)
            text(th_name(i),rem1_name_a(i), num2str(a_count_1rem(i)))
        end

        plot(th_name,rem1_name_b,'Color','r','LineStyle',':')
        for i=1:length(rem1_name_b)
            text(th_name(i),rem1_name_b(i), num2str(b_count_1rem(i)))
        end

        % Plot lines of chosen threshold values in local/process:
        ylimits=ylim;
        line([chosen_th(actual_image_num) chosen_th(actual_image_num)],[ylimits(1) ylimits(2)],...
            'LineWidth',2,'Color','r','LineStyle',':')

        title({['Plot of \bf' method_name{method_num} '\rm scaling exponents \ita \rmand \itb \rmvs. threshold value'];...
        ['\itSubject ' num2str(actual_image_num)]})
        xlabel('Threshold Value'); ylabel('Scaling Exponents \ita \rmand \itb')
        legend('\ita \rm(ALL)','\itb \rm(ALL)','\ita \rm(\beta < 1)','\itb \rm(\gamma < 1)',...
            'Chosen Threshold','Location','NorthWest')  
        hold off
        
        % Increase figure height for PDF printing:
        set(gcf,'PaperPositionMode','manual');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperPosition', [0.25 1 8 9]);
    end
    
    for method_num=5
        
        figure; plot(th_name,a_q_avg_th,'Color','b','LineStyle',':'); hold on
        for i=1:length(a_q_avg_th)
            text(th_name(i),a_q_avg_th(i), num2str(a_count_1rem(i)))
        end

        plot(th_name,b_s_avg_th,'Color','r','LineStyle',':')
        for i=1:length(b_s_avg_th)
            text(th_name(i),b_s_avg_th(i), num2str(b_count_1rem(i)))
        end
        
        % Plot lines of chosen threshold values in local/process:
        ylimits=ylim;
        line([chosen_th(actual_image_num) chosen_th(actual_image_num)],[ylimits(1) ylimits(2)],...
            'LineWidth',2,'Color','r','LineStyle',':')

        title({['Plot of \bf' method_name{method_num} '\rm scaling exponents \ita \rmand \itb \rmvs. threshold value'];...
        ['\itSubject ' num2str(actual_image_num)]})
        xlabel('Threshold Value'); ylabel('Scaling Exponents \ita \rmand \itb')
        legend('\ita \rm(\beta < 1)','\itb \rm(\gamma < 1)',...
            'Chosen Threshold','Location','NorthWest')  
        hold off
        
        % Increase figure height for PDF printing:
        set(gcf,'PaperPositionMode','manual');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperPosition', [0.25 1 8 9]);
     
    end
end

%% SAVE FIGURES TO PDF

psname='ThresholdAnalysisFiveMethods.ps';
print('-dpsc2',psname,'-f1' )
ag=findobj;
nf=max(ag(find(ag==fix(ag))));
for i=2:nf
    print ('-dpsc2',psname,'-append',['-f' num2str(i)])
end 
winopen(psname);