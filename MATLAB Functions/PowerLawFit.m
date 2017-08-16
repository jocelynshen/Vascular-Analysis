function [a_PL,b_PL] = PowerLawFit(rad,len,bad_ind,binning_type)
%This function calculates the scaling exponents a and b using the power law
%fit (distribution) method. Select binning type '1' for linear binning
%and binning type '2' for logarithmic binning. For the max

% Remove unwanted values from data:
rad(bad_ind)=[];
len(bad_ind)=[];

% BIN SIZE DETERMINATION

n_rad=ceil(sqrt(length(rad)));        % Square-root choice (apparently used by many programs)
%n_rad=ceil(log2(length(rad))+1)      % Sturges' Formula
%n_rad=ceil(2*length(rad)^(1/3));     % Rice Rule

n_len=ceil(sqrt(length(len)));
%n_len=ceil(log2(length(len))+1);
%n_len=ceil(2*length(len)^(1/3));    

switch binning_type
    
    case 1  % Linear binning
        
        % FOR a

        [count_n,xcenters]=hist(rad,n_rad);
        empty_bins=find(count_n==0);
        count_n(empty_bins)=[];
        xcenters(empty_bins)=[];
        
        relativefreq=count_n/length(rad);
        
        if length(xcenters)<=1
            a_PL=NaN;
        else
            x_log=log(xcenters);        % Transform to log-log space
            y_log=log(relativefreq);    
            p=polyfit(x_log,y_log,1);   % Use linear regression to determine coefficients
            a_PL=-1/p(1);               % Calculate a from power law
        end

        % For b

        [count_n,xcenters]=hist(len,n_len);
        empty_bins=find(count_n==0);
        count_n(empty_bins)=[];
        xcenters(empty_bins)=[];
        
        relativefreq=count_n/length(len);

        if length(xcenters)<=1
            b_PL=NaN;
        else
            x_log=log(xcenters);        % Transform to log-log space
            y_log=log(relativefreq);    
            p=polyfit(x_log,y_log,1);   % Use linear regression to determine coefficients
            b_PL=-1/p(1);               % Calculate b from power law
        end  
        
    case 2  % Log binning
        
        % Transform rad and len to log:
        rad=log(rad);
        len=log(len);
        
        % FOR a
        [count_n,xcenters]=hist(rad,n_rad);
        empty_bins=find(count_n==0);
        count_n(empty_bins)=[];
        xcenters(empty_bins)=[];
        
        xcenters=exp(xcenters);
        
        relativefreq=count_n/length(rad);
        
        if length(xcenters)<=1
            a_PL=NaN;
        else
            x_log=log(xcenters);        % Transform to log-log space
            y_log=log(relativefreq);    
            p=polyfit(x_log,y_log,1);   % Use linear regression to determine coefficients
            a_PL=-1/p(1);               % Calculate a from power law
        end

        % For b

        [count_n,xcenters]=hist(len,n_len);
        empty_bins=find(count_n==0);
        count_n(empty_bins)=[];
        xcenters(empty_bins)=[];
        
        xcenters=exp(xcenters);
        
        relativefreq=count_n/length(len);

        if length(xcenters)<=1
            b_PL=NaN;
        else
            x_log=log(xcenters);        % Transform to log-log space
            y_log=log(relativefreq);    
            p=polyfit(x_log,y_log,1);   % Use linear regression to determine coefficients
            b_PL=-1/p(1);               % Calculate b from power law
        end  
end
end

