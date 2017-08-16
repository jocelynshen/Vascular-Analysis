function [a_HA,b_HA] = HierarchicalAveraging(rad,beta,len,gamma,n_child,ind_NA,bad_ind,binning_type)
%This function calculates the scaling exponents a and b using the
%hierarchical averaging method. Select binning type '1' for linear binning
%and binning type '2' for logarithmic binning.

% Remove unwanted values from data:
beta([ind_NA bad_ind])=[];
gamma([ind_NA bad_ind])=[];
rad([ind_NA bad_ind])=[];
len([ind_NA bad_ind])=[];

if length(rad)<1 || length(len)<1
    
    a_HA=NaN;
    b_HA=NaN;
    
else

    % BIN SIZE DETERMINATION

    n_bins_a=ceil(sqrt(length(rad)));        % Square-root choice (apparently used by many programs)
    %n_bins_a=ceil(log2(length(rad))+1);     % Sturges' Formula
    %n_bins_a=ceil(2*length(rad)^(1/3));     % Rice Rule

    n_bins_b=ceil(sqrt(length(len)));
    %n_bins_b=ceil(log2(length(len))+1);
    %n_bins_b=ceil(2*length(len)^(1/3));  

    switch binning_type

        case 1  % Linear binning

            % FOR a

            rad_a=[rad -log(beta)/log(n_child)];    % Calculate 'a' from radius data
            rad_a=sortrows(rad_a);                  % Sort 'a' by smallest to largest radius size
            count_n=hist(rad_a(:,1),n_bins_a);      % Place 'a' values in bins according to radius

            count_n(find(count_n==0))=[];           % Remove empty bins

            ari_means_a=mean(rad_a(1:count_n(1),2));
            for i=2:length(count_n)
                bin_range=(sum(count_n(1:(i-1)))+1):sum(count_n(1:i));
                ari_means_a=[ari_means_a mean(rad_a(bin_range,2))];
            end
            a_HA=mean(ari_means_a);       % Calculate arith. mean of bin means

            % FOR b

            len_b=[len -log(gamma)/log(n_child)];
            len_b=sortrows(len_b);
            count_n=hist(len_b(:,1),n_bins_b);

            count_n(find(count_n==0))=[];

            ari_means_b=mean(len_b(1:count_n(1),2));
            for i=2:length(count_n)
                bin_range=(sum(count_n(1:(i-1)))+1):sum(count_n(1:i));
                ari_means_b=[ari_means_b mean(len_b(bin_range,2))];
            end
            b_HA=mean(ari_means_b);

        case 2  % Log binning

            % FOR a

            rad_a=[log(rad) -log(beta)/log(n_child)];   % Calculate 'a' from radius data
            rad_a=sortrows(rad_a);                      % Sort 'a' by smallest to largest log(rad)
            count_n=hist(rad_a(:,1),n_bins_a);          % Place 'a' values in bins according to rad

            count_n(find(count_n==0))=[];           % Remove empty bins

            ari_means_a=mean(rad_a(1:count_n(1),2));
            for i=2:length(count_n)
                bin_range=(sum(count_n(1:(i-1)))+1):sum(count_n(1:i));
                ari_means_a=[ari_means_a mean(rad_a(bin_range,2))];
            end
            a_HA=mean(ari_means_a);       % Calculate arith. mean of bin means

            % FOR b

            len_b=[log(len) -log(gamma)/log(n_child)];
            len_b=sortrows(len_b);
            count_n=hist(len_b(:,1),n_bins_b);

            count_n(find(count_n==0))=[];

            ari_means_b=mean(len_b(1:count_n(1),2));
            for i=2:length(count_n)
                bin_range=(sum(count_n(1:(i-1)))+1):sum(count_n(1:i));
                ari_means_b=[ari_means_b mean(len_b(bin_range,2))];
            end
            b_HA=mean(ari_means_b);
    end
end

