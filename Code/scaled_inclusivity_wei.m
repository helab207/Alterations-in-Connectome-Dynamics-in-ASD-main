function [V_wei, V_vec, Deg_J] = scaled_inclusivity_wei(X, ss)
% Calculating the weighted scaled inclusivity based on the partition
% similarity among the given n partitions
% Inputs,
% X, M-by-N, M partitions for network of N nodes
% ss, type for calculate the similarity of community structure for
% two partitions; 1, normalized mutual information (default); 2, Jaccard similarity
% Ref: Steen PRE 84, 016111 (2011)

if nargin <2
    ss = 1;
end

[M, N] = size(X);   % M, partitions, and N nodes

% Step1, partition similarity among different partitions
J = zeros(M);     % similarity among M partitions

switch ss            % similarity type
    case(1)          % normalized mutual information
        for ii = 1:M
            for jj = ii+1:M
                J(ii,jj) = nmi(X(ii,:)', X(jj,:)');      % normalized mutual information
            end
        end
    case(2)
        for ii = 1:M
            for jj = ii+1:M
                J(ii,jj) = Jaccard_similarity(X(ii,:)', X(jj,:)');
            end
        end
    case(3)
        for ii = 1:M
            for jj = ii+1:M
                J(ii,jj) = ami(X(ii,:), X(jj,:));      % adjusted mutual information, Liao 20161009
                %J(ii,jj) = AMI_T(X(ii,:)', X(jj,:)');      % adjusted mutual information, Liao 20161009
            end
        end
end

J = J+J';
Deg_J = sum(J,2);
Deg_J = Deg_J/sum(Deg_J);             % sum(Deg) =1

% Step2, calculate the scaled vector for each given ref partition
V_wei = zeros(N,1);
V_vec = zeros(N,M);
for ith = 1:M
    V = scaled_inclusivity(X, ith);
    V_wei = V_wei + V*Deg_J(ith);
    V_vec(:,ith) = V;
end

V_wei = 1-V_wei/(M-1);      % Normalize the variability, and chaged to modular variability
end

