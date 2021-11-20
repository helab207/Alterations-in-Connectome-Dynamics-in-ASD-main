function V = scaled_inclusivity( X, ith)
% Calcuating the community assignment consistency of each nodes among
% various network partitions X
% Input,
%         X, M-by-N, M paritions for networks with N nodes, with each row
%         stands for a partition
%         ith, the index of the reference partition
% Ref: Steen et al., PRE 84, 016111(2011)

if nargin<2
    error('Please input the reference partition');
end

[M, N] = size(X);
if ith>M
    error('Please input the right index of the reference partition');
end

V = zeros(N,1);            % N, number of nodes

Ref = X(ith,:); 
idx = 1:N;
MQ = sparse(idx, Ref, 1, N, max(Ref), N);

for ii = 1:M
    if (ii~=ith) 
       %disp(['ii = ', num2str(ii)]);
       A = X(ii,:);  
       MA = sparse(idx, A, 1, N, max(A), N);
       S = comm_comm_similarity(A, Ref);
       Ind = find(S);             % find p and q, with S(p,q)>0
       num =length(Ind);
       
       [P, Q] = ind2sub(size(S), Ind);
       for jj = 1:num
           Ind_com = MA(:, P(jj)).*MQ(:, Q(jj))>0;
           V(Ind_com) = V(Ind_com) + S(Ind(jj));
       end
    end
end
       
end

