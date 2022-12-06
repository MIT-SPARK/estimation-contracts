function bound = computeAposterioriBoundWahba(theta, At, barc, sizeJ)

if nargin < 4
    sizeJ = 3; % size of minimal set
end

min_sigma_At = inf;
I_MC = find(theta>0); % set of selected inliers
possible_J = nchoosek(I_MC,sizeJ);
A = zeros(3*sizeJ,9);
for i=1:size(possible_J,1) % for each possible subset of size sizeJ
    J = possible_J(i,:);
    for k=1:length(J)
        A(blkIndices(k,3),:) = At( :,:,J(k) );
    end
    min_sigma_J = min(svd(A));
    min_sigma_At = min(min_sigma_At , min_sigma_J);
end

bound = 2 * barc * sqrt(sizeJ) / min_sigma_At;

