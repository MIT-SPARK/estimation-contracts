function bound = computeAposterioriBoundWahba_v2(J, At, barc)

A = zeros(3*length(J),9);
for k=1:length(J)
    A(blkIndices(k,3),:) = At( :,:,J(k) );
end
min_sigma_J = min(svd(A));

dim_J = length(J);
bound = 2 * barc * sqrt(dim_J) / min_sigma_J;

