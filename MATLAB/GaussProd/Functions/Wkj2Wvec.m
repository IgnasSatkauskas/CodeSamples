function W_vec = Wkj2Wvec(W_kj, k_start, j_start)

% Make vector representation W_kj to W_vec (dropp coefs below eps)
lin_ind = find(abs(W_kj) > eps);
[j_i, k_i] = ind2sub(size(W_kj),lin_ind);

W_vec = zeros(length(lin_ind),3);
W_vec(:,1) = j_start + (j_i-1);
W_vec(:,2) = k_start + (k_i -1);
W_vec(:,3) = W_kj(lin_ind);