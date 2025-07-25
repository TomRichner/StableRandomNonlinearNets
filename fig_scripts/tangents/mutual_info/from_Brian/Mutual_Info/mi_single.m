function mi4 = mi_single(single_isi)

% INT CODE
% mutual information calculations
% 09/05
%
% based on np050911.m

% Find mutual information Im(del|sigma) that a single interval gives about
% the stimulus variance

% see mi.m to obtain single_isi
mn = -3;
step = 0.05;
mx = 0;
dist = 'log10';
for i = 1:length(single_isi)
    [sx{i} nsx{i} bins] = isi_to_px(single_isi{i},mn,step,mx,dist,1);
end

% calculate MI after a switch between sigma_est and sigma

% find maximum P for each delta
% here: px2 = P(del|sigma) and tmp = P(sig_est|del)
for m = 1:length(sx)
    clear tmp tmp2
    [Y I] = max(sx{m},[],2);
    tmp = zeros(size(sx{m},1),size(sx{m},2));
    for i = 1:size(sx{m},1)
        tmp(i,I(i)) = 1;
    end
    %tmp' = Psigest_delta
    %sx{m} = Pdelta_sig
    Psigest_sig = tmp'*sx{m};
    
    Pdel = sum(nsx{m},2)/sum(sum(nsx{m}));
    PxSigma = sum(nsx{m}) / sum(sum(nsx{m}));
    Psigest = Psigest_sig * PxSigma';
    
    clear tmp3
    % MI (sigest | sig)
    for i = 1:size(Psigest_sig,2)
        tmp3(:,i) = PxSigma(i)*Psigest_sig(:,i).*log2(Psigest_sig(:,i)/Psigest(i));
    end
    tmp3(isnan(tmp3))=0;
    mi3(m) = sum(sum(tmp3));
    
    clear tmp4
    % MI (del | sig)
    for i = 1:size(sx{m},2)
        tmp4(:,i) = PxSigma(i)*sx{m}(:,i).*log2( sx{m}(:,i) ./ Pdel );
    end
    tmp4(isnan(tmp4))=0;
    mi4(m) = sum(sum(tmp4));
end
