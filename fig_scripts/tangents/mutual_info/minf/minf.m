function y=minf(pab,papb)

I=find(papb(:)>1e-12 & pab(:)>1e-12); % function support 
y=pab(I).*log2(pab(I)./papb(I));
