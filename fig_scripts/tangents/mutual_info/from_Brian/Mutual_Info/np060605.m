rand('state',1);    %sum(100*clock));
(length(Stim)-1)/10;
ups = floor(.5*SampleRate);
e = ceil(StimLength./.5);
tmp = randn(1,e+1);
%tmp = [1.9164 -0.1053 -0.2167 0.2979 -1.4508 1.4492 0.5734 -1.6486 ... 
    %-0.3790 0.7585 -0.8167 0.1587 0.2015 0.3809];
E = resample(tmp,ups,1);
E = E/std(E);
%E = E*sqrt(max(stim_var))/max(abs(E)); %scale E by max stim std
E = 35*E+65;
%}