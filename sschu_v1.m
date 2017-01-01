function [bff] = sschu_v1(tim, tp, tpp, bff, vsatk, sav, bm)
%****************
%* use Newton's method on chu's equation
%***************
%*
%***************
%*  set up problem
%***************
iter = 0;
accpt = 0.1e-04; % threshold
hh = bff - bm * sav * log (1.0 + (bff/(bm*sav))) - vsatk * (tim-tp+tpp);
error = 1;

while error > accpt
    iter = iter+1;
	if iter > 200
       disp (['it = ', num2str(iter), ' dhdf = ', num2str(dhdf), ' hh = ', ...
            num2str(hh) ,' bff = ', num2str(bff), ' error = ', num2str(error)]);
       break;
    end
    dhdf = 1.0 - ((bm * sav) /(bm*sav + bff));
    bbfnew = bff - hh/dhdf ;
    bff = bbfnew;
    hh = bff - bm * sav * log (1.0 + (bff/(bm*sav))) - vsatk * (tim-tp+tpp);
    error = abs(hh);
end