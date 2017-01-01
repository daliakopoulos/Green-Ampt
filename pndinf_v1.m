function [vsatk, sav, wcsat, wcini, bm, deltim, stmax, ...
                   ntimes, ttp,ttpp,fpp, ...
                   bf, f, stor, ro, prec, rint, times] = ...
          pndinf_v1(tstart,tend,rfi,tp,tpp,fp,ntimes, ...
                   vsatk, sav, wcsat, wcini, bm, deltim, stmax, ...
                   ttp,ttpp,fpp, ...
                   bf, f, stor, ro, prec, rint, times)

%**********************
%* find the number of time steps...
%*********************
nsteps = floor((tend-tstart)/deltim);
 
 


for kk = 1:nsteps
    ntimes=ntimes+1;
	times(ntimes)=times(ntimes-1)+deltim;
	water= rfi * deltim + stor(ntimes-1);
    
%            ***********************
%            * make a guess for bigf
%            ************************
    bbf = bf(ntimes-1)+water;
    
    
    
	ctime=times(ntimes);
    
	[bbf] = ...
        sschu_v1(ctime,tp,tpp,bbf,vsatk, sav, bm);


    
    delinf = bbf - bf(ntimes-1);
    
	bf(ntimes)= bbf;
	fp=vsatk + (vsatk*bm*sav/bf(ntimes));
	f(ntimes)=fp;
	
    if water > delinf 
        stor(ntimes)=water-delinf;
        if stor(ntimes) > stmax
            ro(ntimes)=ro(ntimes-1)+(stor(ntimes)-stmax);
            stor(ntimes)=stmax;
        else
            ro(ntimes)=ro(ntimes-1);
        end
    else
	    ro(ntimes)=ro(ntimes-1);
	    stor(ntimes)=0.0;
    end
        prec(ntimes) = prec(ntimes-1)+rfi*deltim;
	    rint(ntimes) = rfi;
	    fpp(ntimes) = fp;
        ttp(ntimes) = tp;
	    ttpp(ntimes) = tpp;
end
%*********************
%*  check that we really are at the end of the period
%*********************
dterr = tend - times(ntimes);
if dterr > 0
	ntimes=ntimes+1;
	times(ntimes)=tend;
	water= rfi * dterr+ stor(ntimes-1);
    
    %***********************
    %* make a guess for bigf
    %************************
    bbf = bf(ntimes-1)+water;
	ctime=times(ntimes);
	
    [bbf] = ...
        sschu_v1(ctime, tp, tpp, bbf, vsatk, sav, bm);
    
    delinf = bbf - bf(ntimes-1);
	bf(ntimes)= bbf;
	fp=vsatk + (vsatk*bm*sav/bf(ntimes));
	f(ntimes)=fp;
	
    if water > delinf
        stor(ntimes) = water-delinf;
        if stor(ntimes)> stmax
            ro(ntimes) = ro(ntimes-1)+(stor(ntimes)-stmax);
            stor(ntimes) = stmax;
        else
            ro(ntimes) = ro(ntimes-1);
        end
	else
	    ro(ntimes) = ro(ntimes-1);
	    stor(ntimes) = 0.0;
    end
    prec(ntimes) = prec(ntimes-1)+rfi*dterr;
    rint(ntimes) = rfi;
    fpp(ntimes) = fp;
    ttp(ntimes) = tp;
    ttpp(ntimes) = tpp;
end
 
	
