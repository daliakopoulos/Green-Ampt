function [vsatk, sav, bm, deltim, ...
                   ntimes, tp,tpp,fp,...
                   ttp, ttpp,fpp,...
                   bf,f, stor, ro,prec,rint,times] = ...
          nopond_v1(tstart,tend,rfi,ntimes, ...
                   vsatk, sav, bm, deltim, ...
                   tp, tpp, fp, ...
                   ttp, ttpp,fpp, ...
                   bf,f, stor, ro,prec,rint,times)

%**********************
%c* find the number of time steps...
%*********************
 
nsteps = floor((tend-tstart)/deltim);
for kk=1:nsteps
    ntimes=ntimes+1;
    times(ntimes)=times(ntimes-1)+deltim;
    delinf = rfi * deltim;
    bf(ntimes)=bf(ntimes-1)+delinf;
    fp=vsatk + (vsatk*bm*sav/bf(ntimes));
    f(ntimes)=rfi;
    prec(ntimes)=prec(ntimes-1)+delinf;
    rint(ntimes)=rfi;
    fpp(ntimes)=fp;
    ttp(ntimes)=tp;
    ttpp(ntimes)=tpp;
    stor(ntimes)=0.0;
    ro(ntimes)=ro(ntimes-1);
end
 
 
%*********************
%*  check that we really are at the end of the period
%*********************
dterr = tend - times(ntimes);
if dterr > 0
    ntimes = ntimes+1;
    delinf = rfi * dterr;
    bf(ntimes)=bf(ntimes-1)+delinf;
    fp=vsatk + (vsatk*bm*sav/bf(ntimes));
    f(ntimes)=rfi;
    prec(ntimes)=prec(ntimes-1)+delinf;
    rint(ntimes)=rfi;
    fpp(ntimes)=fp;
    ttp(ntimes)=tp;
    ttpp(ntimes)=tpp;
    stor(ntimes)=0.0;
    ro(ntimes)=ro(ntimes-1);
end
	times(ntimes)=tend;
 

