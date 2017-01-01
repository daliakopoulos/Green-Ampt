function [ntimes, times, ttp, ttpp, rint, prec, bf, fpp, f, stor, ro] = ...
    green_ampt_v2(deltim, vsatk, sav, wcsat, wcini, stmax, rain_vector)

% The program follows the method of Chu (1978), Mein and Larson (1971, 1973)
% and Skaggs and Kaheel (1982) to calculate infiltration for unsteady 
% rainfall using the Green-Ampt equation.
% source: http://abe.ufl.edu/carpena/software/wingampt.shtml

% Green and  Ampt

% deltim (h): saturation time step
% vsatk (cm/hour): saturated hydraulic conductivity 
% sav (cm): Average sunction at the wetting front
% wcsat (cm^3/cm^3): saturated water content
% wcini (cm^3/cm^3): initial water content at start of the storm
% stmax (cm): Soil surface storage (detention storage)
% rain_vector: 3XN start_time (h), end_time, rainfall (mm)


nrain = size(rain_vector,1);
sttime = rain_vector(:,1);
endtim = rain_vector(:,2);
rawrfi = rain_vector(:,3);

bm = wcsat - wcini;

ntimes=1;
times(ntimes)=sttime(1);
bf(ntimes)=0.0;
f(ntimes)=0.0;
stor(ntimes)=0.0;
ro(ntimes)=0.0;
prec(ntimes)=0.0;
rint(ntimes)=rawrfi(1);
ttp(ntimes) = 0.0;      
ttpp(ntimes)=0.0;
fpp(ntimes)=0.0;
ipond=0;
fp = 0;
tp = 0;                 %time to ponding
tpp = 0; 
tnp = 0;

for	jj=1:nrain
      disp(['Timestep:', num2str(jj)]);
	  tstart = sttime(jj); 
	  tend = endtim(jj);
	  %dper = tend - tstart;
	  rfi = rawrfi(jj);
	  %train = rfi * dper;

%*******************************************************
%*   Is there ponding in this period, basically,
%* at the start of the period we need to check two
%* conditions:  
%*   1)  Not ponded at the start of the period, then
%*       a) continue not ponded 
%*       b) become ponded during the period
%*   2)  if the period starts with ponding, then 
%*       a) ponded condition for the entire period.
%*       b) ponding ceases during the period
%*******************************************************
%
%*************
%* Condition 1
%*************
    if ipond <= 0
        
        %*********************************
        %*  find time to ponding         *
        %*  bfp: potential infiltration rate
        %*********************************
        if fp > 0 && rfi >= fp  
            tp = tstart;
            bfp = sav * bm*vsatk/(fp-vsatk);
        elseif rfi > vsatk
            bfp = sav * bm /((rfi/vsatk)-1.0);
            tp = bfp / rfi;
            tp = tp + tstart;
        else
            bfp = 0;
            tp = 9999;
        end
        
        %*********************************
        %*   find tpp                    *
        %*   tpp: time shift
        %*********************************
 
        tpp = (bfp - sav * bm * log(1.0+bfp/(sav*bm)))/vsatk;
  
        if fp > 0 && rfi < fp 
            tpp = tpp + tstart;
        end
        
        if tp > tend
        %************************************
        %* a. no ponding during this period *
        %************************************


            tp = 9999;	
            tpp = 9999;
            
            [vsatk, sav, bm, deltim, ...
                   ntimes, tp,tpp,fp,...
                   ttp, ttpp,fpp,...
                   bf,f, stor, ro,prec,rint,times] = ...
            nopond_v1(tstart,tend,rfi,ntimes, ...
                   vsatk, sav, bm, deltim, ...
                   tp, tpp, fp, ...
                   ttp, ttpp,fpp, ...
                   bf,f, stor, ro,prec,rint,times);
            ipond=0; 
        else
        
            %**********************************
            %* b. ponding at tp               *
            %**********************************
 
            [vsatk, sav, bm, deltim, ...
                   ntimes, tp,tpp,fp,...
                   ttp, ttpp,fpp,...
                   bf,f, stor, ro,prec,rint,times] = ...
             nopond_v1(tstart,tp,rfi,ntimes, ...
                   vsatk, sav, bm, deltim, ...
                   tp, tpp, fp, ...
                   ttp, ttpp,fpp, ...
                   bf,f, stor, ro,prec,rint,times);
 
            %********************************
            %* ponded from tp on to tend    *
            %********************************
 
            [vsatk, sav, wcsat, wcini, bm, deltim, stmax, ...
                   ntimes, ttp,ttpp,fpp, ...
                   bf, f, stor, ...
                   ro, prec, rint, times] = ...
            pndinf_v1(tp,tend,rfi,tp,tpp,fp,ntimes, ...
                   vsatk, sav, wcsat, wcini, bm, deltim, stmax, ...
                   ttp,ttpp,fpp, ...
                   bf, f, stor, ...
                   ro, prec, rint, times);
            ipond=1;
            
           
        end
    else
        %***************
        %* Condition 2 *
        %***************
        %***********************************
        %*  find time to infiltrate fnp    *
        %***********************************
        frate=f(ntimes);
        if rfi < frate
            amtinf= bf(ntimes)+stor(ntimes);
            
            [tnp] = ...
                newtnp_v1(tstart,tend,tnp,tp,tpp,rfi,amtinf, ...
            vsatk, sav, wcsat, wcini, bm, deltim, stmax);
        else 
            %*****************************************
            %*  will not loose ponding, set tnp>tend *
            %*****************************************
            tnp=tend+1.0;
        end
        
        if tnp > tend
            %********************************
            %* 2a. ponding for whole period *
            %********************************
            [vsatk, sav, wcsat, wcini, bm, deltim, stmax, ...
                   ntimes, ttp,ttpp,fpp, ...
                   bf, f, stor, ...
                   ro, prec, rint, times] = ...
            pndinf_v1(tstart,tend,rfi,tp,tpp,fp,ntimes, ...
                   vsatk, sav, wcsat, wcini, bm, deltim, stmax, ...
                   ttp,ttpp,fpp, ...
                   bf, f, stor, ...
                   ro, prec, rint, times);
            ipond=1;
        else 
            %********************************
            %* 2b. ponding ends at tnp      *
            %********************************

            %********************
            %*   ponded portion *
            %********************
            [vsatk, sav, wcsat, wcini, bm, deltim, stmax, ...
                   ntimes, ttp,ttpp,fpp, ...
                   bf, f, stor, ...
                   ro, prec, rint, times] = ...
            pndinf_v1(tstart,tnp,rfi,tp,tpp,fp,ntimes, ...
                   vsatk, sav, wcsat, wcini, bm, deltim, stmax, ...
                   ttp,ttpp,fpp, ...
                   bf, f, stor, ...
                   ro, prec, rint, times);
             %********************
             %*  no pond portion *
             %********************
	         [vsatk, sav, bm, deltim,  ...
                   ntimes, tp,tpp,fp,...
                   ttp, ttpp,fpp,...
                   bf,f, stor, ro,prec,rint,times] = ...
             nopond_v1(tnp,tend,rfi,ntimes, ...
                   vsatk, sav, bm, deltim, ...
                   tp, tpp, fp, ...
                   ttp, ttpp,fpp, ...
                   bf,f, stor, ro,prec,rint,times);

	          ipond=0;
        end
    end
end

%**** infiltrate any water remaining in storage
if ipond > 0
   tstart = tend;
   tend = 1000.0;
   amtinf= bf(ntimes)+stor(ntimes);
   rfi=0.0;
 
   [tnp]= ...
      newtnp_v1(tstart,tend,tnp,tp,tpp,rfi,amtinf, ...
   vsatk, sav, wcsat, wcini, bm, deltim, stmax);

   [vsatk, sav, wcsat, wcini, bm, deltim, stmax, ...
                   ntimes, ttp,ttpp,fpp, ...
                   bf, f, stor, ...
                   ro, prec, rint, times] = ...
    pndinf_v1(tstart,tnp,rfi,tp,tpp,fp,ntimes, ...
                   vsatk, sav, wcsat, wcini, bm, deltim, stmax, ...
                   ttp,ttpp,fpp, ...
                   bf, f, stor, ...
                   ro, prec, rint, times);
   
   
end
