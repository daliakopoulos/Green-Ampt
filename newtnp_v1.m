function [tnp]= ...
          newtnp_v1(tst,tend,tnp,tp,tpp,rfi,amtinf, ...
                vsatk, sav, wcsat, wcini, bm, deltim, stmax)
accpt = 0.1e-10;

iter = 0;
tnp = tend;
bftry = (tnp-tst) * rfi + amtinf;
arglog = 1.0 + bftry/(bm*sav);
hnp = bftry - bm*sav*log(arglog) -  vsatk*(tnp -tp +tpp); 
error = 1;

while error > accpt
    iter = iter+1;
	tnpold = tnp;
	if iter > 1000
	   disp (['it = ', num2str(iter), ' dhnp = ', num2str(dhnp), ' hnp = ', num2str(hnp) ,' tnp = ',num2str(tnp), ' error = ', num2str(error)]);
       break;
    end

	dhnp= rfi - rfi*(1.0/arglog) - vsatk;

	tnp = tnpold - hnp/dhnp;
       
%** if tnp is negative  - fix added jep, 6/16/05
    if tnp < 0 
        tnp = 0;
    end

	bftry = (tnp-tst)*rfi + amtinf;
	arglog = 1.0 + bftry/(bm*sav);

    disp (['tnp = ', num2str(tnp), ' tst = ', num2str(tst), ' rfi = ', ...
        num2str(rfi) ,' amtinf = ', num2str(amtinf), ' bftry = ', ...
        num2str(bftry),' bm = ', num2str(bm),' sav = ', num2str(sav)]);
	
    hnp = bftry - bm*sav*log(arglog) -  vsatk*(tnp -tp +tpp);
	error = abs(hnp);
end