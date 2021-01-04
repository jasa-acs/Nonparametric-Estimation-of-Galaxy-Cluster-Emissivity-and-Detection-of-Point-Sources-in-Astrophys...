function f=generatef(M,sf,type)

%FUNCTION TO GENERATE SIMULATED PROFILE
	if strcmp(type,'truef')  %,[1000 10 0.55]
        load truef.mat
        cutf=900;
        truef=truef(1:cutf);
        x=linspace(0,1,(M/2));
        sf3=sf(3);
        sf2=sf(2);
        sf=sf(1);
        if sf2 ==1 
            xf=linspace(0,1,cutf);
        else
            xf=log(linspace(1,sf2,cutf))/log(sf2);
        end
        f2=interp1(xf,truef,x,'linear','extrap');
        f2=f2.^sf3;
        %minf2=min(f2)
        %f2=f2-min(f2);
        mu0=[fliplr(f2) f2]/sf;
        %mu0=(mu0-min(mu0))/sf3+min(mu0);
        f=mu0;
    elseif strcmp(type,'f4') %[6000 10 0.5]
        load f4.mat
        xf=floor(linspace(1,512,M));
        f2=f4(xf).^sf(3);
        f=f2/sf(1);
        
    elseif strcmp(type,'f3') %[6000 10 0.5]
        load f3.mat
        xf=floor(linspace(1,512,M));
        f2=f3(xf).^sf(3);
        f=f2/sf(1);
    elseif strcmp(type,'blocks')
        x=linspace(0,1,M);
        mu0=waveblocks(x,1);
        mu0=mu0-(min(mu0));
        f=mu0*sf;
    elseif strcmp(type,'f3sym') %[6000 10 0.5]
        load f3.mat
        xf=floor(linspace(1,512,M));
        f2=f3(xf).^sf(3);
        f=f2/sf(1);
        f=f((M/2+1):end);
        f=[fliplr(f) f];
    elseif strcmp(type,'f4sym') %[6000 10 0.5]
        load f4.mat
        xf=floor(linspace(1,512,M));
        f2=f4(xf).^sf(3);
        f=f2/sf(1);
        f=f((M/2+1):end);
        f=[fliplr(f) f];
    elseif strcmp(type,'blockssym')
        x=linspace(0,1,M);
        mu0=waveblocks(x,1);
        mu0=mu0-(min(mu0));
        f=mu0*sf;
        f=f((M/2+1):end);
        f=[fliplr(f) f];
    else
        if strcmp(type,'x')
            x=linspace(0,1,M);
            mu0=100*(500000*(x+0.1).^10.*exp(-20*(x+0.1))+(x+0.1).*exp(-10*(x+0.1)));
        elseif strcmp(type,'ex')
            x=linspace(0,1,M);
            mu0=100*(500000*x.^10.*exp(-20*x)+x.*exp(-10*x));
            mu0=exp(mu0);
        elseif strcmp(type,'rho')
            r=linspace(1/(M/2),1,(M/2));
            rho_s=1e-3;
            mu0=rho_s./(r.*(1+r).^2);
            mu0=[fliplr(mu0) mu0];
        elseif strcmp(type,'null')
            mu0=ones(1,M);
        end

        mu0=mu0./sum(mu0.^2); 
        
        if sf==Inf
            f=mu0*0;
        else
            f=mu0*sf;
        end
	end
end



function y = waveblocks(x, snr)

    pos = [0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81];
    hgt = [4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2];
    y = zeros(1, length(x));
    for(j = 1:length(pos))
      y = y + ((1 + sign(x - pos(j))) * hgt(j))/2;
    end
    y = (snr * y)/1.914;

end