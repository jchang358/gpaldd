function [PSS yobs]= td_gen_y(TLL, ASS, k,e)
    DLL = 800.*(exp(k)*TLL+1).^-1;
    DSS = ASS;
    PSS = (1+exp(-e*(DLL-DSS))).^-1;
    yobs = binornd(1,PSS);
end