function p=ellfit(p0,maschera,psi_in,delta_in,theta0)
    
    % argumments:
    %
    %   1) parameters hint values
    %   2) model mask with NaN for variable parameters and fixed values
    %      everywhere else. Two numbers (real and imag) each index.
    %   3) experimental psi (rad)
    %   4) experimental delta (rad)
    %   5) experiment angles (rad)

    global iterazione
    iterazione=0;
    p=fminsearch('elld',p0,[],maschera,psi_in,delta_in,theta0);
