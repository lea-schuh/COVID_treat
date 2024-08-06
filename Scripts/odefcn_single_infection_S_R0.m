function dydt = odefcn_single_infection_S_R0(t,y,b0,dI,pV,dV,pN,dN,pB,dB,B_thres,alpha,treat,S0)

    dydt = zeros(4,1);
    
    dydt(1) =  pN - dN*y(1) - b0*y(3)*y(1)*(1-y(4)) ; %S_n
    dydt(2) = b0*y(3)*y(1)*(1-y(4)) - dI*y(2); %I_n

    if t < treat
        dydt(3) = pV*y(2) - dV*y(3) - b0*y(3)*y(1)*(1-y(4)); %V_n
        dydt(4) = pB*y(3)*(1-y(4)) - dB*(y(4)-B_thres)*y(4); %B - immune response
    elseif treat <= t && t <= treat+5
        % B_thres2 = 1-dI*dV/(b0*S0*(alpha*pV-dI));
        dydt(3) = pV*alpha*y(2) - dV*y(3) - b0*y(3)*y(1)*(1-y(4)); %V_n
        dydt(4) = pB*y(3)*(1-y(4)) - dB*(y(4)-B_thres)*y(4); %B - immune response
        % dydt(4) = pB*y(3)*(1-y(4)) - dB*(y(4)-B_thres2)*y(4); %B - immune response
    else
        dydt(3) = pV*y(2) - dV*y(3) - b0*y(3)*y(1)*(1-y(4)); %V_n
        dydt(4) = pB*y(3)*(1-y(4)) - dB*(y(4)-B_thres)*y(4); %B - immune response
        % B_thres3 = 1-dI*dV/(b0*S0*((1-(1-alpha)*exp(-2*(t-(treat+5))))*pV-dI));
        % dydt(3) = (1-(1-alpha)*exp(-2*(t-(treat+5))))*pV*y(2) - dV*y(3) - b0*y(3)*y(1)*(1-y(4)); %V_n
        % dydt(4) = pB*y(3)*(1-y(4)) - dB*(y(4)-B_thres3)*y(4); %B - immune response
    end
    
end