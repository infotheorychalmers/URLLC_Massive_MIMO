function R = functionRlocalscatteringApproxUCA(M,theta,ASDdeg)

%This implemenation is based on Section III.B in:
%"Simplified Spatial Correlation Models for
%  Clustered MIMO Channels With Different Array Configurations"

%M              = Number of antennas
%theta          = Nominal angle [radians]
%ASDdeg         = power azimuth spectrum standard deviation [degrees]

ASDdeg = deg2rad(ASDdeg);
type = 'UCA';

if M==1
    R=1;
    return;
end

if strcmp(type, 'UCA')
    
    B=nan(M,M);
    phis = deg2rad(360/M);
    beta = 1-(1-exp(-sqrt(2*pi)/ASDdeg));
    krho = pi/sqrt(2*(1-cos(phis)));
    for m=1:M
        for n=1:M
            den = 1+ (ASDdeg^2/2) * (krho*(sin(theta-(m-1)*phis)-sin(theta-(n-1)*phis)))^2;
            B(m,n) = beta / den;
        end
    end
    a_cua = exp(1i*krho*cos(theta-(0:(M-1))*phis )).';
    
    R = (a_cua*a_cua') .* B;

    % toy example of radius
    %lambda = 3*1e8/(60e9);
    %rho = (lambda/2)/sqrt(2*(1-cos(phis)))
    
elseif strcmp(type,'ULA')
    
    B=nan(M,M);
    beta = 1-(1-exp(-sqrt(2*pi)/ASDdeg));
    kd = pi;
    for m=1:M
        for n=1:M
            den = 1 + (ASDdeg^2/2) * (kd*((m-1)-(n-1)) * cos(theta))^2;
            B(m,n) = beta / den;
        end
    end  
    a_ula = exp(1i*kd*sin(theta)*(0:M-1)).';
    
    R = (a_ula*a_ula') .*B;
end
R= M * (R./trace(R));

end