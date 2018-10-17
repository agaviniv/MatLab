function err = etafitJ(Jw,t,y,handle)

    %J1=3.8e-9; J2=3.4e-9; 
    %______________________________________
    %%%%%%%%%%%%%  Redfield  %%%%%%%%%%%%%%
    Dro=zeros(4,4); 
    
    Iz = 0.5*diag([3 1 -1 -3]);
    C = 1.2e10;
    roeq = (2*Iz+3*eye(4))/12;
    ro = 0.25*eye(4);    

    R01 = trace(ro - roeq);
    R02 = sum(diag(ro - roeq)'*[-1; 1; 1; -1]);
    R03 = sum(diag(ro - roeq)'*[1; 1; -1; -1]);
    R04 = sum(diag(ro - roeq)'*[1; -1; 1; -1]);
    
    np = size(t,2);
    etaf=zeros(1,np);
    etaf(1)=trace((ro-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
    dt = 1e-3;
    for kk = 1:np-1
        Rv = [R01 R02*exp(-2*C*(Jw(1)+Jw(2))*kk*dt) R03*exp(-2*C*Jw(2)*kk*dt) R04*exp(-2*C*Jw(1)*kk*dt)];
        drot(1) = roeq(1,1)+(1/4)*(Rv*[1;-1;1;1]);
        drot(2) = roeq(2,2)+(1/4)*(Rv*[1;1;1;-1]);
        drot(3) = roeq(3,3)+(1/4)*(Rv*[1;1;-1;1]);
        drot(4) = roeq(4,4)+(1/4)*(Rv*[1;-1;-1;-1]);
        rot = diag(drot);

        etaf(kk+1)=trace((rot-0.25*eye(4))^2)/trace((0.25*eye(4))^2);
    end

    err = norm(etaf-y);

    set(gcf,'DoubleBuffer','on');
    set(handle,'ydata',etaf)
    drawnow
    pause(.04)
end