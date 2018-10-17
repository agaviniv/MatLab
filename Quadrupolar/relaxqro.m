%clear all

% Populacoes do estado |00>
load DensidadePPS00W.mat;
mm = max(diag(UAB{end}));
for kk=1:40; 
    pp00(kk,:)=(diag((3/2)*UAB{kk}/mm));
end;

% Populacoes do estado |01>
load DensidadePPS01W.mat;
mm = max(diag(UAB{end}));
for kk=1:40;
    pp01(kk,:)=(diag((3/2)*UAB{kk}/mm));
end;

% Populacoes do estado |10>
load DensidadePPS10W.mat;
mm = max(diag(UAB{end}));
for kk=1:40; 
    pp10(kk,:)=(diag((3/2)*UAB{kk}/mm));
end;

% Populacoes do estado |11>
load DensidadePPS11W.mat;
mm = max(diag(UAB{end}));
for kk=1:40;    
    pp11(kk,:)=(diag((3/2)*UAB{kk}/mm));
end;

tempo = 1000*tempo;
plot(tempo,pp00,'o');hold on;
plot(tempo,pp01,'.'); 
ylim([-1.6 1.6]); xlabel('Time (ms)'); ylabel('Amplitude (a.u.)');
figure;
plot(tempo,pp10,'o');hold on
plot(tempo,pp11,'.');
ylim([-1.6 1.6]); xlabel('Time (ms)'); ylabel('Amplitude (a.u.)');


% Populacoes da matriz de superposicao
load DensidadeSUP.mat;
mm = max(diag(UAB{end}));
for kk=1:40    
    ppsu(kk,:)=diag((3/2)*UAB{kk}/mm);
end
% Coerencias da matriz de superposicao
load DensidadeSUPERW.mat;
coe = tab2(1:40,2:2:8);
tempo2 = 1000*tab2(:,1);
coeds = tab3(1:40,3:2:5);
tempods = 1000*tab3(:,2);
tempop = 1000*tab1(:,1);


figure;
plot(tempo2,coe,'.',tempods,coeds,'o');
xlabel('Time (ms)'); ylabel('Amplitude (a.u.)');
figure;
plot(tempop,ppsu(:,1),'^',tempop,ppsu(:,2),'o',tempop,ppsu(:,3),'d',tempop,ppsu(:,4),'v');
xlabel('Time (ms)'); ylabel('Amplitude (a.u.)');