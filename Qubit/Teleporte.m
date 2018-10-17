% Programa para Estudo do Protocolo de Teleporte
%           3 spins 1/2
%

clear all

Ha=sqrt(0.5)*kron([1 1;1 -1],eye(4));
Hb=sqrt(0.5)*kron(eye(2),kron([1 1;1 -1],eye(2)));
Hc=sqrt(0.5)*kron(eye(4),[1 1;1 -1]);
Xab=kron([1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0],eye(2));
Xbc=kron(eye(2),[1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0]);
Xac=eye(8); w=Xac(5,:);Xac(5,:)=Xac(6,:); Xac(6,:)=w; 
            w=Xac(7,:);Xac(7,:)=Xac(8,:); Xac(8,:)=w;

ro = cell(1,12);            
ro{1}=kron((0.1)*[1 0;0 0],diag([1 0 0 0]));
% Teleporte
%UT = Hc * Xac * Hc * Xbc * Ha * Xab * Xbc * Hb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% APLICACAO DE U - 1a Parte do Protocolo

% Hadamard em b
ro{2} = Hb * ro{1} * Hb';

% Controled-NOT com bit de Controle em b e Alvo em c
ro{3} = Xbc * ro{2} * Xbc';

% Controled-NOT com bit de Controle em a e Alvo em b
ro{4} = Xab * ro{3} * Xab';

% Hadamard em a
ro{5} = Ha * ro{4} * Ha';


% APLICACAO DE U - 2a Parte do Protocolo

% Controled-NOT com bit de controle em b e Alvo em c
ro{6} = Xbc * ro{5} * Xbc';

% Hadamard em c
ro{7} = Hc * ro{6} * Hc';

% Controled-NOT com bit de controle em a e Alvo em c
ro{8} = Xac * ro{7} * Xac';

% Hadamard em c
ro{9} = Hc * ro{8} * Hc';

for k = 1:9
    subplot(3,3,k);
    ro{10}=wig2ro(figwig(ro{k}));
    axis off
end

figure; subplot(1,3,1);ro{11}=wig2ro(figwig(atenufase(0.5,'a',ro{10})));
subplot(1,3,2); ro{12}=wig2ro(figwig(atenufase(0.5,'b',ro{11})));
subplot(1,3,3); rof=wig2ro(figwig((0.5)*kron(0.25*eye(4),[1 1;1 1])));

trace(sqrt(sqrt(ro{12}')*rof*sqrt(ro{12})))
