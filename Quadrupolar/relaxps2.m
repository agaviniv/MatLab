clear all

load DensidadeSUP.mat;
rosu= UAB; %rosua= UA; rosub= UB; 
mm = max(max(rosu{end}));
for kk=1:40
    rosu{kk}=.25*eye(4)+0.25*rosu{kk}/mm;
    pp(kk,:)=diag(rosu{kk});
    eta(kk)=4*trace((rosu{kk}-.25*eye(4))^2);
end
figure;
%subplot(1,2,1);
plot(1000*tempo,pp,'.'); ylim([-0.1 0.6]);
%subplot(1,2,2);plot(1000*tempo,eta,'.'); ylim([0 0.6]);
xlabel('Time (ms)'); ylabel('Amplitude (a.u.)');

clear UA UAB UB WA WB colu fila klm rho tab1 tab2 tab3 Tiempo
clear r11 r14 r22 r23 r33 r44

%%
% Fitting com canal de atenuacao de amplitude
t = tempo; 
figure;
plot(tempo,eta,'ro'); hold on; h = plot(tempo,eta,'b'); hold off;
title('Input data'); ylim([0 0.7])

p(1)=5/6;
p(2)=2/3;
start = [0.08;0.08]; %gama's
options = optimset('Simplex','on','TolX',0.001);
estimated_gama = fminsearch(@(x)etafit(x,t,eta,p,h),start,options)


% Plot com parametros do fitting
ga = estimated_gama(1); gb = estimated_gama(2);
roi=.25*eye(4);
roi=atenuamp(5/6,0.25*ga,'a',roi); %
    roi=atenuamp(2/3,0.25*gb,'b',roi);
for k=1:40
    
    roi=atenuamp(5/6,ga,'a',roi); %
    roi=atenuamp(2/3,gb,'b',roi); %
    eta(k)=4*trace((roi - 0.25*eye(4))^2); 
end
hold on;
plot(tempo,eta,'r');

%%
% Coerencias da matriz de superposicao
load DensidadeSUPERW.mat;

coe = tab2(1:40,2:2:8);
tempo2 = 1000*tab2(:,1);
coeds = tab3(1:40,3:2:5);
tempods = 1000*tab3(:,2);

figure;
plot(tempo2,coe,'.',tempods,coeds,'o')
legend('\rho_{01}','\rho_{23}','\rho_{02}','\rho_{13}','\rho_{12}','\rho_{03}')
xlabel('Time (ms)'); ylabel('Amplitude (a.u.)');



% ######################################################
% Estado do Gato - elementos da diagonal

load DensidadeGATOW.mat;
rogato= UAB; %rogatoa= UA; rogatob= UB;
mm = max(max(rogato{end}));
for kk=1:40
    rogato{kk}=.25*eye(4)+0.25*rogato{kk}/mm; 
    pp(kk,:)=(diag(rogato{kk}));
end
% figure;
% plot(pp,'.'); axis([0 40 -0.1 0.6]); legend('\rho00','\rho01','\rho10','\rho11')

%clear UA UAB UB WA WB colu fila klm rho tab1 tab2 tab3 tempo Tiempo
%clear r11 r14 r22 r23 r33 r44