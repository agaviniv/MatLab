clear all
close all

%%
% Atenuancao de amplitude de um estado maximamente misturado
subplot(2,3,1)
ro = (1/32)*diag(ones(1,32));
figwig(ro);

subplot(2,3,2)
ro1=atenuamp(1,.5,'a',ro);
figwig(ro1);

subplot(2,3,3)
ro2=atenuamp(1,.5,'b',ro1);
figwig(ro2);

subplot(2,3,4)
ro3=atenuamp(1,.5,'c',ro2);
figwig(ro3);

subplot(2,3,5)
ro4=atenuamp(1,.5,'d',ro3);
figwig(ro4);

subplot(2,3,6)
ro5=atenuamp(1,.5,'e',ro4);
figwig(ro5);

%%
% Atenuancao de fase de um estado maximamente misturado
cpn=0.95; 
figure; subplot(2,3,1)
ro=(1/32)*([ones(32,1)]*[ones(1,32)]);
figwig(ro);

subplot(2,3,2)
ro1=atenufase(cpn,'a',ro);
figwig(ro1);

subplot(2,3,3)
ro2=atenufase(cpn,'b',ro1);
figwig(ro2);

subplot(2,3,4)
ro3=atenufase(cpn,'c',ro2);
figwig(ro3);

subplot(2,3,5)
ro4=atenufase(cpn,'d',ro3);
figwig(ro4);

subplot(2,3,6)
ro5=atenufase(cpn,'e',ro4);
figwig(ro5);

%%
% Porta de Hadamard - comecando pelo qbit mais significativo
figure; subplot(2,3,1)
ro=diag([1 zeros(1,31)]);
figwig(ro);

subplot(2,3,2)
ro1=hgate('a',ro);
figwig(ro1);

subplot(2,3,3)
ro2=hgate('b',ro1);
figwig(ro2);

subplot(2,3,4)
ro3=hgate('c',ro2);
figwig(ro3);

subplot(2,3,5)
ro4=hgate('d',ro3);
figwig(ro4);

subplot(2,3,6)
ro5=hgate('e',ro4);
figwig(ro5);

%%
% Porta de Hadamard - comecando pelo qbit menos significativo
figure; subplot(2,3,1)
ro=diag([1 zeros(1,31)]);
figwig(ro);

subplot(2,3,2)
ro1=hgate('e',ro);
figwig(ro1);

subplot(2,3,3)
ro2=hgate('d',ro1);
figwig(ro2);

subplot(2,3,4)
ro3=hgate('c',ro2);
figwig(ro3);

subplot(2,3,5)
ro4=hgate('b',ro3);
figwig(ro4);

subplot(2,3,6)
ro5=hgate('a',ro4);
figwig(ro5);

%%
% Porta de Fase S - comecando pelo qbit mais significativo
figure; subplot(2,3,1)
ro=(1/32)*([ones(32,1)]*[ones(1,32)]);
figwig(ro);

subplot(2,3,2)
ro1=sgate('a',ro);
figwig(ro1);

subplot(2,3,3)
ro2=sgate('b',ro1);
figwig(ro2);

subplot(2,3,4)
ro3=sgate('c',ro2);
figwig(ro3);

subplot(2,3,5)
ro4=sgate('d',ro3);
figwig(ro4);

subplot(2,3,6)
ro5=sgate('e',ro4);
figwig(ro5);


%%
% Porta de Fase S - comecando pelo qbit menos significativo
figure; subplot(2,3,1)
ro=(1/32)*([ones(32,1)]*[ones(1,32)]);
figwig(ro);

subplot(2,3,2)
ro1=sgate('e',ro);
figwig(ro1);

subplot(2,3,3)
ro2=sgate('d',ro1);
figwig(ro2);

subplot(2,3,4)
ro3=sgate('c',ro2);
figwig(ro3);

subplot(2,3,5)
ro4=sgate('b',ro3);
figwig(ro4);

subplot(2,3,6)
ro5=sgate('a',ro4);
figwig(ro5);
