
clear all; clc;

syms a b c d
rom=diag([a b c d]);

P{1}=eye(4);
P{2}=[1 0 0 0;0 0 0 1;0 1 0 0;0 0 1 0];
P{3}=[1 0 0 0;0 0 1 0;0 0 0 1;0 1 0 0];

sro = (P{1}*rom*P{1}'+P{2}*rom*P{2}'+P{3}*rom*P{3}') %nao normalizada
sroa=subs(sro,'b+d+c','1-a')
ropu = (sroa-(1-a)*eye(4)); roid = (1-a)*eye(4);

syms A  %A=hc*wL/kB*T
ropuA=subs(ropu,'a','1-A'); roidA=subs(roid,'a','1-A');
roA=roidA+ropuA

syms E
ropps=simplify(subs(roA,'A','1-E'))
