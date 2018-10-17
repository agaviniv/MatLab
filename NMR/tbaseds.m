% Operadores para dois spins 1/2
% Base de transição em funcao dos operadore de Pauli (Cartesianos)

[x y z]=mangqbts(2);

% Base de Transicao
t{1}=eye(4);
t{2}=2*z{1};
t{3}=2*z{2};
t{4}=4*z{1}*z{2};
t{5}=4*(x{1}*x{2}+y{1}*y{2});
t{6}=4*(x{1}*y{2}-y{1}*x{2});
t{7}=2*x{1};
t{8}=2*y{1};
t{9}=2*x{2};
t{10}=2*y{2};
t{11}=4*x{1}*z{2};
t{12}=4*y{1}*z{2};
t{13}=4*z{1}*x{2};
t{14}=4*z{1}*y{2};
t{15}=4*(x{1}*x{2}-y{1}*y{2});
t{16}=4*(x{1}*y{2}+y{1}*x{2});


%ro = diag([0 1 0 1]/2);
ro = (1/4)*[1;1;1;1]*[1 1 1 1];

tro = zeros(4); c = zeros(1,16);
for k=1:16
    c(k) = 0.5*trace(ro*t{k});
    tro = tro + 0.5*c(k)*t{k};
end    

cc = reshape(c,4,4);