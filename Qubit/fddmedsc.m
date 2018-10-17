function fdd = fddmedsc(Uth,Up,nit)

A = Uth^nit;
B = Up^nit;

fdd = scattcirc(.25*eye(4),A'*B,'z');