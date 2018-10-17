function fdlty = fidelity(ro1,ro2)

numer = trace(ro1*ro2);

fdlty = numer / sqrt(numer);