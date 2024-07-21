function E = LatticeError_hexa_new(x, Qmean, Label_row, Label_col)
    init = x(1:2);
    a1 = x(3:4);
    a2= x(5:6);

    Q_lattice = Label_row * a2 + Label_col * a1 + init;
    E = Q_lattice - Qmean;
    E = sum(E.^2, 'all') / size(Qmean,1);
end