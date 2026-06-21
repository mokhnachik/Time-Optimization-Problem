
restart;
with(LinearAlgebra);
with(Optimization);
with(ComputationalGeometry);
with(plots);
with(plottools);
PolyPlot_3d := proc(X, Col) local plt; display(map(x -> plottools:-polygon(Transpose(X)[x], color = Col), ConvexHull(Transpose(X))), axes = none, style = LINE); end proc;
MinkowskiSum_3d := proc(X, Y) local temp, i, j; temp := <<0, 0, 0>>; for i to ColumnDimension(X) do for j to ColumnDimension(Y) do temp := <temp | Column(X, i) + Column(Y, j)>; end do; end do; SubMatrix(temp, 1 .. 3, convert(convert(convert(ConvexHull(Transpose(temp)), Vector), set), list)); end proc;
UnionPolyhedra_3d := proc(X, Y) local temp; temp := <X | Y>; SubMatrix(temp, 1 .. 3, convert(convert(convert(ConvexHull(Transpose(temp)), Vector), set), list)); end proc;
MinkowskiFunctional_3d := proc(x, X) local i, lambda, temp; lambda := Vector(ColumnDimension(X)); for i to ColumnDimension(X) do lambda[i] := cat('lambda', i); end do; temp := (X . lambda) - x; Minimize(add(lambda[i], i = 1 .. ColumnDimension(X)), {temp[1] = 0, temp[2] = 0, temp[3] = 0}, assume = nonnegative)[1]; end proc;
supPointC := proc(p, C) evalf((1/C) . p/sqrt(((Transpose(p)) . (1/C)) . p)); end proc;
HasPlane := proc(N::Matrix, A::list, n::Vector, a) local j; for j to nops(A) do if Norm(Column(N, j) - n, 2) < tol and abs(A[j] - a) < tol then return true; end if; end do; return false; end proc;
CanonicalizePlanes := proc(N::Matrix, Alpha::list) local i, s, n, a, Nnew, Anew; if ColumnDimension(N) <> nops(Alpha) then error "CanonicalizePlanes: number of normals and alpha values must match"; end if; Nnew := Matrix(3, 0); Anew := []; for i to ColumnDimension(N) do n := Column(N, i); a := Alpha[i]; s := Norm(n, 2); if s < tol then error "CanonicalizePlanes: zero normal at index %1", i; end if; n := evalf(n/s); a := evalf(a/s); if not HasPlane(Nnew, Anew, n, a) then Nnew := <Nnew | n>; Anew := [op(Anew), a]; end if; end do; return Nnew, Anew; end proc;
VertToNorm := proc(Vert::Matrix) local pts, hull, center, i, f, v1, v2, v3, n, a, N, Alpha, s; if RowDimension(Vert) <> 3 then error "VertToNorm: Vert must be a 3 x m matrix"; end if; if ColumnDimension(Vert) < 4 then error "VertToNorm: need at least 4 vertices for a 3D polyhedron"; end if; pts := [seq([evalf(Vert[1, i]), evalf(Vert[2, i]), evalf(Vert[3, i])], i = 1 .. ColumnDimension(Vert))]; hull := ConvexHull(pts); center := add(Column(Vert, i), i = 1 .. ColumnDimension(Vert))/ColumnDimension(Vert); N := Matrix(3, 0); Alpha := []; for i to numelems(hull) do f := hull[i]; if nops(f) < 3 then next; end if; v1 := Column(Vert, f[1]); v2 := Column(Vert, f[2]); v3 := Column(Vert, f[3]); n := CrossProduct(v2 - v1, v3 - v1); s := Norm(n, 2); if s < tol then next; end if; n := evalf(n/s); a := evalf(DotProduct(n, v1)); if a < DotProduct(n, center) then n := -n; a := -a; end if; if not HasPlane(N, Alpha, n, a) then N := <N | n>; Alpha := [op(Alpha), a]; end if; end do; return CanonicalizePlanes(N, Alpha); end proc;
NormToVert := proc(P::Matrix, Alpha::list) local P2, A2, i, j, k, l, M, rhs, x, vals, V, n, t; if RowDimension(P) <> 3 then error "NormToVert: P must be 3 x n"; end if; if ColumnDimension(P) <> nops(Alpha) then error "NormToVert: mismatch between P and Alpha"; end if; P2, A2 := CanonicalizePlanes(P, Alpha); n := ColumnDimension(P2); V := Matrix(3, 0); for i to n - 2 do for j from i + 1 to n - 1 do for k from j + 1 to n do M := Matrix([[P2[1, i], P2[2, i], P2[3, i]], [P2[1, j], P2[2, j], P2[3, j]], [P2[1, k], P2[2, k], P2[3, k]]]); if abs(Determinant(M)) < tol then next; end if; x := LinearSolve(M, Vector([A2[i], A2[j], A2[k]])); x := evalf(x); vals := [seq(evalf(DotProduct(Column(P2, l), x) - A2[l]), l = 1 .. n)]; if max(op(vals)) <= tol then if ColumnDimension(V) = 0 then V := <V | x>; else if tol < max(seq(Norm(Column(V, t) - x, 2), t = 1 .. ColumnDimension(V))) then V := <V | x>; end if; end if; end if; end do; end do; end do; if ColumnDimension(V) = 0 then error "NormToVert: no vertices found"; end if; return V; end proc;
NULL;
SupportValue_3d := proc(u::Vector, Q::Matrix) local i, val, best; best := -infinity; for i to ColumnDimension(Q) do val := evalf(DotProduct(u, Column(Q, i))); if best < val then best := val; end if; end do; return best; end proc;
PolyhedralAlgo := proc(P::Matrix, NormalQ::Matrix, AlphaQ::list, C::Matrix) local NormalP0, AlphaP0, Q0, gaps, i, idx, u, pnew, alphaNew, P1, NormalP1, AlphaP1, NormalQ1, AlphaQ1, Q1; NormalP0, AlphaP0 := VertToNorm(P); Q0 := NormToVert(NormalQ, AlphaQ); gaps := [seq(SupportValue_3d(Column(NormalP0, i), Q0) - AlphaP0[i], i = 1 .. ColumnDimension(NormalP0))]; idx := 1; for i to nops(gaps) do if gaps[idx] < gaps[i] then idx := i; end if; end do; u := Column(NormalP0, idx); pnew := supPointC(u, C); alphaNew := evalf(DotProduct(u, pnew)); P1 := <P | pnew>; NormalP1, AlphaP1 := VertToNorm(P1); NormalQ1 := <NormalQ | u>; AlphaQ1 := [op(AlphaQ), alphaNew]; Q1 := NormToVert(NormalQ1, AlphaQ1); return P1, Q1, NormalP1, AlphaP1, NormalQ1, AlphaQ1, Q0; end proc;
tol := 0.10e-7;
NULL;
NULL;
dark_grey := ColorTools:-Color([0.3, 0.3, 0.4]);
dark_green := ColorTools:-Color([1, 121, 111]);
Color := [dark_grey, black, yellow, orange, pink, purple, blue, yellow, orange, pink, purple, blue, yellow, orange, pink, purple, blue, yellow, orange, pink, purple, blue];
NULL;
NULL;
t := 3;
U := <<2, -1, 1> | <0, -1, 0> | <-2, 0, 1> | <0, 1, 2>>;
U := UnionPolyhedra_3d(U, -U);
V := <<6, 1, 1> | <0, -6, 4> | <-2, 0, 1> | <-1, -2, 1>>;
A := 0.99*IdentityMatrix(3);
x0 := <1.5, 0, -3.5>;
X1 := [-((1/A) . U)];
X2 := [-((1/A) . V)];
plX1 := [PolyPlot_3d(X1[1], pink)];
plX2 := [PolyPlot_3d(t*X2[1], yellow)];
for N while 1 < MinkowskiFunctional_3d(x0, X1[N]) or t < MinkowskiFunctional_3d(x0, X2[N]) do
    X1 := [op(X1), MinkowskiSum_3d((1/A) . (X1[N]), X1[1])];
    X2 := [op(X2), UnionPolyhedra_3d((1/A) . (X2[N]), X2[1])];
    plX2 := [op(plX2), PolyPlot_3d(t*X2[N + 1], yellow)];
    plX1 := [op(plX1), PolyPlot_3d(X1[N + 1], pink)];
    print([N, MinkowskiFunctional_3d(x0, X1[N]), MinkowskiFunctional_3d(x0, X2[N])]);
end do;
print([N, MinkowskiFunctional_3d(x0, X1[N]), MinkowskiFunctional_3d(x0, X2[N])]);
Nmin := N;
M1 := ColumnDimension(U);
M2 := ColumnDimension(V);
T := [t];
X3 := [T[1]*X2[Nmin]];
plX3 := [PolyPlot_3d(X3[1], magenta)];
pic := [];
eta := Vector(M1);
nu := Vector(M2);
for i to M1 do
    eta[i] := cat('eta', i);
end do;
for i to M2 do
    nu[i] := cat('nu', i);
end do;
eta_known := Vector(M1);
nu_known := Vector(M2);
Trajectory := <x0>;
Control := <<0, 0, 0>>;
k := 0;
while k <= Nmin - 2 do
    L1 := ColumnDimension(X1[Nmin - k - 1]);
    L2 := ColumnDimension(X2[Nmin - k - 1]);
    lambda := Vector(L1);
    mu := Vector(L2);
    initp := [];
    for i to L1 do
        lambda[i] := cat('lambda', i);
        initp := [op(initp), lambda[i] = 0];
    end do;
    for i to L2 do
        mu[i] := cat('mu', i);
        initp := [op(initp), mu[i] = 0];
    end do;
    for i to M1 do
        initp := [op(initp), eta[i] = 0];
    end do;
    for i to M2 do
        initp := [op(initp), nu[i] = 0];
    end do;
    xk := Column(Trajectory, k + 1);
    eq1 := (A . xk) + (U . eta) - ((X1[Nmin - k - 1]) . lambda);
    eq2 := (A . xk) + (V . nu) - ((X2[Nmin - k - 1]) . mu);
    eq3 := (U . eta) - (V . nu);
    constr := {eq1[1] = 0, eq1[2] = 0, eq1[3] = 0, eq2[1] = 0, eq2[2] = 0, eq2[3] = 0, eq3[1] = 0, eq3[2] = 0, eq3[3] = 0, add(lambda[i], i = 1 .. L1) = 1 - 0.001*MinkowskiFunctional_3d(x0, X1[Nmin])*(Nmin - k)/(k + 1), add(mu[i], i = 1 .. L2) + add(nu[i], i = 1 .. M2) <= T[k + 1], add(eta[i], i = 1 .. M1) <= 1};
    try
        res := Minimize(add(nu[i], i = 1 .. M2), constr, assume = nonnegative, variables = [op(convert(lambda, list)), op(convert(mu, list)), op(convert(eta, list)), op(convert(nu, list))], initialpoint = initp);
    catch:
        X1 := [op(X1), MinkowskiSum_3d((1/A) . (X1[Nmin]), X1[1])];
        X2 := [op(X2), UnionPolyhedra_3d((1/A) . (X2[Nmin]), X2[1])];
        plX2 := [op(plX2), PolyPlot_3d(t*X2[Nmin + 1], yellow)];
        plX1 := [op(plX1), PolyPlot_3d(X1[Nmin + 1], pink)];
        Nmin := Nmin + 1;
        print('Nmin' = Nmin);
        Trajectory := <x0>;
        Control := <<0, 0, 0>>;
        T := [t];
        X3 := [T[1]*X2[Nmin]];
        plX3 := [PolyPlot_3d(X3[1], magenta)];
        pic := [];
        k := 0;
        next;
    end try;
    lambda_known := Vector(L1);
    mu_known := Vector(L2);
    for i to L1 do
        lambda_known[i] := rhs(res[2][i]);
    end do;
    for i to L2 do
        mu_known[i] := rhs(res[2][L1 + i]);
    end do;
    for i to M1 do
        eta_known[i] := rhs(res[2][L1 + L2 + i]);
    end do;
    for i to M2 do
        nu_known[i] := rhs(res[2][L1 + L2 + M1 + i]);
    end do;
    Control := <Control | U . eta_known>;
    Trajectory := <Trajectory | (A . xk) + Column(Control, k + 2)>;
    T := [op(T), T[k + 1] - add(nu_known[i], i = 1 .. M2)];
    X3 := [op(X3), T[k + 2]*X2[Nmin - k - 1]];
    plX3 := [op(plX3), PolyPlot_3d(T[k + 2]*X2[Nmin - k - 1], magenta)];
    x_new := Column(Trajectory, k + 2);
    m1_old := add(lambda_known[i], i = 1 .. L1);
    m2_old := add(mu_known[i], i = 1 .. L2);
    m1 := MinkowskiFunctional_3d(x_new, X1[Nmin - k - 1]);
    m2 := [MinkowskiFunctional_3d(x_new, T[k + 2]*X2[Nmin - k - 1]), MinkowskiFunctional_3d(x_new, X2[Nmin - k - 1])];
    print(Nmin - k, Nmin - k - 1, x_new, m1_old, m2_old, m1, m2, T[k + 2], res[1]);
    pic := [op(pic), display(PolyPlot_3d(X1[Nmin - k], pink), PolyPlot_3d(T[k + 1]*X2[Nmin - k], yellow), PolyPlot_3d(X1[Nmin - k - 1], pink), PolyPlot_3d(T[k + 2]*X2[Nmin - k - 1], magenta), pointplot3d(xk), pointplot3d(x_new))];
    print(display(pic[k + 1]));
    k := k + 1;
end do;
print('Nmin' = Nmin);
Control := <Control | -(A . (Column(Trajectory, Nmin)))>;
Control := SubMatrix(Control, 1 .. 2, 2 .. Nmin + 1);
Trajectory := <Trajectory | <0, 0, 0>>;
T;
display(plX1, plX3, spacecurve(Transpose(Trajectory), style = pointline, color = blue));
NULL;
matrixA := proc(k) local p, m, phi, S; S := <<1, 0, 1> | <1, 1, 0> | <0, 1, 1>>; phi := 1/4*Pi; p := 1 - 1/(k + 2); m := <<p, 0, 0> | <0, p*cos(phi), p*sin(phi)> | <0, -p*sin(phi), p*cos(phi)>>; MatrixInverse((S . m) . (MatrixInverse(S))); end proc;
A := evalf(matrixA(0, 1));


t := 3;
Cu := DiagonalMatrix(<1, 2, 3>);
Cv := DiagonalMatrix(<1, 3.5, 0.5>);
A := Matrix(3, 3, [[0.99, 0., 0.], [0., 0.99, 0.], [0., 0., 0.99]]);
x0 := <1, 1, 1>;
NULL;
NULL;
NULL;
display(implicitplot3d((<x | y | z> . Cu) . <x | y | z> = 1, x = -1 .. 1, y = -1 .. 1, z = -1 .. 1, color = gray), implicitplot3d((<x | y | z> . Cv) . <x | y | z> = 1, x = -1 .. 1, y = -1 .. 1, z = -2 .. 2, color = dark_green), style = line);
x0;
p0 := <supPointC(<1, 0, 0>, Cu) | supPointC(<-1, 0, 0>, Cu) | supPointC(<0, 1, 0>, Cu) | supPointC(<0, -1, 0>, Cu) | supPointC(<0, 0, 1>, Cu) | supPointC(<0, 0, -1>, Cu)>;
normalq0 := <<1, 0, 0> | <-1, 0, 0> | <0, 1, 0> | <0, -1, 0> | <0, 0, 1> | <0, 0, -1>>;
alphaq0 := [];
for i to ColumnDimension(normalq0) do
    alphaq0 := [op(alphaq0), (Transpose(Column(normalq0, i))) . (Column(p0, i))];
end do;
alphaq0;
q0 := NormToVert(normalq0, alphaq0);
NULL;
NULL;
P0 := <supPointC(<1, 0, 0>, Cv) | supPointC(<-1, 0, 0>, Cv) | supPointC(<0, 1, 0>, Cv) | supPointC(<0, -1, 0>, Cv) | supPointC(<0, 0, 1>, Cv) | supPointC(<0, 0, -1>, Cv)>;
NormalQ0 := <<1, 0, 0> | <-1, 0, 0> | <0, 1, 0> | <0, -1, 0> | <0, 0, 1> | <0, 0, -1>>;
AlphaQ0 := [];
for i to ColumnDimension(NormalQ0) do
    AlphaQ0 := [op(AlphaQ0), (Transpose(Column(NormalQ0, i))) . (Column(P0, i))];
end do;
AlphaQ0;
Q0 := NormToVert(NormalQ0, AlphaQ0);
p := p0;
norml := normalq0;
alphaq := alphaq0;
n_in := 0;
n_out := 2;
P := P0;
Norml := NormalQ0;
AlphaQ := AlphaQ0;
N_in := 1;
N_out := 2;
Ku_in := 6;
Ku_out := 8;
Kv_in := 6;
Kv_out := 8;
while n_in <> n_out or N_in <> N_out do if n_in <> n_out then res := PolyhedralAlgo(p, norml, alphaq, Cu); p := res[1]; q := res[2]; norml := res[5]; alphaq := res[6]; Ku_in := Ku_in + 1; Ku_out := Ku_out + 1; print(display(implicitplot3d((<x | y | z> . Cu) . <x | y | z> = 1, x = -1 .. 1, y = -1 .. 1, z = -1 .. 1, color = grey, style = line), PolyPlot_3d(p, Color[2]), PolyPlot_3d(q, Color[1]))); U := p; X1 := [-((1/A) . U)]; plX1_in := [PolyPlot_3d(X1[1], black)]; for N while 1 < MinkowskiFunctional_3d(x0, X1[N]) do X1 := [op(X1), MinkowskiSum_3d((1/A) . (X1[N]), X1[1])]; plX1_in := [op(plX1_in), PolyPlot_3d(X1[N + 1], black)]; end do; n_in := N; U := q; X1_out := [-((1/A) . U)]; plX1_out := [PolyPlot_3d(X1_out[1], dark_grey)]; for N while 1 < MinkowskiFunctional_3d(x0, X1_out[N]) do X1_out := [op(X1_out), MinkowskiSum_3d((1/A) . (X1_out[N]), X1_out[1])]; plX1_out := [op(plX_out), PolyPlot_3d(X1_out[N + 1], Color[N])]; end do; n_out := N; end if; print(n_in, n_out); if N_in <> N_out then res := PolyhedralAlgo(P, Norml, AlphaQ, Cv); P := res[1]; Q := res[2]; Norml := res[5]; AlphaQ := res[6]; Kv_in := Kv_in + 1; Kv_out := Kv_out + 1; print(display(implicitplot3d((<x | y | z> . Cv) . <x | y | z> = 1, x = -2 .. 2, y = -2 .. 2, z = -2 .. 2, color = dark_grey, style = line), PolyPlot_3d(P, Color[2]), PolyPlot_3d(Q, Color[1]))); V := P; X2 := [-((1/A) . V)]; plX2_in := [PolyPlot_3d(X2[1], black)]; for N while t < MinkowskiFunctional_3d(x0, X2[N]) do X2 := [op(X2), UnionPolyhedra_3d((1/A) . (X2[N]), X2[1])]; plX2_in := [op(plX2_in), PolyPlot_3d(X2[N + 1], black)]; end do; N_in := N; V := Q; X2_out := [-((1/A) . V)]; plX2_out := [PolyPlot_3d(X2_out[1], dark_grey)]; for N while t < MinkowskiFunctional_3d(x0, X2_out[N]) do X2_out := [op(X2_out), UnionPolyhedra_3d((1/A) . (X2_out[N]), X2_out[1])]; plX2_out := [op(plX2_out), PolyPlot_3d(X2_out[N + 1], Color[N])]; end do; N_out := N; end if; print(N_in, N_out); end do;
print(Ku_in, Ku_out);
print(Kv_in, Kv_out);
U := p;
V := P;
X1 := [-((1/A) . U)];
X2 := [-((1/A) . V)];
plX1 := [PolyPlot_3d(X1[1], pink)];
plX2 := [PolyPlot_3d(t*X2[1], yellow)];
for N while 1 < MinkowskiFunctional_3d(x0, X1[N]) or t < MinkowskiFunctional_3d(x0, X2[N]) do
    X1 := [op(X1), MinkowskiSum_3d((1/A) . (X1[N]), X1[1])];
    X2 := [op(X2), UnionPolyhedra_3d((1/A) . (X2[N]), X2[1])];
    plX2 := [op(plX2), PolyPlot_3d(t*X2[N + 1], yellow)];
    plX1 := [op(plX1), PolyPlot_3d(X1[N + 1], pink)];
    print([N, MinkowskiFunctional_3d(x0, X1[N]), MinkowskiFunctional_3d(x0, X2[N])]);
end do;
print([N, MinkowskiFunctional_3d(x0, X1[N]), MinkowskiFunctional_3d(x0, X2[N])]);
Nmin := N;
M1 := ColumnDimension(U);
M2 := ColumnDimension(V);
T := [t];
X3 := [T[1]*X2[Nmin]];
plX3 := [PolyPlot_3d(X3[1], magenta)];
pic := [];
eta := Vector(M1);
nu := Vector(M2);
for i to M1 do
    eta[i] := cat('eta', i);
end do;
for i to M2 do
    nu[i] := cat('nu', i);
end do;
eta_known := Vector(M1);
nu_known := Vector(M2);
Trajectory := <x0>;
Control := <<0, 0, 0>>;
k := 0;
while k <= Nmin - 2 do
    L1 := ColumnDimension(X1[Nmin - k - 1]);
    L2 := ColumnDimension(X2[Nmin - k - 1]);
    lambda := Vector(L1);
    mu := Vector(L2);
    initp := [];
    for i to L1 do
        lambda[i] := cat('lambda', i);
        initp := [op(initp), lambda[i] = 0];
    end do;
    for i to L2 do
        mu[i] := cat('mu', i);
        initp := [op(initp), mu[i] = 0];
    end do;
    for i to M1 do
        initp := [op(initp), eta[i] = 0];
    end do;
    for i to M2 do
        initp := [op(initp), nu[i] = 0];
    end do;
    xk := Column(Trajectory, k + 1);
    eq1 := (A . xk) + (U . eta) - ((X1[Nmin - k - 1]) . lambda);
    eq2 := (A . xk) + (V . nu) - ((X2[Nmin - k - 1]) . mu);
    eq3 := (U . eta) - (V . nu);
    constr := {eq1[1] = 0, eq1[2] = 0, eq1[3] = 0, eq2[1] = 0, eq2[2] = 0, eq2[3] = 0, eq3[1] = 0, eq3[2] = 0, eq3[3] = 0, add(lambda[i], i = 1 .. L1) = 1 - 0.001*MinkowskiFunctional_3d(x0, X1[Nmin])*(Nmin - k)/(k + 1), add(mu[i], i = 1 .. L2) + add(nu[i], i = 1 .. M2) <= T[k + 1], add(eta[i], i = 1 .. M1) <= 1};
    try
        res := Minimize(add(nu[i], i = 1 .. M2), constr, assume = nonnegative, variables = [op(convert(lambda, list)), op(convert(mu, list)), op(convert(eta, list)), op(convert(nu, list))], initialpoint = initp);
    catch:
        X1 := [op(X1), MinkowskiSum_3d((1/A) . (X1[Nmin]), X1[1])];
        X2 := [op(X2), UnionPolyhedra_3d((1/A) . (X2[Nmin]), X2[1])];
        plX2 := [op(plX2), PolyPlot_3d(t*X2[Nmin + 1], yellow)];
        plX1 := [op(plX1), PolyPlot_3d(X1[Nmin + 1], pink)];
        Nmin := Nmin + 1;
        print('Nmin' = Nmin);
        Trajectory := <x0>;
        Control := <<0, 0, 0>>;
        T := [t];
        X3 := [T[1]*X2[Nmin]];
        plX3 := [PolyPlot_3d(X3[1], magenta)];
        pic := [];
        k := 0;
        next;
    end try;
    lambda_known := Vector(L1);
    mu_known := Vector(L2);
    for i to L1 do
        lambda_known[i] := rhs(res[2][i]);
    end do;
    for i to L2 do
        mu_known[i] := rhs(res[2][L1 + i]);
    end do;
    for i to M1 do
        eta_known[i] := rhs(res[2][L1 + L2 + i]);
    end do;
    for i to M2 do
        nu_known[i] := rhs(res[2][L1 + L2 + M1 + i]);
    end do;
    Control := <Control | U . eta_known>;
    Trajectory := <Trajectory | (A . xk) + Column(Control, k + 2)>;
    T := [op(T), T[k + 1] - add(nu_known[i], i = 1 .. M2)];
    X3 := [op(X3), T[k + 2]*X2[Nmin - k - 1]];
    plX3 := [op(plX3), PolyPlot_3d(T[k + 2]*X2[Nmin - k - 1], magenta)];
    x_new := Column(Trajectory, k + 2);
    m1_old := add(lambda_known[i], i = 1 .. L1);
    m2_old := add(mu_known[i], i = 1 .. L2);
    m1 := MinkowskiFunctional_3d(x_new, X1[Nmin - k - 1]);
    m2 := [MinkowskiFunctional_3d(x_new, T[k + 2]*X2[Nmin - k - 1]), MinkowskiFunctional_3d(x_new, X2[Nmin - k - 1])];
    print(Nmin - k, Nmin - k - 1, x_new, m1_old, m2_old, m1, m2, T[k + 2], res[1]);
    pic := [op(pic), display(PolyPlot_3d(X1[Nmin - k], pink), PolyPlot_3d(T[k + 1]*X2[Nmin - k], yellow), PolyPlot_3d(X1[Nmin - k - 1], pink), PolyPlot_3d(T[k + 2]*X2[Nmin - k - 1], magenta), pointplot3d(xk), pointplot3d(x_new))];
    print(display(pic[k + 1]));
    k := k + 1;
end do;
print('Nmin' = Nmin);
Control := <Control | -(A . (Column(Trajectory, Nmin)))>;
Control := SubMatrix(Control, 1 .. 2, 2 .. Nmin + 1);
Trajectory := <Trajectory | <0, 0, 0>>;
T := [op(T), T[k + 1] - MinkowskiFunctional_3d(Column(Trajectory, Nmin), X2[1])];

T[1] - T[Nmin + 1];
ColumnDimension(Trajectory);
display(plX1, spacecurve(Transpose(Trajectory), style = pointline, color = blue));
display(plX3, spacecurve(Transpose(Trajectory), style = pointline, color = blue));
NULL;
NULL;
NULL;
matrixA := proc(k) local p, m, phi, S; S := <<1, 0, 1> | <1, 1, 0> | <0, 1, 1>>; phi := 1/4*Pi; p := 1 - 1/(k + 2); m := <<p, 0, 0> | <0, p*cos(phi), p*sin(phi)> | <0, -p*sin(phi), p*cos(phi)>>; MatrixInverse((S . m) . (MatrixInverse(S))); end proc;
matrixA(0);
t := 1;
x0 := <2, 4, 2>;
A := Matrix(3, 3, [[-cos(0.25) + 2, sin(0.25), -2*cos(0.25) + 2], [sin(0.25), cos(0.25), 2*sin(0.25)], [cos(0.25) - 1, -sin(0.25), 2*cos(0.25) - 1]]);
U := Matrix(3, 4, [[0.0362796597, -0.0258954969, 0.0258954969, -0.0362796597], [0.3095791163, -0.1852288023, 0.1852288023, -0.3095791163], [0.2137203403, 0.2758954969, -0.2758954969, -0.2137203403]]);
Acont := <<0 | 1 | 0>, <1 | 0 | 2>, <0 | -1 | 0>>;
NULL;
Lambda_cont, S_cont := Eigenvectors(Acont);
phi1 := exp(Re(Lambda_cont[2])*tt)*(Re(Column(S_cont, 2))*cos(Im(Lambda_cont[2])*tt) - Im(Column(S_cont, 2))*sin(Im(Lambda_cont[2])*tt));
phi2 := exp(Re(Lambda_cont[2])*tt)*(Im(Column(S_cont, 2))*cos(Im(Lambda_cont[2])*tt) + Re(Column(S_cont, 2))*sin(Im(Lambda_cont[2])*tt));
phi3 := exp(Lambda_cont[1]*tt)*Column(S_cont, 1);
Phi := <phi1 | phi2 | phi3>;
diff(Phi . <C1, C2, C3>, tt) = (Acont . Phi) . <C1, C2, C3>;
ReS := <Re(Column(S_cont, 2)) | Im(Column(S_cont, 2)) | Column(S_cont, 1)>;
J1 := <<Re(Lambda_cont[2]) | Im(Lambda_cont[2])>, <-Im(Lambda_cont[2]) | Re(Lambda_cont[2])>>;
Phi1 := -((ReS . <<MatrixInverse(J1) | <0, 0>>, <0 | 0 | 0>>) . (MatrixInverse(ReS))) + (((Column(S_cont, 1)) . <0 | 0 | 1>) . (MatrixInverse(ReS)))*tt;
Adiscrete := Phi . (MatrixInverse(subs([tt = 0], Phi)));
Bdiscrete := Phi1 - (Adiscrete . (subs([tt = 0], Phi1)));
Solution := (Adiscrete . <y01, y02, y03>) + (Bdiscrete . <u1, u2, u3>);
simplify(diff(Solution, tt) - (Acont . Solution) - <u1, u2, u3>);
simplify(subs([tt = 0], Solution));
Delta := 0.25;
U := (subs([tt = Delta], Bdiscrete)) . <<0, 1, 1> | <0, -1, 1> | <0, 1, -1> | <0, -1, -1>>;
A := evalf(subs([tt = Delta], Adiscrete));
V := (subs([tt = Delta], Bdiscrete)) . <<0, 1, 1> | <0, -1, 1> | <0, 1, -1> | <0, -1, -1>>;
<v1, v2, v3>;
(MatrixInverse(Bdiscrete)) . (<v1, v2, v3>[1]);
supPointC := proc(p, C) evalf((1/C) . p/sqrt(((Transpose(p)) . (1/C)) . p)); end proc;

NULL;
