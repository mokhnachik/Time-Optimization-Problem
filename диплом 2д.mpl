
restart;
with(LinearAlgebra);
with(Optimization);
with(ComputationalGeometry);
with(plots);
PolyPlot := proc(X, col) local T; T := SubMatrix(X, 1 .. 2, ConvexHull(Transpose(X))); plot(Transpose(<T | Column(T, 1)>), color = col, thickness = 4); end proc;
MinkowskiSum := proc(X, Y) local temp, i, j; temp := <<0, 0>>; for i to ColumnDimension(X) do for j to ColumnDimension(Y) do temp := <temp | Column(X, i) + Column(Y, j)>; end do; end do; SubMatrix(temp, 1 .. 2, ConvexHull(Transpose(temp))); end proc;
UnionPolyhedra := proc(X, Y) local temp; temp := <X | Y>; SubMatrix(temp, 1 .. 2, ConvexHull(Transpose(temp))); end proc;
MinkowskiFunctional := proc(x, X) local i, lambda, temp; lambda := Vector(ColumnDimension(X)); for i to ColumnDimension(X) do lambda[i] := cat('lambda', i); end do; temp := (X . lambda) - x; Minimize(add(lambda[i], i = 1 .. ColumnDimension(X)), {temp[1] = 0, temp[2] = 0}, assume = nonnegative)[1]; end proc;
IntersectPolygons := proc(X1, X2) local K1, K2, i, vert, v, r, Matr; K1 := SubMatrix(X1, 1 .. 2, ConvexHull(Transpose(X1))); K2 := SubMatrix(X2, 1 .. 2, ConvexHull(Transpose(X2))); vert := []; for i to ColumnDimension(K1) do v := Transpose(K1)[i]; if member(PointInPolygon(v, Transpose(K2)), {"inside", "boundary"}) then vert := [op(vert), v]; end if; end do; for i to ColumnDimension(K2) do v := Transpose(K2)[i]; if member(PointInPolygon(v, Transpose(K1)), {"inside", "boundary"}) then vert := [op(vert), v]; end if; end do; r := MultiSegmentIntersect([Transpose(K1), Transpose(K2)], mode = ["polygon", "polygon"], output = "coordinate"); vert := [op(vert), op(r)]; Matr := Matrix(2, nops(vert), (i, j) -> vert[j][i]); return SubMatrix(Matr, 1 .. 2, ConvexHull(Transpose(Matr))); end proc;
VertToNorm := proc(P0) local num, P, temp, NormalP0, AlphaP0, i, area, v1, v2; num := ConvexHull(Transpose(P0)); P := SubMatrix(P0, 1 .. 2, num); area := 0; for i to ColumnDimension(P) do v1 := Column(P, i); v2 := Column(P, 1 + (i mod ColumnDimension(P))); area := area + v1[1]*v2[2] - v2[1]*v1[2]; end do; if area < 0 then P := SubMatrix(P, 1 .. 2, [seq(ColumnDimension(P) - i + 1, i = 1 .. ColumnDimension(P))]); end if; temp := <<0, -1> | <1, 0>> . (SubMatrix(P, 1 .. 2, [2 .. ColumnDimension(P), 1]) - P); NormalP0 := <<0, 0>>; for i to ColumnDimension(temp) do NormalP0 := <NormalP0 | Column(temp, i)/Norm(Column(temp, i), 2)>; end do; NormalP0 := SubMatrix(NormalP0, 1 .. 2, 2 .. ColumnDimension(NormalP0)); AlphaP0 := []; for i to ColumnDimension(P) do AlphaP0 := [op(AlphaP0), (Transpose(Column(NormalP0, i))) . (Column(P, i))]; end do; NormalP0, AlphaP0; end proc;
NormToVert := proc(Normal, AlphaQ) local Q0, temp, tempAlpha, i, angles, perm, pairs, NormalQ, AlphaQ0; NormalQ, AlphaQ0 := SortByAngle(Normal, AlphaQ); Q0 := <<0, 0>>; temp := <NormalQ | Column(NormalQ, 1)>; tempAlpha := [op(AlphaQ0), AlphaQ0[1]]; for i to ColumnDimension(NormalQ) do Q0 := <Q0 | LinearSolve(Transpose(<Column(temp, i) | Column(temp, i + 1)>), <AlphaQ0[i], tempAlpha[i + 1]>)>; end do; Q0 := SubMatrix(Q0, 1 .. 2, 2 .. ColumnDimension(Q0)); end proc;
SortByAngle := proc(N, Alpha) local i, k, theta, perm; k := ColumnDimension(N); theta := [seq(evalf(arctan(Column(N, i)[2], Column(N, i)[1])), i = 1 .. k)]; perm := sort([seq(i, i = 1 .. k)], (i, j) -> evalb(theta[i] < theta[j])); return SubMatrix(N, 1 .. 2, perm), [seq(Alpha[perm[i]], i = 1 .. k)]; end proc;
PolyhedralAlgo := proc(P, NormalQ, AlphaQ, C) local temp, Q0, NormalP0, AlphaP0, AlphaQ0, Q0constr, n, NormalQ1, AlphaQ1, P1, Q1, i; NormalP0, AlphaP0 := VertToNorm(P); Q0 := NormToVert(NormalQ, AlphaQ); Q0constr := {}; for i to ColumnDimension(NormalQ) do Q0constr := Q0constr union {<x | y> . (Column(NormalQ, i)) <= AlphaQ[i]}; end do; temp := []; for i to ColumnDimension(NormalP0) do temp := [op(temp), Maximize(<x | y> . (Column(NormalP0, i)), Q0constr)[1] - AlphaP0[i]]; end do; n := max[index](temp); P1 := <P | supPointC(Column(NormalP0, n), C)>; NormalQ1 := <NormalQ | Column(NormalP0, n)>; AlphaQ1 := [op(AlphaQ), (Transpose(Column(NormalP0, n))) . (supPointC(Column(NormalP0, n), C))]; Q1 := NormToVert(NormalQ1, AlphaQ1); return P1, Q1, NormalQ1, AlphaQ1, Q0; end proc;
supPointC := proc(p, C) evalf((1/C) . p/sqrt(((Transpose(p)) . (1/C)) . p)); end proc;
NULL;
t := 15;
U := <<2, 0> | <0, -1> | <-2, 0> | <0, 1>>;
V := <<1, 1> | <0, 1> | <-1, -1> | <0, -1>>;
S := RandomMatrix(2, 2);
A := evalf(((0.9*S) . <<cos(1), -sin(1)> | <sin(1), cos(1)>>) . (1/S));
A := Matrix(2, 2, [[0.2289125271, 0.3782702239], [-1.691313155, 0.7436316236]]);
x0 := <13.5, 13.5>;
X1 := [-((1/A) . U)];
X2 := [-((1/A) . V)];
plX1 := [PolyPlot(X1[1], "orchid")];
plX2 := [PolyPlot(t*X2[1], yellow)];
for N while 1 < MinkowskiFunctional(x0, X1[N]) or t < MinkowskiFunctional(x0, X2[N]) do
    X1 := [op(X1), MinkowskiSum((1/A) . (X1[N]), X1[1])];
    X2 := [op(X2), UnionPolyhedra((1/A) . (X2[N]), X2[1])];
    plX2 := [op(plX2), PolyPlot(t*X2[N + 1], yellow)];
    plX1 := [op(plX1), PolyPlot(X1[N + 1], orchid)];
    print([N, MinkowskiFunctional(x0, X1[N]), MinkowskiFunctional(x0, X2[N])]);
end do;
print([N, MinkowskiFunctional(x0, X1[N]), MinkowskiFunctional(x0, X2[N])]);
Nmin := N;
display(plX1, plX2, pointplot(x0));
M1 := ColumnDimension(U);
M2 := ColumnDimension(V);
T := [t];
X3 := [T[1]*X2[Nmin]];
plX3 := [PolyPlot(X3[1], magenta)];
inter := [PolyPlot(IntersectPolygons(X1[Nmin], X3[1]), "DeepPink")];
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
Control := <<0, 0>>;
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
    constr := {eq1[1] = 0, eq1[2] = 0, eq2[1] = 0, eq2[2] = 0, eq3[1] = 0, eq3[2] = 0, add(lambda[i], i = 1 .. L1) = 1 - 0.001*MinkowskiFunctional(x0, X1[Nmin])*(Nmin - k)/(k + 1), add(mu[i], i = 1 .. L2) + add(nu[i], i = 1 .. M2) <= T[k + 1], add(eta[i], i = 1 .. M1) <= 1};
    try
        res := Minimize(add(nu[i], i = 1 .. M2), constr, assume = nonnegative, variables = [op(convert(lambda, list)), op(convert(mu, list)), op(convert(eta, list)), op(convert(nu, list))], initialpoint = initp);
    catch:
        X1 := [op(X1), MinkowskiSum((1/A) . (X1[Nmin]), X1[1])];
        X2 := [op(X2), UnionPolyhedra((1/A) . (X2[Nmin]), X2[1])];
        plX2 := [op(plX2), PolyPlot(t*X2[Nmin + 1], yellow)];
        plX1 := [op(plX1), PolyPlot(X1[Nmin + 1], orchid)];
        Nmin := Nmin + 1;
        print('Nmin' = Nmin);
        Trajectory := <x0>;
        Control := <<0, 0>>;
        T := [t];
        X3 := [T[1]*X2[Nmin]];
        plX3 := [PolyPlot(X3[1], magenta)];
        inter := [PolyPlot(IntersectPolygons(X1[Nmin], X3[1]), "DeepPink")];
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
    plX3 := [op(plX3), PolyPlot(T[k + 2]*X2[Nmin - k - 1], magenta)];
    inter := [op(inter), PolyPlot(IntersectPolygons(X1[Nmin - k - 1], X3[k + 2]), "DeepPink")];
    x_new := Column(Trajectory, k + 2);
    m1_old := add(lambda_known[i], i = 1 .. L1);
    m2_old := add(mu_known[i], i = 1 .. L2);
    m1 := MinkowskiFunctional(x_new, X1[Nmin - k - 1]);
    m2 := [MinkowskiFunctional(x_new, T[k + 2]*X2[Nmin - k - 1]), MinkowskiFunctional(x_new, X2[Nmin - k - 1])];
    print(Nmin - k, Nmin - k - 1, x_new, m1_old, m2_old, m1, m2, T[k + 2], res[1]);
    pic := [op(pic), display(PolyPlot(X1[Nmin - k], pink), PolyPlot(T[k + 1]*X2[Nmin - k], yellow), PolyPlot(X1[Nmin - k - 1], pink), PolyPlot(T[k + 2]*X2[Nmin - k - 1], magenta), pointplot(xk), pointplot(x_new), inter[k + 2])];
    print(display(pic[k + 1]));
    k := k + 1;
end do;
print('Nmin' = Nmin);
Control := <Control | -(A . (Column(Trajectory, Nmin)))>;
Control := SubMatrix(Control, 1 .. 2, 2 .. Nmin + 1);
Trajectory := <Trajectory | <0, 0>>;
T := [op(T), T[k + 1] - MinkowskiFunctional(Column(Trajectory, Nmin), X2[1])];
T[1] - T[Nmin + 1];
display(plX1, plot(Transpose(Trajectory), style = pointline, color = black, thickness = 3));
display(inter, plot(Transpose(Trajectory), style = pointline, color = black, thickness = 3));
display(plX3, plot(Transpose(Trajectory), style = pointline, color = black, thickness = 3));
NULL;



NULL;
dark_grey := ColorTools:-Color([0.3, 0.3, 0.4]);
dark_green := ColorTools:-Color([1, 121, 111]);
Color := [dark_grey, black, yellow, orange, pink, purple, blue, yellow, orange, pink, purple, blue, yellow, orange, pink, purple, blue, yellow, orange, pink, purple, blue];
t := 5;
Cu := <<1, 1.5> | <1.5, 4>>;
Cv := <<2, 1> | <1, 2>>;
A := 0.9*IdentityMatrix(2);
A := Matrix(2, 2, [[0.2289125271, 0.3782702239], [-1.691313155, 0.7436316236]]);
x0 := <5, 10>;
normalq0 := evalf(<<1, 1>/Norm(<1, 1>, 2) | <-5, -4>/Norm(<-5, -4>, 2) | <1, -3>/Norm(<1, -3>, 2)>);
alphaq0 := [seq(DotProduct(Column(normalq0, i), supPointC(Column(normalq0, i), Cu)), i = 1 .. ColumnDimension(normalq0))];
p0 := <supPointC(Column(normalq0, 1), Cu) | supPointC(Column(normalq0, 2), Cu) | supPointC(Column(normalq0, 3), Cu)>;
q0 := NormToVert(normalq0, alphaq0);
NormalQ0 := evalf(<<1, 1>/Norm(<1, 1>, 2) | <-5, -4>/Norm(<-5, -4>, 2) | <1, -3>/Norm(<1, -3>, 2)>);
AlphaQ0 := [seq(DotProduct(Column(NormalQ0, i), supPointC(Column(NormalQ0, i), Cv)), i = 1 .. ColumnDimension(NormalQ0))];
P0 := <supPointC(Column(NormalQ0, 1), Cv) | supPointC(Column(NormalQ0, 2), Cv) | supPointC(Column(NormalQ0, 3), Cv)>;
Q0 := NormToVert(NormalQ0, AlphaQ0);
display(implicitplot((<x | y> . Cu) . <x | y> = 1, x = -2 .. 2, y = -2 .. 2, color = dark_grey, thickness = 3), implicitplot((<x | y> . Cv) . <x | y> = 1, x = -2 .. 2, y = -2 .. 2, color = dark_green, thickness = 3), PolyPlot(p0, dark_grey), PolyPlot(P0, dark_green));
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
Ku := 3;
Kv := 3;
while n_in <> n_out or N_in <> N_out do if n_in <> n_out then res := PolyhedralAlgo(p, norml, alphaq, Cu); p := res[1]; q := res[2]; norml := res[3]; alphaq := res[4]; Ku := Ku + 1; print(display(implicitplot((<x | y> . Cu) . <x | y> = 1, x = -2 .. 2, y = -2 .. 2, color = dark_grey, thickness = 3), PolyPlot(p, Color[2]), PolyPlot(q, Color[1]))); U := p; X1 := [-((1/A) . U)]; plX1_in := [PolyPlot(X1[1], black)]; for N while 1 < MinkowskiFunctional(x0, X1[N]) do X1 := [op(X1), MinkowskiSum((1/A) . (X1[N]), X1[1])]; plX1_in := [op(plX1_in), PolyPlot(X1[N + 1], black)]; end do; n_in := N; U := q; X1_out := [-((1/A) . U)]; plX1_out := [PolyPlot(X1_out[1], dark_grey)]; for N while 1 < MinkowskiFunctional(x0, X1_out[N]) do X1_out := [op(X1_out), MinkowskiSum((1/A) . (X1_out[N]), X1_out[1])]; plX1_out := [op(plX_out), PolyPlot(X1_out[N + 1], Color[N])]; end do; n_out := N; end if; print(n_in, n_out); if N_in <> N_out then res := PolyhedralAlgo(P, Norml, AlphaQ, Cv); P := res[1]; Q := res[2]; Norml := res[3]; AlphaQ := res[4]; Kv := Kv + 1; print(display(implicitplot((<x | y> . Cv) . <x | y> = 1, x = -2 .. 2, y = -2 .. 2, color = dark_green, thickness = 3), PolyPlot(P, Color[2]), PolyPlot(Q, Color[1]))); U := P; X2 := [-((1/A) . V)]; plX2_in := [PolyPlot(X2[1], black)]; for N while t < MinkowskiFunctional(x0, X2[N]) do X2 := [op(X2), UnionPolyhedra((1/A) . (X2[N]), X2[1])]; plX2_in := [op(plX2_in), PolyPlot(X2[N + 1], black)]; end do; N_in := N; U := Q; X2_out := [-((1/A) . V)]; plX2_out := [PolyPlot(X2_out[1], dark_grey)]; for N while t < MinkowskiFunctional(x0, X2_out[N]) do X2_out := [op(X2_out), UnionPolyhedra((1/A) . (X2_out[N]), X2_out[1])]; plX2_out := [op(plX2_out), PolyPlot(X2_out[N + 1], Color[N])]; end do; N_out := N; end if; print(N_in, N_out); end do;
print(Ku, Kv);
U := p;
V := P;
X1 := [-((1/A) . U)];
X2 := [-((1/A) . V)];
plX1 := [PolyPlot(X1[1], pink)];
plX2 := [PolyPlot(t*X2[1], yellow)];
for N while 1 < MinkowskiFunctional(x0, X1[N]) or t < MinkowskiFunctional(x0, X2[N]) do
    X1 := [op(X1), MinkowskiSum((1/A) . (X1[N]), X1[1])];
    X2 := [op(X2), UnionPolyhedra((1/A) . (X2[N]), X2[1])];
    plX2 := [op(plX2), PolyPlot(t*X2[N + 1], yellow)];
    plX1 := [op(plX1), PolyPlot(X1[N + 1], pink)];
    print([N, MinkowskiFunctional(x0, X1[N]), MinkowskiFunctional(x0, X2[N])]);
end do;
print([N, MinkowskiFunctional(x0, X1[N]), MinkowskiFunctional(x0, X2[N])]);
Nmin := N;
display(plX1, plX2, pointplot(x0));
M1 := ColumnDimension(U);
M2 := ColumnDimension(V);
T := [t];
X3 := [T[1]*X2[Nmin]];
plX3 := [PolyPlot(X3[1], magenta)];
inter := [PolyPlot(IntersectPolygons(X1[Nmin], X3[1]), green)];
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
Control := <<0, 0>>;
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
    constr := {eq1[1] = 0, eq1[2] = 0, eq2[1] = 0, eq2[2] = 0, eq3[1] = 0, eq3[2] = 0, add(lambda[i], i = 1 .. L1) = 1 - 0.001*MinkowskiFunctional(x0, X1[Nmin])*(Nmin - k)/(k + 1), add(mu[i], i = 1 .. L2) + add(nu[i], i = 1 .. M2) <= T[k + 1], add(eta[i], i = 1 .. M1) <= 1};
    try
        res := Minimize(add(nu[i], i = 1 .. M2), constr, assume = nonnegative, variables = [op(convert(lambda, list)), op(convert(mu, list)), op(convert(eta, list)), op(convert(nu, list))], initialpoint = initp);
    catch:
        X1 := [op(X1), MinkowskiSum((1/A) . (X1[Nmin]), X1[1])];
        X2 := [op(X2), UnionPolyhedra((1/A) . (X2[Nmin]), X2[1])];
        plX2 := [op(plX2), PolyPlot(t*X2[Nmin + 1], yellow)];
        plX1 := [op(plX1), PolyPlot(X1[Nmin + 1], pink)];
        Nmin := Nmin + 1;
        print('Nmin' = Nmin);
        Trajectory := <x0>;
        Control := <<0, 0>>;
        T := [t];
        X3 := [T[1]*X2[Nmin]];
        plX3 := [PolyPlot(X3[1], magenta)];
        inter := [PolyPlot(IntersectPolygons(X1[Nmin], X3[1]), green)];
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
    plX3 := [op(plX3), PolyPlot(T[k + 2]*X2[Nmin - k - 1], magenta)];
    inter := [op(inter), PolyPlot(IntersectPolygons(X1[Nmin - k - 1], X3[k + 2]), green)];
    x_new := Column(Trajectory, k + 2);
    m1_old := add(lambda_known[i], i = 1 .. L1);
    m2_old := add(mu_known[i], i = 1 .. L2);
    m1 := MinkowskiFunctional(x_new, X1[Nmin - k - 1]);
    m2 := [MinkowskiFunctional(x_new, T[k + 2]*X2[Nmin - k - 1]), MinkowskiFunctional(x_new, X2[Nmin - k - 1])];
    print(Nmin - k, Nmin - k - 1, x_new, m1_old, m2_old, m1, m2, T[k + 2], res[1]);
    pic := [op(pic), display(PolyPlot(X1[Nmin - k], pink), PolyPlot(T[k + 1]*X2[Nmin - k], yellow), PolyPlot(X1[Nmin - k - 1], pink), PolyPlot(T[k + 2]*X2[Nmin - k - 1], magenta), pointplot(xk), pointplot(x_new), inter[k + 2])];
    print(display(pic[k + 1]));
    k := k + 1;
end do;
print('Nmin' = Nmin);
Control := <Control | -(A . (Column(Trajectory, Nmin)))>;
Control := SubMatrix(Control, 1 .. 2, 2 .. Nmin + 1);
Trajectory := <Trajectory | <0, 0>>;
NULL;
display(plX1, plX3, plot(Transpose(Trajectory), style = pointline));
display(plot(Transpose(Trajectory), style = pointline), inter);
T;

NULL;
