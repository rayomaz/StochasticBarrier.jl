%% SOS 3D Polynomial System
clc; clear all;

deg_list = [2, 4, 6, 8, 12];
results = struct('degree', [], 'Bx', [], 'betaval', [], 'gam', [], 'Ps', [], 'time', []);

for kk = 1:length(deg_list)
    deg = deg_list(kk);
    
    % Start timer for this experiment
    tic
    [Bxpolys, betaval, gam, Ps] = runSOS3D(deg);
    elapsedTime = toc;
    
    results(kk).degree = deg;
    results(kk).Bx = Bxpolys;
    results(kk).betaval = betaval;
    results(kk).gam = gam;
    results(kk).Ps = Ps;
    results(kk).time = elapsedTime;
    
    fprintf('Degree: %d, gam = %.4f, betaval = %.4f, Ps = %.4f, time = %.2f s\n', ...
        deg, gam, betaval, Ps, elapsedTime);
end

%% Function definition
function [Bxpolys, betaval, gam, Ps] = runSOS3D(deg)
    alpha = 1; gx = 0.1; N = 10; x0 = [18;18;18];
    
    % Symbols
    syms z1 z2 z3 x1 x2 x3 betasym gamsym real
    EXP = 0;
    
    solver_opt.solver = 'sdpt3';
    
    % System dynamics
    fx = [0.4 + 0.031*x2 + 0.929*x1;
          0.4 + 0.031*x3 + 0.898*x2 + 0.031*x1;
          0.4 + 0.031*x2 + 0.929*x3];
    
    % Setup SOS program
    prog = sosprogram([x1, x2, x3], [betasym, gamsym]);
    Zmon = monomials([x1, x2, x3], 0:deg);
    
    [prog, B]       = sospolyvar(prog, Zmon, 'wscoeff');
    [prog, sig_u1]  = sospolyvar(prog, Zmon);
    [prog, sig_u2]  = sospolyvar(prog, Zmon);
    [prog, sig_u3]  = sospolyvar(prog, Zmon);
    [prog, sig_x1]  = sospolyvar(prog, Zmon);
    [prog, sig_x2]  = sospolyvar(prog, Zmon);
    [prog, sig_x3]  = sospolyvar(prog, Zmon);
    [prog, sig_o1]  = sospolyvar(prog, Zmon);
    [prog, sig_o2]  = sospolyvar(prog, Zmon);
    [prog, sig_o3]  = sospolyvar(prog, Zmon);
    
    % Nonnegativity
    prog = sosineq(prog, betasym);
    prog = sosineq(prog, sig_u1); prog = sosineq(prog, sig_u2); prog = sosineq(prog, sig_u3);
    prog = sosineq(prog, sig_x1); prog = sosineq(prog, sig_x2); prog = sosineq(prog, sig_x3);
    prog = sosineq(prog, sig_o1); prog = sosineq(prog, sig_o2); prog = sosineq(prog, sig_o3);
    prog = sosineq(prog, B);
    
    % Unsafe set: x notin [17,29]^3
    prog = sosineq(prog, B - sig_u1*(x1 - 17)*(29 - x1) - 1);
    prog = sosineq(prog, B - sig_u2*(x2 - 17)*(29 - x2) - 1);
    prog = sosineq(prog, B - sig_u3*(x3 - 17)*(29 - x3) - 1);
    
    % Initial region: x in [18,19]^3
    prog = sosineq(prog, -B - sig_o1*(x1 - 18)*(19 - x1) + gamsym);
    prog = sosineq(prog, -B - sig_o2*(x2 - 18)*(19 - x2) + gamsym);
    prog = sosineq(prog, -B - sig_o3*(x3 - 18)*(19 - x3) + gamsym);
    
    prog = sosineq(prog, gamsym);
    prog = sosineq(prog, 1 - gamsym - 1e-6);
    
    % Noise substitution
    stdvar1 = gx; stdvar2 = gx; stdvar3 = gx;
    x1 = fx(1) + z1; x2 = fx(2) + z2; x3 = fx(3) + z3;
    
    Bsub = expand(subs(B));
    clear x1 x2 x3;
    syms x1 x2 x3 real;
    termlist = children(Bsub);
    
    % Compute expected value term-by-term
    for ii = 1:length(termlist)
        z1count = 0; z2count = 0; z3count = 0;
        EXPz = 0;
        factored = cell2sym(termlist(ii));
        factoredterm = factor(factored);
        
        for jj = 1:length(factoredterm)
            if isequaln(factoredterm(jj), z1), z1count = z1count + 1; end
            if isequaln(factoredterm(jj), z2), z2count = z2count + 1; end
            if isequaln(factoredterm(jj), z3), z3count = z3count + 1; end
        end
        
        if mod(z1count,2)==1 || mod(z2count,2)==1 || mod(z3count,2)==1
            EXPz = 0;
        else
            nonNoise = factoredterm(~ismember(factoredterm,[z1 z2 z3]));
            EXPz = prod(nonNoise) * prod(1:2:z1count)*stdvar1^z1count ...
                                  * prod(1:2:z2count)*stdvar2^z2count ...
                                  * prod(1:2:z3count)*stdvar3^z3count;
        end
        EXP = EXP + EXPz;
    end
    
    % Hypercube expectation constraint X = [17,29]^3
    prog = sosineq(prog, -EXP + B/alpha + betasym ...
        - sig_x1*(x1 - 17)*(29 - x1) ...
        - sig_x2*(x2 - 17)*(29 - x2) ...
        - sig_x3*(x3 - 17)*(29 - x3));
    
    % Objective
    objfunc = gamsym + betasym;
    prog = sossetobj(prog, objfunc);
    
    % Solve
    prog = sossolve(prog, solver_opt);
    Bxpolys = sosgetsol(prog, B);
    betaval = double(sosgetsol(prog, betasym));
    gam = double(sosgetsol(prog, gamsym));
    
    x1 = x0(1); x2 = x0(2); x3 = x0(3);
    if alpha == 1
        probvalue = double(subs(Bxpolys) + betaval*N);
    end
    
    Ps = 1 - gam - betaval*N;
end
