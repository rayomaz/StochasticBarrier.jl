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
    
    syms z x1 x2 x3 betasym gamsym real
    EXP = 0;
    
    solver_opt.solver = 'sdpt3';
    
    % System dynamics
    fx = [0.4 + 0.031*x2 + 0.929*x1;
          0.4 + 0.031*x3 + 0.898*x2 + 0.031*x1;
          0.4 + 0.031*x2 + 0.929*x3];
    
    % Setup SOS program
    prog = sosprogram([x1, x2, x3], [betasym, gamsym]);
    Zmon = monomials([x1, x2, x3], 0:deg);
    
    [prog, B] = sospolyvar(prog, Zmon, 'wscoeff');
    [prog, sig_u1] = sospolyvar(prog, Zmon);
    [prog, sig_u2] = sospolyvar(prog, Zmon);
    [prog, sig_u3] = sospolyvar(prog, Zmon);
    [prog, sig_x]  = sospolyvar(prog, Zmon);
    [prog, sig_o1] = sospolyvar(prog, Zmon);
    [prog, sig_o2] = sospolyvar(prog, Zmon);
    [prog, sig_o3] = sospolyvar(prog, Zmon);
    
    % Basic inequalities
    prog = sosineq(prog, betasym);
    prog = sosineq(prog, sig_u1);
    prog = sosineq(prog, sig_u2);
    prog = sosineq(prog, sig_u3);
    prog = sosineq(prog, sig_x);
    prog = sosineq(prog, sig_o1);
    prog = sosineq(prog, sig_o2);
    prog = sosineq(prog, sig_o3);
    prog = sosineq(prog, B);
    
    % Unsafe set constraints: x in [17,29]^3
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
    stdvar = gx;
    x1 = fx(1) + z;
    x2 = fx(2) + z;
    x3 = fx(3) + z;
    
    Bsub = expand(subs(B));
    clear x1 x2 x3;
    syms x1 x2 x3 real;
    termlist = children(Bsub);
    
    % Compute expected value
    for ii = 1:length(termlist)
        zcount = 0; x1count = 0; x2count = 0; x3count = 0; EXPz = 0;
        factored = cell2sym(termlist(ii));
        factoredterm = factor(factored);
        
        for jj = 1:length(factoredterm)
            if isequaln(factoredterm(jj),z), zcount = zcount + 1; end
            if isequaln(factoredterm(jj),x1), x1count = x1count + 1; end
            if isequaln(factoredterm(jj),x2), x2count = x2count + 1; end
            if isequaln(factoredterm(jj),x3), x3count = x3count + 1; end
        end
        
        if zcount == 0
            EXPz = factored;
        elseif mod(zcount,2) == 0
            EXPz = prod(factoredterm(factoredterm~=z)) * prod(1:2:zcount) * stdvar^zcount;
        end
        
        EXP = EXP + EXPz;
    end
    
    % Expectation constraint
    prog = sosineq(prog, -EXP + B/alpha + betasym - sig_x*(29^2 - x1^2 - x2^2 - x3^2));
    
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
