%% SOS Two Tank
clc; clear all;

deg_list = [4, 6, 8];
results = struct('degree', [], 'Bx', [], 'betaval', [], 'gam', [], 'Ps', [], 'time', []);

for kk = 1:length(deg_list)
    deg = deg_list(kk);
    
    % Start timer for this experiment
    tic
    [Bxpolys, betaval, gam, Ps] = runSOS2D(deg);
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
function [Bxpolys, betaval, gam, Ps] = runSOS2D(deg)
    alpha = 1; gx = 0.01; N = 10; x0 = [0;0];
    
    syms z x1 x2 betasym gamsym real
    EXP = 0;
    
    solver_opt.solver = 'sdpt3';
    fx = [0.90*x1 + 0.1*x2 + 0.45; 0.90*x2 - 0.3];
    
    prog = sosprogram([x1, x2], [betasym, gamsym]);
    Zmon = monomials([x1, x2], 0:deg);
    [prog, B] = sospolyvar(prog, Zmon, 'wscoeff');
    [prog, sig_u1] = sospolyvar(prog, Zmon);
    [prog, sig_u2] = sospolyvar(prog, Zmon);
    [prog, sig_x]  = sospolyvar(prog, Zmon);
    [prog, sig_o1] = sospolyvar(prog, Zmon);
    [prog, sig_o2] = sospolyvar(prog, Zmon);
    
    prog = sosineq(prog, betasym);
    prog = sosineq(prog, sig_u1);
    prog = sosineq(prog, sig_u2);
    prog = sosineq(prog, sig_x);
    prog = sosineq(prog, sig_o1);
    prog = sosineq(prog, B);
    
    prog = sosineq(prog, B - sig_u1*(x1 - 1)*(9 - x1) - 1);
    prog = sosineq(prog, B - sig_u2*(x2 - 1)*(9 - x2) - 1);
    
    prog = sosineq(prog, -B - sig_o1*(x1 - 2.75)*(3.25 - x1) + gamsym);
    prog = sosineq(prog, -B - sig_o2*(x2 - 2.75)*(3.25 - x2) + gamsym);
    
    prog = sosineq(prog, gamsym);
    prog = sosineq(prog, 1 - gamsym - 1e-6);
    
    stdvar = gx;
    x1 = fx(1) + z;
    x2 = fx(2) + z;
    
    Bsub = expand(subs(B));
    clear x1 x2;
    syms x1 x2 real;
    termlist = children(Bsub);
    
    for ii = 1:length(termlist)
        zcount = 0; x1count = 0; x2count = 0; EXPz = 0;
        factored = cell2sym(termlist(ii));
        factoredterm = factor(factored);
        
        for jj = 1:length(factoredterm)
            if isequaln(factoredterm(jj),z), zcount = zcount + 1; end
            if isequaln(factoredterm(jj),x1), x1count = x1count + 1; end
            if isequaln(factoredterm(jj),x2), x2count = x2count + 1; end
        end
        
        if zcount == 0
            EXPz = factored;
        elseif mod(zcount,2) == 0
            EXPz = prod(factoredterm(factoredterm~=z)) * prod(1:2:zcount) * stdvar^zcount;
        end
        
        EXP = EXP + EXPz;
    end
    
    prog = sosineq(prog, -EXP + B/alpha + betasym - sig_x*(2^2 - x1^2 - x2^2));
    
    objfunc = gamsym + betasym;
    prog = sossetobj(prog, objfunc);
    
    prog = sossolve(prog, solver_opt);
    Bxpolys = sosgetsol(prog, B);
    betaval = double(sosgetsol(prog, betasym));
    gam = double(sosgetsol(prog, gamsym));
    
    x1 = x0(1); x2 = x0(2);
    if alpha == 1
        probvalue = double(subs(Bxpolys) + betaval*N);
    end
    
    Ps = 1 - gam - betaval*N;
end
