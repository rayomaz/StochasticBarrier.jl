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
    
    % declare separate noise variables z1,z2 and x1,x2
    syms z1 z2 x1 x2 betasym gamsym real
    EXP = 0;
    
    solver_opt.solver = 'sdpt3';
    fx = [0.90*x1 + 0.1*x2 + 0.45; 0.90*x2 - 0.3];
    
    prog = sosprogram([x1, x2], [betasym, gamsym]);
    Zmon = monomials([x1, x2], 0:deg);
    [prog, B]       = sospolyvar(prog, Zmon, 'wscoeff');
    [prog, sig_u1]  = sospolyvar(prog, Zmon);
    [prog, sig_u2]  = sospolyvar(prog, Zmon);
    [prog, sig_x1]  = sospolyvar(prog, Zmon);   % separate multipliers for box
    [prog, sig_x2]  = sospolyvar(prog, Zmon);
    [prog, sig_o1]  = sospolyvar(prog, Zmon);
    [prog, sig_o2]  = sospolyvar(prog, Zmon);
    
    % SOS nonnegativity
    prog = sosineq(prog, betasym);
    prog = sosineq(prog, sig_u1);
    prog = sosineq(prog, sig_u2);
    prog = sosineq(prog, sig_x1);
    prog = sosineq(prog, sig_x2);
    prog = sosineq(prog, sig_o1);
    prog = sosineq(prog, B);
    
    prog = sosineq(prog, B - sig_u1*(x1 - 1)*(9 - x1) - 1);
    prog = sosineq(prog, B - sig_u2*(x2 - 1)*(9 - x2) - 1);
    
    prog = sosineq(prog, -B - sig_o1*(x1 - 2.75)*(3.25 - x1) + gamsym);
    prog = sosineq(prog, -B - sig_o2*(x2 - 2.75)*(3.25 - x2) + gamsym);
    
    prog = sosineq(prog, gamsym);
    prog = sosineq(prog, 1 - gamsym - 1e-6);
    
    % separate noise magnitudes (allow different values if desired)
    stdvar1 = gx;
    stdvar2 = gx;
    
    % substitute independent noises into each axis
    x1 = fx(1) + z1;
    x2 = fx(2) + z2;
    
    Bsub = expand(subs(B));
    clear x1 x2;
    syms x1 x2 real;
    termlist = children(Bsub);
    
    for ii = 1:length(termlist)
        z1count = 0; z2count = 0; x1count = 0; x2count = 0; EXPz = 0;
        factored = cell2sym(termlist(ii));
        factoredterm = factor(factored);
        
        for jj = 1:length(factoredterm)
            if isequaln(factoredterm(jj), z1), z1count = z1count + 1; end
            if isequaln(factoredterm(jj), z2), z2count = z2count + 1; end
            if isequaln(factoredterm(jj), x1), x1count = x1count + 1; end
            if isequaln(factoredterm(jj), x2), x2count = x2count + 1; end
        end
        
        % expectation logic for independent Gaussian-like noises:
        % odd power on any axis -> zero; otherwise multiply even moments
        if (mod(z1count,2) == 1) || (mod(z2count,2) == 1)
            EXPz = 0;
        elseif (z1count == 0) && (z2count == 0)
            EXPz = factored;
        elseif z1count == 0
            % only z2 present (even)
            nonNoise = factoredterm(~(factoredterm == z2));
            EXPz = prod(nonNoise) * prod(1:2:z2count) * stdvar2^z2count;
        elseif z2count == 0
            % only z1 present (even)
            nonNoise = factoredterm(~(factoredterm == z1));
            EXPz = prod(nonNoise) * prod(1:2:z1count) * stdvar1^z1count;
        else
            % both present with even powers
            nonNoise = factoredterm(~( (factoredterm == z1) | (factoredterm == z2) ));
            EXPz = prod(nonNoise) * ...
                   prod(1:2:z1count) * stdvar1^z1count * ...
                   prod(1:2:z2count) * stdvar2^z2count;
        end
        
        EXP = EXP + EXPz;
    end
    
    prog = sosineq(prog, -EXP + B/alpha + betasym ...
        - sig_x1*(x1 - 1)*(9 - x1) ...
        - sig_x2*(x2 - 1)*(9 - x2));
    
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
