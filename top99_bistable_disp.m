function top99_bistable_disp(nelx,nely,volfrac,penal,rmin)
%%%%%%%Parameter
nelx    = 150; 
nely    = 50;  
volfrac = 0.4; 
penal   = 3.0; 
rmin    = 3.0;
%%%%%%%%%%%%%%
lx = 150.0;  
ly = 25.0;   % full design domain
dx = lx / nelx;
dy = ly / nely;
fprintf('Full domain: %.1f x %.1f mm, mesh: %d x %d -> element: %.4f x %.4f mm\n', ...
        lx, ly, nelx, nely, dx, dy);

%%%%%%%%% Mechanical Property
E0   = 26.0;        
Emin = 0.001 * E0;     
nu   = 0.40;
%%%%%%%%%%%%%%%%

%%%%%%%%%%% displacement of two stable stage
s1 = 0.45 * ly;
s2 = 1.2  * lx;
fprintf('Using s1 = %.4f mm, s2 = %.4f mm\n', s1, s2);

x = volfrac * ones(nely, nelx);   
loop = 0; change = 1.0;
%%%%%%%%%%

while change > 1e-2&&loop<100
    loop = loop + 1;
    xold = x;

    [K, KE] = assembleK(nelx, nely, x, penal, nu, E0, Emin);

    U1 = solve_disp_case(K, nelx, nely, s1);  
    U2 = solve_disp_case(K, nelx, nely, s2);

    ce1 = zeros(nely, nelx); ce2 = zeros(nely, nelx);
    c = 0.0;
    for ely = 1:nely
        for elx = 1:nelx
            n1 = (nely+1)*(elx-1) + ely;
            n2 = (nely+1)* elx      + ely;
            edofs = [2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2];
            Ue1 = U1(edofs); Ue2 = U2(edofs);
            ce1(ely,elx) = Ue1' * KE * Ue1;
            ce2(ely,elx) = Ue2' * KE * Ue2;
            Ex = Emin + x(ely,elx)^penal * (E0 - Emin);
            c = c + Ex * (ce2(ely,elx)/s2 - ce1(ely,elx)/s1);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%% OC
    dc = zeros(nely, nelx);
    for ely = 1:nely
        for elx = 1:nelx
            sraw = (E0 - Emin) * penal * x(ely,elx)^(penal-1) * ( ce2(ely,elx)/s2 - ce1(ely,elx)/s1 );
            dc(ely,elx) = - sraw;    % negative sign to match OC usage
        end
    end

    % filter
    dc = check(nelx,nely,rmin,x,dc);

    % OC update
    x = OC(nelx,nely,x,volfrac,dc);

    % report and plot
    change = max(abs(x(:) - xold(:)));
    fprintf('It: %3d  Obj(f2-f1): %12.6f  Vol: %6.3f  ch: %7.5f\n', loop, c, mean(x(:)), change);
    colormap(gray); imagesc(-x); axis equal tight off; drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[K, KE] = assembleK(nelx, nely, x, penal, nu, E0, Emin);

s_min = max(1e-3, 0.01 * Ly);  
s_max = 1.5 * s2;
nsamp = 101;
s_vals = linspace(s_min, s_max, nsamp);
f_vals = zeros(size(s_vals));

for i = 1:nsamp
    sval = s_vals(i);
    U = solve_disp_case(K, nelx, nely, sval);
    f_vals(i) = (U' * (K * U)) / sval;
end

% Plot
figure;
plot(s_vals, f_vals, 'b-', 'LineWidth', 1.5); hold on;
plot([s1 s1], ylim, 'k--', 'LineWidth', 1); 
plot([s2 s2], ylim, 'k--', 'LineWidth', 1);  
plot(s1, interp1(s_vals,f_vals,s1), 'ro', 'MarkerFaceColor','r');
plot(s2, interp1(s_vals,f_vals,s2), 'go', 'MarkerFaceColor','g');
xlabel('Prescribed displacement s (mm)');
ylabel('Reaction force f(s) (units consistent with E)');
title('Force-Displacement curve f(s) for optimized design');
legend('f(s)','s1','s2','f(s1)','f(s2)','Location','Best');
grid on;
hold off;

function [K, KE] = assembleK(nelx,nely,x,penal,nu,E0,Emin)
    KE = lk(nu);                    
    ndof = 2 * (nelx+1) * (nely+1);
    K = sparse(ndof, ndof);
    for ex = 1:nelx
        for ey = 1:nely
            n1 = (nely+1)*(ex-1) + ey;
            n2 = (nely+1)* ex      + ey;
            edofs = [2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2];
            Ex = Emin + x(ey,ex)^penal * (E0 - Emin);
            K(edofs, edofs) = K(edofs, edofs) + Ex * KE;
        end
    end
end

function U = solve_disp_case(K, nelx, nely, sval)
    ndof = size(K,1);
    U = zeros(ndof,1);
    Ny = nely+1; Nx = nelx+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BC
    left_nodes = 1:Ny;
    left_dofs = [2*left_nodes-1, 2*left_nodes]'; left_dofs = left_dofs(:);
    right_nodes = Ny * nelx + (1:Ny);
    right_dofs = [2*right_nodes-1, 2*right_nodes]'; right_dofs = right_dofs(:);
    fixeddofs = unique([left_dofs; right_dofs]);

    mid_col = ceil(Nx/2);
    mid_node = Ny * (mid_col - 1) + Ny;  
    mid_Uy = 2 * mid_node;              
    U(mid_Uy) = -abs(sval);                

    fixeddofs = unique([fixeddofs; mid_Uy]);

    alldofs = (1:ndof).';
    freedofs = setdiff(alldofs, fixeddofs);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U(freedofs) = K(freedofs,freedofs) \ ( - K(freedofs, fixeddofs) * U(fixeddofs) );
    U(fixeddofs) = U(fixeddofs);  % keep prescribed values (others zero)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xnew = OC(nelx,nely,x,volfrac,dc)
    l1 = 0; l2 = 1e5; move = 0.2;
    while (l2 - l1 > 1e-4)
        lmid  = 0.5 * (l2 + l1);
        xcand = x .* sqrt(-dc ./ lmid);
        xnew = max(0.001, max(x - move, min(1.0, min(x + move, xcand))));
        if sum(xnew(:)) - volfrac * nelx * nely > 0
            l1 = lmid;
        else
            l2 = lmid;
        end
    end
end

function dcn = check(nelx,nely,rmin,x,dc)
    dcn = zeros(nely,nelx);
    for i = 1:nelx
        for j = 1:nely
            sumw = 0.0; val = 0.0;
            i1 = max(i - floor(rmin), 1);
            i2 = min(i + floor(rmin), nelx);
            j1 = max(j - floor(rmin), 1);
            j2 = min(j + floor(rmin), nely);
            for k = i1:i2
                for l = j1:j2
                    fac = rmin - sqrt((i-k)^2 + (j-l)^2);
                    if fac > 0
                        w = fac;
                        sumw = sumw + w;
                        val  = val  + w * x(l,k) * dc(l,k);
                    end
                end
            end
            dcn(j,i) = val / max(1e-9, x(j,i)*sumw);
        end
    end
end

function KE = lk(nu)
    k = [ 1/2 - nu/6,  1/8 + nu/8,  -1/4 - nu/12, -1/8 + 3*nu/8, ...
         -1/4 + nu/12, -1/8 - nu/8,  nu/6,          1/8 - 3*nu/8 ];
    KE = 1/(1 - nu^2) * [ ...
        k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8);
        k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3);
        k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2);
        k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5);
        k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4);
        k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7);
        k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6);
        k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1) ];
end

end
