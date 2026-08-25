%TEST1 - forward accuracy test, changing condition number

clc; close all; clear; rng(1);

m = 1000;
n = 800;
epsln = eps('double')/2;
kkappa = logspace(3,15,20);
f1 = zeros(length(m), 5);
f2 = zeros(length(m), 5);
f3 = zeros(length(m), 5);
f4 = zeros(length(m), 5);
bound1 = zeros(length(m),5);
bound2 = zeros(length(m),5);

for mode = 1:5

    for i = 1:length(kkappa)

        mm = m;
        nn = n;
        kappa = kkappa(i); 

        A = gallery('randsvd', [mm,nn], kappa, mode);
        sref = reference_singular_values(A);

        [U1,S1,V1,nos1,scnd] = mposj(A, 3, true);
        [f1(i,mode),~,~,~] = compute_error(A, U1, S1, V1, sref);

        [U2,S2,V2,sva2,work2,info2] = dgesvj_mex(A,'G','U','V',nn,eye(nn),max(6,mm+nn));
        if info2 ~= 0
            fprintf("Error: DGESVJ does not converge.\n");
            break;
        end
        [f2(i,mode),~,~,~] = compute_error(A, U2, S2, V2, sref);

        [U3,S3,V3,sva3,work3,iwork3,info3] = dgejsv_mex(A,'C','U','V','R','N','N');
        if info3 ~= 0
            fprintf("Error: DGEJSV does not converge.\n");
            break;
        end
        [f3(i,mode),~,~,~] = compute_error(A, U3, S3, V3, sref);

        [U4,S4,V4] = svd(A,'econ');
        [f4(i,mode),~,~,~] = compute_error(A, U4, S4, V4, sref);

        bound1(i,mode) = scond(A, 'C') * sqrt(mm * nn)* epsln;
        bound2(i,mode) = scnd * sqrt(mm * nn) * epsln;

        fprintf("Finished MODE = %d, %d of %d\n", mode, i, length(kkappa));

    end
end

savedata = 1;
if savedata == 1
    save("./data/fwderr_diff_cond.mat");
end

%% 
close all; 
C1 = "#1171BE";
C2 = "#DD5400";
C3 = "#EDB120";
C4 = "#3BAA32";
    

for mode = 1:5
    figure(mode)
    
    % %%%%%%%%%%%%%%%%%%
    loglog(kkappa,f1(:,mode),'LineStyle','none','Marker','*', 'Color', C1); hold on;
    loglog(kkappa,f2(:,mode),'LineStyle','none','Marker','pentagram', 'Color', C2);
    loglog(kkappa,f3(:,mode),'LineStyle','none','Marker','square', 'Color', C3);
    loglog(kkappa,f4(:,mode),'LineStyle','none','Marker','diamond', 'Color', C4);
    loglog(kkappa,bound1(:,mode),'LineStyle','--','Marker', 'none', 'Color','k');
    loglog(kkappa,bound2(:,mode),'LineStyle',':','Marker', 'none', 'Color','k');
    set(findall(gcf, 'Type', 'Line'), 'LineWidth', 1);
    % %%%%%%%%%%%%%%%%%%

    axis square
    
    xlabel('Condition number $\kappa_2(A)$','FontSize',10);
    xlim([1e3,1e15]);
    xticks([1e3,1e6,1e9,1e12,1e15]);
    ylim([1e-16,1e2]);
    yticks([1e-16, 1e-13, 1e-10,1e-7,1e-4, 1e-1, 1e2]);
    
    % Only get ylabel if the subplot is on the left
    if mod(mode,2) == 1
        ylabel('$\mathrm{max}_k {\varepsilon}^{(k)}_{fwd}$','FontSize',10);
    end

    % Set the title 
    alphabet = ['a','b','c','d','e'];
    t = sprintf('(%s) MODE = %d', alphabet(mode), mode); 
    title(t, 'FontWeight', 'normal', 'FontSize', 10); 
end
