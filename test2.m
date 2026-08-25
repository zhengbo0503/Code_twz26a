%TEST2 -- forward accuracy test, changing number of columns

clc; close all; clear; rng(1);

m = ones(15,1)*1000;
n = round(logspace(1,3,15));
epsln = eps('double')/2;
kappa = 1e8;
f1 = zeros(length(m), 5);
f2 = zeros(length(m), 5);
f3 = zeros(length(m), 5);
f4 = zeros(length(m), 5);
bound1 = zeros(length(m),5);
bound2 = zeros(length(m),5);

for mode = 1:5

    for i = 1:length(m)

        mm = m(i);
        nn = n(i);
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

        fprintf("Finished MODE = %d, %d of %d\n", mode, i, length(m));

    end
end

savedata = 1;
if savedata == 1
    save("./data/fwderr_diff_cols.mat")
end


%% 
close all; 

C1 = "#1171BE";
C2 = "#DD5400";
C3 = "#EDB120";
C4 = "#3BAA32";

for mode = 1:5
    figure(mode)
    
    loglog(n,f1(:,mode),'LineStyle','none','Marker','*','Color',C1);
    hold on;
    loglog(n,f2(:,mode),'LineStyle','none','Marker','pentagram','Color',C2);
    loglog(n,f3(:,mode),'LineStyle','none','Marker','square','Color',C3);
    loglog(n,f4(:,mode),'LineStyle','none','Marker','diamond','Color',C4);
    loglog(n,bound2(:,mode),'LineStyle',':','Marker','none','Color','k');
    set(findall(gcf, 'Type', 'Line'), 'LineWidth', 1);
    %

    axis square
    xlim([10,1000]);
    xlabel('$n$', 'FontSize', 10);
    yticks([1e-16,1e-13,1e-10,1e-7,1e-4])
    ylim([1e-16, 1e-4]);
    % Only get ylabel if the subplot is on the left
    if mod(mode,2) == 1
        ylabel('$\mathrm{max}_k {\varepsilon}^{(k)}_{fwd}$', 'FontSize', 10);
    end

    % Set the title 
    alphabet = ['a','b','c','d','e'];
    t = sprintf('(%s) MODE = %d', alphabet(mode), mode); 
    title(t, 'FontWeight', 'normal', 'FontSize', 10); 
end
