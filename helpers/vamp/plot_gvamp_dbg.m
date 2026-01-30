function plot_gvamp_dbg(VAMP_dbg_mode, EXIT_dbg_mode, n_trig, IqYX_itvec, IqYX, i, ...
        dbg_alpha1_alg, dbg_alpha1_mc, ...
        dbg_alpha2_alg, dbg_alpha2_mc, ...
        dbg_nu_U1_alg, dbg_nu_U1_mc, ...
        dbg_nu_U2_alg, dbg_nu_U2_mc, ...
        dbg_beta1_alg, dbg_beta1_mc, ...
        dbg_beta2_alg, dbg_beta2_mc, ...
        dbg_nu_W1_alg, dbg_nu_W1_mc, ...
        dbg_nu_W2_alg, dbg_nu_W2_mc)

    % Plot variance trajectory over EXIT chart
    if EXIT_dbg_mode == 1
        figure(843);
        iteroff = 1;
        hold on;

        nu_W1_up = upsample(dbg_nu_W1_mc, 2);
        nu_W2_up = upsample(dbg_nu_W2_mc, 2);

        tmp_w1 = nu_W1_up + circshift(nu_W1_up, 1);
        traj_w1 = tmp_w1(2:end); % Remove first index

        tmp_w2 = nu_W2_up + circshift(nu_W2_up, 1);
        traj_w2 = tmp_w2(1:end - 1); % Remove last index

        loglog(traj_w1(iteroff:end), traj_w2(iteroff:end), 'k.-', 'MarkerSize', 8);
        loglog(traj_w1(n_trig * 2), traj_w2(n_trig * 2), 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', 'black'); %iteration where annealing is turned off
        loglog(traj_w1(end), traj_w2(end), 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', 'black');

        text(0.8 * traj_w1(n_trig * 2), 1.3 * traj_w2(n_trig * 2), ...
            strcat("Rate with annealing: ", num2str(IqYX_itvec(n_trig)), " bpcu"), ...
            'Color', 'black', 'FontSize', 6, 'HorizontalAlignment', 'right');

        text(0.85 * traj_w1(end), 1.4 * traj_w2(end), ...
            strcat("Actual Rate: ", num2str(IqYX), " bpcu"), ...
            'Color', 'black', 'FontSize', 6, 'HorizontalAlignment', 'right');
    end

    % Plot rate trajectory
    if VAMP_dbg_mode >= 1
        figure(873488);
        plot(1:i, IqYX_itvec(1:i), 'bo-','MarkerFaceColor','auto', 'MarkerSize',3);
        xlabel('Iteration');
        ylabel('Rate')
        xlim([1, i]);
        title('Information Rates')
        grid on;
        hold on;
        drawnow;
    end

    % Detailed debugging plots
    if VAMP_dbg_mode == 2
        figure(5356);
        clf();

        subplot(4, 2, 1);
        semilogy(1:i, dbg_alpha1_alg, 'r.-'); hold on;
        semilogy(1:i, dbg_alpha1_mc, 'b.-');
        legend('Algorithm', 'Ground Truth');
        title('alpha1');
        xlabel('Iterations');
        grid on;

        subplot(4, 2, 2);
        title('alpha2');
        %semilogy(1:i, dbg_alpha2_alg, 'r.-'); hold on; semilogy(1:i,dbg_alpha2_mc, 'b.-'); %Not relevant
        xlabel('Iterations');
        grid on;

        subplot(4, 2, 3);
        semilogy(1:i, dbg_nu_U1_alg, 'r.-'); hold on;
        semilogy(1:i, dbg_nu_U1_mc, 'b.-');
        title('nuU1');
        xlabel('Iterations');
        grid on;

        subplot(4, 2, 4);
        semilogy(1:i, dbg_nu_U2_alg, 'r.-'); hold on;
        semilogy(1:i, dbg_nu_U2_mc, 'b.-');
        title('nuU2');
        xlabel('Iterations');
        grid on;

        subplot(4, 2, 5);
        semilogy(1:i, dbg_beta1_alg, 'r.-'); hold on;
        semilogy(1:i, dbg_beta1_mc, 'b.-');
        title('beta1');
        xlabel('Iterations');
        grid on;

        subplot(4, 2, 6);
        semilogy(1:i, dbg_beta2_alg, 'r.-'); hold on;
        semilogy(1:i, dbg_beta2_mc, 'b.-');
        title('beta2');
        xlabel('Iterations');
        grid on;

        subplot(4, 2, 7);
        semilogy(1:i, dbg_nu_W1_alg, 'r.-'); hold on;
        semilogy(1:i, dbg_nu_W1_mc, 'b.-');
        title('nuW1');
        xlabel('Iterations');
        grid on;

        subplot(4, 2, 8);
        semilogy(1:i, dbg_nu_W2_alg, 'r.-'); hold on;
        semilogy(1:i, dbg_nu_W2_mc, 'b.-');
        title('nuW2');
        xlabel('Iterations');
        grid on;

        drawnow;
    end
end
