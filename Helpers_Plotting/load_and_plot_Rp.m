function ReturnStruct = load_and_plot_Rp(basePath, Target_dir, opts)
    if nargin < 3
        opts = struct('fontsize', 16);
        opts.scalar = 1;
        opts.Pe_range = linspace(0.0001, 2, 50);
        opts.Da_range = linspace(0.0001, 2, 50);
        opts.ShowPlots = true;
    end
    FONTSIZE = opts.fontsize;
    scalar = opts.scalar;
    Pe_range = opts.Pe_range .* scalar;
    Da_range = opts.Da_range .* scalar;

    % Construct the expected full‐path to ANumeric_Rp.mat
    folderPath = fullfile(basePath, Target_dir);
    matFile    = fullfile(folderPath, "ANumeric_Rp.mat");

    % If the file does not exist, print a message, return empty
    if ~exist(matFile, "file")
        fprintf(">> load_and_plot_Rp: File not found:\n   %s\n", matFile);
        Rp_array = [];
        return;
    end

    % Attempt to load the .mat file
    S = load(matFile);
    Burnout = 2
    % Check that saveStruct.Rp_array actually exists
    if isfield(S, "saveStruct") && isfield(S.saveStruct, "Rp_array")
        Rp_array = S.saveStruct.Rp_array;

        Rp_inner = Rp_array((1+Burnout):end-1, (1+Burnout):end-1);
        [~, idx_lin] = max(Rp_inner(:));
        [i_local, j_local] = ind2sub(size(Rp_inner), idx_lin);
        i_max = i_local + Burnout;    % actual row in Rp_array
        j_max = j_local + Burnout;    % actual column in Rp_array
        Pe_max = Pe_range(i_max);
        Da_max = Da_range(j_max);
        ReturnStruct = struct('Pe_max', Pe_max, "Da_max", Da_max);

        if opts.ShowPlots
            % Plot Rp_array as a heatmap
            fig = figure('position', [1 1, 1460, 840]);
            hold on;
            colormap('hot');
            %imagesc(Da_range((1+Burnout):end-1), Pe_range((1+Burnout):end-1), Rp_inner ./ max(max(Rp_inner)));
            contourf(Da_range((1+Burnout):end-1), Pe_range((1+Burnout):end-1), Rp_inner./max(max(Rp_inner)), 'linestyle', ':');
            line([Da_max, Da_max], [Pe_range(Burnout) Pe_range(end-1)],'linewidth', 2.5, 'linestyle', ':', 'color', 'g');
            line([Da_range(Burnout) Da_range(end-1)], [Pe_max Pe_max],'linewidth', 2.5, 'linestyle', ':', 'color', 'g');
            plot(Da_max, Pe_max, 'marker' ,'*', 'linewidth', 2, 'linestyle', ':', 'color', 'k')
            hold off;
            set(gca, "YDir", "normal");     % so that row‐1 is at bottom
            for ax = findobj(gcf,'Type','Axes')'
                set(ax,'fontsize', FONTSIZE-4);
            end
            xlabel("Da", 'Fontsize', FONTSIZE, 'interpreter', 'latex');
            ylabel("Pe",'rotation',0,'Fontsize', FONTSIZE, 'interpreter', 'latex');
            %title("Reaction production $R_p$", 'Fontsize', FONTSIZE, 'interpreter', 'latex');
            xlim([Da_range(1+Burnout) Da_range(end-1)]);
            ylim([Pe_range(1+Burnout) Pe_range(end-1)]);
            subR = Rp_inner((1+Burnout):(end-1), (1+Burnout):(end-1));
            %clim = [ min(subR(:)),  max(subR(:)) ];
            caxis([0.4 1]);
            cb = colorbar();
            for ax = findobj(gcf,'Type','Axes')'
                set(ax,'fontsize', FONTSIZE-4);
            end
            legend('$R_p$', 'Optimal Pe/Da','fontsize',FONTSIZE, 'interpreter', 'latex')
            %set(get(cb, "label"), "string", "$R_p$", "Fontsize", FONTSIZE, "rotation", 0, "interpreter", "latex")
            waitforbuttonpress;
        end
    else
        error("load_and_plot_Rp: the file %s does not contain saveStruct.Rp_array.", matFile);
    end
end
    