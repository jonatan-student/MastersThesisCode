clc; clear;
addpath("\Helpers_Plotting");
addpath("\Helpers_analytical");

Data_dir = "Final_results";
Test_dir = "SuperLargeScanSlipsAndSlopes";

Target_dir = "Slopes=10";


basePath = fullfile(Data_dir, Test_dir);

if ~isfolder(basePath)
    error('Folder does not exist: %s\n', basePath);
end

% get list of sub‐folders
dirInfo    = dir(basePath);
isSubdir   = [dirInfo.isdir];
allFolders = {dirInfo(isSubdir).name};
allFolders = allFolders(~ismember(allFolders, {'.','..'}));

vals = cellfun(@(s) str2double( s(strfind(s,'=')+1:end) ), allFolders);

[~, order] = sort(vals);
Folders = allFolders(order);

FONTSIZE = 20


% plot_types can be any or all of these -->
% * {"Velocity Field", "Velocity Profile", 
%    * "Pressure in x", "S of x", 
%    * "Concentrations",  "Concentration a Profile"
%    * "Velocity validation"}

Target_Da = 2;
Target_Pe = 2;
Target_test = "conc_Pe2_Da2.mat";

StrIn = {basePath, {Target_dir}, Target_test};
plot_types = {"none"};
[these_stats, Rp, ReScaler, xScalar] = plot_specific(StrIn, plot_types, struct('fontsize', FONTSIZE));


fprintf('\n___________\n------------\nPe: %d, Da: %d\n------------\n', Target_Pe*ReScaler, Target_Da*ReScaler);

Pe_range = linspace(0.0001, 2, 50);
Da_range = linspace(0.0001, 2, 50);

options = struct('Pe_range', Pe_range, 'Da_range', Da_range, 'fontsize', FONTSIZE, 'scalar', ReScaler, 'ShowPlots', true);
m_array = [0,2,4,6,8,10];
Ls_array = [0,4,8,12,16,20];
RpCalcFromLoad = false;
if RpCalcFromLoad

    Rp_max_Pe_m = zeros(1, numel(m_array));
    Rp_max_Da_m = zeros(1, numel(m_array));
    Rp_max_Pe_Ls = zeros(1, numel(m_array));
    Rp_max_Da_Ls = zeros(1, numel(m_array));

    for i = 1:numel(m_array)
        Target_dir = sprintf("Slopes=%d", m_array(i))
        rS = load_and_plot_Rp(basePath, Target_dir, options);
        fprintf('\nslope %d, Pe_max = %d, Da_max = %d\n', m_array(i)*0.02, rS.Pe_max, rS.Da_max);
        Rp_max_Pe_m(i) = rS.Pe_max;
        Rp_max_Da_m(i) = rS.Da_max;
    end
    for i = 1:numel(Ls_array)
        Target_dir = sprintf("SlipLengths=%d", Ls_array(i))
        rS = load_and_plot_Rp(basePath, Target_dir, options);
        fprintf('\nLs %d, Pe_max = %d, Da_max = %d\n', Ls_array(i)*0.02, rS.Pe_max, rS.Da_max);
        Rp_max_Pe_Ls(i) = rS.Pe_max;
        Rp_max_Da_Ls(i) = rS.Da_max;
    end

    fig = figure('position', [1 1, 1460, 840]);
    hold on;
    subplot(1,2,1);
    m_sp = m_array .* 0.02;
    hold on;
    plot(m_sp, Rp_max_Pe_m, 'color', 'b', 'linestyle', ':', "marker",'o','LineWidth', 1.5);
    plot(m_sp, Rp_max_Da_m, 'color', 'r','linestyle', ':', "marker",'o','LineWidth', 1.5);
    hold off;
    xlabel("Slope ($n_x$/$n_y$)", 'Fontsize', FONTSIZE, 'interpreter', 'latex');
    ylabel("Pe/Da",'Fontsize', FONTSIZE, 'interpreter', 'latex');
    legend('Pe', 'Da', 'Fontsize', FONTSIZE, 'interpreter', 'latex');
    grid on;

    subplot(1,2,2);
    Ls_sp = Ls_array *xScalar;
    hold on;
    plot(Ls_sp, Rp_max_Pe_Ls, 'color', 'b','linestyle', ':', "marker",'o','LineWidth', 1.5)
    plot(Ls_sp, Rp_max_Da_Ls,'color', 'r', 'linestyle', ':', 'marker', 'o', 'LineWidth', 1.5)
    hold off;
    xlabel("Slip Length (nm)", 'Fontsize', FONTSIZE, 'interpreter', 'latex');
    ylabel("Pe/Da",'Fontsize', FONTSIZE, 'interpreter', 'latex');
    legend('Pe', 'Da', 'Fontsize', FONTSIZE, 'interpreter', 'latex');
    grid on;
    waitforbuttonpress;
end


parSweep = false;
if parSweep
    saveStruct = struct();
    plot_types = {"none"};
    Rp_array = zeros(numel(Pe_range), numel(Da_range));
    for i = 1:numel(Pe_range)
        for j = 1:numel(Da_range)
            Pe = Pe_range(i);
            Da = Da_range(j);
            Pe_string = strrep(num2str(Pe, '%.5g'), '.', '_');
            Da_string = strrep(num2str(Da, '%.5g'), '.', '_');
            Target_test = sprintf('conc_Pe%s_Da%s.mat', Pe_string, Da_string);
            StrIn = {basePath, {Target_dir}, Target_test};
            [these_stats, Rp, RS] = plot_specific(StrIn, plot_types);
            Rp_array(i,j) = Rp;
        end
    end
    saveStruct.Rp_array = Rp_array;
    fld = fullfile(basePath, Target_dir);
    save(fullfile(fld, "ANumeric_Rp.mat"),'saveStruct');
    load_struct = load(fullfile(fld, "ANumeric_Rp.mat"));
    disp(load_struct)
    RpLoaded = load_struct.saveStruct.Rp_array;

    figure;
    colormap('jet');
    imagesc(Da_range, Pe_range, Rp_array);
    set(gca,'YDir','normal');   % so Pe increases up
    xlabel('Da');
    ylabel('Pe');
    title('Reaction production R_p');
    waitforbuttonpress;
end


bigplotBool = true;
if bigplotBool
    
    specifier = "SlipLengths=";  

    isMatch = cellfun(@(x) ~isempty(strfind(x, specifier)), Folders);
    Folders_filtered = Folders(isMatch);

    stat_struct= struct();

    for i = 1:length(Folders_filtered)
        Folder = {Folders_filtered{i}};
        keyname = genvarname(Folder{1})
        plot_types = {"None"};

        StrIn = {basePath, Folder, Target_test};
        [these_stats, Rp]= plot_specific(StrIn, plot_types);
        stat_struct.(keyname) = these_stats;
    end




    % ==============================================================
    nF          = numel(Folders_filtered);
    Ux_MAPE       = zeros(1,nF);   Ux_MAR  = zeros(1,nF);
    S_MAPE        = zeros(1,nF);   S_MAR   = zeros(1,nF);
    dp_MAPE       = zeros(1,nF);   dp_MAR  = zeros(1,nF);
    a_MAPE        = zeros(1,nF);   a_MAR   = zeros(1,nF);

    for k = 1:nF
        Fkey          = genvarname(Folders_filtered{k});
        st            = stat_struct.(Fkey);

        Ux_MAPE(k)      = st.MAPE_ux;     Ux_MAR(k)  = st.MAR_ux;
        S_MAPE(k)       = st.MAPE_S;      S_MAR(k)   = st.MAR_S;
        dp_MAPE(k)      = st.MAPE_dp;   dp_MAR(k)  = st.MAR_dp;
        a_MAPE(k)       = st.MAPE_a;      a_MAR(k)   = st.MAR_a;
    end

    folderLabels = strrep(Folders, '_', '\_');   % escape underscores for axes

    % --- Choose what to plot ---
    param_to_plot = "a";   % options: "ux", "s_x", "dpdx", "a"

    switch lower(param_to_plot)
        case {"ux", "u_x", "velocity"}
            mape_array = Ux_MAPE;  mar_array  = Ux_MAR;  label = "mid-U_x";
        case {"s", "s_x"}
            mape_array = S_MAPE;   mar_array = S_MAR;   label = "S(x)";
        case {"dpdx", "dp_dx", "pressure"}
            mape_array = dp_MAPE;  mar_array = dp_MAR;  label = "dp/dx";
        case {"a", "concentration", "a_x"}
            mape_array = a_MAPE;   mar_array = a_MAR;   label = "a(x)";
        otherwise
            error("Unknown parameter name: %s", param_to_plot);
    end

    % --- Plot the selected parameter's metrics ---
    figure('position', [1, 1, 1460, 840]);

    subplot(1,2,1);
    mape_sp = linspace(0, numel(mape_array), numel(mape_array))
    hold on;
    plot(Ls_array*xScalar, mape_array, 'color','r','linestyle', ':', 'LineWidth', 1.5, 'marker', 'o');
    for ax = findobj(gcf,'Type','Axes')
        set(ax,'fontsize', FONTSIZE-6);
    end
    hold off;
    title('Max Abs $\%$ Error', 'Fontsize', FONTSIZE, "interpreter", "latex"); 
    xlabel("$L_s$ (nm)", 'Fontsize', FONTSIZE, 'interpreter', 'latex')
    ylabel('MAPE ($\%$)','Fontsize',FONTSIZE , 'margin', 1, 'interpreter', 'latex');
    grid on;

    subplot(1,2,2);
    hold on;
    plot(Ls_array*xScalar, mar_array,'color','r','linestyle', ':', 'LineWidth', 1.5, 'marker', 'o');
    for ax = findobj(gcf,'Type','Axes')
        set(ax,'fontsize', FONTSIZE-6);
    end
    hold off;

    title('Max Abs Rel Error ','Fontsize' ,FONTSIZE, "interpreter", "latex")
    xlabel("$L_s$ (nm)", 'Fontsize', FONTSIZE, 'interpreter', 'latex')
    ylabel("MARE ($Rel_{err}$)",'Fontsize',FONTSIZE , 'margin', 1, 'interpreter', 'latex')
    grid on;
    waitforbuttonpress;
%{

    pan = @(r,c) subplot(4,2, (r-1)*2 + c);      % tiny helper

    % ---- row 1 : mid-channel Ux ----------------------------------
    pan(1,1);  plot(Ux_MAPE ,'-o');
    title('Mean Absolute Percent Error: mid-U_x');   grid on;

    pan(1,2);  plot(Ux_MAR,'-o');
    title('Max Absolute Relative Error: mid-U_x');   grid on;

    % ---- row 2 : average S(x) ------------------------------------
    pan(2,1);  plot(S_MAPE ,'-o');
    title('Mean Absolute Percent Error:  S(x)');     grid on;

    pan(2,2);  plot(S_MAR,'-o');
    title('Max Absolute Relative Error: S(x)');     grid on;

    % ---- row 3 : pressure gradient -------------------------------
    pan(3,1);  plot(dp_MAPE ,'-o');
    title('Mean Absolute Percent Error:  dp/dx');    grid on;

    pan(3,2);  plot(dp_MAR,'-o');
    title('Max Absolute Relative Error: dp/dx');    grid on;

    % ---- row 4 : concentration profile ---------------------------
    pan(4,1);  plot(a_MAPE ,'-o'); 
    title('Mean Absolute Percent Error:  a(x)');     grid on; 

    pan(4,2);  plot(a_MAR,'-o');
    title('Max Absolute Relative Error: a(x)');     grid on;

    % ---- label all x-axes with folder names ----------------------
    axesList = findobj(gcf,'Type','Axes');
    for ax = axesList.'
        set(ax,'XTick',1:numel(Folders), ...
            'XTickLabel',strrep(Folders,'_','\_'), ...
            'XTickLabelRotation',45);
    end
    xlabel(axesList( end ), 'Folder');

    waitforbuttonpress;

%}
end