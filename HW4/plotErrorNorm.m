function plotErrorNorm(x_num, x_planned, t)
    % Calc norm
    Err_x2 = (x_num - x_planned).^2;
    x_nrm = sqrt(sum(Err_x2, 2));
    full_path_len = 0;
    for i=1:length(t)-1
        dx_planned = x_planned(i+1,:) - x_planned(i,:);
        full_path_len = full_path_len + sqrt(sum(dx_planned.^2)); % root of sum squared distance between every two points
    end
    x_nrm_percent = 100 * x_nrm ./ full_path_len;
    
    figure("Name","Position Error Norm", "Color",'w', 'NumberTitle', 'off'); 
    set(gcf, "Color", 'w'); box on;
    plot(t, x_nrm_percent, "LineWidth",1.5, "DisplayName","Error Norm (\%)"); grid on; box on;
    title("Position Error Norm", "FontSize", 18);
    xlabel("$t\ (s)$","Interpreter","latex")
    ylabel("$Error\ (\%)$","Interpreter","latex")
    yline(0.15, 'Color','r', 'LineWidth',2, 'LineStyle',':', "DisplayName","0.15 \%");
    legend('show')
end

