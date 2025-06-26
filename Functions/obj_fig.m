classdef obj_fig
    %% FIG_OBJ
    %% How to use
    % This file acts as a holder for all information relating to rendering
    % the figures, using data generated from simulation files. Manually
    % update these parameters and feed the file through to the
    % "render_figure.m" function.
    %
    % 7/25/2024, JRW

    %% Properties ---------------------------------------------------------
    properties
        render_for_paper = false;
        save_figure = true;
        fig_num = 1;
        fig_title = "";
        fig_data_type = "BER";
        xlabel_vec;
        ylim_vec;
        % legend_color_vec;
        legend_name_vec;
        var_list;
        var1_spec = "";
        var2_spec = "";
        gen_spec = "";
        legend_loc = "southwest";
    end

    properties (Dependent)
        legend_color_vec
    end

    methods
        %% DEPENDENT VARIABLES --------------------------------------------

        function result = get.legend_color_vec(obj)
            % variants = combinations(obj.var1_spec,obj.var2_spec,obj.pt_spec);


            % Create grids
            [A_grid, B_grid, C_grid] = ndgrid(obj.var2_spec, obj.var1_spec, obj.gen_spec);

            % Convert the grids to linear arrays
            A_list = A_grid(:);
            B_list = B_grid(:);
            C_list = C_grid(:);

            % Concatenate the strings
            result = strcat(A_list,B_list,C_list);
        end

        %% INTERNAL FUNCTIONS ---------------------------------------------

        
    end
end