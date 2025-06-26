classdef obj_results
    %% How to use
    % This file is updated with the relevent data from a simulation setup,
    % then the needed statistics update automatically for retreival. Save
    % this updated object to be pulled from later, in order to render
    % figures whenever.
    %
    % 7/25/2024, JRW

    %% Properties ---------------------------------------------------------
    properties
        % System properties
        max_frames;
        max_trials;
        M_ary;
        syms_per_frame;

        % Significance level
        alpha = 0.05;

        % Trial information
        bit_errors;
        sym_errors
        frame_errors;
    end

    properties (Dependent)
        % Trial flags
        l2_norm;
        var_norm;
        ci_norm;

        % Overall data
        BER_total;
        SER_total;
        FER_total;
        Thr_total;

        % Average data
        BER;
        SER;
        FER;
        Thr;
    end

    methods
        %% DEPENDENT VARIABLES --------------------------------------------

        function result = get.l2_norm(obj)
            if obj.max_trials == 1
                result = NaN;
            else
                BER_mean = mean(obj.BER_total(1:obj.max_trials));
                norm_BER = (obj.BER_total(1:obj.max_trials) - BER_mean) / BER_mean;
                result = norm(norm_BER);
            end
        end

        function result = get.var_norm(obj)
            if obj.max_trials == 1
                result = NaN;
            else
                BER_mean = mean(obj.BER_total(1:obj.max_trials));
                BER_var = var(obj.BER_total(1:obj.max_trials));
                result = BER_var / BER_mean;
            end
        end

        function result = get.ci_norm(obj)
            SEM = std(obj.BER_total(1:obj.max_trials))/sqrt(length(obj.BER_total(1:obj.max_trials)));
            ts = tinv([(obj.alpha/2)  (1-obj.alpha/2)],length(obj.BER_total(1:obj.max_trials))-1);
            CI = mean(obj.BER_total(1:obj.max_trials)) + ts*SEM;
            result = (CI(2) - CI(1)) / mean(obj.BER_total(1:obj.max_trials));
        end

        function result = get.BER_total(obj)
            result = sum(obj.bit_errors(1:obj.max_frames,1:obj.max_trials),1) ./ (obj.syms_per_frame .* size(obj.sym_errors(1:obj.max_frames,1:obj.max_trials),1) .* log2(obj.M_ary));
        end

        function result = get.SER_total(obj)
            result = sum(obj.sym_errors(1:obj.max_frames,1:obj.max_trials),1) ./ (obj.syms_per_frame .* size(obj.sym_errors(1:obj.max_frames,1:obj.max_trials),1));
        end

        function result = get.FER_total(obj)
            result = sum(obj.frame_errors(1:obj.max_frames,1:obj.max_trials),1) ./ size(obj.sym_errors(1:obj.max_frames,1:obj.max_trials),1);
        end

        function result = get.Thr_total(obj)
            result = log2(obj.M_ary) .* (1-obj.FER_total(1:obj.max_trials));
        end

        function result = get.BER(obj)
            result = mean(obj.BER_total);
        end

        function result = get.SER(obj)
            result = mean(obj.SER_total);
        end

        function result = get.FER(obj)
            result = mean(obj.FER_total);
        end

        function result = get.Thr(obj)
            result = mean(obj.Thr_total);
        end

        %% INTERNAL FUNCTIONS ---------------------------------------------

        
    end
end