function Ambiguity_View

% Build the GUI
init_ui;

% This funciton sets up the figure and slider.
% note the CreateFcn and the Callback.
% The CreatFcn runs the first time while the Callback
% runs each time you move the slider.
    function init_ui()
        fig = figure;
        alpha_slider = uicontrol('Style', 'slider',...
            'Min',0.001,'Max',1,'Value',1,...
            'Position', [90 0 400 20],...
            'Callback',  @solve_and_plot,...
            'String',"Alpha",...
            "BackgroundColor",[1,1,1],...
            "SliderStep",[1/10 2/10]);
        q_enable = uicontrol('Style', 'checkbox',...
            'Position', [490 380 70 20],...
            'Callback',  @solve_and_plot,...
            'String',"Truncation",...
            "BackgroundColor",[1,1,1]);
        q_drop = uicontrol('Style', 'popupmenu',...
            'Position', [520 350 40 20],...
            'Callback',  @solve_and_plot,...
            'String',["2","4","6","8"],...
            "BackgroundColor",[1,1,1]);
        fun_drop = uicontrol('Style', "popupmenu",...
            'Position', [20 350 40 20],...
            'Callback',  @solve_and_plot,...
            'String',["rect","sinc","rrc"],...
            "BackgroundColor",[1,1,1]);

        % Store data in figure
        fig.UserData = struct("alpha",alpha_slider,"truncation",q_enable,"q",q_drop,"shape",fun_drop);
        solve_and_plot(fig);
    end

% This funciton it what is run when you first run the
% code and each time you move the slider.
    function solve_and_plot(src,~)

        % Get the current slider value
        fig = ancestor(src,"figure","toplevel");
        data = fig.UserData;
        alpha = data.alpha.Value;
        truncation = data.truncation.Value;
        q = 2*data.q.Value;
        shape = data.shape.Value;
        switch shape
            case 1
                shape = "rect";
            case 2
                shape = "sinc";
            case 3
                shape = "rrc";
        end


        % Update parameters
        T = 1;
        res_range = 100;
        res_int = 10;

        % Set time and frequency ranges
        t_range = linspace(-4*T,4*T,res_range);
        f_range = linspace(-4/T,4/T,res_range);

        % Convert to meshes
        t_range_mat = t_range.' .* ones(res_range,length(f_range));
        f_range_mat = f_range .* ones(res_range,length(f_range));

        % Render meshgrid
        if shape ~= "rect"
            if truncation
                result_3d = zeros(res_range);
                for m = 1:length(t_range)
                    for n = 1:length(f_range)
                        result_3d(m,n) = ambig_direct(t_range(m),f_range(n),T,shape,alpha,q,res_int);
                    end
                end
            else
                result_3d = ambig(t_range_mat,f_range_mat,T,shape,alpha);
            end
        else
            result_3d = ambig(t_range_mat,f_range_mat,T,shape,alpha);
        end

        % Update plot
        mesh(t_range,f_range,abs(result_3d));
        title(sprintf("Amplitude of ambiguity function using %s pulses",shape))
        ylabel("t (in T)")
        xlabel("f (in 1/T)")
        zlabel(sprintf("|A_{%s}(t,f)|",shape))
    end
end

function result = ambig(t,f,Ts,shape,alpha)
% This function results the result of the ambiguity function of two
% of the same shaped pulse filters. Supported shapes are "rect", "sinc" and
% "rrc". (alpha is only used for "rrc" and can usually be left out)
%
% Instructions:
% 1. t:         enter an m-by-n matrix of time values (in seconds)
% 2. f:         enter an m-by-n matrix of frequency values (in Hertz)
% 3. Ts:        enter symbol period (in seconds)
% 4. shape:     enter shape of pulse filters ("rect"/"sinc"/"rrc")
% 5. alpha:     enter value of roll-off factor
%                   (unused if "rrc" not selected)
%
% Output: an m-by-n sized matrix of ambiguity values, with each value
%   calculated element-wise from the given t and f matrices
%
% Coded by Jeremiah Rhys Wimer, 4/23/2024

% Tolerance for exceptions
tol = .001 * Ts;

% Premake result, reset t and f, set tolerance
result = t * 0;

% Make result for each shape
switch shape
    case "rect"
        % Define f=0 condition
        f_zero = abs(f) < tol;

        % Typical condition definitions
        cond1n = not(f_zero) & (t>0 & t<=Ts);
        cond2 = not(f_zero) & (t<=0 & t>=-Ts);

        % Limit condition definitions
        cond3n = f_zero & (t>0 & t<=Ts);
        cond4n = f_zero & (t<=0 & t>=-Ts);

        % Positive and negative definitions
        result(cond1n) = (exp(1j.*2.*pi.*f(cond1n).*Ts) - exp(1j.*2.*pi.*f(cond1n).*t(cond1n))) ./ (1j.*2.*pi.*f(cond1n).*Ts);
        result(cond2) = (exp(1j.*2.*pi.*f(cond2).*(Ts+t(cond2))) - 1) ./ (1j.*2.*pi.*f(cond2).*Ts);

        % Positive and negative limits as f goes to 0
        result(cond3n) = 1 - t(cond3n) / Ts;
        result(cond4n) = 1 + t(cond4n) / Ts;

    case "sinc"
        % Define f=0 condition
        t_zero = abs(t) < tol;

        % Typical condition definitions
        cond1n = not(t_zero) & (f>0 & f<=1/Ts);
        cond2 = not(t_zero) & (f<=0 & f>=-1/Ts);

        % Limit condition definitions
        cond3n = t_zero & (f>0 & f<=1/Ts);
        cond4n = t_zero & (f<=0 & f>=-1/Ts);

        % Positive and negative definitions
        result(cond1n) = (Ts ./ (pi.*t(cond1n))) .* sin(((pi.*t(cond1n)) ./ Ts) - pi.*f(cond1n).*t(cond1n)) .* exp(1j.*pi.*f(cond1n).*t(cond1n));
        result(cond2) = (Ts ./ (pi.*t(cond2))) .* sin(((pi.*t(cond2)) ./ Ts) + pi.*f(cond2).*t(cond2)) .* exp(1j.*pi.*f(cond2).*t(cond2));

        % Positive and negative limits as t goes to 0
        result(cond3n) = (1 - f(cond3n).*Ts);
        result(cond4n) = (1 + f(cond4n).*Ts);

    case "rrc"
        % Swap f for definition
        f = -f;
        t = t / Ts;
        f = f * Ts;
        Ts = 1;

        % RRC definitions
        W0 = 1/(2*Ts);
        W = (1+alpha)*W0;

        % Internal definitions for RRC function
        x1 = 2*W0-W;
        x2 = W;
        a = (pi/4) / (W-W0);
        phi1 = pi*(W-2*W0)/(4*(W-W0));
        phi2 = pi*(-f+W-2*W0)/(4*(W-W0));
        phi3 = pi*(f+W-2*W0)/(4*(W-W0));
        r = zeros(4,1);
        if alpha <= 0.5
            r(1) = x2-x1;
            r(2) = 2*x1;
        else
            r(1) = 2*x1;
            r(2) = x2-x1;
        end
        r(3) = x2+x1;
        r(4) = 2*x2;

        % Typical condition definitions
        cond1n = -r(1) <= f & f < 0;
        cond21n = -r(2) <= f & f < -r(1) & alpha >=0.5;
        cond22n = -r(2) <= f & f < -r(1) & alpha <0.5;
        cond3n = -r(3) <= f & f < -r(2);
        cond4n = -r(4) <= f & f <= -r(3);
        cond1p = 0 <= f & f < r(1);
        cond21p = r(1) <= f & f < r(2) & alpha >=0.5;
        cond22p = r(1) <= f & f < r(2) & alpha <0.5;
        cond3p = r(2) <= f & f < r(3);
        cond4p = r(3) <= f & f <= r(4);

        % Negative range definitions
        result(cond1n) = h2(t(cond1n),-x1+f(cond1n),-x2,-a,-a,phi1,phi3(cond1n))...
            + h1(t(cond1n),-x1,-x1+f(cond1n),-a,phi1)...
            + h0(t(cond1n),x1+f(cond1n),-x1)...
            + h1(t(cond1n),x1,x1+f(cond1n),a,phi2(cond1n))...
            + h2(t(cond1n),x2+f(cond1n),x1,a,a,phi1,phi2(cond1n));
        result(cond21n) = h2(t(cond21n),-x1+f(cond21n),-x2,-a,-a,phi1,phi3(cond21n))...
            + h1(t(cond21n),x1+f(cond21n),-x1+f(cond21n),-a,phi1)...
            + h2(t(cond21n),-x1,x1+f(cond21n),-a,a,phi1,phi2(cond21n))...
            + h1(t(cond21n),x1,-x1,a,phi2(cond21n))...
            + h2(t(cond21n),x2+f(cond21n),x1,a,a,phi1,phi2(cond21n));
        result(cond22n) = h1(t(cond22n),-x1,-x2,-a,phi1)...
            + h0(t(cond22n),x1+f(cond22n),-x1)...
            + h1(t(cond22n),x2+f(cond22n),x1+f(cond22n),a,phi2(cond22n));
        result(cond3n) = h1(t(cond3n),x1+f(cond3n),-x2,-a,phi1)...
            + h2(t(cond3n),-x1,x1+f(cond3n),-a,a,phi1,phi2(cond3n))...
            + h1(t(cond3n),x2+f(cond3n),-x1,a,phi2(cond3n));
        result(cond4n) = h2(t(cond4n),x2+f(cond4n),-x2,-a,a,phi1,phi2(cond4n));

        % Positive range definitions
        result(cond1p) = h2(t(cond1p),-x1,-x2+f(cond1p),-a,a,phi1,-phi3(cond1p))...
            + h1(t(cond1p),-x1+f(cond1p),-x1,a,-phi3(cond1p))...
            + h0(t(cond1p),x1,-x1+f(cond1p))...
            + h1(t(cond1p),x1+f(cond1p),x1,a,phi1)...
            + h2(t(cond1p),x2,x1+f(cond1p),a,-a,phi1,-phi2(cond1p));
        result(cond21p) = h2(t(cond21p),-x1,-x2+f(cond21p),-a,a,phi1,-phi3(cond21p))...
            + h1(t(cond21p),x1,-x1,a,-phi3(cond21p))...
            + h2(t(cond21p),-x1+f(cond21p),x1,a,a,phi1,-phi3(cond21p))...
            + h1(t(cond21p),x1+f(cond21p),-x1+f(cond21p),a,phi1)...
            + h2(t(cond21p),x2,x1+f(cond21p),a,-a,phi1,-phi2(cond21p));
        result(cond22p) = h1(t(cond22p),-x1+f(cond22p),-x2+f(cond22p),a,-phi3(cond22p))...
            + h0(t(cond22p),x1,-x1+f(cond22p))...
            + h1(t(cond22p),x2,x1,a,phi1);
        result(cond3p) = h1(t(cond3p),x1,-x2+f(cond3p),a,-phi3(cond3p))...
            + h2(t(cond3p),-x1+f(cond3p),x1,a,a,phi1,-phi3(cond3p))...
            + h1(t(cond3p),x2,-x1+f(cond3p),a,phi1);
        result(cond4p) = h2(t(cond4p),x2,-x2+f(cond4p),a,a,phi1,-phi3(cond4p));

    otherwise
        error("Unsupported shape selected! Please choose: 'rect' / 'sinc' / 'rc'")
end

end

function result = h0(t,f_pos,f_neg)
% Integrates exp(-1j*2*pi*f*t), from f1 to f2

% Set tolerance
tol = 1e-6;

% Define conditions
cond1 = abs(t) > tol;
cond2 = not(cond1);

% Define expression for usual cases
result_pos_vec = (1 ./ (-1j.*2.*pi.*t)) .* exp(-1j.*2.*pi.*f_pos.*t);
result_neg_vec = (1 ./ (-1j.*2.*pi.*t)) .* exp(-1j.*2.*pi.*f_neg.*t);
result_vec = result_pos_vec - result_neg_vec;

% Set results based on conditions
result(cond1) = result_vec(cond1);

% Define expression for abnormal cases
result_vec = f_pos - f_neg;
result(cond2) = result_vec(cond2);

end

function result = h1(t,f_pos,f_neg,w,phi)
% Integrates exp(-1j*2*pi*f*t)*cos(w*f+phi), from f1 to f2

% Set tolerance
tol = 1e-6;

% Define conditions
cond1 = abs(w.^2 - (2.*pi.*t).^2) > tol;
cond2 = not(cond1);

% Define expression for usual cases
result_pos_vec = exp(-1j.*2.*pi.*f_pos.*t) ./ (w.^2 - (2.*pi.*t).^2) .* (w.*sin(w.*f_pos+phi) - 1j.*2.*pi.*t.*cos(w.*f_pos+phi));
result_neg_vec = exp(-1j.*2.*pi.*f_neg.*t) ./ (w.^2 - (2.*pi.*t).^2) .* (w.*sin(w.*f_neg+phi) - 1j.*2.*pi.*t.*cos(w.*f_neg+phi));
result_vec = result_pos_vec - result_neg_vec;

% Set results based on usual conditions
result(cond1) = result_vec(cond1);

% Define expression for abnormal cases
result_neg_vec = (4.*1i.*pi.^2.*f_pos.*(2.*1i.*pi.*f_pos.*1j.*cos(f_pos.*w+phi).*t + 1i.*f_pos.*w.*sin(f_pos.*w+phi) - 2.*1j.*cos(f_pos.*w+phi)).*exp(-2.*1i.*pi.*f_pos.*t)) ./ (-8.*pi.^2);
result_pos_vec = (4.*1i.*pi.^2.*f_neg.*(2.*1i.*pi.*f_neg.*1j.*cos(f_neg.*w+phi).*t + 1i.*f_neg.*w.*sin(f_neg.*w+phi) - 2.*1j.*cos(f_neg.*w+phi)).*exp(-2.*1i.*pi.*f_neg.*t)) ./ (-8.*pi.^2);
result_vec = result_pos_vec - result_neg_vec;

% Set results based on abnormal conditions
result(cond2) = result_vec(cond2);

end

function result = h2(t,f_pos,f_neg,w1,w2,phi1,phi2)
% Integrates exp(-1j*2*pi*f*t)*cos(w1*f+phi1)*cos(w2*f+phi2), from f1 to f2

% Use other h1 function for this (cos(a)cos(b) splits into 2 parts)
part1 = h1(t,f_pos,f_neg,w1+w2,phi1+phi2);
part2 = h1(t,f_pos,f_neg,w1-w2,phi1-phi2);
result = 0.5 * (part1 + part2);

end

function result = ambig_direct(t,f,Ts,shape,alpha,q,res)
% This function returns the result of the ambiguity function of two
% of the same shaped pulse filters, truncated to [0 q*Ts]. Supported shapes
% are "rect", "sinc" and "rrc". (alpha is only used for "rrc")
%
% Instructions:
% 1. t:         enter a time value (in seconds)
% 2. f:         enter a frequency value (in Hertz)
% 3. Ts:        enter period (in seconds)
% 4. shape:     enter shape of pulse filters ("rect"/"sinc"/"rrc")
% 5. alpha:     enter value of roll-off factor
%                   (unused if "rrc" not selected)
% 6. q:         enter number of sample periods before time-domain cutoff
% 7. res:       enter the resolution of each integration, where dt=Ts/res
%
% Output: an ambiguity value, calculated from the given t and f matrices
% Note: unlike the closed form function, this can only find values one at a
% time
%
% Coded by Jeremiah Rhys Wimer, 6/4/2024

% Define resolution of "integration"
dt = Ts / res;
t_range = 0:dt:q*Ts;

% Define TX/RX filters
switch shape
    case "rect"
        filter_1 = gen_rect_filter(t_range-t,Ts,q/2);
        filter_2 = gen_rect_filter(t_range,Ts,q/2);
    case "sinc"
        filter_1 = sinc_trunc((t_range-t - q*Ts/2),Ts,q/2);
        filter_2 = sinc_trunc((t_range - q*Ts/2),Ts,q/2);
    case "rrc"
        filter_1 = RRCt_trunc((t_range-t - q*Ts/2),alpha,Ts,q/2);
        filter_2 = RRCt_trunc((t_range - q*Ts/2),alpha,Ts,q/2);
end

% Defind integration function
fun_vec = filter_1 .* filter_2 .* exp(1j.*2.*pi.*t_range.*f);
norm_val = sum((filter_2 .* filter_2).*dt,"all");
result = sum(fun_vec.*dt,"all") / norm_val;

end

function result = gen_rect_filter(t,Ts,~)
% Generates a rectangular pulse of amplitude sqrt(1/Ts) from 0 to Ts
% (the third input can be taken off, I'm just too lazy to change all my
% code)

tol = 0.0001 * Ts;

cond1 = t < tol | (t - Ts) >= tol;
cond2 = not(cond1);

result(cond1) = 0;
result(cond2) = sqrt(1/Ts);

end

function result = sinc_trunc(t,Ts,q)
% Generates a truncated sinc pulse from -qTs to qTs of amplitude sqrt(1/Ts)

tol = 0.0001 * Ts;

cond1 = abs(t) - q*Ts > tol;
cond2 = not(cond1);

result = 0 * t;
result(cond1) = 0;
result(cond2) = sqrt(1/Ts) * sinc(t(cond2)./Ts);

end

function result = RRCt_trunc(t,alpha,Ts,q)
% Generates a VERSION OF the RRC pulse (there are other definitions but
% this is valid)

tol = 0.0001 * Ts;

result = 0 * t;
R = 0.5/Ts;

cond1 = abs(t) < tol;
cond2 = abs(abs(t) - (64 * alpha^2 * R^2)^(-1/2)) < tol;
cond3 = (abs(t) - q*Ts) > tol;
cond4 = not(cond1 | cond2 | cond3);

result(cond1) = sqrt(2.*R) .* (1 - alpha + 4.*alpha./pi); % check this is correct
result(cond2) = (sqrt(2.*R) ./ (2.*pi.*R.*(1 - 192.*(alpha.^2).*(R.^2).*(t(cond2).^2)))).*(2.*pi.*R.*(1-alpha).*cos(2.*pi.*R.*(1-alpha).*t(cond2)) + 8.*R.*alpha.*(cos(2.*pi.*R.*(1+alpha).*t(cond2)) - 2.*pi.*R.*(1+alpha).*t(cond2).*sin(2.*pi.*R.*(1+alpha).*t(cond2))));
result(cond3) = 0;
result(cond4) = (sqrt(2.*R) ./ (1 - 64.*alpha.^2.*R.^2.*t(cond4).^2)) .* (sin(2.*pi.*R.*(1-alpha).*t(cond4))./(2.*pi.*R.*t(cond4)) + 4.*alpha./pi .* cos(2.*pi.*R.*(1+alpha).*t(cond4)));

end
