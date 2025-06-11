function Y = ode_RKn_deci_bounded(odefun,tspan,y0, RK_method_num, dispMethodOrder, deci, min_max_range, varargin)
  %   Solve differential equations with fixed step RK methods from Wikipedia's explicit "List of Runge–Kutta methods"
  %   Initially based on a historic version of ode5 from Mathworks 
  %   I've added min_max_range enforces the min and max range of the output state 
  %   min_max_range is nan for dimensions with no range limits.
  %   set min_max_range = []; if no states are bounded.
  %   I've added more ButcherTables for orders 1, 2, and 4 in a subfunction below

  %   Example 
  %         tspan = 0:0.1:20;
  %         y = ode5(@vdp1,tspan,[2 0]);  
  %         plot(tspan,y(:,1));
  %     solves the system y' = vdp1(t,y) with a constant step size of 0.1, 
  %     and plots the first component of the solution.   

  persistent warning_shown;

  if ~isnumeric(tspan)
    error('TSPAN should be a vector of integration steps.');
  end

  if ~isnumeric(y0)
    error('Y0 should be a vector of initial conditions.');
  end

  % Validate deci parameter early
  if ~isnumeric(deci) || ~isscalar(deci) || deci < 1 || floor(deci) ~= deci
      error('DECI must be a positive integer scalar.');
  end

  if ~iscolumn(y0)
    error('Y0 must be a column vector.');
  end

  n_states = length(y0); % number of dynamical state variables

  % Validate RK_method_num parameter
  if ~isnumeric(RK_method_num) || ~isscalar(RK_method_num) || floor(RK_method_num) ~= RK_method_num || ~ismember(RK_method_num, [1, 2, 3, 4, 5, 6])
    error('RK_method_num must be an integer and one of the following values: 1, 2, 3, 4, 5, 6');
  end

  h = diff(tspan);
  if any(sign(h(1))*h <= 0)
    error('Entries of TSPAN are not in order.') 
  end  

  % checks for min_max_range
  if ~isnumeric(min_max_range)
      error('MIN_MAX_RANGE must be numeric.');
  end

  if ~isempty(min_max_range)
      % Check dimensions
      if size(min_max_range, 2) ~= 2
          error('MIN_MAX_RANGE must have exactly 2 columns [min_vals, max_vals].');
      end
      
      if size(min_max_range, 1) ~= n_states
          if isempty(warning_shown) || ~warning_shown
              warning('ODE_RKN_DECI_BOUNDED:DimensionMismatch', 'MIN_MAX_RANGE row count (%d) does not match state vector length (%d). Bounds will not be applied.', size(min_max_range, 1), n_states);
              warning_shown = true;
          end
          min_logical = false(n_states,1);
          max_logical = false(n_states,1);
          min_max_range = NaN(n_states,2); % Set to non-bounding state to prevent errors in later logic
      else
          % This is the case where bounds can be applied.
          min_logical = ~isnan(min_max_range(:,1));
          max_logical = ~isnan(min_max_range(:,2));
      
          % Check that min <= max where both are not NaN
          both_finite_and_specified = min_logical & max_logical;
      
          if any(both_finite_and_specified & (min_max_range(:,1) > min_max_range(:,2)))
              error('MIN_MAX_RANGE: minimum values must be less than or equal to maximum values for rows where both are specified.');
          end
      
          % Check that initial conditions are within bounds
          y0_col = y0(:);  % Ensure column vector for consistent indexing
      
          % Check minimum bounds
          if any(min_logical & (y0_col < min_max_range(:,1)))
              violating_indices = find(min_logical & (y0_col < min_max_range(:,1)));
              error('Initial condition Y0(%d) = %g violates minimum bound %g.', ...
                    violating_indices(1), y0_col(violating_indices(1)), min_max_range(violating_indices(1),1));
          end
      
          % Check maximum bounds
          if any(max_logical & (y0_col > min_max_range(:,2)))
              violating_indices = find(max_logical & (y0_col > min_max_range(:,2)));
              error('Initial condition Y0(%d) = %g violates maximum bound %g.', ...
                    violating_indices(1), y0_col(violating_indices(1)), min_max_range(violating_indices(1),2));
          end
      end
      
  else
      min_max_range = NaN(n_states,2); % Ensure min_max_range has a defined structure
      min_logical = false(n_states,1);
      max_logical = false(n_states,1);
  end

  f0 = feval(odefun,tspan(1),y0,varargin{:});

  y0 = y0(:);   % Make a column vector.
  if ~isequal(size(y0),size(f0))
    error('Inconsistent sizes of Y0 and f(t0,y0).');
  end  

  nt= length(tspan);
  Y = zeros(n_states,nt,'double');

  % Method coefficients -- Butcher's tableau
  %  
  %   C | A
  %   --+---
  %     | B


  %% RK5
  % C = double([1/5; 3/10; 4/5; 8/9; 1]);
  % A = double([ 1/5,          0,           0,            0,         0
  %               3/40,         9/40,        0,            0,         0
  %               44/45        -56/15,       32/9,         0,         0
  %               19372/6561,  -25360/2187,  64448/6561,  -212/729,   0
  %               9017/3168,   -355/33,      46732/5247,   49/176,   -5103/18656]);
  % B = double([35/384, 0, 500/1113, 125/192, -2187/6784, 11/84]);

  [A, B, C, order, MethodName] = getRKnButcherTableau(RK_method_num);

  if dispMethodOrder % for assurance or debug
    display([MethodName ', order = ' num2str(order)])
  end

  % Determine the number of stages from the A matrix (Butcher A)
  nstages = size(A, 1);

  % C is the c vector (nodes) from Butcher tableau.
  % K_stages will store the k_i values for each stage
  K_stages = zeros(n_states,nstages,'double');

  yi = y0;

  Y(:,1) = y0;
  
  for i = 2:nt
      
    ti = tspan(i-1);
    hi = h(i-1);

    % Standard explicit Runge-Kutta framework
    for s_idx = 1:nstages % Iterate through each stage s_idx = 1, ..., nstages
        
        % Calculate y_intermediate = yi + hi * sum_{j=1}^{s_idx-1} A(s_idx, j) * K_stages(:,j)
        y_sum_prev_K = zeros(n_states,1,'double');
        if s_idx > 1
            % Vectorized sum: K_stages(:, 1:s_idx-1) is n_states x (s_idx-1)
            % A(s_idx, 1:s_idx-1).' is (s_idx-1) x 1 (using Butcher A matrix)
            y_sum_prev_K = K_stages(:, 1:s_idx-1) * A(s_idx, 1:s_idx-1).';
        end
        
        y_arg_for_f = yi + hi * y_sum_prev_K;
        t_arg_for_f = ti + C(s_idx) * hi; % C(s_idx) is c_s from Butcher tableau
        
        K_stages(:,s_idx) = feval(odefun, t_arg_for_f, y_arg_for_f, varargin{:});
    end
    
    % Update yi: yi_new = yi + hi * sum_{j=1}^{nstages} B(j) * K_stages(:,j)
    % K_stages is n_states x nstages
    % B is 1 x nstages (row vector)
    % K_stages * B.' results in an n_states x 1 vector
    yi = yi + hi * (K_stages * B.');
    
    % Apply bounds if specified
    % min_logical and max_logical correctly indicate which states have bounds.
    % min_max_range is guaranteed to be n_states-by-2 (possibly with NaNs).
    if any(min_logical)
      actually_violating_min = min_logical & (yi < min_max_range(:,1));
      if any(actually_violating_min)
          yi(actually_violating_min) = min_max_range(actually_violating_min, 1);
      end
    end
    if any(max_logical)
      actually_violating_max = max_logical & (yi > min_max_range(:,2));
      if any(actually_violating_max)
          yi(actually_violating_max) = min_max_range(actually_violating_max, 2);
      end
    end

    Y(:,i) = yi;
  end

  Y = Y(:,1:deci:end); % decimate at the end. This uses more memory but it is simple, and simple is preferred.

  Y = Y.'; % make it nt_decimated x n_states

end

function [A, b, c, order, MethodName] = getRKnButcherTableau(methodNumber)
  % getRKNButcherTableau retrieves Butcher tableau coefficients for Runge-Kutta methods.
  %
  %   Inputs:
  %       methodNumber - An integer identifying the desired Runge-Kutta method.
  %
  %   Outputs:
  %       A          - The A matrix of the Butcher tableau.
  %       b          - The b vector (weights) of the Butcher tableau.
  %       c          - The c vector (nodes) of the Butcher tableau.
  %       order      - The order of accuracy of the method.
  %       MethodName - The name of the Runge-Kutta method.
  %
  %   Available methods:
  %       1: Forward Euler (order 1)
  %       2: Explicit Midpoint Method (order 2)
  %       3: Ralston's 2nd-Order Method (order 2)
  %       4: Heun's Method (order 2)
  %       5: Classic RK4 (order 4)
  %       6: Ralston's 4th-Order Method (order 4)
  %       X: Dormand-Prince 5th Order (6 stages) (order 5) %  not good, commented out. some details possibly wrong.

  % Input validation
  if ~isnumeric(methodNumber) || ~isscalar(methodNumber) || floor(methodNumber) ~= methodNumber || methodNumber < 1 || methodNumber > 6
      error('methodNumber must be an integer between 1 and 6.');
  end

  BT = struct('MethodName', {}, 'A', {}, 'b', {}, 'c', {}, 'order', {});

  % 1. Forward Euler (order 1)
  BT(1).MethodName = 'Forward Euler';
  BT(1).A          = 0;
  BT(1).b          = 1;
  BT(1).c          = 0;
  BT(1).order      = 1;

  % 2. Explicit Midpoint Method (order 2)
  BT(2).MethodName = 'Explicit Midpoint';
  BT(2).A          = [0   0;
                      0.5 0];
  BT(2).b          = [0   1];
  BT(2).c          = [0; 0.5];
  BT(2).order      = 2;

  % 3. Ralston's 2nd-Order Method (order 2)
  BT(3).MethodName = 'Ralston 2nd Order';
  BT(3).A          = [0     0;
                      2/3   0];
  BT(3).b          = [1/4   3/4];
  BT(3).c          = [0; 2/3];
  BT(3).order      = 2;

  % 4. Heun's Method (order 2)
  % Note: There are multiple methods named Heun's method.
  % This typically refers to the improved Euler method or trapezoidal rule for ODEs.
  BT(4).MethodName = 'Heun''s Method'; % Using two apostrophes for a single quote in a string
  BT(4).A          = [0   0;
                      1   0];
  BT(4).b          = [0.5   0.5];
  BT(4).c          = [0; 1];
  BT(4).order      = 2;

  % 5. Classic RK4 (order 4)
  BT(5).MethodName = 'Classic RK4';
  BT(5).A          = [0   0   0   0;
                      0.5 0   0   0;
                      0   0.5 0   0;
                      0   0   1   0];
  BT(5).b          = [1/6 1/3 1/3 1/6];
  BT(5).c          = [0; 0.5; 0.5; 1];
  BT(5).order      = 4;

  % 6. Ralston's 4th-Order Method (order 4)
  % This method aims for high accuracy with minimal function evaluations for a 4th order method.
  % It has a smaller error coefficient compared to the classic RK4 for some problems.
  BT(6).MethodName = 'Ralston 4th Order';
  BT(6).A          = [0,                               0,                                   0,                                    0;
                      2/5,                             0,                                   0,                                    0;
                      (-2889 + 1428*sqrt(5))/1024,     (3785 - 1620*sqrt(5))/1024,          0,                                    0;
                      (-3365 + 2094*sqrt(5))/6040,     (-975 - 3046*sqrt(5))/2552,          (467040 + 203968*sqrt(5))/240845,     0];
  BT(6).b          = [(263 + 24*sqrt(5))/1812, (125 - 1000*sqrt(5))/3828, (3426304 + 1661952*sqrt(5))/5924787, (30 - 4*sqrt(5))/123];
  BT(6).c          = [0; 2/5; (14 - 3*sqrt(5))/16; 1];
  BT(6).order      = 4;

  % % 7. Dormand-Prince 5th Order (6 stages)
  % % Coefficients for a 6-stage, 5th order Dormand-Prince method.
  % % c_nodes are c_2 to c_s (s=6)
  % % A_coeffs are A_{i+1, j} for i=1..s-1, j=1..i
  % BT(7).MethodName = 'Dormand-Prince 5th Order';
  % BT(7).c          = [0;             % c_1
  %                     1/5;           % c_2
  %                     3/10;          % c_3
  %                     4/5;           % c_4
  %                     8/9;           % c_5
  %                     1];            % c_6
  % BT(7).A          = [0,            0,            0,            0,            0,            0;
  %                     1/5,          0,            0,            0,            0,            0;
  %                     3/40,         9/40,         0,            0,            0,            0;
  %                     44/45,       -56/15,       32/9,         0,            0,            0;
  %                     19372/6561,  -25360/2187,  64448/6561,  -212/729,     0,            0;
  %                     9017/3168,   -355/33,      46732/5247,   49/176,      -5103/18656,  0];
  % BT(7).b          = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
  % BT(7).order      = 5;

  % Assign outputs
  A = BT(methodNumber).A;
  b = BT(methodNumber).b;
  c = BT(methodNumber).c;
  order = BT(methodNumber).order;
  MethodName = BT(methodNumber).MethodName;

end