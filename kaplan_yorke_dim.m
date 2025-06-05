% --- Function to calculate Kaplan-Yorke Dimension ---
function KY_dim = kaplan_yorke_dim(LE_spectrum_sorted)
    % Calculates the Kaplan-Yorke dimension from a sorted Lyapunov spectrum.
    % Assumes LE_spectrum_sorted is already sorted in descending order.
    LE_spectrum_s = sort(LE_spectrum_sorted, 'descend'); % Ensure sorted
    
    KY_dim = 0;
    sum_LE = 0;
    j = 0;
    for i = 1:length(LE_spectrum_s)
        if sum_LE + LE_spectrum_s(i) >= 0
            sum_LE = sum_LE + LE_spectrum_s(i);
            j = j + 1;
        else
            % Last exponent made sum negative
            if abs(LE_spectrum_s(i)) < eps % Avoid division by zero if LE is tiny
                 KY_dim = j; % If next LE is zero, dimension is integer j
                 return;
            end
            KY_dim = j + sum_LE / abs(LE_spectrum_s(i));
            return;
        end
    end
    % If sum of all exponents is positive (e.g. expanding system or only positive LEs considered)
    KY_dim = j; 
end