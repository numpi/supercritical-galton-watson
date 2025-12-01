function empirical_distr_W = distr_W_from_simulation(p, T, maxtime)
    % Computes T (approximate) samples of the distribution W corresponding
    % to probability vector p, stopping the simulation at time = maxtime
    % Warning: gets really slow if m is large...
    empirical_distr_W = zeros(1, T);
    m = sum(p'.*(0:length(p)-1));
    f = waitbar(0, 'Starting');
    for i = 1:T
        empirical_distr_W(i) = simulate(p, m, maxtime);
        waitbar(i/T, f, sprintf('Progress: %d %%', floor(i/T*100)));
    end
    close(f)
end

function w = simulate(p, m, maxtime)
    % Computes one sample of the evolution of the population corresponding
    % to probability vector p, stopping time = maxtime, and normalized by
    % m^maxtime.
    w = 1;
    for t = 1:maxtime
        if (w == 0) 
            return;
        end
        w = sum(randsample(0:length(p)-1, w, true, p));
    end
    w = w / m^maxtime;
end
