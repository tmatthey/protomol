function result = WHAM(T,Phi,PhiE,P,V,referenceTemp)
    % T - temperature
    % Phi - dihedrals
    % PhiE - dihedral energies
    % P - Potential Energy
    % V - biasing potentials (each potential should be a column in a
    % matrix)
    
    % The number of simulations
    R = length(T);
    
    % Divide the dihedral values into bins
        dihedrals = [];
        for i = 1:R
            dihedrals = [dihedrals; Phi(1:end,i)];
        end
        % 1) Determine the number of bins to use (Sturgis Rule)
        numBins = 1 + round(log2(length(dihedrals)));
        numBins = 5 * numBins;
        % 2) Divide into the specified number of bins
        [binCounts binCenters] = hist(dihedrals,numBins);
        % 3) Determine the bin edges
        for i = 1:(length(binCenters) - 1)
            binEdges(i) = (binCenters(i) + binCenters(i+1)) / 2;
        end
        
    % Sum the biasing potentials at each timestep
    if length(V(1,1:end)) > 1
        V_sum = transpose(sum(transpose(V)));
    else
        V_sum = V;
    end
    
    % Compute the free energies
    k_B = 0.001987191;
    B = 1 ./ (k_B * T);
    Beta = 1 / (k_B * referenceTemp);
    %exp_f = ones(1,R);
    %oldexp_f = exp_f;
    %converged = 0;
    %while converged == 0
    %    for i = 1:R
    %        sum = 0;
    %        for k = 1:R
    %            thisPhi = Phi(1:end,k);
    %            for t = 1:length(thisPhi)
    %                numer = exp(-B(i) * P(t,k));
    %                denom = 0;
    %                for m = 1:length(T)
    %                    denom = denom + length(thisPhi) * oldexp_f(m) * exp(-B(m) * P(t,k));
    %                end
    %                sum = sum + (numer / denom);
    %            end
    %        end
    %        exp_f(i) = 1 / sum;
    %    end
        
    %    % Check for convergence
    %    delta = abs(exp_f - oldexp_f) ./ oldexp_f;
    %    if max(delta) <= 0.01
    %        converged = 1;
    %    end
    %    
    %    oldexp_f = exp_f;
    %end
    
    %f = log(exp_f);
    %for i = 1:length(binCenters)
    %    bc(i) = mean(binCenters(i,1:end));
    %end
    

    % Compute the probabilities
        numerVal = 0;
        numer = zeros(numBins,1);
        count = zeros(numBins,1);
        %xwham = zeros(numBins,1);
        for r = 1:R
            dihedral = Phi(1:end,r);
            dihedralEnergy = PhiE(1:end,r);
            potential = P(1:end,r);
            %value = X(1:end,r);
            for i = 1:length(dihedral)
                j = -1;
                for k = 1:(numBins-1)
                    if j == -1
                        if binEdges(k) > dihedral(i)
                            j = k;
                        end
                    end
                end
                if j == -1
                    j = numBins;
                end
                
                %xwham(j) = xwham(j) + exp((B(r) - Beta) * potential(i) + (B(r) * V_sum(i))) * value(i);
                numer(j) = numer(j) + exp((B(r) - Beta) * potential(i) + (B(r) * V_sum(i)));
                numerVal = numerVal + dihedralEnergy(i) * exp((B(r) - Beta) * potential(i) + (B(r) * V_sum(i)));
                count(j) = count(j) + 1;
            end
        end
        
        % Normalize the probabilities
        denom=trapezoidRule(binCenters,numer)
        denomVal=sum(numer)
        P = numer ./ denom;
        Val = numerVal ./ denomVal
        binCenters = 180 / pi * transpose(binCenters);
        
        %% Create figure
        figure1 = figure;
 
        %% Create axes
        axes1 = axes('YTick',[0 0.25 0.5 0.75 1 1.25 1.5],'Parent',figure1);
        xlabel(axes1,'Dihedral Angle','FontSize',12);
        hold(axes1,'all');
        
        plot(binCenters,P);
        res = parabolaFit(binCenters,P);
        %plot(res(1:end,1),res(1:end,2),'k');
        result = res;
        %plot(binCenters,P,'k');