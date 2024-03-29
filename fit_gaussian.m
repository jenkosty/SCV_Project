function results = fit_gaussian(qc_ts,tag_no,i,choice)

        if choice == 1
            qc_ts(tag_no).ds.spice_anom(:,i) = -qc_ts(tag_no).ds.spice_anom(:,i);
        end
    
        % Grab amplitude and depth of max spice anomaly
        spike.A  = max(qc_ts(tag_no).ds.spice_anom(:,i));
        spike.P = qc_ts(tag_no).ds.pres(find(qc_ts(tag_no).ds.spice_anom(:,i) == spike.A),i);
        
        % Get range of allowable parameters
        prng = [-0.2:0.05:0.2];  % allow pressure peak to vary between +- 20% of height
        arng  = [0.8:0.05:1.2];  % allow amplitude range of +- 20% of spice anomaly peak
        hrng  = [50:10:500];     % allow height to vary  between 50 and 850m
        
        % Set up matrix for least-squared error calculation
        lse = nan(length(prng),length(arng),length(hrng));
        
        % Go through all possible combinations
        hcnt = 0; % reset h counter
        for h = hrng
            hcnt = hcnt + 1; % increase 'h' counter
            acnt = 0;        % reset 'a' counter
            for a = arng
                acnt = acnt + 1; % increase 'a' counter
                pcnt = 0;        % reset 'p' counter
                for p = prng
                    pcnt = pcnt + 1; % increase 'p'
                    
                    % Center Gaussian model around spike.P + p*h
                    zo = [];
                    zo = double(qc_ts(tag_no).ds.pres(:,i) - [spike.P + p*(4)*sqrt(h^2/2)]);
                    sa = double(qc_ts(tag_no).ds.spice_anom(:,i));
                    
                    % Reduce to where data exists
                    dat      = sa + zo;
                    sa       = sa(~isnan(dat));
                    zo       = zo(~isnan(dat));
                    
                    % Generate gaussian model using updated amplitude, center, and height
                    gauss = (spike.A*a)*exp((-(zo.^2))/(h.^2));
                    
                    % Get gaussian limits for testing
                    pl = [spike.P + p*(4)*sqrt(h^2/2)] - 2*sqrt((h^2)/2); pl  = round(pl/10)*10;
                    ph = [spike.P + p*(4)*sqrt(h^2/2)] + 2*sqrt((h^2)/2); ph  = round(ph/10)*10;
                    
                    % Grab results
                    zp     = [zo + spike.P + p*(4)*sqrt(h^2/2)];
                    dataX  = qc_ts(tag_no).ds.spice_anom(pl <= qc_ts(tag_no).ds.pres(:,i) & qc_ts(tag_no).ds.pres(:,i) <= ph, i);
                    dataY  = qc_ts(tag_no).ds.pres(pl <= qc_ts(tag_no).ds.pres(:,i) & qc_ts(tag_no).ds.pres(:,i) <= ph, i);
                    dataY  = round(dataY, 6);
                    modelX = gauss(pl <= zp & zp <= ph);
                    modelY = zp(pl <= zp & zp <= ph);
                    modelY = round(modelY, 6);
                    
                    % Check that depths of model and data intersect
                    if length(dataX) < length(modelX) | length(modelX) < length(dataX)
                        [c,~,~] = intersect(dataY,modelY);
                        ind     = find(min(c) <= dataY & dataY <= max(c));
                        dataX   = dataX(ind);   dataY = dataY(ind);
                        ind     = find(min(c) <= modelY & modelY <= max(c));
                        modelX  = modelX(ind); modelY = modelY(ind);
                    end
                    
                    % Calculate R^2
                    R2(pcnt,acnt,hcnt) = corr2(dataX,modelX).^2;
                    
                    % Save least-squared error results (ignore if bad R2 value (i.e. < 0.5))
                    if R2(pcnt,acnt,hcnt) < 0.5
                        lse(pcnt,acnt,hcnt) = NaN;
                    else
                        lse(pcnt,acnt,hcnt) = sum([dataX-modelX].^2);
                    end
                end
            end
        end
        
        % Find best zo,A,H combo according to lse
        [minlse,idxlse] = min(lse(:));
        if isnan(minlse)
            results.rejected = 1;
            results.reason = "Bad R^2 Value (Gaussian)";
            return
        end
        [a,b,c] = ind2sub(size(lse),idxlse);
        
        % Update parameters
        results.A    = spike.A*arng(b);
        results.H    = hrng(c);
        results.P    = spike.P + prng(a)*(4)*sqrt(results.H^2/2);
        results.Plow = spike.P - 2*sqrt((results.H^2)/2);
        results.Phih = spike.P + 2*sqrt((results.H^2)/2);
        results.Plow = round(results.Plow/10)*10;
        results.Phih = round(results.Phih/10)*10;
        
        % Update zo,zp,gauss for final model
        zo    = double(qc_ts(tag_no).ds.pres(:,i) - [results.P]);
        zp    = zo + results.P;
        gauss = results.A*exp((-(zo.^2))/(results.H.^2));
        
        % Save final model
        results.X = gauss;
        results.Y = zp;
        
        % Determining whether fit passes
        if minlse < 0.05 && results.A > 0.08
            results.rejected = 0;
        elseif minlse >=0.05 && results.A <= 0.08
            results.rejected = 1;
            results.reason = "Failed Gaussian Fit (Error + Amplitude)";
        elseif minlse >= 0.05
            results.rejected = 1;
            results.reason = "Failed Gaussian Fit (Error)";
        elseif results.A <= 0.08
            results.rejected = 1;
            results.reason = "Failed Gaussian Fit (Amplitude)";
        end
        
        %// Finally, fix orientation of model if minty
        if choice == 1
            results.X = -results.X;
            qc_ts(tag_no).ds.spice_anom(:,i) = -qc_ts(tag_no).ds.spice_anom(:,i);
        end
        
%                 figure()
%         plot(qc_ts(tag_no).ds.spice_anom(:,i),qc_ts(tag_no).ds.pres(:,i),'k','linewidth',2)
%         hold on; grid on; set(gca,'YDir','Reverse')
%         plot(results.X,results.Y,'Color','r','LineWidth',3,'LineStyle','-.')
%         xlabel('m')
%         ylabel('dbar');
        
   
