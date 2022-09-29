function prof = dha_check(qc_ts, tag_no, i)

        opts1 = optimset('display','off','UseParallel',false);

        %%% Extract vertical velocity and horizontal structure modes of climatology
        [meop_profile(i).ref.wmodes, meop_profile(i).ref.pmodes, ~, ~] = dynmodes(qc_ts(tag_no).ps.ref_N2(~isnan(qc_ts(tag_no).ps.ref_N2(:,i)),i), qc_ts(tag_no).ps.pres(~isnan(qc_ts(tag_no).ps.ref_N2(:,i)),i),1);

        %%% Grab pressure levels of mode decomposition, add zero level (surface)
        meop_profile(i).ref.mode_pres = qc_ts(tag_no).ps.pres(~isnan(qc_ts(tag_no).ps.ref_N2(:,i)),i);
        meop_profile(i).ref.mode_pres = [0;meop_profile(i).ref.mode_pres];

        % Interpolate 1st baroclinic mode to pressure of SCV cast
        meop_profile(i).dyn_pres = qc_ts(tag_no).ps.pres(~isnan(qc_ts(tag_no).ps.dyn_height_anom(:,i)),i);
        meop_profile(i).dyn_height_anom = qc_ts(tag_no).ps.dyn_height_anom(~isnan(qc_ts(tag_no).ps.dyn_height_anom(:,i)),i);
        meop_profile(i).ref.BC1_data = meop_profile(i).ref.pmodes(:,1);
        meop_profile(i).ref.BC1_pres = meop_profile(i).ref.mode_pres;
        meop_profile(i).ref.BC1_data = interp1(meop_profile(i).ref.mode_pres(~isnan(meop_profile(i).ref.BC1_data)),meop_profile(i).ref.BC1_data(~isnan(meop_profile(i).ref.BC1_data)),meop_profile(i).dyn_pres);
        meop_profile(i).ref.BC1_pres = meop_profile(i).dyn_pres;

        % Create function that describes residuals between projected BC1 and dyn_height_anom
        % Exclude data inbetween SCV limits for better fit to first mode
        dat = [];
        dat = [meop_profile(i).ref.BC1_data + meop_profile(i).dyn_height_anom];
        x_o = [];
        x_o = meop_profile(i).ref.BC1_data(~isnan(dat));
        x_p = [];
        x_p = meop_profile(i).ref.BC1_pres(~isnan(dat));
        x_f = [];
        x_f = meop_profile(i).dyn_height_anom(~isnan(dat));

        % Get limits
        pl = qc_ts(tag_no).gauss_fit{1,i}.Plow;
        ph = qc_ts(tag_no).gauss_fit{1,i}.Phih;
        
        % Remove values between upper/lower limits of SCV to avoid bad fit
        ind = [];
        ind = find(pl < x_p & x_p < ph);
        if isempty(ind)
            meop_profile(i).rejected = 1;
            meop_profile(i).reason = "DHA Preprocessing";
            prof = meop_profile(i);
            return;
        end
        if ind(end) == length(x_p)
            ind = ind(1:end-1);
        end
        x_o(ind) = [];
        x_f(ind) = [];
        %x_p(ind) = [];

        % Remove mixed layer depths (Lynne Talley method, first density greater than 0.03 from sfc value
        ind      = [];
        mld_dens = qc_ts(tag_no).ps.sigma0(~isnan(qc_ts(tag_no).ps.sigma0(:,i)),i);
        mld_pres = qc_ts(tag_no).ps.pres(~isnan(qc_ts(tag_no).ps.sigma0(:,i)),i);
        ind      = find(mld_dens > mld_dens(1)+0.03);
        mld_pres = mld_pres(ind(1));
        ind      = find(x_p < mld_pres);
        if length(ind) >= length(x_o)
            meop_profile(i).rejected = 1;
            meop_profile(i).reason = "DHA Preprocessing";
            prof = meop_profile(i);
            return;
        end
        x_o(ind) = [];
        x_f(ind) = [];

        % f simply evaluates a given alpha (modal amplitude) and returns the
        % difference between the input DHanom profile and the projected 1st mode
        % We want to restrict our solutions such that the bottom of the projected
        % profile is equal to the bottom of the DHanom profile
        % SO let alpha2 = DHanom(end) - alpha*BT1(end)
        f = [];
        f = @(alpha) (alpha*x_o - x_f + (x_f(end) - alpha*x_o(end)));
        x0  = 0.05; % First guess

        % Solve for best modal amplitude
        alpha = [];
        alpha = lsqnonlin(f,x0,[-1],[1],opts1);

        % Redfine x_o and x_f with full profile
        x_o = meop_profile(i).ref.BC1_data(~isnan(dat));
        x_p = meop_profile(i).ref.BC1_pres(~isnan(dat)); %%% NOTE: Changed from meop_profile(i).ref.mode_pres - need to check with Danny
        x_f = meop_profile(i).dyn_height_anom(~isnan(dat));

        % Fix dynamic height anomaly by removing projected 1st mode, add back in barotopic mode
        meop_profile(i).dyn_height_anom_BC1 = [x_f] - [x_o*alpha + (x_f(end) - alpha*x_o(end))];
        meop_profile(i).dyn_height_pres_BC1 = meop_profile(i).dyn_pres(~isnan(dat));

        % Save VMD results
        meop_profile(i).ref.VMD.x_f      = x_f;
        meop_profile(i).ref.VMD.x_o      = x_o;
        meop_profile(i).ref.VMD.x_p      = x_p;
        meop_profile(i).ref.VMD.alpha    = alpha;

        % Get mode decomposition results
        BC1 = meop_profile(i).ref.VMD.x_o*meop_profile(i).ref.VMD.alpha;
        BC1 = BC1 - BC1(end); %// Set bottom to zero
        
        %%% Checking if positive maximum is achieved within anomaly height
        local_max = find(islocalmax(meop_profile(i).dyn_height_anom_BC1));
        
        ind = local_max(meop_profile(i).dyn_height_pres_BC1(local_max) > pl & meop_profile(i).dyn_height_pres_BC1(local_max) < ph);
        
        if isempty(ind)
            meop_profile(i).rejected = 1;
            meop_profile(i).reason = "Failed DHA Test (Max Not Within SCV Limits)";
        end
        
        for j = ind
            meop_profile(i).rejected = [];
            if meop_profile(i).dyn_height_anom_BC1(j) > 0
                meop_profile(i).rejected = 0;
                meop_profile(i).reason = strings;
                break
            end
        end
        
        if isempty(meop_profile(i).rejected)
            meop_profile(i).rejected = 1;
            meop_profile(i).reason = "Failed DHA Test (Negative Max)";
        end
        
        prof = meop_profile(i);
        
        %            
%         % Plot results
% %         figure();
% %         subplot(121)
% %         plot(-meop_profile(i).ref.pmodes(:,1),meop_profile(i).ref.mode_pres,'r','linewidth',2)
% %         hold on; grid on; set(gca,'YDir','Reverse')
% %         plot(meop_profile(i).ref.pmodes(:,2),meop_profile(i).ref.mode_pres,'b','linewidth',2)
% %         plot(meop_profile(i).ref.pmodes(:,3),meop_profile(i).ref.mode_pres,'g','linewidth',2)
% %         plot(meop_profile(i).ref.pmodes(:,4),meop_profile(i).ref.mode_pres,'y','linewidth',2)
% %         plot(meop_profile(i).ref.pmodes(:,5),meop_profile(i).ref.mode_pres,'color',[0.5 0 0.5],'linewidth',2)
% %         title({'\it\bf\fontsize{8}\fontname{Helvetica}Horizontal Velocity','Modes'})
% %         set(gca,'XTick',[0])
% %         ylabel('Pressure (dbar)')
% %         [l,~] = legend('Mode-1','Mode-2','Mode-3','Mode-4','Mode-5','location','southeast');
% %         l.Box = 'off';
% %         ylim([0 400]);
% %         
% %         subplot(122)
% %         plot(meop_profile(i).dyn_height_anom,meop_profile(i).dyn_pres,'k','linewidth',2)
% %         hold on; grid on; set(gca,'YDir','Reverse')
% %         plot(BC1,meop_profile(i).ref.VMD.x_p,':r','linewidth',2)
% %         plot(meop_profile(i).dyn_height_anom_BC1,meop_profile(i).dyn_height_pres_BC1,':k','linewidth',2)
% %         legend('DH''_{orig}','BC1_{fit}','DH''_{adj}','location','southeast');
% %         xlabel('m^2/s^2')
% %         title({'\it\bf\fontsize{8}\fontname{Helvetica}Dynamic Height','Anomaly'})
% %         set(gca,'YTickLabel',[])
% %         ylim([0 400]);
