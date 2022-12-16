function OW = Okubo_Weiss(dxg, dyg, dxc, dyc, u, v)
    
    dxg = double(dxg);
    dyg = double(dyg);
    dxc = double(dxc);
    dyc = double(dyc);
    u = double(u);
    v = double(v);

    dudx = NaN(size(dxg));
    dudy = NaN(size(dxg));
    dvdx = NaN(size(dxg));
    dvdy = NaN(size(dxg));

    for i = 2:size(dudx,2)-1
        dudx(i,:) = (u(i+1,:)-u(i,:)) ./ dxg(i,:); % Cell center
        dvdx(i,:) = (v(i,:)-v(i-1,:)) ./ dxc(i,:); % Cell corner
    end

    for j = 2:size(dudx,2)-1
        dudy(:,j) = (u(:,j)-u(:,j-1)) ./ dyc(:,j); % Cell corner
        dvdy(:,j) = (v(:,j+1)-v(:,j)) ./ dyg(:,j); % Cell center
    end

    %%% Calculating normal strain
    s_n_center = dudx - dvdy; % Cell center

    %%% Averaging to get normal strain on the cell corner
    s_n = NaN(size(s_n_center));
    for i = 2:size(s_n_center,1)
        for j = 2:size(s_n_center,2)
            s_n(i,j) = mean([s_n_center(i,j) s_n_center(i-1,j) s_n_center(i-1,j-1) s_n_center(i,j-1)], 'omitnan');
        end
    end
    
    %%% Calculating shear strain
    s_s = dvdx + dudy; % Cell corner
    
    %%% Calculating relative vorticity
    omega = dvdx - dudy; % Cell corner
    
    %%% Calculating Okubo Weiss
    OW = s_n.^2 + s_s.^2 - omega.^2; % Cell corner
    
end