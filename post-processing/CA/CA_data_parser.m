clear all
close all

tic
enable_plotting = false;

for replica=0:2

    filename_dump = sprintf('../../lammps/CA_test/%d-replica/wettability_%d.dump', replica,replica);
    filename_txt = sprintf('../../lammps/CA_test/%d-replica/wettability_%d.txt', replica,replica);
    movefile(filename_dump, filename_txt)
    
    data_tot = readtable(filename_txt,'ReadVariableNames',false);
    movefile(filename_txt, filename_dump)

    data_tot(:,"Var15") = [];
    data_tot(:,"Var16") = [];
    data_tot = data_tot(:, {'Var1', 'Var2', 'Var3','Var14', 'Var5', 'Var6', 'Var7', 'Var8', 'Var9', 'Var10', 'Var11', 'Var12', 'Var13'});
    data_tot(:,"Var8") = [];
    data_tot(:,"Var9") = [];
    data_tot(:,"Var10") = [];
    
    if ~exist("indexes.mat","file")
    	kk = 0;
    	index = zeros(500,1);
    	num_ref = data_tot{2,1};
    	for ii = 2:(size(data_tot,1)-1)
        	    if data_tot{ii,1} == num_ref && isnan(data_tot{ii+1,1}) && isnan(data_tot{ii-1,1})
                	kk = kk+1; 
                	index(kk) = ii;
        	    end
    	end
    	index(kk+1:end) = [];
    
    	save indexes index
    end
    load indexes.mat
    
    tot = struct('data',[],'CA_tan',[],'CA_min',[],'CA_max',[],'CA_dev',[],'x_solv',[],'y_solv',[],'z_solv',[],'x_pol',[],'y_pol',[],'z_pol',[]);
    tot(length(index),1) = tot;
    
    for ii = 1:length(index)
        if ii == length(index)
            tot(ii).data = data_tot{(index(ii)+6):end,:};
        else
            tot(ii).data = data_tot{(index(ii)+6):(index(ii+1)-4),:};
        end
    end
    
    for rr = 1:length(index)
        
        clearvars -except rr data_tot index tot replica enable_plotting
    
        qq = linspace(0.07,0.1,20)';
        contact_angle_global = zeros(size(qq,1),3); 
        
        for pp = 1:size(qq,1)
        
            clearvars -except contact_angle_global pp qq rr data_tot index tot replica enable_plotting
        
            data = tot(rr).data;
            num_atoms = size(data,1);
            lim_solv = 3;
            jj = 1;
            kk = 1;
            data_solv = zeros(num_atoms,10);
            data_pol = zeros(num_atoms,10);
    
            for ii = 1:num_atoms
                if data(ii,3)<lim_solv
                    data_solv(jj,:) = data(ii,:);
                    jj = jj+1;
                else
                    data_pol(kk,:) = data(ii,:);
                    kk = kk+1;
                end
            end
            
            data_solv(jj:end,:) = [];
            data_pol(kk:end,:) = [];
            
            data_solv = sortrows(data_solv,7);
            data_pol = sortrows(data_pol,7);
            
            num_at_solv = size(data_solv,1);
            num_at_pol = size(data_pol,1);
            
            coord_mean_1 = data_pol(end,7);
            coord_mean_2 = data_solv(1,7);
            ii = 1;
            flag = 1;
            
            while abs(coord_mean_1-coord_mean_2)>0.01
                ii = ii+1;
                if ii>size(data_pol,1) || ii>size(data_solv,1)
                    % flag = 0;
                    coord_mean_1 = data_solv(1,7);
                    coord_mean_2 = data_solv(1,7);
                else
                    coord_mean_1 = mean(data_pol((num_at_pol-ii+1):end,7));
                    coord_mean_2 = mean(data_solv(1:ii,7));
                end
            end
            
            if flag == 0
                break
            end
            
            coord_mean = mean([coord_mean_1,coord_mean_2]);
            
            center_x = mean(data_solv(:,5));
            center_y = mean(data_solv(:,6));
            
            temp = zeros(num_at_solv,13);
            temp(:,1:10) = data_solv;
            for ii = 1:num_at_solv
                temp(ii,11) = temp(ii,5)-center_x;
                temp(ii,12) = temp(ii,6)-center_y;
                temp(ii,13) = sqrt((temp(ii,11))^2+(temp(ii,12))^2);
            end
            
            rad_temp = mean(temp(:,13))*qq(pp)*2;
            
            cil = zeros(num_at_solv,13);
            kk = 1;
            for ii = 1:num_at_solv
                if temp(ii,13)<(rad_temp) && temp(ii,7)>max(data_pol(:,7))
                    cil(kk,:) = temp(ii,:);
                    kk = kk+1;
                end
            end
            cil(kk:end,:) = [];
            cil = sortrows(cil,7);
            
            delta_3 = rad_temp/10;
            interface_3 = zeros(size(cil,1),6);
            ii = 1;
            jj = 1;
            first = max(data_pol(:,7));
            last = first+delta_3;
            
            interface_3(jj,3) = first;
            while ii <= size(cil,1)
                if cil(ii,7)>=first && cil(ii,7)<last
                    interface_3(jj,1) = interface_3(jj,1)+cil(ii,7);
                    interface_3(jj,2) = interface_3(jj,2)+1;
                    ii = ii+1;
                else
                    if interface_3(jj,2) > 0
                        interface_3(jj,4) = last;
                        jj = jj+1;
                        first = last;
                        last = last+delta_3;
                        interface_3(jj,3) = first;
                    else
                        first = last;
                        last = last+delta_3;
                    end
                end
            end
            interface_3(jj,4) = last;
            
            interface_3((jj+1):end,:) = [];
            interface_3(:,1) = interface_3(:,1)./(interface_3(:,2));
            interface_3(:,5) = pi*(rad_temp)^2*(interface_3(:,4)-interface_3(:,3));
            interface_3(:,6) = interface_3(:,2)./interface_3(:,5);
            
            rounded = round(size(interface_3,1)/10);
            if rounded == 0
                rounded = 1;
            end
            filtered_3 = movmean(interface_3(:,6),rounded);
            
            num_el = round(size(filtered_3,1)/2);
            ref_value_3 = mean(interface_3(1:num_el,6));
            
            temp_2 = zeros(size(interface_3,1),6);
            kk = 1;
            for ii = 1:size(interface_3,1)
                if interface_3(ii,6)>(0.2*ref_value_3)
                    temp_2(kk,:) = interface_3(ii,:);
                    kk = kk+1;
                end
            end
            temp_2(kk:end,:) = [];
            
            h_d = temp_2(end,1)-coord_mean;
            
            delta = h_d*qq(pp);
    
            lim_1 = coord_mean+delta;
            lim_2 = lim_1+delta;
            
            plane_1 = zeros(num_at_solv,13);
            jj = 1;
            
            for ii = 1:num_at_solv
                if data_solv(ii,7)>lim_1 && data_solv(ii,7)<lim_2
                    plane_1(jj,1:10) = data_solv(ii,:);
                    jj = jj+1;
                end
            end
            
            plane_1(jj:end,:) = [];
            
            for ii = 1:size(plane_1,1)
                plane_1(ii,11) = plane_1(ii,5)-center_x;
                plane_1(ii,12) = plane_1(ii,6)-center_y;
                plane_1(ii,13) = sqrt((plane_1(ii,11))^2+(plane_1(ii,12))^2);
            end
            
            delta_2 = delta/5;
            plane_1 = sortrows(plane_1,13);
            z_1 = mean(plane_1(:,7));
            interface_1 = zeros(size(plane_1,1),6);
            ii = 1;
            jj = 1;
            first = 0;
            last = first+delta_2;
            
            interface_1(jj,3) = first;
            while ii <= size(plane_1,1)
                if plane_1(ii,13)>=first && plane_1(ii,13)<last
                    interface_1(jj,1) = interface_1(jj,1)+plane_1(ii,13);
                    interface_1(jj,2) = interface_1(jj,2)+1;
                    ii = ii+1;
                else
                    if interface_1(jj,2) > 0
                        interface_1(jj,4) = last;
                        jj = jj+1;
                        first = last;
                        last = last+delta_2;
                        interface_1(jj,3) = first;
                    else
                        first = last;
                        last = last+delta_2;
                    end
                end
            end
            interface_1(jj,4) = last;
            
            interface_1((jj+1):end,:) = [];
            interface_1(:,1) = interface_1(:,1)./(interface_1(:,2));
            interface_1(:,5) = pi*(interface_1(:,4).^2-interface_1(:,3).^2).*delta;
            interface_1(:,6) = interface_1(:,2)./interface_1(:,5);
            
            rounded = round(size(interface_1,1)/10);
            if rounded == 0
                rounded = 1;
            end
            filtered_1 = movmean(interface_1(:,6),rounded);
            
            num_el = round(size(filtered_1,1)/2);
            ref_value_1 = mean(interface_1(1:num_el,6));
            
            ref_value_1_1 = 0.5*ref_value_1;
            ii = size(filtered_1,1)-1;
            radius_1_1 = 0;
            while ii > 0
                if filtered_1(ii)>ref_value_1_1 && filtered_1(ii+1)<ref_value_1_1
                    radius_1_1 = ((ref_value_1_1-filtered_1(ii))/(filtered_1(ii+1)-filtered_1(ii)))*(interface_1(ii+1,1)-interface_1(ii,1))+interface_1(ii,1);
                    ii = 0;
                elseif filtered_1(ii) == ref_value_1_1
                    radius_1_1 = interface_1(ii,1);
                    ii = 0;
                else
                    ii = ii-1;
                end
            end
            if radius_1_1 == 0
                radius_1 = max(interface_1(:,1));
            else
                ref_value_1_2 = 0.2*ref_value_1;
                ii = size(filtered_1,1)-1;
                radius_1_2 = 0;
                while ii > 0
                    if filtered_1(ii)>ref_value_1_2 && filtered_1(ii+1)<ref_value_1_2
                        radius_1_2 = ((ref_value_1_2-filtered_1(ii))/(filtered_1(ii+1)-filtered_1(ii)))*(interface_1(ii+1,1)-interface_1(ii,1))+interface_1(ii,1);
                        ii = 0;
                    elseif filtered_1(ii) == ref_value_1_2
                        radius_1_2 = interface_1(ii,1);
                        ii = 0;
                    else
                        ii = ii-1;
                    end
                end
                radius_1 = 0;
                jj = 0;
                if radius_1_2 == 0
                    radius_1_2 = max(interface_1(:,1));
                end
                for ii = 1:size(plane_1,1)
                    if plane_1(ii,13)>=radius_1_1 && plane_1(ii,13)<=radius_1_2
                        radius_1 = radius_1+plane_1(ii,13);
                        jj = jj+1;
                    end
                end
                radius_1 = radius_1/jj;
            end
            
            lim_3 = lim_2;
            lim_4 = lim_3+delta;
            
            plane_2 = zeros(num_at_solv,13);
            jj = 1;
            
            for ii = 1:num_at_solv
                if data_solv(ii,7)>lim_3 && data_solv(ii,7)<lim_4
                    plane_2(jj,1:10) = data_solv(ii,:);
                    jj = jj+1;
                end
            end
            
            plane_2(jj:end,:) = [];
            
            for ii = 1:size(plane_2,1)
                plane_2(ii,11) = plane_2(ii,5)-center_x;
                plane_2(ii,12) = plane_2(ii,6)-center_y;
                plane_2(ii,13) = sqrt((plane_2(ii,11))^2+(plane_2(ii,12))^2);
            end
            
            plane_2 = sortrows(plane_2,13);
            z_2 = mean(plane_2(:,7));
            interface_2 = zeros(size(plane_2,1),6);
            ii = 1;
            jj = 1;
            first = 0;
            last = first+delta_2;
            
            interface_2(jj,3) = first;
            while ii <= size(plane_2,1)
                if plane_2(ii,13)>=first && plane_2(ii,13)<last
                    interface_2(jj,1) = interface_2(jj,1)+plane_2(ii,13);
                    interface_2(jj,2) = interface_2(jj,2)+1;
                    ii = ii+1;
                else
                    if interface_2(jj,2) > 0
                        interface_2(jj,4) = last;
                        jj = jj+1;
                        first = last;
                        last = last+delta_2;
                        interface_2(jj,3) = first;
                    else
                        first = last;
                        last = last+delta_2;
                    end
                end
            end
            interface_2(jj,4) = last;
            
            interface_2((jj+1):end,:) = [];
            interface_2(:,1) = interface_2(:,1)./(interface_2(:,2));
            interface_2(:,5) = pi*(interface_2(:,4).^2-interface_2(:,3).^2).*delta;
            interface_2(:,6) = interface_2(:,2)./interface_2(:,5);
            
            rounded = round(size(interface_2,1)/10);
            if rounded == 0
                rounded = 1;
            end
            filtered_2 = movmean(interface_2(:,6),rounded);
    
            num_el = round(size(filtered_2,1)/2);
            ref_value_2 = mean(interface_2(1:num_el,6));
            
            ref_value_2_1 = 0.5*ref_value_2;
            ii = size(filtered_2,1)-1;
            radius_2_1 = 0;
            while ii > 0
                if filtered_2(ii)>ref_value_2_1 && filtered_2(ii+1)<ref_value_2_1
                    radius_2_1 = ((ref_value_2_1-filtered_2(ii))/(filtered_2(ii+1)-filtered_2(ii)))*(interface_2(ii+1,1)-interface_2(ii,1))+interface_2(ii,1);
                    ii = 0;
                elseif filtered_2(ii) == ref_value_2_1
                    radius_2_1 = interface_2(ii,1);
                    ii = 0;
                else
                    ii = ii-1;
                end
            end
            if radius_2_1 == 0
                radius_2 = max(interface_2(:,1));
            else
                ref_value_2_2 = 0.2*ref_value_2;
                ii = size(filtered_2,1)-1;
                radius_2_2 = 0;
                while ii > 0
                    if filtered_2(ii)>ref_value_2_2 && filtered_2(ii+1)<ref_value_2_2
                        radius_2_2 = ((ref_value_2_2-filtered_2(ii))/(filtered_2(ii+1)-filtered_2(ii)))*(interface_2(ii+1,1)-interface_2(ii,1))+interface_2(ii,1);
                        ii = 0;
                    elseif filtered_2(ii) == ref_value_2_2
                        radius_2_2 = interface_2(ii,1);
                        ii = 0;
                    else
                        ii = ii-1;
                    end
                end
                radius_2 = 0;
                jj = 0;
                if radius_2_2 == 0
                    radius_2_2 = max(interface_2(:,1));
                end
                for ii = 1:size(plane_2,1)
                    if plane_2(ii,13)>=radius_2_1 && plane_2(ii,13)<=radius_2_2
                        radius_2 = radius_2+plane_2(ii,13);
                        jj = jj+1;
                    end
                end
                radius_2 = radius_2/jj;
            end
            
            dz = z_2-z_1;
            dr = radius_2-radius_1;
            
            if dr<0
                CA_tan = rad2deg(atan(abs(dz/dr))); 
            else
                CA_tan = 180-rad2deg(atan(dz/dr));
            end
            
            contact_angle_global(pp,1) = CA_tan;
            contact_angle_global(pp,2) = qq(pp);
            contact_angle_global(pp,3) = delta;
            
        end
        
        if flag == 1 && size(contact_angle_global,1) ~= 1
            CA_tan = contact_angle_global(:,1);
            CA_tan = sort(CA_tan);
            CA_tan([1,2,length(CA_tan)-1,length(CA_tan)]) = [];
            CA_tan_fin = mean(CA_tan);
        elseif flag == 1 && size(contact_angle_global,1) == 1 
            CA_tan_fin = contact_angle_global(1,1);
        elseif flag == 0 
            CA_tan_fin = 90;
        end
    
        if flag == 1
            tot(rr).CA_min = min(CA_tan);
            tot(rr).CA_max = max(CA_tan);
            tot(rr).CA_dev = std(CA_tan);
        else
            tot(rr).CA_min = 90;
            tot(rr).CA_max = 90;
            tot(rr).CA_dev = 0;
        end
        
        tot(rr).CA_tan = CA_tan_fin;
        tot(rr).x_solv = data_solv(:,5);
        tot(rr).y_solv = data_solv(:,6);
        tot(rr).z_solv = data_solv(:,7);
        tot(rr).x_pol = data_pol(:,5);
        tot(rr).y_pol = data_pol(:,6);
        tot(rr).z_pol = data_pol(:,7);
        
        if enable_plotting
            if flag == 0
                rr
                figure (7*(rr-1)+1)
                plot3(data_pol(:,5),data_pol(:,6),data_pol(:,7),'b.','MarkerSize',10)
                hold all
                plot3(data_solv(:,5),data_solv(:,6),data_solv(:,7),'r.','MarkerSize',10)
            else
                if (rr == 10 || rr == 20 || rr == 40)
                    rr
                    figure (7*(rr-1)+1)
                    semilogy(interface_1(:,1),interface_1(:,6))
                    grid on
                    hold all
                    semilogy(interface_1(:,1),filtered_1)
                    hold all
                    semilogy(interface_1(:,1),ref_value_1*ones(size(interface_1,1),1))
                    figure (7*(rr-1)+2)
                    plot(interface_1(:,1),interface_1(:,6))
                    grid on
                    hold all
                    plot(interface_1(:,1),filtered_1)
                    hold all
                    plot(interface_1(:,1),ref_value_1*ones(size(interface_1,1),1))
                    hold all
                    plot(interface_1(:,1),ref_value_1_1*ones(size(interface_1,1),1))
                    hold all
                    plot(interface_1(:,1),ref_value_1_2*ones(size(interface_1,1),1))
                    hold all
                    if radius_1_1 ~= 0
                        plot([radius_1_1,radius_1_1],[0.01,0.06])
                        hold all
                        plot([radius_1_2,radius_1_2],[0.01,0.06])
                    end
                    
                    figure (7*(rr-1)+3)
                    semilogy(interface_2(:,1),interface_2(:,6))
                    grid on
                    hold all
                    semilogy(interface_2(:,1),filtered_2)
                    hold all
                    semilogy(interface_2(:,1),ref_value_2*ones(size(interface_2,1),1))
                    figure (7*(rr-1)+4)
                    plot(interface_2(:,1),interface_2(:,6))
                    grid on
                    hold all
                    plot(interface_2(:,1),filtered_2)
                    hold all
                    plot(interface_2(:,1),ref_value_2*ones(size(interface_2,1),1))
                    hold all
                    plot(interface_2(:,1),ref_value_2_1*ones(size(interface_2,1),1))
                    hold all
                    plot(interface_2(:,1),ref_value_2_2*ones(size(interface_2,1),1))
                    hold all
                    if radius_2_1 ~= 0
                        plot([radius_2_1,radius_2_1],[0.01,0.06])
                        hold all
                        plot([radius_2_2,radius_2_2],[0.01,0.06])
                    end
                    
                    figure (7*(rr-1)+5)
                    semilogy(interface_3(:,1),interface_3(:,6))
                    grid on
                    hold all
                    semilogy(interface_3(:,1),filtered_3)
                    hold all
                    semilogy(interface_3(:,1),ref_value_3*ones(size(interface_3,1),1))
                    figure (7*(rr-1)+6)
                    plot(interface_3(:,1),interface_3(:,6))
                    grid on
                    hold all
                    plot(interface_3(:,1),filtered_3)
                    hold all
                    plot(interface_3(:,1),ref_value_3*ones(size(interface_3,1),1))
                    
                    x_min = min(data(:,5));
                    x_max = max(data(:,5));
                    y_min = min(data(:,6));
                    y_max = max(data(:,6));
                    z_min = min(data(:,7));
                    z_max = max(data(:,7));
                    
                    figure (7*(rr-1)+7)
                    plot3(data_pol(:,5),data_pol(:,6),data_pol(:,7),'b.','MarkerSize',10)
                    hold all
                    plot3(data_solv(:,5),data_solv(:,6),data_solv(:,7),'r.','MarkerSize',10)
                    hold all
                    plot3(ones(2,1)*center_x,ones(2,1)*center_y,[coord_mean/2,z_max*1.1],'b-*')
                    hold all
                    plot3([x_min,x_max],[y_min,y_max],lim_1*ones(2,1),'m-*')
                    hold all
                    plot3([x_min,x_max],[y_max,y_min],lim_1*ones(2,1),'m-*')
                    hold all
                    plot3([x_min,x_max],[y_min,y_max],lim_2*ones(2,1),'c-*')
                    hold all
                    plot3([x_min,x_max],[y_max,y_min],lim_2*ones(2,1),'c-*')
                    hold all
                    plot3([x_min,x_max],[y_min,y_max],lim_3*ones(2,1),'c-*')
                    hold all
                    plot3([x_min,x_max],[y_max,y_min],lim_3*ones(2,1),'c-*')
                    hold all
                    plot3([x_min,x_max],[y_min,y_max],lim_4*ones(2,1),'y-*')
                    hold all
                    plot3([x_min,x_max],[y_max,y_min],lim_4*ones(2,1),'y-*')
                    hold all
                    plot3([x_min,x_max],[y_min,y_max],(temp_2(end,1))*ones(2,1),'b-*')
                    hold all
                    plot3([x_min,x_max],[y_max,y_min],(temp_2(end,1))*ones(2,1),'b-*')
                    hold all
                    plot3([x_min,x_max],[y_min,y_max],(coord_mean)*ones(2,1),'g-*')
                    hold all
                    plot3([x_min,x_max],[y_max,y_min],(coord_mean)*ones(2,1),'g-*')
                    hold all
                    
                    circle_x = linspace((-radius_1),(radius_1),100);
                    circle_y_1 = sqrt(radius_1^2-circle_x.^2);
                    circle_y_2 = -sqrt(radius_1^2-circle_x.^2);
                    circle_x = circle_x+center_x;
                    circle_y_1 = circle_y_1+center_y;
                    circle_y_2 = circle_y_2+center_y;
                    
                    plot3(circle_x,circle_y_1,(lim_1+lim_2)*0.5*ones(100,1),'b-')
                    hold all
                    plot3(circle_x,circle_y_2,(lim_1+lim_2)*0.5*ones(100,1),'b-')
                    hold all
                    
                    circle_x = linspace((-radius_2),(radius_2),100);
                    circle_y_1 = sqrt(radius_2^2-circle_x.^2);
                    circle_y_2 = -sqrt(radius_2^2-circle_x.^2);
                    circle_x = circle_x+center_x;
                    circle_y_1 = circle_y_1+center_y;
                    circle_y_2 = circle_y_2+center_y;
                    
                    plot3(circle_x,circle_y_1,(lim_3+lim_4)*0.5*ones(100,1),'b-')
                    hold all
                    plot3(circle_x,circle_y_2,(lim_3+lim_4)*0.5*ones(100,1),'b-')
                    hold all
                
                    axis equal
                    
                end
            end
        end
        
    end
    
    savefile = sprintf("wet_%d.mat", replica);
    save(savefile,'tot')
end

fprintf("Completed in %.1f min\n",toc/60)
