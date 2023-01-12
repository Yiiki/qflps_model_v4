function y=path_intn_par(pa,vgs_c,vds_c,nid,AFE)
    y=vgs_c;
    for i=1:length(vgs_c)
        y(i)=-sign(vds_c(i)).*path_intn(pa,vgs_c(i)+max(0,-vds_c(i)),abs(vds_c(i)),nid(i))./AFE;
    end
%     y=-sign(vds_c).*path_intn_par(pa,vgs_c+max(0,-vds_c),abs(vds_c),nid)./AFE;
end