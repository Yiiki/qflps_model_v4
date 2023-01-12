function y=vec2cell(x,args_pts)

if length(args_pts)~=7
    error('invalid args_pts length')
end

if sum(args_pts)~=length(x)
    error('invalid x met')
end

for i=1:7
    %     y.ue=x(1:args_pts(1));
    %     y.vthn=x(args_pts(1)+1:sum(args_pts(1:2)));
    %     y.uh=x(sum(args_pts(1:2))+1:sum(args_pts(1:3)));
    %     y.vthp=x(sum(args_pts(1:3))+1:sum(args_pts(1:4)));
    %     y.nid=x(sum(args_pts(1:4))+1:sum(args_pts(1:5)));
    %     y.phi=x(sum(args_pts(1:5))+1:sum(args_pts(1:6)));
    %     y.phi0=x(sum(args_pts(1:6))+1:sum(args_pts(1:7)));
    y{1}=x(1:args_pts(1));
    y{2}=x(args_pts(1)+1:sum(args_pts(1:2)));
    y{3}=x(sum(args_pts(1:2))+1:sum(args_pts(1:3)));
    y{4}=x(sum(args_pts(1:3))+1:sum(args_pts(1:4)));
    y{5}=x(sum(args_pts(1:4))+1:sum(args_pts(1:5)));
    y{6}=x(sum(args_pts(1:5))+1:sum(args_pts(1:6)));
    y{7}=x(sum(args_pts(1:6))+1:sum(args_pts(1:7)));
end

end