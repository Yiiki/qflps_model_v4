function idx=find_vec(vec,val)
    % it returns the closest loc in vec
    idx=find(abs(vec-val)==min(abs(vec-val)));
end