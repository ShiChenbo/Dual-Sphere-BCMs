function indexMatrix = indexMatrix(degree)
% This function is to calculate the mapping between 'alpha' index and 'tau sigma l m index' index for SVWFs
% In this file, l index -> id; m index -> io; tau index -> it; sigma index -> is
Nsize = 2 * degree * (degree + 2);
n_index = zeros(Nsize, 4);
for it = 1:2
    for id = 1:degree
        for is = 0:1
            for io = 0:id
                if io == 0 && is == 1
                    continue;
                end
                n = 2 * (id * (id + 1) + io * (-1)^is + -1) + it;
                n_index(n, 1) = id;
                n_index(n, 2) = io;
                n_index(n, 3) = is;
                n_index(n, 4) = it;
            end
        end
    end
end
indexMatrix=n_index';
end

