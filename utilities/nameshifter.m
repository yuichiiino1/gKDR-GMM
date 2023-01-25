function shiftednames = nameshifter(names)
   
m = length(names);
shiftednames = cell(m,1);
shifts = [1 9; 1 4; 0 0];
for shifti = 1:3
    shift = shifts(shifti,:);
    for row = shifti:3:m
        name = names{row};
        shiftedname = [name repmat(' ',1,shift(1)) repmat('-',1, shift(2))];
        shiftednames{row} = shiftedname;
    end
end

end