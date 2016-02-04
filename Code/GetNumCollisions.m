function numCol = GetNumCollisions(inputX, vols)

numdebris = inputX(2:2:end);
numsp = inputX(1:2:end);

dens = (numdebris(:)) ./ vols(:);

c1 = 4.7842 / (36500);
c2 = 1.7114e+04 / (36500)*1.35;

numCol(:,1) = dens .* numdebris(:) * c1;
numCol(:,2) = dens .* numsp(:) * c2;

%%%%%%%%%%%%%

end