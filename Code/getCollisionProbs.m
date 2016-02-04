function collisionProbs = getCollisionProbs(inputX, vols)

%basealt = 200;
maxalt = 1600;
earthR = 6371;

a1 = 10^(-10);
a2 = 2*10^(-6);

ca11 = (sqrt(a1)+sqrt(a1))^2;

ca21 = (sqrt(a1)+sqrt(a2))^2;


numdebris = inputX(1:2:end);
numsp = inputX(1:2:end);

%%%%%%%spatial density%%%%%%%

%diff = (maxalt - basealt) / length(numdebris);
%vols = ((altitudes(2:end)+earthR).^3 - (altitudes(1:end-1)+earthR).^3);

beta = (numdebris(:)) ./ vols(:);

collisionProbs(:,1) = 7 * beta * ca11 .* (numdebris(:).^2);
collisionProbs(:,2) = 7 * beta * ca21 .* numdebris(:) .* numsp(:);

collisionProbs = collisionProbs*3600*24;

%%%%%%%%%%%%%

end