 function [G, E, update] = calcgradest(block_prev, block_curr)
   [gx, gy] = gradient(block_prev);
   E = block_curr - block_prev;
   G = [gx(:) gy(:)]';
   %GtG = [sum(sum(gx .* gx)) sum(sum(gx .* gy)); sum(sum(gx .* gy)) sum(sum(gy .* gy))];
   GtG = [sum(gx(:).^2) sum(gx(:).*gy(:)); sum(gx(:).*gy(:)) sum(gy(:).^2)];

   %update = -inv(GtG) * G * E(:); 
   update = -pinv(GtG) * G * E(:);



 end
