function lu  = compute_Lu(ubar, u, dubar, du, ddubar, ddu, Y, f, nu)
	[Ns Nr] = size(f);
	N  = size(Y,2);
	lu = zeros(size(f));
	
	% terms without Y_i's
	lu = repmat( -ubar.*dubar + nu*ddubar, 1, Nr );	

	% terms with Y_i's
	for i=1:N
		tmp = -du(:,i).*ubar - u(:,i).*dubar + nu*ddu(:,i);
		lu = lu + tmp*Y(:,i)';
	end

% 	terms with Y_i Y_j's 
	for i=1:N
		for j=1:N
			lu = lu - (u(:,i).*du(:,j)) * (Y(:,i).*Y(:,j))';
		end
	end
	
	lu = lu+f;
end

