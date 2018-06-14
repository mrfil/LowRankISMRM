 function [xs, info] = qpwls_pcg(x, G, W, yi, nder1, C, M, niter, mask, chat)
%function [xs, info] = qpwls_pcg(x, G, W, yi, nder1, C, M, niter, mask, chat)
%
%	quadratic penalized weighted least squares via
%	preconditioned conjugate gradients algorithm
%	cost:	Psi(x) = (y-Gx)'W(y-Gx)/2 - n'(y-Gx) + x'C'Cx/2
%	Input
%		G	[nd,np]		system matrix
%		W	[nd,nd]		data weighting matrix
%		yi	[nd,1]		noisy data
%		nder1	[nd,1]		linear term (or zero)
%		x	[np,1] or [nx,ny]	initial estimate
%		C	[nc,np]		penalty 'derivatives' (R = \Half C'*C)
%		M	[np,np]		preconditioner (or object) (or 1)
%				use 'circ0' for center-based circulant precon!
%		niter			# total iterations
%		mask	[nx,ny]		which pixels are updated
%	Output
%		xs	[np,niter]	estimates each iteration
%		info	[niter, 4]	gamma, step, time, flops
%
%	Copyright Jan 1998	Jeff Fessler, The University of Michigan

if nargin < 8 | nargin > 10, help(mfilename), error(mfilename), end

if ~isvar('chat') | isempty(chat), chat = logical(0); end

%
%	if requested, build circulant preconditioner based on center pixel
%
if ischar(M)
	if ~streq(M, 'circ0'), error 'unknown preconditioner string', end
	if ~isvar('mask'), error 'mask required for circulant precon', end
	e = zeros(size(mask));
	e(end/2+1,end/2+1) = 1;	% impulse at center
	e = e(mask(:));
	gwg = G' * (W * (G * e));
	ccc = C' * (C * e);
	gwg = embed(gwg, mask);
	ccc = embed(ccc, mask);
	f.gwg = fft2(fftshift(gwg));
	f.ccc = fft2(fftshift(ccc));
	f.ccc = Real(f.ccc, 'warn');		% these should be nearly real
	f.gwg = Real(f.gwg, 'warn');
	if min(f.gwg(:)) < 0
		Sprintf('setting %g%% to zero', ...
			min(f.gwg(:)) / max(f.gwg(:)) * 100)
		f.gwg = max(f.gwg, 0);
	end
	f.h = f.gwg + f.ccc;	% approximate hessian in freq. domain
	if chat
		Sprintf('approximate condition number: %g', max(f.h) / min(f.h))
	end

	M = 1 ./ f.h;
	iscirc = logical(1);
else
	iscirc = logical(0);
end

np = ncol(G);
if length(x(:)) ~= np & isvar('mask')
	x = x(mask(:));
end

yi = yi(:);
xs = zeros(np, niter);
xs(:,1) = x;

if chat
	clf, im(121, embed(x,mask), 'xinit')
end

info	= zeros(niter,4);

%
%	initialize projections
%
Gx = G * x;
Cx = C * x;

%tic
%flops0 = flops;

%
%	iterate
%
for ii=2:niter

	if chat & ~rem(ii,10), Sprintf('CG iteration %d', ii), end

	%
	%	(negative) gradient
	
	grad = G' * (W * (yi-Gx) - nder1) - C' * Cx;

	%
	%	preconditioned gradient
	%
	if ~iscirc
		pregrad = M * grad;
	else
		pregrad = ifft2(M .* fft2(embed(grad,mask)));
		pregrad = pregrad(mask(:));
                if isreal(grad)
  		  pregrad = Real(pregrad, 'warn');
 		end
	end

	%
	%	direction
	%
	newinprod = grad' * pregrad;
	if ii == 2
		dir = pregrad;
		gamma = 0;
	else
	%	gamma = (newinprod - oldgrad' * pregrad) / oldinprod;
		gamma = newinprod / oldinprod;
		dir = pregrad + gamma * dir;
	end
	oldgrad = grad;
	oldinprod = newinprod;

	Gdir = G * dir;
	Cdir = C * dir;

	%	check if descent direction
	if dir' * grad < 0
		warning 'wrong direction'
		keyboard
	end

	%
	%	step size in search direction
	%
	%step = real((dir' * grad) / (Gdir'*(W*Gdir) + Cdir'*Cdir));
	step = real(dir' * grad) / real(Gdir'*(W*Gdir) + Cdir'*Cdir);
	if step < 0
		warning 'downhill?'
		keyboard
	end

	%
	%	update
	%
	Gx	= Gx  + step * Gdir;
	Cx	= Cx  + step * Cdir;
	x	= x + step * dir;
	xs(:,ii) = x;
	
	if chat & (ii < 10 | rem(ii,10)==0)
		im(121, embed(x,mask), 'x')
		title(sprintf('%d step=%g gamma=%g', ii, step, gamma))
		if isobject(M)
			ptype = M.type;
		else
			ptype = 'unknown';
		end
		im(122, embed(pregrad,mask), ['PreGrad ' ptype])
		drawnow
	end

	info(ii,1) = gamma;
	info(ii,2) = step;
	%info(ii,3) = toc;		% accum. time
%	info(ii,4) = flops;		% accum. flops
end

%	info(:,4) = info(:,4) - flops0;
