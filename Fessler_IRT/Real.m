 function y = Real(x, arg2)
%function y = Real(x)
%function y = Real(x, tol)
%function y = Real(x, 'warn')
%function y = Real(x, 'error')
%function y = Real(x, 'report')
%	return real part of complex matrix
%	checking that imaginary part is negligible

com = 'error';
tol = 1e-13;
if nargin > 1
	if ischar(arg2)
		com = arg2;
	elseif isnumeric(arg2)
		tol = arg2;
	end
end

	if strcmp(com, 'warn')
		onlywarn = 1;
	elseif strcmp(com, 'error')
		onlywarn = 0;
	elseif strcmp(com, 'report')
		;
	else
		error(sprintf('bad argument %s', com))
	end

	if max(abs(x(:))) == 0
		if any(imag(x(:)) ~= 0)
			error 'max real 0, but imaginary!'
		else
			y = real(x);
			return
		end
	end

	frac = max(abs(imag(x(:)))) / max(abs(x(:)));
	if strcmp(com, 'report')
		disp(sprintf('imaginary part %g%%', frac * 100))
		return
	end

	if frac > tol
		t = sprintf('imaginary fraction %g', frac);
		if onlywarn
			warning(t)
		else
			error(t)
		end
	end
	y = real(x);
