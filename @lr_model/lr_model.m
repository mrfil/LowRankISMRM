function ob = lr_model(time_basis,ky_ind,N)
%function ob = lr_model()
%	Construct a low rank MRI object, which can do Ax and A'y operations

%	default object
ob.time_basis = 0;
ob.ky_ind = 0;
ob.N = 0;
ob.N_tpt = 0;
ob.L = 0;
ob.is.empty	= logical(1);
ob.is.transpose = logical(0);
%ob.version = 1.0;

if nargin == 0
	ob = class(ob, 'lr_model');
	return
end


if nargin ~= 3  %7
	help lr_model
	error nargin
end

	%	fill object
    ob.time_basis = time_basis;
	ob.ky_ind = ky_ind;
    ob.N_tpt = size(time_basis,1);
    ob.L = size(time_basis,2);
    ob.N = N;

	ob.is.empty	= logical(0);

	ob = class(ob, 'lr_model');



