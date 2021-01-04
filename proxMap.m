function  x = proxMap(z, coef, method, blocksize) 
% Proxmap function
% 
% Copyright 2011, Hatef Monajemi (monajemi@stanford.edu)
% Date: 15 Nov. 2011
% Feb 21 2014: Changed to include blocksize for thresholding 
%
% input 
% z : input vector (can be complex)
% coef: proxMap coefficient
% method: 
%		- 'SOFT'
%       - 'HARD'  
%       - 'JS' 
% blocksize: block size for thresholding. Set 4 for quaternion
%            if you are using complex variable blocksize = 1 will do;
% 			 A quatrenion can also be represented with complex variables (simplex, perplex)
%            and blocksize = 2; It all depends on how you define your data
%            we threshold 'blocksize' elements together. So if you are using 
%            real, choose blocksize= 4 for quaternion.
if(nargin ==1)
  error('number of args must be 2 or 3');
elseif (nargin == 2)
  method = 'SOFT';
  blocksize = 1  ;
elseif (nargin == 3);
  blocksize = 1;
end 


x = zeros(length(z),1);


nblocks = floor(length(z)/blocksize);
if(nblocks*blocksize~=length(z)); error('bad data'); end;
% Arrange data into a matrix of size (blocksize)*(nblocks);
z = reshape(z,[blocksize,nblocks]);
% threshold each column separately
absZ = sqrt(sum(abs(z).^2,1));  % this is a row vector of size "nblocks" -- abs to include complex variables

% compute for each column the thresholding factor
	switch method
	case 'JS'
		error('Not yet implemented');
	case 'SOFT'
		thrshFactor = max(1- coef./absZ, 0);
		x  			= z * sparse(1:length(thrshFactor),1:length(thrshFactor),thrshFactor);
	case 'PosSOFT'
		thrshFactor = max(1- coef./absZ, 0);
		x 			= (z.*(z>0)) * sparse(1:length(thrshFactor),1:length(thrshFactor),thrshFactor);
	case 'HARD'
		thrshFactor = (absZ > coef);
		x		    = z * sparse(1:length(thrshFactor),1:length(thrshFactor),thrshFactor);	
	otherwise
	warning('Unexpected option');
	end

% reshape it back to a vector
x = x(:);


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original definition (older version)
%if (blocksize==1)
%	switch method
%	case 'JS'
%		error('Not yet implemented');
%	case 'SOFT'
%		x = z .* max(1- coef./abs(z), 0) ; 
%	case 'PosSOFT'
%		x = z .* max(1- coef./abs(z), 0) .* (z > 0)	;
%	case 'HARD'
%		x = z .* ( abs(z) > coef);	
%	otherwise
%	warning('Unexpected option');
%	end
% end
