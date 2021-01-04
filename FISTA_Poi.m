function sol = FISTA(A, u, b, lambda, Lf, options,qmf,S,divX,lambdamax)
% This routine solves the l1 minimization
% problem for poisson distributed noise
% using Fast Iterative Shrinkage
% Tresholding Algorithm (FISTA) 
% 
%
%
% Date: 22 Feb 2012
% Last Updated: 2 Nov 2014
%
% input : 
% A     : an n-by-m explicit matrix or an operator. An operator 
% (refer to ASP code of Michael Saunders and see class '@operator')
% example:
% A  = operator(@(x,mode) myfunc(x, mode),m,n);
% If mode == 1, returns  v1 =   A * x;
% if mode == 2, returns  v2 =   A'* x; 
%
% b	: data 
% lambda:  regularization parameter
% x0	: initial guess (complex column vector)
% crit	: relative variation of objectove function (for insatnce 10^-6, 10^-7) 
% Lf (Lipschitz constant)   ( [] : equivalent to backtracking )
%
% options: option structure (check  amp_setparams)
%
% \frac{|F(x_k) - mu_F(x_k)|}{mu_F(x_k)} < crit,  where 'mu_F' 
% is the 'average' of min{k,10} prior values of F(x)
% This termination criterion is based on NESTA paper (Becker, Bobin, Candes)
%
% output
% sol.x : argument that minimizes F(x) = f(x) + g(x) when 
% f(x) =  sum_i { (Ax)_i - b_i * log(Ax_i)} and g = lambda * || x ||_1
% sol.k : number of fista iterations 
%
% Copyright 2012, Hatef Monajemi (monajemi@stanford.edu)
% http://www.stanford.edu/~monajemi

% Check number of arguments
if(nargin < 3 || nargin> 10 )
error('incorrect number of arguments')
end

% check Lambda
if(nargin > 3 && ~isempty(lambda) )
lambda_auto = false;
else
lambda_auto = true ;
end



% Check Lf
if(nargin >4 && ~isempty(Lf) )
    back_track = false;
else
    back_track = true ;
end

% Set default parameters or grab user-defined params
if (nargin < 6 || isempty(options) )
   opts = fista_setparams();
else
   opts = fista_setparams(options);
end

%-----------------------------------
%        Extract info from opts
%------------------------------------
x0 		= opts.fistax0;
ind1    = opts.ind1;
ind0    = opts.ind0;
k_max	= opts.fistaitr;
qPrint  = opts.fistalog;
qPeriod = opts.fistalogperiod;
crit	= opts.fistatol;
blocksize = opts.fistablocksize;
indpos  = opts.indpos;

if (isempty(x0)) 
    x0 = zeros(size(A,2),1);
end

if isempty(ind0)
    ind0=1:length(x0);
    ind0(ind1)=[];
end

if(lambda_auto)
    tau = opts.tau;
end

N = size(A,2)/blocksize;
n = size(A,1)/blocksize;


delta = n/N; % undersampling ratio
%---------------------------------
% define objective function F
% and Quadratic approximation of it
% at a given point Y
%---------------------------------

%func =  @(A,x,b)  sum(abs(A*x) - b .* log(abs(A*x)));
%gradfunc = @(A,x,b) A' * ( diag(1./(abs(A*x))) * (abs(A * x) - b) ) ;

func =  @(A,x,u,b,w)  sum(w.*(A*x+u - b .* log(A*x+u)));
gradfunc = @(A,x,u,b,w) A' * ( w.*(A * x+u - b) ./ (A*x+u) ) ;

w=ones(1,length(b))';
%load weights.mat;w=w';w=w.^0.5;

%F_Func = @(A,x,b, lambda)  func(A,x,u,b) + lambda * norm(x,1);
F_Func = @(A,x,b,w, lambda)  func(A,x,u,b,w) + lambda * norm(x(ind0),1);
Q_Func = @(A,x,b,w,lambda,y,Lp) func(A,y,u,b,w) + (x-y)' *  gradfunc(A,y,u,b,w) ...
                          +0.5*Lp*(norm(x-y,2))^2  + lambda* norm(x(ind0),1);
%+ 0.5*Lp*(norm(x-y,2))^2  + lambda* norm(x,1);



if(~isempty(opts.truesol))
convergence_track = true;
else
convergence_track = false;
end



tol               = Inf;
k 				  = 0;
s0  			  = 1;
y   			  = x0;
prior_FuncVals    = zeros(1,10);

%-----------------------------------------------------------------------
% 		Start solution field.
%-----------------------------------------------------------------------
sol.x            = [];
sol.k            = [];
sol.nnz          = [];
sol.sigma        = [];
if(convergence_track)
%sol.convergence     = [];
sol.err             =[];
sol.funcval        = [];
end


%-----------------------------------------------------------------------
% 		Print log header.
%-----------------------------------------------------------------------
    if (qPrint==1) 
      fprintf('\n');
      fprintf(' %s\n',repmat('=',1,90));
      fprintf(' FISTA  (%s)\n', date);
      fprintf(' Copyright 2012, Hatef Monajemi (monajemi@stanford.edu)\n');
      fprintf(' %s\n',repmat('=',1,90));
      if(~lambda_auto)
      fprintf(' %-20s: %8.2e\n'  	   	,'lambda'                ,lambda		);
      else
      fprintf(' %-20s: %8s\n'  	   	,'lambda'                ,'auto selection'	);
      end
      fprintf(' %-20s: %8i %5s'    		,'No. rows'          	 ,size(A,1) ,''	);
      fprintf(' %-20s: %8i \n'    		,'No. columns'       	 ,size(A,2)    	);
      fprintf(' %-20s: %8.2e %5s'    	,'FISTA Optim. tol'    	 ,crit 		,''	);
      fprintf(' %-20s: %8i \n'    		,'FISTA Max itr'		 ,k_max			);
      fprintf(' %s\n',	repmat('=',1,90));
      fprintf('%-10s %-20s %-20s %-10s %-20s %10s\n', 'itr #','rel. err','func. value', 'Lf', 'nnz', 'lambda')
     end


F_k=0;
if (back_track==false)
% FISTA WITHOUT BACKTRACKING  
	while (abs(tol) > crit && k< k_max && imag(F_k)==0)
  
 		k 		= k+1;
 		tmp1 	= gradfunc(A,y,u,b);
        tmp 	= y - tmp1./Lf;

        sol.sigma(k) = EstimateNoise(tmp1, delta, blocksize);

        if(lambda_auto)
        lambda    = tau * sol.sigma(k) ;
        sol.lambda(k) = lambda;
        end

        if(k == 1)
        prior_FuncVals(1) = F_Func(A,x0,b, lambda);
        end
 		% Tresholding
        
        tmptmp=tmp(2:size(tmp,1));
 		x1tmp  		= proxMap(tmptmp, lambda/Lf, 'SOFT', blocksize);
        x1=[tmp(1);x1tmp];
        %x1  		= proxMap(tmp, lambda/Lf, 'SOFT', blocksize);
        
        % Keep only the biggest n entries (HARD)
        %[~,Idx]      = sort(abs(x1), 'descend');
        %x1(Idx(n+1:end))     = 0;

        %=============== UPDATES =================
        % stepsize
 		%s1 		= (1+ sqrt(1+ 4*s0^2))/2;
        s1          = 1;
		 % Psudo point
 		y 		= x1 + ((s0-1)/s1) .* (x1 - x0) ;

 		% tolerance
 		F_mu 	= sum(prior_FuncVals) / nnz(prior_FuncVals);
 		F_k  	= F_Func(A,x1,b,lambda)
 		tol 	= abs( F_k - F_mu) / F_mu;
 		TOL(k) = tol;
 		sol.nnz(k) = findnnz(x1,blocksize);
 		if( prior_FuncVals( mod(k,10) + 10*(mod(k,10)==0)) == F_k )
 		% we detected A'*A = I (orthonormal transform)
 		    tol = 0.0;
 		    fprintf('%-10i %-20e %-20e %-10f %-20i %-10e\n', k, tol, F_k, Lf, sol.nnz(k), lambda); 
 		  	fprintf('Warning: It seems the sensing matrix is orthogonal.\n');
 		  	break;
 		else
 		% prior values of objective function for termination purpose
 			prior_FuncVals(1 + mod(k,10)) = F_k;
 		end
 
 		% update varaibles for the next iteration
		x0 	= x1;
 		s0 	= s1;  
 	 
       
        
        if(convergence_track)
        sol.err(k)      = norm(opts.truesol - x0)/norm(opts.truesol);
        %sol.convergence(k) = F_k  - F_Func(A,opts.truesol,b,opts.truelambda);
        sol.funcval(k) = F_k;
        end
 	    if(mod(k,qPeriod)==0 & qPrint==1)
 	    fprintf('%-10i %-20e %-20e %-10f %-20i %-10e\n', k, tol, F_k, Lf, sol.nnz(k), lambda);
        end

        

	end
	
else
% FISTA WITH BACKTRACKING
   	    Lf   = 1;
    	eta  = 2.0;

        lambdaopt=lambda;
        %if(~isempty(lambdamax)) lambda=lambdamax; end
%         if lambda==0 
%             slambda=1;
%         else
%             slambda=(lambdamax/lambdaopt)^(1/100);
%         end
		while ( (abs(tol) > crit && k< k_max ))

 		k 		= k+1;
        tmp1 	= gradfunc(A,y,u,b,w);
        sol.sigma(k) = EstimateNoise(tmp1, delta, blocksize);
%         if(~isempty(lambdamax))
%             lambda=lambda/slambda;
%             if(lambda<=lambdaopt) slambda=1; end
%         end       
        if(lambda_auto)
        lambda    = tau * sol.sigma(k) ;
        sol.lambda(k) = lambda;
        end
    
        
        
        
        if(k == 1)
        prior_FuncVals(1) = F_Func(A,x0,b,w, lambda);
        end
        
 			% backtracking for step k
 			ik   = 0;
 			beta = Inf;
            
 			%while(beta > 0||imag(FktempBT)~=0||controllog~=0)
            
            while(beta > 0||imag(FktempBT)~=0||controllog~=0)

                L_bar   = (eta^ik ) * Lf; 
                if L_bar==Inf
                    sol.conv=0;
                    return
                end
                tmp 	= y - tmp1./L_bar;
                %tmp(1)=y(1)-tmp1(1);
                % Tresholding 

                %tmp(1)=115;

                %tmptmp=tmp(2:size(tmp,1));

                tmptmp=tmp(ind0);

                %tmptmp=tmp;
                x1tmp 		= proxMap(tmptmp, lambda/L_bar, 'SOFT', blocksize);
                x1=[tmp(1);x1tmp];

                x1=tmp;
                x1(ind1)=tmp(ind1);
                x1(ind0)=x1tmp;
                %if(x1(1)<0) x1(1)=0; end;


%                 if length(x1)>length(b)
%                     dps=x1((length(x1)-length(b)+1):length(x1));
%                     dps(dps<0)=0;
%                     %dps(dps<max(dps)/100)=0;
%                     x1((length(x1)-length(b)+1):length(x1))=dps;
%                 end
                x11pos=x1(indpos);
                x11pos(x11pos<0)=0;
                x1(indpos)=x11pos;
                %x1(1)= 2.2642e-12*divX(1);
                %x1(abs(x1)<1e-4)=0;

                % Keep only the biggest n entries (hard)
                %[~,Idx]      = sort(abs(x1), 'descend');
                %x1(Idx(380+1:end))     = 0;


                FktempBT=F_Func (A, x1, b,w, lambda);
                % update backtracking parameters
                beta 	= real(FktempBT - Q_Func(A, x1, b,w, lambda, y, L_bar));
                controllog=sum(imag(log(A*x1+u))~=0);

                %fff=IWT_PO(x1./divX',0,qmf);
                ik 		= ik+1;	
            end
            
            if(length(x1)>length(b))
                ps=1;
            else
                ps=0;
            end
            
%             subplot(1,2,1);
% 
%             XX=x1./divX';
%             XX=XX(1:(length(XX)-ps*length(b)));
%             yy1=0;yy2=0;
%             
%             plot(find(XX~=0),XX(XX~=0),'.')
%             hold on
%             plot(find(XX==0),XX(XX==0),'*r');
%             hold off
% 
%             %plot current estimation
%             aaa=(~isempty(qmf))+(~isempty(S));
%             if(lambda~=0)
%                 if(~isempty(qmf))
%                     if(~isempty(S))
%                         yy1=IWT_PO(XX(1:(length(XX)/2)),0,qmf);
%                     else
%                         yy1=IWT_PO(XX(1:(length(XX))),0,qmf);
%                     end
%                 end
%                 if(~isempty(S))
%                     if(~isempty(qmf))
%                         yy2=S*XX((length(XX)/2+1):end);
%                     else
%                         yy2=S*XX;
%                     end
%                 end
%                 subplot(1,2,2)
%                 yy=yy1+yy2;
%                 plot(yy)
%                 drawnow
%             end
          
		% update the initial guess for the next step
        Lf=L_bar;
        if L_bar>0 
            Lf   = Lf/2;
        end
        if mod(k,100)==0
            if Lf>=1
                Lf=1;
            end
        end
            
			
	    %=============== UPDATES =================
		 % stepsize
 		%s1 		= (1+ sqrt(1+ 4*s0^2))/2;
        s1 = 1;
         % Psudo point
 		y 		= x1 + ((s0-1)/s1) .* (x1 - x0) ;

 		% tolerance
 		F_mu 	= sum(prior_FuncVals) / nnz(prior_FuncVals);
 		F_k  	= F_Func(A,x1,b,w,lambda);
 		tol 	= abs( F_k - F_mu) / F_mu;
 		TOL(k) = tol;
 		
        sol.nnz(k) = findnnz(x1,blocksize);

 		if( prior_FuncVals( mod(k,10) + 10*(mod(k,10)==0)) == F_k)
 		% we detected A'*A = I (orthonormal transform)
 		    tol = 0.0;
 		    fprintf('%-10i %-20e %-20e %-10f %-20i %-10e\n', k, tol, F_k, Lf, sol.nnz(k), lambda); 
 		    fprintf('Warning: It seems the sensing matrix is orthogonal.\n');
 		  	break;
 		else
 		% prior values of objective function for termination purpose
 			prior_FuncVals(1 + mod(k,10)) = F_k;
 		end
 		% update varaibles for the next iteration
		x0 	= x1;
 		s0 	= s1;  
        

        
        if(convergence_track)
        sol.err(k)      = norm(opts.truesol - x0)/norm(opts.truesol);
        %sol.convergence(k) = F_k  - F_Func(A,opts.truesol,b,opts.truelambda);
        sol.funcval(k) = F_k;
        end

        
        if(mod(k,qPeriod)==0 & qPrint==1)
 	    fprintf('%-10i %-20e %-20e %-10f %-20i %-10e\n', k, tol, F_k, Lf, findnnz(x1, blocksize), lambda);
        end




       end

end

sol.conv=1;
sol.x  = x1;
sol.k  = k ;

if( k <= k_max && tol < crit)
  fprintf('FISTA: converged in %i iterations (nnz(x)= %i)\n', k, findnnz(x1,blocksize));
  fprintf('FISTA: Relative error: %-.2e\n', tol); 
else  
  warning(['FISTA: not converged to desired tolerance, relative error achieved: ', num2str(tol),'( nnz(x)= ', num2str(findnnz(x1,blocksize)) ,')']);
end


%figure
%semilogy(TOL, 'linewidth', 3)
%set(gca, 'fontsize', 17)
%ylabel('rel. err','fontsize', 18)
%xlabel('iteration','fontsize', 18)
%title('Convergence - FISTA','fontsize', 18)
end

%%% Private Function
function sigma_hat = EstimateNoise(x, delta, blocksize)

nblocks 	 = floor(length(x)/blocksize);
if(nblocks*blocksize~=length(x)); error('bad data'); end;
% Arrange data into a matrix of size (blocksize)*(nblocks);
x 			 = reshape(x,[blocksize,nblocks]);
absX		 = sqrt(sum(abs(x).^2,1));
sigma_hat 		 = (1/0.6745) * median(absX);
end

function  k  = findnnz(x, blocksize)
nblocks 	 = floor(length(x)/blocksize);
x 			 = reshape(x,[blocksize,nblocks]);
absX		 = sqrt(sum(abs(x).^2,1));
k            = length(find(absX>0));
end
