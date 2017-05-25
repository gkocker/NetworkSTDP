%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% bisection for self-consistent steady-state of EIF neuron with voltage-activated conductance


function [P0,p0,J0,r0,x0] = theory0(mu_in, sigma2, params, xi)   

gx = params(13);
x0_in = [0;.5;1];
x0_out = zeros(3,1);
tol = 1e-6;
it = 0;

if gx > 0 % bisection for self-consistent x0
    while x0_in(3)-x0_in(1)>tol
        it = it+1;
        x0_in(2) = x0_in(1)+.5*(x0_in(3)-x0_in(1));
        
        if it == 1
            for xloop = 1:3
                %            [~,~,~,~,x0_out(xloop)] = thin_x_mcode(x0_in(xloop),mu_in,gx);
                [~,~,~,~,x0_out(xloop)] = thin_x(params,x0_in(xloop),mu_in,sigma2,xi);
            end
        else
            %        [~,~,~,~,x0_out(2)] = thin_x_mcode(x0_in(2),mu_in,gx);
            [~,~,~,~,x0_out(2)] = thin_x(params,x0_in(2),mu_in,sigma2,xi);
        end
        
        if sign(x0_out(1)-x0_in(1)) ~= sign(x0_out(2)-x0_in(2))
            x0_in(3) = x0_in(2);
            x0_out(3) = x0_out(2);
        elseif sign(x0_out(3)-x0_in(3)) ~= sign(x0_out(2)-x0_in(2))
            x0_in(1) = x0_in(2);
            x0_out(1) = x0_out(2);
        else
            error('Root not in bracket');
        end
        
    end
end
    
%%% self-consistent solution
% [P0 J0 p0 r0 x0] = thin_x_mcode(x0_out(2),mu_in,gx);
[P0,p0,J0,r0,x0] = thin_x(params,x0_out(2),mu_in,sigma2,xi);


end