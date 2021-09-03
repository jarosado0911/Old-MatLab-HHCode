classdef sim_functions
    methods(Static)
        %% Pure Forward Euler Step Function
        function out = FE(X1,X2,K)
            out = X1 + K.*X2;
        end

        %% RK4 Step
        function out = RK4(ss,vv,k,f)
            y1 = ss;
            y2 = y1 + 0.5.*k.*f(y1,vv);
            y3 = y1 + 0.5.*k.*f(y2,vv);
            y4 = y1 + 1.0.*k.*f(y3,vv);
            out = y1 + (k/6.0).*(f(y1,vv) + 2.0.*f(y2,vv) + 2.0.*f(y3,vv)+f(y4,vv));
        end
        
        %% RK4 for reaction with diffusion, set A = 0 if only reaction
        function out = RK4_react_diff(nn,mm,hh,vv,k,f,A)
            y1 = vv;
            y2 = y1 + 0.5.*k.*(A*y1+f(nn,mm,hh,y1));
            y3 = y1 + 0.5.*k.*(A*y2+f(nn,mm,hh,y2));
            y4 = y1 + 1.0.*k.*(A*y3+f(nn,mm,hh,y3));
            out = y1+(k/6.0).*(A*y1+f(nn,mm,hh,y1) + ...
                               2.0.*(A*y2+f(nn,mm,hh,y2))+ ...
                               2.0.*(A*y3+f(nn,mm,hh,y3))+ ...
                               1.0.*(A*y4+f(nn,mm,hh,y4)));
        end
    end
end

