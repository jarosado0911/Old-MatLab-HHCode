classdef sim_functions
    methods(Static)
        %% Pure Forward Euler Step Function
        function out = FE(X1,X2,K)
            out = X1 + K.*X2;
        end

        %% RK Step
        function out = RK4(ss,vv,k,f)
            y1 = ss;
            y2 = y1 + 0.5.*k*f(y1,vv);
            y3 = y1 + 0.5.*k*f(y2,vv);
            y4 = y1 + 1.0.*k*f(y3,vv);
            out = y1 + (k/6.0).*(f(y1,vv) + 2.0.*f(y2,vv) + 2.0.*f(y3,vv)+f(y4,vv));
        end
        
        function make_plot(x,u,t)
            clf()
            hold on
            plot(x,u,'k','Linewidth',1.5)
            ylim([-20 120])
            xlim([x(1) x(end)])
            title(sprintf('time = %0.2f [ms]',t))

            scatter(x,0.*x-15,80,u,'filled','s');
            xlim([x(1) x(end)])
            axis off
            colormap jet
            caxis([-15 100])
            colorbar('southoutside') 
            drawnow
        end
    end
end

