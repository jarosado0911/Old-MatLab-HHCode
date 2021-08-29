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
        
        %% for making plots
        function make_plot(x,u,nn,mm,hh,t)
            s1 = subplot(1,4,1);
            cla(s1)
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
            
            s2 = subplot(1,4,2);
            cla(s2)
            plot(u,nn,'k');
            ylim([0 1])
            xlim([-20 120])
            xlabel('voltage [mV]')
            ylabel('n-state []')
            
            s3 = subplot(1,4,3);
            cla(s3)
            plot(u,mm,'k');
            ylim([0 1])
            xlim([-20 120])
            xlabel('voltage [mV]')
            ylabel('m-state []')
            
            s4 = subplot(1,4,4);
            cla(s4)
            plot(u,hh,'k');
            ylim([0 1])
            xlim([-20 120])
            xlabel('voltage [mV]')
            ylabel('h-state []')
            
            drawnow
        end
    end
end

