classdef gating_functions
       methods(Static)
        %% Gating functions
        function out=an(vin)
        out=(0.01).*(10-vin)./(exp((10-vin)./10)-1);
        end

        function out=bn(vin)
        out=(0.125).*exp(-vin./80);
        end

        function out=am(vin)
        out=(0.1).*(25-vin)./(exp((25-vin)./10)-1);
        end

        function out=bm(vin)
        out=4.*exp(-vin./18);
        end

        function out=ah(vin)
        out=(0.07).*exp(-vin./20);
        end

        function out=bh(vin)
        out=1./(exp((30-vin)./10)+1);
        end
    end
end

