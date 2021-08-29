classdef gating_functions
       methods(Static)
        %% Gating functions
        function out=an(vin)
            vin = vin.*1e3;
            out=(0.01).*(10-65-vin)./(exp((10-65-vin)./10)-1);
            out = out*1e3;
        end

        function out=bn(vin)
            vin = vin.*1e3;
            out=(0.125).*exp((-65-vin)./80);
            out = out*1e3;
        end

        function out=am(vin)
            vin = vin.*1e3;
            out=(0.1).*(25-65-vin)./(exp((25-65-vin)./10)-1);
            out = out*1e3;
        end

        function out=bm(vin)
            vin = vin.*1e3;
            out=4.*exp((-65-vin)./18);
            out = out*1e3;
        end

        function out=ah(vin)
            vin = vin.*1e3;
            out=(0.07).*exp((-65-vin)./20);
            out = out*1e3;
        end

        function out=bh(vin)
            vin = vin.*1e3;
            out=1./(exp((30-65-vin)./10)+1);
            out = out*1e3;
        end
    end
end

