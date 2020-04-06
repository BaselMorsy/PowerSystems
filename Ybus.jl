function Y_Bus(LineData,N)

    Y_bus = zeros(N,N)*0*im
    M = size(LineData,1);
    otherBus = 0;
    for i in 1:N
        BusIndex = i;
        for j in 1:M
            b1 = LineData[j,1];
            b2= LineData[j,2];
            b = [b1 b2];

            if b1 == BusIndex || b2 == BusIndex
                for k in 1:2
                    if b[k] != BusIndex
                        otherBus = b[k];
                    end
                end

                r = LineData[j,3];
                x = LineData[j,4];

                Z = r + x*im;
                
                Y_sh = LineData[j,5]*im;
                Y_bus[BusIndex,BusIndex] = Y_bus[BusIndex,BusIndex] + 1/Z + 0.5*Y_sh;
                Y_bus[BusIndex,otherBus] = Y_bus[BusIndex,otherBus] - 1/Z
             end
         end
     end

     b = zeros(N,N)
    for l in 1:M
        b[LineData[l,:fbus],LineData[l,:tbus]] = LineData[l,:b]
        b[LineData[l,:tbus],LineData[l,:fbus]] = LineData[l,:b]
    end

     return Y_bus, b
 end
