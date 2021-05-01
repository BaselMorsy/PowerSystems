function print_report(System ::System_Struct)

    N_Loads = length(unique!(Array(System.BusData[:Pd])))-1 ;
    P_Loads = sum(System.BusData[:Pd])
    Q_Loads = sum(System.BusData[:Qd])
    Cost = System.Operating_Cost;
    LineLoading = System.LineLoading
    BusData = System.BusData_output
    Pg = BusData[:Pg]
    Qg = BusData[:Qg]
    Pd = BusData[:Pd]
    Qd = BusData[:Qd]
    V = BusData[:V]
    PLoss = LineLoading[:PLoss]
    QLoss = LineLoading[:QLoss]

    println("==================================================================================================")
    println("|          System Summary                                                                        |")
    println("==================================================================================================")
    println("Objective Function Value = ", round(Cost,digits=3), " USD/hr")
    println("How many?                     How much?                    P(MW)                   Q(MVAr)")
    println("------------------           ----------------              ---------------         -------------   ")
    println("Buses          ",System.N_bus,"           Total Gen Capacity               ",
        sum(System.Gen_Data[:,5]),"                  ",
         sum(System.Gen_Data[:,8])," to ",
         sum(System.Gen_Data[:,7]))
    println("Generators     ",size(System.Gen_Data,1),"            On-Line  Capacity               ",sum(System.Gen_Data[:,5]),"                  ",
        sum(System.Gen_Data[:,8])," to ",
        sum(System.Gen_Data[:,7]))
    println("Commited Gens  ",size(System.Gen_Data,1),"          Generation (Actual)               ",round(sum(Pg),digits=3),"               ",round(sum(Qg),digits=3))
    println("Loads          ",N_Loads,"            Load                            ",P_Loads,"                    ",Q_Loads  )
    println()
    println("                               Minimmum                             Maximum ")
    println("                            --------------                       --------------")
    println("Voltage Magnitude            ",round(minimum(V),digits = 3), " @ bus ",findall(x -> x==minimum(V),V),"                       ",round(maximum(V),digits = 3), " @ bus ",findall(x -> x==maximum(V),V) )

    i = findall(x -> x==maximum(PLoss),PLoss);
    println("P Losses                           -                             ",round(maximum(PLoss),digits = 3)," @ Line ",Int64(LineLoading[i,:FromBus][1,1])," - ",Int64(LineLoading[i,:ToBus][1,1]) )
    i = findall(x -> x==maximum(QLoss),QLoss);
    println("Q Losses                           -                             ",round(maximum(QLoss),digits = 3)," @ Line ",Int64(LineLoading[i,:FromBus][1,1])," - ",Int64(LineLoading[i,:ToBus][1,1]) )
    println()
    println("Bus Data")
    println(BusData)
    println("                           ----------------------------------------")
    println("                   Total:  " ,round(sum(Pg),digits=3),"      ",round(sum(Qg),digits=3),"       ",sum(Pd),"       ",sum(Qd))
    println()
    println("Line Loading")
    println(LineLoading)
    println("                                       ---------------------------------------------------------------")
    println("                                       Total:           ", round(sum(PLoss),digits=3), "                                 ",round(sum(QLoss),digits=3))
    println("=========================== END OF REPORT ========================================")
end

function print_latex_data(System ::System_Struct)

    bus_cols = [:bus_i,:Pd]
    ind = findall(x->x!=0,System.BusData[:Pd])
    bus_data = System.BusData[ind,bus_cols]
    show(stdout, MIME("text/latex"),bus_data)

    line_cols = [:fbus,:tbus,:x,:rate]
    line_data = System.LineData[line_cols]
    show(stdout, MIME("text/latex"),line_data)

    gen_cols = [:bus,:C1,:C0,:Pmax,:Pmin]
    gen_data = System.Gen_Data[gen_cols]
    show(stdout, MIME("text/latex"),gen_data)
end
