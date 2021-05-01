include("Ybus.jl")

import JuMP, Ipopt
using JuMP, Ipopt

module System_Definition
    using Parameters
    using DataFrames
    export System_Struct,SubStation

    @with_kw mutable struct System_Struct
        y_bus
        b_line
        LineData
        BusData
        N_bus
        Sbase
        load_multiplier = 1
        generation_multiplier = 1

        Operating_Cost = Float64[]
        state_error=Float64[]
        BusData_output = DataFrame(Bus = Int64[], V = Float64[], δ = Float64[],
         Pg = Float64[],Qg = Float64[],Pd = Float64[], Qd = Float64[])

        LineLoading = DataFrame(FromBus = Int64[],ToBus = Int64[]
            ,PL_1 = Float64[],PL_2 = Float64[],
            PLoss = Float64[],QL_1 = Float64[]
            ,QL_2 = Float64[],QLoss = Float64[])

        Line_Constraints= DataFrame(fbus = Union{Float64,Missing}[],
            tbus = Union{Float64,Missing}[], SL_Rating = Union{Float64,Missing}[])

        Gen_Data =  DataFrame(bus = Int64[],
                Qmax = Union{Float64,Missing}[],
                Qmin = Union{Float64,Missing}[],
                Pmax = Union{Float64,Missing}[],
                Pmin = Union{Float64,Missing}[],
                C2 = Union{Float64,Missing}[],
                C1 = Union{Float64,Missing}[],
                C0 = Union{Float64,Missing}[])
        Pd_total = Float64[];

        static_line_index = []
        non_static_line_index = []

        Line_Duals = []
        Bus_Duals = []

        SubStations = []

        Active_buses = []
        Inactive_buses = []
        A =[]#Incidence matrix
        A_ex = [] #extended Incidence Matrix
        Y_s=[] #Sorted admittance vector
        Y_pr=[] #primitive admittance matrix
        I=[] #branch index matrix
        Δ=[] #duplicat buses array
        Y_bus_inc=[] #Ybus from Incidence matrix
        D = [] #No. of duplicates for each bus
        z = [] #Branch/duplicates relatation array (for tbus duplicat)
        ρ = [] #Branch/duplicates relational array (for fbus duplicat)
    end

    @with_kw mutable struct SubStation
        SubStationID = []
        flag = false
        Num_BusBars = []
        SubStation_Buses = []
        Aux_Buses = []
        Non_Aux_Buses =[]
        SubStation_AuxLines = []
        SubStation_NonAuxLines = DataFrame(fbus=[],tbus=[],r=[],x=[],b=[],rate=[],ID=[],is_aux=[])
        SubStation_Gen = []
        SubStation_Loads = []
        Lines_at_bus = []



    end

end  # module


using Main.System_Definition

function SystemAssembly(LineData,Sbase,BusData_sheet=nothing,
     GenData=nothing)


    Lines = [];
    busData = [];
    YBUS  = [];
    bLine = [];
    genData = [];
    gen_constraints=[];
    S_rating_lines = [];
    N_bus =[];



    if LineData != []
        #Lines = DataFrame(CSV.read(LineData;copycols = true))
        Lines = LineData
        N_bus = length(unique!(vcat(unique!(Vector(Lines[:fbus])),unique!(Vector(Lines[:tbus])))))
        YBUS,bLine = Y_Bus(Lines,N_bus)
        S_rating_lines = DataFrame(fbus = Lines[!,:fbus],tbus = Lines[!,:tbus]
        ,SL_Rating = Lines[!,:rate])
    end

    substations = Dict();

    if BusData_sheet != nothing
        #busData = DataFrame(CSV.read(BusData_sheet;copycols = true))
        busData = BusData_sheet
        for i in 1:length(busData[:,1])
            substations[string("SubStation_",i)] = System_Definition.SubStation(SubStationID=busData[i,1],SubStation_Buses = i,
                SubStation_Gen = GenData[findall(x->x==i,GenData[:,:bus]),:],SubStation_Loads=busData[findall(x->x==i,busData[:,:bus_i]),:],
                SubStation_NonAuxLines = LineData[vcat(findall(x->x==i,LineData[:,1]),findall(x->x==i,LineData[:,2])),:])
        end
    end

    if GenData != nothing
        # gendata = DataFrame(CSV.read(GenData,copycols = true))
        # gen_data = DataFrame(bus = gendata[:bus],C2 = gendata[:c2],
        #     C1 = gendata[:c1],C0 = gendata[:c0],
        #     Pmax = convert.(Float64,gendata[:Pmax]),
        #     Pmin = convert.(Float64,gendata[:Pmin]),
        #     Qmax = convert.(Float64,gendata[:Qmax]),
        #     Qmin = convert.(Float64,gendata[:Qmin]));
        gendata = GenData
    end
    Sys = System_Struct(y_bus = YBUS,b_line = bLine, LineData = Lines,BusData =busData
        ,N_bus=N_bus,Sbase = Sbase,SubStations = substations)
    append!(Sys.Line_Constraints,S_rating_lines);
    Sys.Gen_Data=GenData;
    Sys.Y_bus_inc = YBUS
    Sys.Pd_total = sum(busData[:Pd]);
    return Sys
end

function show_cases(verbose = true)

    lst = readdir(string(pwd(),"\\Systems\\csv"))
    if verbose
        for i in 1:length(lst)
            println("[INFO] ", lst[i]," - ", i)
        end
    end
    return lst
end

function load_case(case_id,Sbase)

    lst = show_cases(false)
    CASE_DIR = string(string(pwd(),"\\Systems\\csv"),"/",lst[case_id])

    line_data_raw = DataFrame(CSV.read(string(CASE_DIR,"/","branch.csv");copycols = true))
    LineData = DataFrame(fbus = line_data_raw[:,1],tbus = line_data_raw[:,2],r = line_data_raw[:,3],x = line_data_raw[:,4],b = line_data_raw[:,5],rate = line_data_raw[:,6]./Sbase,ID = collect(1:length(line_data_raw[:,6])), is_aux = zeros(length(line_data_raw[:,1]),))

    bus_data_raw = DataFrame(CSV.read(string(CASE_DIR,"/","bus.csv");copycols = true))
    BusData = DataFrame(bus_i = bus_data_raw[:,1],type = bus_data_raw[:,2],Pd = bus_data_raw[:,3],Qd = bus_data_raw[:,4],Vmax = bus_data_raw[:,12],Vmin = bus_data_raw[:,13],is_aux = zeros(length(bus_data_raw[:,1]),))

    gen_raw =  DataFrame(CSV.read(string(CASE_DIR,"/","gen.csv");copycols = true))
    gen_cost_raw =  DataFrame(CSV.read(string(CASE_DIR,"/","gencost.csv");copycols = true))
    GenData = DataFrame(bus = gen_raw[:,1],C2 = gen_cost_raw[:,5],C1 = gen_cost_raw[:,6],C0 = gen_cost_raw[:,7],Pmax = gen_raw[:,9],Pmin = gen_raw[:,10], Qmax = gen_raw[:,4], Qmin = gen_raw[:,5])


    return SystemAssembly(LineData,Sbase,BusData,GenData)
end

function set_load_multiplier!(System ::System_Struct,value)

    System.BusData[!,3] = convert.(Float64,System.BusData[!,3])
    System.BusData[!,4] = convert.(Float64,System.BusData[!,4])

    System.BusData[!,3] = System.BusData[!,3]*(value)
    System.BusData[!,4] = System.BusData[!,4]*(value)
    SubStations = deepcopy(System.SubStations)
    for substation in SubStations
        substation[2].SubStation_Loads[!,:Pd] = substation[2].SubStation_Loads[!,:Pd]*(value)
        substation[2].SubStation_Loads[!,:Qd] = substation[2].SubStation_Loads[!,:Qd]*(value)
    end
    System.load_multiplier = value
    System.SubStations = SubStations
end

function set_generation_multiplier!(System ::System_Struct,value)

    System.Gen_Data[!,5] = convert.(Float64,System.Gen_Data[!,5])
    System.Gen_Data[!,7] = convert.(Float64,System.Gen_Data[!,7])
    System.Gen_Data[!,8] = convert.(Float64,System.Gen_Data[!,8])

    System.Gen_Data[!,5] = System.Gen_Data[!,5]*(value/System.generation_multiplier)
    System.Gen_Data[!,7] = System.Gen_Data[!,7]*(value/System.generation_multiplier)
    System.Gen_Data[!,8] = System.Gen_Data[!,8]*(value/System.generation_multiplier)
    System.generation_multiplier = value
end

function set_line_capacity_multiplier!(System ::System_Struct,value)
    System.Line_Constraints[!,3] = System.Line_Constraints[!,3]*value
    System.LineData[!,:rate] = System.LineData[!,:rate]*value
end

function Solve_OPF!(System ::System_Struct, method = "ACOPF", print_status = true)

    #OPF solver
    THETAMAX = 0.6
    G = real(System.y_bus);
    B = imag(System.y_bus);
    b = System.b_line;

    N = System.N_bus
    N_Lines = size(System.LineData[:,1],1)
    Nodes_set = System.BusData[:,:bus_i];
    Gen_set = 1:length(System.Gen_Data[!,:bus]);

    Sbase = System.Sbase;

    BusData = System.BusData;
    GenData = System.Gen_Data;
    BranchData = System.LineData;

    #equivalence array
    eq_array = zeros(N,2)
    for i in 1:N
        eq_array[i,1] = i
        eq_array[i,2] = System.BusData[i,1]
    end
    S = zeros(N,N)
    eq_array_lines = zeros(N_Lines,2)

        for l in 1:N_Lines
            #fbus_real_index = findfirst(x->x==Int64(BranchData[l,:fbus]),System.Δ[:,3])
            #tbus_real_index = findfirst(x->x==Int64(BranchData[l,:tbus]),System.Δ[:,3])

            #fbus_real = System.Δ[fbus_real_index,1]
            #tbus_real = System.Δ[tbus_real_index,1]
            fbus = Int64(eq_array[findfirst(x->x==Int64(BranchData[l,:fbus]),eq_array[:,2]),1])
            tbus = Int64(eq_array[findfirst(x->x==Int64(BranchData[l,:tbus]),eq_array[:,2]),1])
            S[fbus,tbus] = BranchData[l,:rate]
            S[tbus,fbus] = BranchData[l,:rate]
            eq_array_lines[l,1] = fbus
            eq_array_lines[l,2] = tbus
        end

    if method == "ACOPF"


            m = JuMP.Model(JuMP.with_optimizer(Ipopt.Optimizer, print_level =0))

            # 2.1.Variables
            JuMP.@variable(m, BusData[i, :Vmin] ≤ v[i in Nodes_set] ≤ BusData[i, :Vmax])
            JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set] ≤ THETAMAX)

            JuMP.@variable(m, GenData[(GenData[:bus].==g), :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[(GenData[:bus].==g), :Pmax][1,1])
            JuMP.@variable(m, GenData[(GenData[:bus].==g), :Qmin][1,1] ≤ q[g in Gen_set] ≤ GenData[(GenData[:bus].==g), :Qmax][1,1])

            JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set])
            JuMP.@variable(m, qij[i = Nodes_set, j = Nodes_set])

            # 2.2. Constraints
            #JuMP.@constraint(m, ReferenceAngle,
            #    (δ[1] ==  0.0))
            JuMP.@constraint(m, ReferenceAngle,
                (δ[Gen_set[1]] ==  0.0))

            # ACTIVE POWER THROUGH LINE N-M
            JuMP.@NLconstraint(m, p_line[i in Nodes_set, j in Nodes_set],
                 (pij[i,j] ==  Sbase*(v[i]*v[j]*(G[i,j]*cos(δ[i]-δ[j])+B[i,j]*sin(δ[i]-δ[j])) -(v[i]^2)*G[i,j] )))

            # REACTIVE POWER THROUGH LINE N-M
            JuMP.@NLconstraint(m, q_line[i in Nodes_set, j in Nodes_set],
                 (qij[i,j] ==  Sbase*(v[i]*v[j]*(G[i,j]*sin(δ[i]-δ[j]) - B[i,j]*cos(δ[i]-δ[j])) +(v[i]^2)*(B[i,j]-0) )));

            # ACTIVE NODAL BALANCE
            JuMP.@constraint(m, Pnodal[i in Nodes_set],
                sum(pij[i,j] for j = Nodes_set) == sum(p[g] for g in Gen_set if GenData[(GenData[:bus].==g),:bus][1,1]==i) - BusData[i,:Pd])

            # REACTIVE NODAL BALANCE
            JuMP.@constraint(m, Qnodal[i in Nodes_set],
                sum(qij[i,j] for j = Nodes_set) == sum(q[g] for g in Gen_set if GenData[(GenData[:bus].==g),:bus][1,1]==i) - BusData[i,:Qd])

            # LINE CAPACITY

            JuMP.@NLconstraint(m, Smax[i in Nodes_set, j in Nodes_set],
                pij[i,j]^2 + qij[i,j]^2 ≤ ((Sbase)^2)*S[i,j]^2)

            #OBJECTIVE

            JuMP.@objective(m,Min,sum(GenData[(GenData[:bus].==g),:C2][1,1]*p[g]^2+GenData[(GenData[:bus].==g),:C1][1,1]*p[g]+GenData[(GenData[:bus].==g),:C0][1,1] for g in Gen_set))

            JuMP.optimize!(m)
            solution_status = raw_status(m)

            if solution_status != "Infeasible_Problem_Detected"

                save_duals = false
                if save_duals == true
                    Line_Duals = []
                    for i in 1:N_Lines
                        fbus = System.LineData[i,1]
                        tbus = System.LineData[i,2]
                        current_dual = JuMP.dual(p_line[fbus,tbus])
                        Line_Duals = append!(Line_Duals,current_dual)
                    end

                    Bus_Duals = []
                    for i in 1:N
                        current_dual_b = JuMP.dual(Pnodal[i])
                        Bus_Duals = append!(Bus_Duals,current_dual_b)
                    end

                    System.Bus_Duals = Bus_Duals
                    System.Line_Duals = Line_Duals
                end




                    Pg = JuMP.value.(p)
                    Qg = JuMP.value.(q)
                    Pij = JuMP.value.(pij)
                    Qij = JuMP.value.(qij)
                    V = JuMP.value.(v)
                    δ = JuMP.value.(δ)
                    Pd = System.BusData[!,3]
                    Qd = System.BusData[!,4]

                    System.Operating_Cost = JuMP.objective_value(m)
                    println(Pg)
                    Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1])]) : 0 for i in 1:System.N_bus]
                    Qg = [i in System.Gen_Data[!,1] ? Qg.data[findall(x -> x==i,System.Gen_Data[!,1])][1,1] : 0 for i in 1:System.N_bus]
                    df = DataFrame(Bus = 1:System.N_bus, V = V.data, δ = δ.data,
                     Pg = Pg ,Qg = Qg, Pd = Pd, Qd = Qd)
                    System.BusData_output = df

                    Pij_1 = [ Pij[eq_array_lines[i,1],eq_array_lines[i,2]] for i in 1:size(System.Line_Constraints,1) ]
                    Pij_2 = [ Pij[eq_array_lines[i,2],eq_array_lines[i,1]] for i in 1:size(System.Line_Constraints,1) ]
                    P_Loss = Pij_1 + Pij_2

                    Qij_1 = [ Qij[eq_array_lines[i,1],eq_array_lines[i,2]] for i in 1:size(System.Line_Constraints,1) ]
                    Qij_2 = [ Qij[eq_array_lines[i,2],eq_array_lines[i,1]] for i in 1:size(System.Line_Constraints,1) ]
                    Q_Loss = Qij_1 + Qij_2
                    df2 = DataFrame(FromBus =System.Line_Constraints[!,1] ,ToBus =System.Line_Constraints[!,2]
                        ,PL_1 = Pij_1,PL_2 = Pij_2,
                        PLoss = P_Loss,QL_1 = Qij_1
                        ,QL_2 = Qij_2,QLoss = Q_Loss)
                    System.LineLoading = df2

                    if print_status
                        print_report(System)
                    end
                end

                return solution_status

            elseif method =="DCOPF"

                B_ = -B[1:N,1:N];
                #JuMP.with_optimizer(Gurobi.Optimizer,OutputFlag=0),quiet=true

                m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,OutputFlag=0))


                    # 2.1.Variables
                    JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set] ≤ THETAMAX)

                    JuMP.@variable(m, GenData[g, :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[g, :Pmax][1,1])

                    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set])

                    # 2.2. Constraints
                    JuMP.@constraint(m, ReferenceAngle,
                        (δ[Nodes_set[1]] ==  0.0))

                    JuMP.@constraint(m,Nodal_balance[i in Nodes_set],
                            sum(pij[i,j] for j = Nodes_set) == sum(p[g] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])
                    JuMP.@constraint(m,pl[i in Nodes_set,j in Nodes_set],
                            pij[i,j] == Sbase*B_[i,j]*(δ[i]-δ[j]))
                    JuMP.@constraint(m,pl_rate[i in Nodes_set,j in Nodes_set],
                            -Sbase*S[i,j] ≤ pij[i,j] ≤ Sbase*S[i,j])

                    JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g]+GenData[g,:C0][1,1] for g in Gen_set))
                    JuMP.optimize!(m)

                    solution_status = raw_status(m)
                    println(solution_status)
                    if solution_status != "Model was proven to be infeasible."
                        Pg = JuMP.value.(p)
                        Pij = JuMP.value.(pij)
                        δ = JuMP.value.(δ)
                        V = ones(N)
                        Pd = System.BusData[!,3]
                        Qd = System.BusData[!,4]
                        System.Operating_Cost = JuMP.objective_value(m)
                        System.Gen_Data[:Pg] = Pg.data
                        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1])]) : 0 for i in 1:System.N_bus]
                        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data,
                         Pg = Pg , Pd = Pd)
                        System.BusData_output = df

                        Pij_1 = [ Pij[eq_array_lines[i,1],eq_array_lines[i,2]] for i in 1:size(System.Line_Constraints,1) ]
                        Pij_2 = [ Pij[eq_array_lines[i,2],eq_array_lines[i,1]] for i in 1:size(System.Line_Constraints,1) ]
                        P_Loss = Pij_1 + Pij_2

                        U = []
                        for i in 1:N_Lines
                            uu = round((abs(Pij_1[i])/(System.Line_Constraints[i,3]*System.Sbase))*100,digits = 2)
                            append!(U,uu)
                        end
                        df2 = DataFrame(FromBus =System.Line_Constraints[!,1] ,ToBus =System.Line_Constraints[!,2]
                            ,PL_1 = Pij_1,PL_2 = Pij_2,
                            PLoss = P_Loss,Utilization = U)
                        System.LineLoading = df2

                        Line_Duals = []
                        for i in 1:N_Lines
                            fbus = System.LineData[i,1]
                            tbus = System.LineData[i,2]
                            current_dual_1 = JuMP.dual(pl_rate[fbus,tbus])
                            current_dual_2 = JuMP.dual(pl_rate[tbus,fbus])
                            if current_dual_1 != 0
                                Line_Duals = append!(Line_Duals,current_dual_1)
                            elseif current_dual_2 != 0
                                Line_Duals = append!(Line_Duals,current_dual_2)
                            else
                                Line_Duals = append!(Line_Duals,0)
                            end


                        end

                        Bus_Duals = []
                        for i in 1:N
                            current_dual_b = JuMP.getdual(Nodal_balance[i])
                            Bus_Duals = append!(Bus_Duals,current_dual_b)
                        end

                        System.Bus_Duals = Bus_Duals
                        System.Line_Duals = Line_Duals

                        if print_status
                            println("==================================================================================================")
                            println("|          System Summary                                                                        |")
                            println("==================================================================================================")
                            println("Objective Function Value = ", round(System.Operating_Cost,digits=3), " USD/hr")
                            println(System.BusData_output)
                            println()
                            println(System.LineLoading)
                            println("=========================== END OF REPORT ========================================")
                        end
                    end

                    return solution_status
    end
end

function Solve_DCOPF!(System ::System_Struct, print_status = true)

    #OPF solver
    THETAMAX = 0.6
    G = real(System.y_bus);
    B = imag(System.y_bus);
    b = System.b_line;

    N = System.N_bus
    N_Lines = size(System.LineData[:,1],1)
    Nodes_set = System.BusData[:,:bus_i];
    Gen_set = 1:length(System.Gen_Data[!,:bus]);

    Sbase = System.Sbase;

    BusData = System.BusData;
    GenData = System.Gen_Data;
    BranchData = System.LineData;

    df_lines = System.LineData[:,[:fbus,:tbus]]
    lines_set = Set.(convert.(Array, row for row in eachrow(df_lines)))
    lines_rates = Dict()
    for i in 1:length(System.LineData[:,1])
        lines_rates[lines_set[i]]=System.LineData[i,:rate]
    end
    Nodes_at_Lines = Dict()
    for i in System.LineData[:ID]
        Nodes_at_Lines[i] = Set([System.LineData[i,:fbus],System.LineData[i,:tbus]])
    end
    B_ = -B;
    #JuMP.with_optimizer(Gurobi.Optimizer,OutputFlag=0),quiet=true

    m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,OutputFlag=0))


        # 2.1.Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set] ≤ THETAMAX)

        JuMP.@variable(m, GenData[g, :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[g, :Pmax][1,1])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set ; Set([i,j]) in lines_set])

        # 2.2. Constraints
        # JuMP.@constraint(m, ReferenceAngle,
        #     (δ[Nodes_set[1]] ==  0.0))

        JuMP.@constraint(m,Nodal_balance[i in Nodes_set],
                sum(pij[i,j] for j = Nodes_set if Set([i,j]) in lines_set) == sum(p[g] for g in Gen_set if System.Gen_Data[g,:bus][1,1] == i) - sum(BusData[findall(x->x==i,BusData[:bus_i]),:Pd]))
        JuMP.@constraint(m,pl[i in Nodes_set,j in Nodes_set; Set([i,j]) in lines_set],
                pij[i,j] == Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])][1],:x]))*(δ[i]-δ[j]))
        # JuMP.@constraint(m,consistency[i in Nodes_set,j in Nodes_set; Set([i,j]) in lines_set],
        #         pij[i,j] == -p[j,i])
        JuMP.@constraint(m,pl_rate[i in Nodes_set,j in Nodes_set; Set([i,j]) in lines_set],
                -Sbase*lines_rates[Set([i,j])] ≤ pij[i,j] ≤ Sbase*lines_rates[Set([i,j])] )

        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g]+GenData[g,:C0][1,1] for g in Gen_set))
        JuMP.optimize!(m)

        solution_status = raw_status(m)
        if solution_status != "Model was proven to be infeasible."
            Pg = JuMP.value.(p)
            Pij = JuMP.value.(pij)
            δ = JuMP.value.(δ)
            V = ones(N)
            Pd = System.BusData[!,3]
            Qd = System.BusData[!,4]
            System.Operating_Cost = JuMP.objective_value(m)
            System.Gen_Data[:Pg] = Pg.data
            Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1])]) : 0 for i in Nodes_set]
            df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data,
             Pg = Pg , Pd = Pd)
            System.BusData_output = df

            Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus]] for i in 1:size(System.LineData,1) ]
            Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus]] for i in 1:size(System.LineData,1) ]
            P_Loss = Pij_1 + Pij_2

            U = []
            for i in 1:N_Lines
                uu = round((abs(Pij_1[i])/(System.Line_Constraints[i,3]*System.Sbase))*100,digits = 2)
                append!(U,uu)
            end
            df2 = DataFrame(FromBus =System.LineData[!,1] ,ToBus =System.LineData[!,2]
                ,PL_1 = Pij_1,PL_2 = Pij_2,
                PLoss = P_Loss,Utilization = U)
            System.LineLoading = df2

            Line_Duals = []
            for i in 1:N_Lines
                fbus = System.LineData[i,1]
                tbus = System.LineData[i,2]
                current_dual_1 = JuMP.dual(pl_rate[fbus,tbus])
                current_dual_2 = JuMP.dual(pl_rate[tbus,fbus])
                if current_dual_1 != 0
                    Line_Duals = append!(Line_Duals,current_dual_1)
                elseif current_dual_2 != 0
                    Line_Duals = append!(Line_Duals,current_dual_2)
                else
                    Line_Duals = append!(Line_Duals,0)
                end


            end

            Bus_Duals = []
            for i in Nodes_set
                current_dual_b = JuMP.dual(Nodal_balance[i])
                Bus_Duals = append!(Bus_Duals,current_dual_b)
            end

            System.Bus_Duals = Bus_Duals
            System.Line_Duals = Line_Duals

            if print_status
                println("==================================================================================================")
                println("|          System Summary                                                                        |")
                println("==================================================================================================")
                println("Objective Function Value = ", round(System.Operating_Cost,digits=3), " USD/hr")
                println(System.BusData_output)
                println()
                println(System.LineLoading)
                println("=========================== END OF REPORT ========================================")
            end
        end

        return JuMP.value.(p).data,solution_status
end

function Solve_OTS!(System ::System_Struct,method = "DCOTS",use_dc_flag = true , print_status = true)
    #This function needs the shitty library called Juniper, that's why I decided to create a new method
    # Optimal transmission switching solver - still under development
    System_Init = deepcopy(System)
    #Solve_OPF!(System_Init,"ACOPF",false)
    G = real(System.y_bus);
    B = imag(System.y_bus);
    b = System.b_line;

    N = System.N_bus
    N_Lines = size(System.Line_Constraints,1)
    Nodes_set = 1:N;
    Gen_set = System.Gen_Data[!,:bus];

    fbus_set = transpose([Int(System.Line_Constraints[i,1]) for i in 1:N_Lines])
    tbus_set = transpose([Int(System.Line_Constraints[i,2]) for i in 1:N_Lines])
    lines_set_1 = vcat(fbus_set,tbus_set)
    lines_set_1 =[lines_set_1[:,i] for i in 1:N_Lines]
    lines_set_2 = vcat(tbus_set,fbus_set)
    lines_set_2 =[lines_set_2[:,i] for i in 1:N_Lines]
    Sbase = System.Sbase;

    BusData = System.BusData;
    GenData = System.Gen_Data;
    BranchData = System.Line_Constraints;
    S = zeros(N,N)

        for l in 1:N_Lines
            S[Int64(BranchData[l,:fbus]),Int64(BranchData[l,:tbus])] = BranchData[l,:SL_Rating]
            S[Int64(BranchData[l,:tbus]),Int64(BranchData[l,:fbus])] = BranchData[l,:SL_Rating]
        end

    if method == "ACOTS"
            if use_dc_flag
                dc_line_status_unremoved,dc_line_status_removed,_ =  Solve_OTS!(System,"DCOTS",false,false)
            else
                dc_line_status_removed = lines_set_1
                dc_line_status_unremoved = []
            end
            #M = maximum(S)^2;
            optimizer = Juniper.Optimizer
            nl_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
            mip_solver = JuMP.optimizer_with_attributes(Cbc.Optimizer,"logLevel" => 0)
            #m = Model(optimizer_with_attributes(optimizer, "nl_solver" => nl_solver,"mip_solver" => mip_solver, "branch_strategy" => :StrongPseudoCost,
            #    "strong_restart" => true, "feasibility_pump" => true,
            #    "traverse_strategy" => :DFS, "incumbent_constr" => true,
            #    "num_resolve_root_relaxation" => 25))
            m = Model(with_optimizer(optimizer, nl_solver = nl_solver,mip_solver = mip_solver, atol = 1e-4, num_resolve_root_relaxation = 20))


            # 2.1.Variables
            JuMP.@variable(m, BusData[i, :Vmin] ≤ v[i in Nodes_set] ≤ BusData[i, :Vmax])
            JuMP.@variable(m, -2*π ≤ δ[i in Nodes_set] ≤ 2*π)

            JuMP.@variable(m, GenData[(GenData[:bus].==g), :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[(GenData[:bus].==g), :Pmax][1,1])
            JuMP.@variable(m, GenData[(GenData[:bus].==g), :Qmin][1,1] ≤ q[g in Gen_set] ≤ GenData[(GenData[:bus].==g), :Qmax][1,1])

            JuMP.@variable(m, pij[i in Nodes_set,j in Nodes_set; [i,j] in lines_set_1 || [i,j] in lines_set_2])
            JuMP.@variable(m, qij[i in Nodes_set,j in Nodes_set; [i,j] in lines_set_1 || [i,j] in lines_set_2])


            if ! isempty(dc_line_status_removed)
                JuMP.@variable(m, a[i in Nodes_set,j in Nodes_set; [i,j] in dc_line_status_removed || [j,i] in dc_line_status_removed],Bin)
                JuMP.@constraint(m,a_sym[i in Nodes_set,j in Nodes_set; [i,j] in dc_line_status_removed || [j,i] in dc_line_status_removed],
                    a[i,j] == a[j,i])

                JuMP.@NLconstraint(m, p_line[i in Nodes_set,j in Nodes_set; [i,j] in dc_line_status_removed || [j,i] in dc_line_status_removed],
                      pij[i,j] ==  a[i,j]*(Sbase*(v[i]*v[j]*(G[i,j]*cos(δ[i]-δ[j])+B[i,j]*sin(δ[i]-δ[j])) -(v[i]^2)*G[i,j]) ))

                 JuMP.@NLconstraint(m, q_line[i in Nodes_set,j in Nodes_set; [i,j] in dc_line_status_removed || [j,i] in dc_line_status_removed],
                      qij[i,j] ==  a[i,j]*(Sbase*(v[i]*v[j]*(G[i,j]*sin(δ[i]-δ[j]) - B[i,j]*cos(δ[i]-δ[j])) +(v[i]^2)*(B[i,j]-b[i,j])) ));
            else
                a = ones(N,N)
            end

            #JuMP.@constraint(m,max_rem,N_Lines-sum(a[i,j] for i in Nodes_set, j = Nodes_set if ([i,j] in lines_set_1)) ≤ max_remove)
            #JuMP.@constraint(m, ReferenceAngle,
            #    (δ[1] ==  0.0))
            JuMP.@constraint(m, ReferenceAngle,
                (δ[Gen_set[1]] ==  0.0))

            if use_dc_flag
                # ACTIVE POWER THROUGH LINE N-M
                JuMP.@NLconstraint(m, p_line_unrem[i in Nodes_set,j in Nodes_set; [i,j] in dc_line_status_unremoved || [j,i] in dc_line_status_unremoved],
                    pij[i,j] ==  (Sbase*(v[i]*v[j]*(G[i,j]*cos(δ[i]-δ[j])+B[i,j]*sin(δ[i]-δ[j])) -(v[i]^2)*G[i,j]) ))

                # REACTIVE POWER THROUGH LINE N-M
                JuMP.@NLconstraint(m, q_line_unrem[i in Nodes_set,j in Nodes_set; [i,j] in dc_line_status_unremoved || [j,i] in dc_line_status_unremoved],
                     qij[i,j] ==  (Sbase*(v[i]*v[j]*(G[i,j]*sin(δ[i]-δ[j]) - B[i,j]*cos(δ[i]-δ[j])) +(v[i]^2)*(B[i,j]-b[i,j])) ));
            end


            # ACTIVE NODAL BALANCE
            JuMP.@constraint(m, Pnodal[i in Nodes_set],
                sum(pij[i,j] for j = Nodes_set if ([i,j] in lines_set_1 || [i,j] in lines_set_2)) == sum(p[g] for g in Gen_set if GenData[(GenData[:bus].==g),:bus][1,1]==i) - BusData[i,:Pd])

                # REACTIVE NODAL BALANCE
            JuMP.@constraint(m, Qnodal[i in Nodes_set],
                sum(qij[i,j] for j = Nodes_set if ([i,j] in lines_set_1 || [i,j] in lines_set_2)) == sum(q[g] for g in Gen_set if GenData[(GenData[:bus].==g),:bus][1,1]==i) - BusData[i,:Qd])

            # LINE CAPACITY

            JuMP.@NLconstraint(m, Smax_1[i in Nodes_set,j in Nodes_set; [i,j] in lines_set_1 || [i,j] in lines_set_2],
                pij[i,j]^2 + qij[i,j]^2 ≤ (Sbase^2)*(S[i,j]^2)  )

            #OBJECTIVE

            JuMP.@NLobjective(m,Min,sum(GenData[(GenData[:bus].==g),:C2][1,1]*p[g]^2+GenData[(GenData[:bus].==g),:C1][1,1]*p[g]+GenData[(GenData[:bus].==g),:C0][1,1] for g in Gen_set))

            JuMP.optimize!(m)


                println()
                Pg = JuMP.value.(p)
                Qg = JuMP.value.(q)
                Pij = JuMP.value.(pij)
                Qij = JuMP.value.(qij)
                V = JuMP.value.(v)
                δ = JuMP.value.(δ)
                Pd = System.BusData[!,3]
                Qd = System.BusData[!,4]

                System.Operating_Cost = JuMP.objective_value(m)
                Pg = [i in System.Gen_Data[!,1] ? Pg.data[findall(x -> x==i,System.Gen_Data[!,1])][1,1] : 0 for i in 1:System.N_bus]
                Qg = [i in System.Gen_Data[!,1] ? Qg.data[findall(x -> x==i,System.Gen_Data[!,1])][1,1] : 0 for i in 1:System.N_bus]
                df = DataFrame(Bus = 1:System.N_bus, V = V.data, δ = δ.data,
                 Pg = Pg ,Qg = Qg, Pd = Pd, Qd = Qd)
                System.BusData_output = df

                Pij_1 = [ Pij[System.Line_Constraints[i,1],System.Line_Constraints[i,2]] for i in 1:size(System.Line_Constraints,1) ]
                Pij_2 = [ Pij[System.Line_Constraints[i,2],System.Line_Constraints[i,1]] for i in 1:size(System.Line_Constraints,1) ]
                P_Loss = Pij_1 + Pij_2

                Qij_1 = [ Qij[System.Line_Constraints[i,1],System.Line_Constraints[i,2]] for i in 1:size(System.Line_Constraints,1) ]
                Qij_2 = [ Qij[System.Line_Constraints[i,2],System.Line_Constraints[i,1]] for i in 1:size(System.Line_Constraints,1) ]
                Q_Loss = Qij_1 + Qij_2
                df2 = DataFrame(FromBus =System.Line_Constraints[!,1] ,ToBus =System.Line_Constraints[!,2]
                    ,PL_1 = Pij_1,PL_2 = Pij_2,
                    PLoss = P_Loss,QL_1 = Qij_1
                    ,QL_2 = Qij_2,QLoss = Q_Loss)
                System.LineLoading = df2

                if print_status
                    print_report(System)
                end
                if ! isempty(dc_line_status_removed)
                    Z = JuMP.value.(a)
                else
                    Z = a
                end

                line_status = zeros(size(System.LineLoading)[1],1)

                for i in 1:size(System.LineLoading)[1]
                    fbus = System.Line_Constraints[i,1]
                    tbus = System.Line_Constraints[i,2]
                    if [fbus,tbus] in dc_line_status_unremoved || [tbus,fbus] in dc_line_status_unremoved
                        line_status[i] = 1
                    elseif [fbus,tbus] in dc_line_status_removed || [tbus,fbus] in dc_line_status_removed
                        line_status[i] = round(Z[fbus,tbus])
                    end
                end

                return line_status

            elseif method  == "DCOTS"

                M = 20*maximum([BusData[i,:Pd] for i in Nodes_set])
                B_ = -B[1:N,1:N];

                #m = JuMP.Model(JuMP.with_optimizer(Cbc.Optimizer, logLevel=1))
                m = JuMP.Model(JuMP.with_optimizer(Cbc.Optimizer))

                # 2.1.Variables
                JuMP.@variable(m, -2*π ≤ δ[i in Nodes_set] ≤ 2*π)

                JuMP.@variable(m, GenData[(GenData[:bus].==g), :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[(GenData[:bus].==g), :Pmax][1,1])

                JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set; [i,j] in lines_set_1 || [j,i] in lines_set_1])

                JuMP.@variable(m, z[i in Nodes_set,j in Nodes_set; [i,j] in lines_set_1 || [j,i] in lines_set_1],Bin)

                #JuMP.@variable(m,max_rem_var)

                # 2.2. Constraints
                JuMP.@constraint(m, ReferenceAngle,
                    (δ[Gen_set[1]] ==  0.0))


                JuMP.@constraint(m,Nodal_balance[i in Nodes_set],
                        sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set_1 || [j,i] in lines_set_1) == sum(p[g] for g in Gen_set if GenData[(GenData[:bus].==g),:bus][1,1]==i) - BusData[i,:Pd])

                JuMP.@constraint(m,balance, sum(p[g] for g in Gen_set) == sum(BusData[i,:Pd] for i in Nodes_set))

                JuMP.@constraint(m,z_sym[i in Nodes_set,j in Nodes_set; [i,j] in lines_set_1 || [j,i] in lines_set_1], z[i,j]==z[j,i])

                JuMP.@constraint(m,pl_1[i in Nodes_set,j in Nodes_set; [i,j] in lines_set_1 || [j,i] in lines_set_1],
                        -pij[i,j] + Sbase*B_[i,j]*(δ[i]-δ[j])≤ (1-z[i,j])*M )

                JuMP.@constraint(m,pl_2[i in Nodes_set,j in Nodes_set; [i,j] in lines_set_1 || [j,i] in lines_set_1],
                        -(1-z[i,j])*M ≤ -pij[i,j] + Sbase*B_[i,j]*(δ[i]-δ[j]) )



                JuMP.@constraint(m,pl_rate_1[i in Nodes_set,j in Nodes_set; [i,j] in lines_set_1 || [j,i] in lines_set_1],
                        -Sbase*S[i,j]*z[i,j]  ≤ pij[i,j])


                JuMP.@constraint(m,pl_rate_2[i in Nodes_set,j in Nodes_set; [i,j] in lines_set_1 || [j,i] in lines_set_1],
                         pij[i,j] ≤ Sbase*S[i,j]*z[i,j] )


                #JuMP.@constraint(m,max_rem,sum(pij[Int(System.Line_Constraints[i,1]),Int(System.Line_Constraints[i,2])] .==0.0 for i in 1:N_Lines) <= max_remove)
                #JuMP.@constraint(m,max_rem,max_rem_var <= max_remove)
                #JuMP.@constraint(m,max_rem_equality,max_rem_var == N_Lines-sum(z[Int(System.Line_Constraints[i,1]),Int(System.Line_Constraints[i,2])] for i in 1:N_Lines))

                JuMP.@objective(m,Min,sum(GenData[(GenData[:bus].==g),:C1][1,1]*p[g]+GenData[(GenData[:bus].==g),:C0][1,1] for g in Gen_set))

                JuMP.optimize!(m)

                Pg = JuMP.value.(p)
                Pij = JuMP.value.(pij)
                δ = JuMP.value.(δ)
                V = ones(N)
                Pd = System.BusData[!,3]
                Qd = System.BusData[!,4]
                System.Operating_Cost = JuMP.objective_value(m)
                Z = JuMP.value.(z)

                Pg = [i in System.Gen_Data[!,1] ? Pg.data[findall(x -> x==i,System.Gen_Data[!,1])][1,1] : 0 for i in 1:System.N_bus]
                df = DataFrame(Bus = 1:System.N_bus, V = V, δ = δ.data,
                 Pg = Pg , Pd = Pd)
                System.BusData_output = df

                Pij_1 = [ Pij[System.Line_Constraints[i,1],System.Line_Constraints[i,2]] for i in 1:size(System.Line_Constraints,1) ]
                Pij_2 = [ Pij[System.Line_Constraints[i,2],System.Line_Constraints[i,1]] for i in 1:size(System.Line_Constraints,1) ]
                P_Loss = Pij_1 + Pij_2

                df2 = DataFrame(FromBus =System.Line_Constraints[!,1] ,ToBus =System.Line_Constraints[!,2]
                    ,PL_1 = Pij_1,PL_2 = Pij_2,
                    PLoss = P_Loss)
                System.LineLoading = df2

                if print_status
                    println(Z)
                    println(System.BusData_output)
                    println()
                    println(System.LineLoading)
                    println(System.Operating_Cost)
                    #println(JuMP.value.(max_rem_var))
                    println(sum(Z[Int(System.Line_Constraints[i,1]),Int(System.Line_Constraints[i,2])] for i in 1:N_Lines))
                end
                k = 1
                for i in 1:size(System.LineLoading)[1]
                    if round(Z[System.Line_Constraints[i,1],System.Line_Constraints[i,2]]) == 1
                        k = k+1
                    end
                end
                m = 1
                l = 1
                unrem_num = k-1
                rem_num = N_Lines-unrem_num
                println(rem_num)
                println(N_Lines)
                println(unrem_num)
                rem_info = zeros(N_Lines,1)
                line_status_unremoved = Vector{Any}(undef,unrem_num)
                line_status_removed = Vector{Any}(undef,rem_num)
                for i in 1:size(System.LineLoading)[1]
                    if ceil(Z[System.Line_Constraints[i,1],System.Line_Constraints[i,2]]) == 1
                        line_status_unremoved[m] = [System.Line_Constraints[i,1],System.Line_Constraints[i,2]]
                        rem_info[i] = 1
                        m = m+1
                    else
                        line_status_removed[l] =[System.Line_Constraints[i,1],System.Line_Constraints[i,2]]
                        l = l+1
                    end
                end

                return line_status_unremoved,line_status_removed,rem_info
    end
end

function remove_line!(System ::System_Struct,line_status)

    for i in 1:length(line_status)
        if line_status[i] == 0
            System.LineData[i,3] = Inf
            System.LineData[i,4] = Inf
            System.LineData[i,5] = 0
        end
    end
    System.LineData = System.LineData[System.LineData[:r] .!= Inf,:]
    YBUS,bLine = Y_Bus(System.LineData,System.N_bus)
    System.y_bus = YBUS
    System.b_line = bLine
end

function get_network_matrices!(System ::System_Struct)
    LineData = System.LineData
    l = size(LineData,1)
    l = l[1]
    n = System.N_bus
    A = zeros(l,n)
    I = zeros(n,n)
    Y_sorted = zeros(l,1)*im
    r_new = zeros(l,1)*im
    x_new = zeros(l,1)*im
    b_new = zeros(l,1)
    rate_new = zeros(l,1)
    fbus_new = zeros(l,1)
    tbus_new = zeros(l,1)
    for i in 1:l
        if i ≤ n
            A[i,i] = 1;
            if i in LineData[!,1]
                j = findfirst(x->x==i,LineData[!,1])
                k = LineData[j[1],2]
            else
                j = findfirst(x->x==i,LineData[!,2])
                k = LineData[j[1],1]
            end
            A[i,k] = -1
            I[i,k] = i
            I[k,i] = i
            r_new[i] = LineData[j[1],3]
            x_new[i] = LineData[j[1],4]
            fbus_new[i] = i
            tbus_new[i] = k
            b_new[i] = LineData[j[1],5]
            rate_new[i] = LineData[j[1],6]
            z = r_new[i]+im*x_new[i]
            Y_sorted[i] = 1/z
            deleterows!(LineData,j[1])

        else
            fbus_ = LineData[1,1]
            tbus_ = LineData[1,2]
            A[i,fbus_] = 1
            A[i,tbus_] = -1
            I[fbus_,tbus_] = i
            I[tbus_,fbus_] = i
            r_new[i] = LineData[1,3]
            x_new[i] = LineData[1,4]
            fbus_new[i] = fbus_
            tbus_new[i] = tbus_
            b_new[i] = LineData[1,5]
            rate_new[i] = LineData[1,6]
            z = r_new[i]+im*x_new[i]
            Y_sorted[i] = 1/z
            deleterows!(LineData,1)
        end
    end

    LineData_new = DataFrame(fbus=Int.(vec(fbus_new)), tbus =Int.(vec(tbus_new)),r =vec(r_new) ,
        x = vec(x_new),b = vec(b_new), rate = vec(rate_new))
    System.LineData = LineData_new
    Y_pr = diagm(vec(Y_sorted))
    System.Y_bus_inc = transpose(A)*Y_pr*A
    System.Y_pr = Y_pr
    System.A = A
    System.I = I

    Y_bus_der, b_der = Y_Bus(LineData_new,n)
    System.y_bus = Y_bus_der
    System.b_line = b_der
end

function get_network_matrices_unsorted!(System ::System_Struct)
    LineData = System.LineData
    l = size(LineData,1)
    n = System.N_bus
    A_0 = zeros(l,n)
    I = zeros(n,n)
    Y_unsorted = zeros(l,1)*im
    r_new = zeros(l,1)*im
    x_new = zeros(l,1)*im
    b_new = zeros(l,1)
    rate_new = zeros(l,1)

    for i in 1:l
        fbus = LineData[i,1]
        tbus = LineData[i,2]
        r = LineData[i,3]
        x = LineData[i,4]
        z = r+1*im*x
        Y_unsorted[i] = 1/z
        I[fbus,tbus] = i
        I[tbus,fbus] = i

        if fbus == i
            A_0[i,fbus] = 1
            A_0[i,tbus] = -1

        elseif tbus == i
            A_0[i,tbus] = 1
            A_0[i,fbus] = -1

        else
            A_0[i,fbus] = 1
            A_0[i,tbus] = -1
        end


    end

    A = deepcopy(A_0)
    static_index = -1

    for i in 1:n
        j = findfirst(x->x==1,A[:,i])
        if j != nothing
            if j[1] == i
                static_index = hcat(static_index,j[1])
                A[i,:] = zeros(1,size(A,2))
            end

        end
    end
    static_index = static_index[2:end]

    non_static_index = setdiff(1:l,static_index)

    Y_pr = diagm(vec(Y_unsorted))
    System.Y_bus_inc = transpose(A_0)*Y_pr*A_0
    System.Y_pr = Y_pr
    System.A = A_0
    System.I = I
    System.static_line_index = static_index
    System.non_static_line_index = non_static_index
end

function create_extended_incidence_matrix!(System ::System_Struct)

    LineData = System.LineData
    A = System.A
    n = System.N_bus
    m , _ = size(LineData)
    # No. of duplicat buses
    N_dup = 2*m-n

    #Build D array
    D = zeros(n,1)
    for i in 1:n
        D[i] = count(x->x==-1,A[:,i])+count(x->x==1,A[:,i])-1
    end

    #Build Δ
    Δ = zeros(N_dup,3)
    k = 0
    for i in 1:n
        if D[i] != 0
            for j in 1:Int(D[i])
                k = k + 1
                Δ[k,1] = i
                Δ[k,2] = j
                Δ[k,3] = n + k
            end
        end
    end

    #Build z and ρ
    z = Vector{Any}(undef, m)
    ρ = Vector{Any}(undef, m-n)
    for i in 1:m
        temp_z = LineData[i,2]
        j = findall(x->x==LineData[i,2],Δ[:,1])
        for k in 1:length(j)
            temp_z = hcat(temp_z,Δ[j[k],3])
        end
        z[i] = temp_z

        if i>n
            temp_ρ = LineData[i,1]
            j = findall(x->x==LineData[i,1],Δ[:,1])
            for k in 1:length(j)
                temp_ρ = hcat(temp_ρ,Δ[j[k],3])
            end
            ρ[i-n] = temp_ρ
        end
    end

    A_ex = hcat(A,zeros(m,N_dup))

    System.D = D
    System.Δ = Δ
    System.z = z
    System.ρ = ρ
    System.A_ex = A_ex
end

function create_extended_matrices_unsorted!(System ::System_Struct)
    LineData = System.LineData
    A = System.A
    n = System.N_bus
    m , _ = size(LineData)
    # No. of duplicat buses
    N_dup = 2*m-n

    #Build D array
    D = zeros(n,1)
    for i in 1:n
        D[i] = count(x->x==-1,A[:,i])+count(x->x==1,A[:,i])-1
    end

    #Build Δ
    Δ = zeros(N_dup,3)
    k = 0
    for i in 1:n
        if D[i] != 0
            for j in 1:Int(D[i])
                k = k + 1
                Δ[k,1] = i
                Δ[k,2] = j
                Δ[k,3] = n + k
            end
        end
    end

    non_static_lines = System.non_static_line_index
    static_lines = System.static_line_index

    #Build z and ρ
    z = Vector{Any}(undef, m)
    ρ = Vector{Any}(undef, length(non_static_lines))

    non_static_counter = 0
    for i in 1:m

        if i in static_lines

            temp_z = LineData[i,findfirst(x->x != i,LineData[i,1:2])]
            j = findall(x->x==temp_z,Δ[:,1])
            for k in 1:length(j)
                temp_z = hcat(temp_z,Int(Δ[j[k],3]))
            end
            z[i] = temp_z

        elseif i in non_static_lines
            non_static_counter = non_static_counter + 1

            temp_z = Int(LineData[i,2])
            j = findall(x->x==LineData[i,2],Δ[:,1])
            for k in 1:length(j)
                temp_z = hcat(temp_z,Int(Δ[j[k],3]))
            end
            z[i] = temp_z

            temp_ρ = Int(LineData[i,1])
            j = findall(x->x==LineData[i,1],Δ[:,1])
            for k in 1:length(j)
                temp_ρ = hcat(temp_ρ,Int(Δ[j[k],3]))
            end
            ρ[non_static_counter] = temp_ρ

        else
            println("Error!!!")
        end
    end



    A_ex = hcat(A,zeros(m,N_dup))

    System.D = D
    System.Δ = Δ
    System.z = z
    System.ρ = ρ
    System.A_ex = A_ex
end

function system_updater(System ::System_Struct,A_ex)

    System_new = deepcopy(System)
    System_new.A_ex = A_ex
    _, N_bus_total = size(A_ex)

    #find buses in service
    Active_buses = []
    Inactive_buses= []
    for i in 1:N_bus_total
        a = findall(x->x==1,A_ex[:,i])
        b = findall(x->x==-1,A_ex[:,i])

        if  iszero(A_ex[:,i])
            append!(Inactive_buses,i)
        else
            append!(Active_buses,i)
        end
    end
    A_ex2 = A_ex[:,setdiff(1:end,Inactive_buses)]
    Ybus_new = transpose(A_ex2)*System.Y_pr*A_ex2
    System_new.Y_bus_inc = Ybus_new
    System_new.y_bus = Ybus_new
    System_new.Active_buses = Active_buses
    System_new.Inactive_buses = Inactive_buses
    N_buses_active = length(Active_buses)

    #Add new buses
    df = DataFrame(bus_i=Int64[],type=Int64[],Pd=Float64[],
        Qd=Float64[],Vmax=Float64[],Vmin=Float64[])

    Current_buses = System_new.BusData[:bus_i]

    for i in Active_buses
        if i ∉ Current_buses
            push!(df,[i -1 0 0 1.1 0.9])
        end
    end
    append!(System_new.BusData,df)

    deleted_buses = []
    for i in 1:length(System_new.Inactive_buses)
        append!(deleted_buses,findall(x->x==System_new.Inactive_buses[i],System_new.BusData[:bus_i]))
    end

    if !isempty(deleted_buses)
        for i in 1:length(deleted_buses)
            #Removing inactive buses
                deleterows!(System_new.BusData, deleted_buses[i])
        end
    end
    #Updating line connections
    N_lines,_ = size(System_new.Line_Constraints)
    for i in 1:N_lines

        fbus = findfirst(x->x==1,A_ex[i,:])
        tbus = findfirst(x->x==-1,A_ex[i,:])
        System_new.Line_Constraints[i,1] = Int64(fbus[1])
        System_new.Line_Constraints[i,2] = Int64(tbus[1])
    end
    System_new.N_bus = N_buses_active
    return System_new
end

function system_perturber(System ::System_Struct,l)

    A_new = deepcopy(System.A_ex)
    z = deepcopy(System.z)
    ρ = deepcopy(System.ρ)
    z_conn,ρ_conn = get_connection_status(System)
    if l<= System.N_bus
        z_i = z[l]
        zp = z_conn[l]
        p = length(z_i)

        j = findfirst(x->x==-1,zp)[1]
        k = Int(z_i[j])

        j_new = rand(setdiff(1:p,j))
        k_new = Int(z_i[j_new])

        A_new[l,k] = 0
        A_new[l,k_new] = -1

        zp[j] = 0
        zp[j_new] = -1
        return l,zp,A_new,1
    else
        coin = rand(1:2)
        if coin==1
            #perturb z

            z_i = z[l]
            zp = z_conn[l]
            p = length(z_i)

            j = findfirst(x->x==-1,zp)[1]
            k = Int(z_i[j])

            j_new = rand(setdiff(1:p,j))
            k_new = Int(z_i[j_new])

            A_new[l,k] = 0
            A_new[l,k_new] = -1
            zp[j] = 0
            zp[j_new] = -1
            return l,zp,A_new,coin
        else
            #perturb ρ

            ρ_p = ρ_conn[l-System.N_bus]
            pp = length(ρ_p)
            ρ_i = ρ[l-System.N_bus]
            println(ρ_i)
            jj= findfirst(x->x==1,ρ_p)[1]
            kk = Int(ρ_i[jj])
            jj_new = rand(setdiff(1:pp,jj))
            kk_new = Int(ρ_i[jj_new])
            A_new[l,kk] = 0
            A_new[l,kk_new] = 1
            ρ_p[jj] = 0
            ρ_p[jj_new] = 1
            return l,ρ_p,A_new,coin
        end
    end

end

function system_perturber_unsorted(System ::System_Struct,l)

    A_new = deepcopy(System.A_ex)
    z = deepcopy(System.z)
    ρ = deepcopy(System.ρ)
    z_conn,ρ_conn = get_connection_status_unsorted(System)
    if l in System.static_line_index
        z_i = z[l]
        zp = z_conn[l]
        p = length(z_i)

        j = findfirst(x->x==-1,zp)[1]
        k = Int(z_i[j])


        if p > 1
            j_new = rand(setdiff(1:p,j))
            k_new = Int(z_i[j_new])

            A_new[l,k] = 0
            A_new[l,k_new] = -1

            zp[j] = 0
            zp[j_new] = -1
        end
        return l,zp,A_new,1
    else
        coin = rand(1:2)
        if coin==1
            #perturb z

            z_i = z[l]
            zp = z_conn[l]
            p = length(z_i)

            j = findfirst(x->x==-1,zp)[1]
            k = Int(z_i[j])


            if p > 1
                j_new = rand(setdiff(1:p,j))
                k_new = Int(z_i[j_new])
                A_new[l,k] = 0
                A_new[l,k_new] = -1
                zp[j] = 0
                zp[j_new] = -1
            end


            return l,zp,A_new,coin
        else
            #perturb ρ
            m = findfirst(x->x==l,System.non_static_line_index)
            ρ_p = ρ_conn[m]
            pp = length(ρ_p)
            ρ_i = ρ[m]
            jj= findfirst(x->x==1,ρ_p)[1]
            kk = Int(ρ_i[jj])
            jj_new = rand(setdiff(1:pp,jj))
            kk_new = Int(ρ_i[jj_new])
            A_new[l,kk] = 0
            A_new[l,kk_new] = 1
            ρ_p[jj] = 0
            ρ_p[jj_new] = 1
            return l,ρ_p,A_new,coin
        end
    end
end

function get_connection_status(System ::System_Struct)

    z = System.z
    ρ = System.ρ
    A = System.A_ex
    z_conn = Vector{Any}(undef, length(z))
    ρ_conn = Vector{Any}(undef, length(ρ))

    for i in 1:length(z)

        z_conn[i] = A[i,Int.(z[i])]

        if i>System.N_bus
            ρ_conn[i-System.N_bus] = A[i,Int.(ρ[i-System.N_bus])]
        end
    end

    return z_conn,ρ_conn
end

function get_connection_status_unsorted(System ::System_Struct)

    z = System.z
    ρ = System.ρ
    A = System.A_ex
    z_conn = Vector{Any}(undef, length(z))
    ρ_conn = Vector{Any}(undef, length(ρ))
    k = 0
    for i in 1:length(z)

        z_conn[i] = A[i,Int.(z[i])]

        if i in System.non_static_line_index
            k = k+1
            ρ_conn[k] = A[i,Int.(ρ[k])]
        end
    end

    return z_conn,ρ_conn
end

function OBS_optimizer!(System ::System_Struct,OPF_method,itr_lim)

    get_network_matrices!(System)
    create_extended_incidence_matrix!(System)
    if Solve_OPF!(System,OPF_method,false) == "Infeasible_Problem_Detected"
        return "Infeasible Problem"
    end
    m,_ = size(System.LineData)
    best_cost = System.Operating_Cost
    best_A = System.A
    cost = Vector{Any}(undef,m)
    A_ex = Vector{Any}(undef,m)
    itr = 0
    sys_new = deepcopy(System)
    while true
        for i in 1:m
                l,s,A_new,_ = system_perturber(System,i)
                sys_new = system_updater(System,A_new)
                ss = Solve_OPF!(sys_new,OPF_method,false)

                if ss == "Infeasible_Problem_Detected"
                    cost[i] = Inf
                    A_ex[i] = deepcopy(A_new)
                else
                    cost[i] = deepcopy(sys_new.Operating_Cost)
                    A_ex[i] = deepcopy(A_new)
                end
        end

        best_P_cost = minimum(cost)
        best_index = findfirst(x->x==best_P_cost,cost)
        println(best_P_cost)
        if best_P_cost < best_cost
            best_cost = best_P_cost
            A_new = A_ex[best_index[1]]
            System = system_updater(System,A_new)
            itr = 0
            continue
        else
            itr = itr + 1
            if itr ≥ itr_lim
                println("Iteration limit")
                break
            else
                continue
            end
        end
    end
    Solve_OPF!(System,OPF_method,true)
    return System
end

function OBS_optimizer_unsorted!(System ::System_Struct,
    OPF_method,itr_lim,max_depth,
    optimizer_strategy,
    OPF_eval_ratio,
    lines_subset_size = size(System.LineData,1))

    get_network_matrices_unsorted!(System)
    create_extended_matrices_unsorted!(System)
    solution_status = Solve_OPF!(System,OPF_method,false)

    if solution_status == "Infeasible_Problem_Detected" || solution_status == "Model was proven to be infeasible."
        return "Infeasible Problem"
    end
    if optimizer_strategy == "rsp"
        msg = @sprintf "Initial Cost: %.3f" System.Operating_Cost
        println(msg)
        System = random_sequential_purterber(System,max_depth,itr_lim,lines_subset_size,OPF_method)
    elseif optimizer_strategy == "dsp"
        msg = @sprintf "Initial Cost: %.3f" System.Operating_Cost
        println(msg)
        System = determinant_search_perturber(System,max_depth,itr_lim,lines_subset_size,OPF_method,OPF_eval_ratio)
    elseif optimizer_strategy == "dsp_s"
        msg = @sprintf "Initial Cost: %.3f" System.Operating_Cost
        println(msg)
        System = determinant_search_perturber_seq(System,max_depth,itr_lim,lines_subset_size,OPF_method,OPF_eval_ratio)
    elseif optimizer_strategy == "dsp_s_a"
        msg = @sprintf "Initial Cost: %.3f" System.Operating_Cost
        println(msg)
        System = approx_determinant_search_perturber_seq(System,max_depth,itr_lim,lines_subset_size,OPF_method,OPF_eval_ratio)
    end

    Solve_OPF!(System,OPF_method,true)
    return System
end

function determinant_search_perturber(System,max_depth,itr_lim,lines_subset_size,OPF_method,OPF_eval_ratio)
    n_OPF = 0
    m,_ = size(System.LineData)
    initial_cost = System.Operating_Cost
    best_cost = System.Operating_Cost
    best_A = System.A
    cost = Vector{Any}(undef,lines_subset_size)
    ybus_det_abs = Vector{Any}(undef,lines_subset_size)
    itr = 0
    global_iteration = 0
    sys_new = deepcopy(System)
    msg_2 = ""
    msg = ""
    flag2 = 0
    while true
        ybus_det_abs = Vector{Any}(undef,lines_subset_size)
        A_ex = Vector{Any}(undef,lines_subset_size)
        global_iteration = global_iteration +1
        lines_spanned = []
        flag = 0

        println("Perturbing System...")
        for i in 1:lines_subset_size
                j = rand(setdiff(1:m,lines_spanned))
                append!(lines_spanned,j)
                l,s,A_new,_ = system_perturber_unsorted(deepcopy(System),j)
                sys_new = system_updater(System,A_new)
                ybus_det_abs[i] = abs(det(sys_new.Y_bus_inc))
                A_ex[i] = A_new
        end

        n_itr_unsolved = 0

        println("Finding Best Perturbation...")
        for counter in 1:lines_subset_size

            max_det = maximum(ybus_det_abs)
            max_det_index = findall(x->x == max_det,ybus_det_abs)
            A_max_det = A_ex[max_det_index[1]]
            sys_new = system_updater(System,A_max_det)
            println("Solving DCOPF...")
            ss = Solve_OPF!(sys_new,OPF_method,false)
            n_OPF = n_OPF +1
            if ss == "Model was proven to be infeasible." || sys_new.Operating_Cost ≥ best_cost
                deleteat!(ybus_det_abs,max_det_index[1])
                deleteat!(A_ex,max_det_index[1])
                if ss == "Model was proven to be infeasible."
                    msg1 = @sprintf "Infeasible Perturbation.."
                else
                    msg1 =  @sprintf "Non-optimal Perturbation.."
                end
                n_itr_unsolved = n_itr_unsolved+1

                if n_itr_unsolved ≥ ceil(OPF_eval_ratio*m)
                        flag2 = 1
                        println(msg1)
                        println("Maximum OPF evaluations reached...")
                        return System
                else
                    println(msg1)
                    continue
                end
            else
                System = deepcopy(sys_new)
                best_cost = System.Operating_Cost
                flag = 1
                itr = 0
                enahncement = ((initial_cost-best_cost)/initial_cost)*100
                msg_2 = @sprintf "--Depth %.0f , Obtained Cost: %.3f , Cost Reduction: %.2f %% , DCOPF Evaluations: %0.f " global_iteration best_cost enahncement n_OPF
                println(msg_2)
                println([best_cost,max_det])

                break
            end
        end



         if global_iteration ≥ max_depth

             println("Max depth reached.")
             break
         end

         if flag == 0
             if itr ≥ itr_lim
                 println("Iteration limit.")
                 break
             else
                 itr = itr + 1
                 continue
             end
         end

         if flag2 == 1
             println("A stable solution has been achieved...")
             break
         end

    end
    return System
end

function determinant_search_perturber_seq(System,max_depth,itr_lim,lines_subset_size,OPF_method,OPF_eval_ratio)
    n_OPF = 0
    m,_ = size(System.LineData)
    initial_cost = System.Operating_Cost
    best_cost = System.Operating_Cost
    best_A = System.A
    cost = Vector{Any}(undef,lines_subset_size)
    #ybus_det_abs = Vector{Any}(undef,lines_subset_size)
    itr = 0
    global_iteration = 0
    sys_new = deepcopy(System)
    msg_2 = ""
    msg = ""
    flag2 = 0
    while true
        A_ex = []
        global_iteration = global_iteration +1
        flag = 0

        println("Perturbing System...")
        for i in 1:lines_subset_size
                append!(A_ex,system_perturber_unsorted_seq(deepcopy(System),i))
        end

        ybus_det_abs = Vector{Any}(undef,length(A_ex))


        for i in 1:length(A_ex)
            ybus_det_abs[i] = abs(det(transpose(A_ex[i])*System.Y_pr*A_ex[i]))
        end

        n_itr_unsolved = 0

        println("Finding Best Perturbation...")
        for counter in 1:length(A_ex)

            max_det = maximum(ybus_det_abs)
            max_det_index = findall(x->x == max_det,ybus_det_abs)
            A_max_det = A_ex[max_det_index[1]]
            sys_new = system_updater(System,A_max_det)
            println("Solving DCOPF...")
            ss = Solve_OPF!(sys_new,OPF_method,false)
            n_OPF = n_OPF +1
            if ss == "Model was proven to be infeasible." || sys_new.Operating_Cost ≥ best_cost
                deleteat!(ybus_det_abs,max_det_index[1])
                deleteat!(A_ex,max_det_index[1])
                if ss == "Model was proven to be infeasible."
                    msg1 = @sprintf "Infeasible Perturbation.."
                else
                    msg1 =  @sprintf "Non-optimal Perturbation.."
                end
                n_itr_unsolved = n_itr_unsolved+1

                if n_itr_unsolved ≥ ceil(OPF_eval_ratio*length(A_ex))
                        flag2 = 1
                        println(msg1)
                        println("Maximum OPF evaluations reached...")
                        return System
                else
                    println(msg1)
                    continue
                end
            else
                System = deepcopy(sys_new)
                best_cost = System.Operating_Cost
                flag = 1
                itr = 0
                enahncement = ((initial_cost-best_cost)/initial_cost)*100
                msg_2 = @sprintf "--Depth %.0f , Obtained Cost: %.3f , Cost Reduction: %.2f %% , DCOPF Evaluations: %0.f " global_iteration best_cost enahncement n_OPF
                println(msg_2)
                println([best_cost,max_det])

                break
            end
        end



         if global_iteration ≥ max_depth

             println("Max depth reached.")
             break
         end

         if flag == 0
             if itr ≥ itr_lim
                 println("Iteration limit.")
                 break
             else
                 itr = itr + 1
                 continue
             end
         end

         if flag2 == 1
             println("A stable solution has been achieved...")
             break
         end

    end
    return System
end

function approx_determinant_search_perturber_seq(System,max_depth,itr_lim,lines_subset_size,OPF_method,OPF_eval_ratio)
    n_OPF = 0
    m,_ = size(System.LineData)
    initial_cost = System.Operating_Cost
    best_cost = System.Operating_Cost
    best_A = System.A
    cost = Vector{Any}(undef,lines_subset_size)
    #ybus_det_abs = Vector{Any}(undef,lines_subset_size)
    itr = 0
    global_iteration = 0
    sys_new = deepcopy(System)
    msg_2 = ""
    msg = ""
    flag2 = 0
    while true
        A_ex = []
        global_iteration = global_iteration +1
        flag = 0

        println("Perturbing System...")
        for i in 1:lines_subset_size
                append!(A_ex,system_perturber_unsorted_seq(deepcopy(System),i))
        end

        ybus_det_abs = Vector{Any}(undef,length(A_ex))
        n_ex = length(A_ex)


        for i in 1:length(A_ex)
            temp = diag(transpose(A_ex[i])*A_ex[i])
            temp2 = [temp[i] for i in 1:length(temp) if temp[i] != 0]
            ybus_det_abs[i] = prod(temp2)
        end

        n_itr_unsolved = 0

        println("Finding Best Perturbation...")
        for counter in 1:length(A_ex)

            max_det = minimum(ybus_det_abs)
            println(max_det)

            max_det_index = findall(x->x == max_det,ybus_det_abs)
            A_max_det = A_ex[max_det_index[1]]
            sys_new = system_updater(System,A_max_det)
            println("Solving DCOPF...")
            ss = Solve_OPF!(sys_new,OPF_method,false)
            n_OPF = n_OPF +1
            if ss == "Model was proven to be infeasible." || sys_new.Operating_Cost ≥ best_cost
                if OPF_eval_ratio < 1
                    deleteat!(ybus_det_abs,max_det_index[1])
                    deleteat!(A_ex,max_det_index[1])
                else
                    deleteat!(ybus_det_abs,max_det_index)
                    deleteat!(A_ex,max_det_index)
                end

                if ss == "Model was proven to be infeasible."
                    msg1 = @sprintf "Infeasible Perturbation.."
                else
                    msg1 =  @sprintf "Non-optimal Perturbation. Cost: %.3f .." sys_new.Operating_Cost
                end
                n_itr_unsolved = n_itr_unsolved+1

                if n_itr_unsolved ≥ ceil(OPF_eval_ratio*n_ex)
                        flag2 = 1
                        println(msg1)
                        println("Maximum OPF evaluations reached...")
                        return System
                else
                    println(msg1)
                    if isempty(ybus_det_abs)
                        println("All Possible Perturbations Have Been Traversed!")
                        println("Aborting Optimizer...")
                        flag2 = 1
                        break
                    end
                    continue
                end
            else
                System = deepcopy(sys_new)
                best_cost = System.Operating_Cost
                flag = 1
                itr = 0
                enahncement = ((initial_cost-best_cost)/initial_cost)*100
                msg_2 = @sprintf "--Depth %.0f , Obtained Cost: %.3f , Cost Reduction: %.2f %% , DCOPF Evaluations: %0.f " global_iteration best_cost enahncement n_OPF
                println(msg_2)
                println([best_cost,max_det])

                break
            end

        end



         if global_iteration ≥ max_depth

             println("Max depth reached.")
             break
         end

         if flag == 0
             if itr ≥ itr_lim
                 println("Iteration limit.")
                 break
             else
                 itr = itr + 1
                 continue
             end
         end

         if flag2 == 1
             println("A stable solution has been achieved...")
             break
         end

    end
    return System
end


function system_perturber_unsorted_seq(System ::System_Struct,l)

    A_ex = deepcopy(System.A_ex)
    z = deepcopy(System.z)
    ρ = deepcopy(System.ρ)
    z_conn,ρ_conn = get_connection_status_unsorted(System)



    if l in System.static_line_index
        z_i = z[l]
        zp = z_conn[l]
        p = length(z_i)

        A_new_vec = Vector{Any}(undef,p)

        j = findfirst(x->x==-1,zp)[1]
        k = Int(z_i[j])

        for i in 1:p
            k_new = Int(z_i[i])
            A_new = deepcopy(A_ex)
            A_new[l,k] = 0
            A_new[l,k_new] = -1
            A_new_vec[i] = A_new
        end

        return A_new_vec
    else

        z_i = z[l]
        zp = z_conn[l]
        p = length(z[l])
        m = findfirst(x->x==l,System.non_static_line_index)
        ρ_i = ρ[m]
        ρ_p = ρ_conn[m]
        pp = length(ρ_p)
        A_new_vec = Vector{Any}(undef,p+pp)
        n_ex_mat = p+pp

        j_z = findfirst(x->x==-1,zp)[1]
        k_z = Int(z_i[j_z])

        for i in 1:n_ex_mat
            if i ≤ p
                A_new = deepcopy(A_ex)
                k_new = Int(z_i[i])
                A_new[l,k_z] = 0
                A_new[l,k_new] = -1
                A_new_vec[i] = A_new

            else
                A_new = deepcopy(A_ex)
                jj= findfirst(x->x==1,ρ_p)[1]
                kk = Int(ρ_i[jj])
                jj_new = i-p
                kk_new = Int(ρ_i[jj_new])
                A_new[l,kk] = 0
                A_new[l,kk_new] = 1
                A_new_vec[i] = A_new

            end
        end

            return A_new_vec
    end

end

function random_sequential_purterber(System,max_depth,itr_lim,lines_subset_size,OPF_method)

    m,_ = size(System.LineData)
    initial_cost = System.Operating_Cost
    best_cost = System.Operating_Cost
    best_A = System.A
    cost = Vector{Any}(undef,lines_subset_size)
    A_ex = Vector{Any}(undef,lines_subset_size)
    itr = 0
    global_iteration = 0
    sys_new = deepcopy(System)
    msg_2 = ""
    msg = ""
    while true
        global_iteration = global_iteration +1
        lines_spanned = []
        for i in 1:lines_subset_size
                j = rand(setdiff(1:m,lines_spanned))
                append!(lines_spanned,j)
                l,s,A_new,_ = system_perturber_unsorted(System,j)
                sys_new = system_updater(System,A_new)
                ss = Solve_OPF!(sys_new,OPF_method,false)

                if ss == "Model was proven to be infeasible." || ss == "Infeasible_Problem_Detected"
                    cost[i] = Inf
                    A_ex[i] = deepcopy(A_new)
                else
                    cost[i] = deepcopy(sys_new.Operating_Cost)
                    A_ex[i] = deepcopy(A_new)
                end
                msg_2 = @sprintf "--Depth %.0f , Perturbed line %.0f , Obtained Cost: %.3f" global_iteration j cost[i]
                println(msg_2)
        end

        best_P_cost = minimum(cost)
        best_index = findfirst(x->x==best_P_cost,cost)
        enahncement = ((initial_cost-best_P_cost)/initial_cost)*100
         msg = @sprintf "+Global iteration number: %.0f. Best Cost: %.3f. Cost Reduction: %.2f %%" global_iteration best_P_cost enahncement
         println(msg)

         if global_iteration ≥ max_depth
             if best_P_cost < best_cost
                 best_cost = best_P_cost
                 A_new = A_ex[best_index[1]]
                 System = system_updater(System,A_new)

             end

             println("Max depth reached.")
             break
         end
        if best_P_cost < best_cost
            best_cost = best_P_cost
            A_new = A_ex[best_index[1]]
            System = system_updater(System,A_new)
            itr = 0
            continue
        else
            itr = itr + 1
            if itr ≥ itr_lim
                println("Iteration limit.")
                break
            else
                continue
            end
        end
    end
    return System
end
