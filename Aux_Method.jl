include("SystemAssembly.jl")

function convert_to_aux!(System ::System_Struct,BusID,NumBusBars=2)

    System_converted = deepcopy(System)


    substation_selected = System_converted.SubStations[string("SubStation_",BusID)]
    substation_selected.Num_BusBars = NumBusBars
    substation_selected.flag = true
    Num_Lines = length(substation_selected.SubStation_NonAuxLines[:,1])
    if sum(substation_selected.SubStation_Loads[:,:Pd]) == 0
        Num_Loads = 0
    else
        Num_Loads = length(substation_selected.SubStation_Loads[:,1])
    end
    Num_Gens = length(substation_selected.SubStation_Gen[:,1])

    N_SubStation_Total = Num_Lines + Num_Loads + Num_Gens + NumBusBars
    N_System_Current = System_converted.N_bus
    N_System_New = N_System_Current + N_SubStation_Total - 1
    System_converted.N_bus = N_System_New
    New_Buses = collect(N_System_Current+1:N_System_New)
    Bus_bars = vcat(BusID,New_Buses[1:NumBusBars-1])
    Load_aux_Buses = New_Buses[NumBusBars:Num_Loads+NumBusBars-1]
    Gen_aux_Buses = New_Buses[NumBusBars+Num_Loads:Num_Loads+Num_Gens+NumBusBars-1]
    Conn_aux_Buses = New_Buses[Num_Gens+Num_Loads+NumBusBars:Num_Loads+Num_Gens+Num_Lines+NumBusBars-1]

    substation_selected.SubStation_Buses = vcat(BusID,New_Buses)
    substation_selected.Aux_Buses = vcat(Load_aux_Buses,Gen_aux_Buses,Conn_aux_Buses)
    substation_selected.Non_Aux_Buses = Bus_bars

    #Clearing root bus
    System_converted.BusData[BusID,:Pd] = 0.0
    System_converted.BusData[BusID,:Qd] = 0.0

    #Adding additional auxilliary buses
    for i in 1:NumBusBars-1
        push!(System_converted.BusData,[Bus_bars[i+1] 2 0.0 0.0 1.1 0.9 0.0])
    end

    #Adding Load Buses
    for i in 1:Num_Loads
        push!(System_converted.BusData,[
            Load_aux_Buses[i] 3 substation_selected.SubStation_Loads[i,:Pd] substation_selected.SubStation_Loads[i,:Qd] 1.1 0.9 1.0
            ])
    end

    #Adding Gen Buses
    Gen_ind = findall(x->x==BusID,System_converted.Gen_Data[:,1])
    for i in 1:Num_Gens
        push!(System_converted.BusData,[Gen_aux_Buses[i] 1 0.0 0.0 1.1 0.9 1.0])
        System_converted.Gen_Data[Gen_ind[i],:bus] = Gen_aux_Buses[i]
    end


    #Adding Connector Buses
    for i in 1:Num_Lines
        push!(System_converted.BusData,[Conn_aux_Buses[i] 2 0.0 0.0 1.1 0.9 1.0])
    end

    #------------------------------------------------------------#
    #Adding and modifing lines
    #------------------------------------------------------------#

    #Connecting Lines to Connector buses
    LineIDs = substation_selected.SubStation_NonAuxLines[:ID]
    k = 1
    for i in 1:length(System_converted.LineData[:,1])
        if System_converted.LineData[i,:fbus] == BusID
            System_converted.LineData[i,:fbus] = Conn_aux_Buses[k]
            k = k+1
        elseif System_converted.LineData[i,:tbus] == BusID
            System_converted.LineData[i,:tbus] = Conn_aux_Buses[k]
            k = k+1
        end
    end

    #Adding auxilliary Lines
    N_aux_lines = NumBusBars*(Num_Lines+Num_Gens+Num_Loads)
    fbus_array = Int16.(BusID*ones(Num_Lines+Num_Gens+Num_Loads,))
    for i in 1:NumBusBars-1
        fbus_array=vcat(fbus_array,Int16.(Bus_bars[i+1]*ones(Num_Lines+Num_Gens+Num_Loads,)))
    end
    tbus_array = vcat(Conn_aux_Buses,Load_aux_Buses,Gen_aux_Buses)
    for i in 1:NumBusBars-1
        tbus_array=vcat(tbus_array,tbus_array)
    end

    N_Lines = length(System_converted.LineData[:,1])
    substation_selected.SubStation_AuxLines = DataFrame(fbus = Int64[],tbus = Int64[],r =Float64[],x=Float64[],b = Float64[],rate = Float64[],ID = Int64[],is_aux=Float64[])
    for i in 1:N_aux_lines
        push!(System_converted.LineData,[
        #fbus tbus r x b rate ID is_aux
        fbus_array[i] tbus_array[i] 0 10e-4 0 50 N_Lines+i 1.0
        ])

        push!(substation_selected.SubStation_AuxLines,[
        #fbus tbus r x b rate ID is_aux
        fbus_array[i] tbus_array[i] 0 10e-4 0 50 N_Lines+i 1.0
        ])
    end

    Lines_at_bus = Dict()
    for i in substation_selected.SubStation_Buses
        ind = vcat(findall(x->x==i,System_converted.LineData[:fbus]),findall(x->x==i,System_converted.LineData[:tbus]))
        Lines_at_bus[i] = [System_converted.LineData[j,:ID] for j in ind]
    end
    substation_selected.Lines_at_bus = Lines_at_bus
    System_converted.SubStations[string("SubStation_",BusID)] = substation_selected
    System_converted.y_bus,System_converted.b_line = Y_Bus(System_converted.LineData,System_converted.N_bus)
    #System = deepcopy(System_converted)
    System.LineData = System_converted.LineData
    System.BusData = System_converted.BusData
    System.N_bus = System_converted.N_bus
    System.SubStations = System_converted.SubStations
    System.Gen_Data = System_converted.Gen_Data
    System.y_bus = System_converted.y_bus
    System.b_line = System_converted.b_line
    System.Line_Constraints = DataFrame(fbus = System.LineData[!,:fbus],tbus = System.LineData[!,:tbus],SL_Rating = System.LineData[!,:rate])
end

function convert_to_aux_coupled!(System ::System_Struct,BusID,NumBusBars=2)

    System_converted = deepcopy(System)


    substation_selected = System_converted.SubStations[string("SubStation_",BusID)]
    substation_selected.Num_BusBars = NumBusBars
    substation_selected.flag = true
    Num_Lines = length(substation_selected.SubStation_NonAuxLines[:,1])
    if sum(substation_selected.SubStation_Loads[:,:Pd]) == 0
        Num_Loads = 0
    else
        Num_Loads = length(substation_selected.SubStation_Loads[:,1])
    end
    Num_Gens = length(substation_selected.SubStation_Gen[:,1])

    N_SubStation_Total = Num_Lines + Num_Loads + Num_Gens + NumBusBars
    N_System_Current = System_converted.N_bus
    N_System_New = N_System_Current + N_SubStation_Total - 1
    System_converted.N_bus = N_System_New
    New_Buses = collect(N_System_Current+1:N_System_New)
    Bus_bars = vcat(BusID,New_Buses[1:NumBusBars-1])
    Load_aux_Buses = New_Buses[NumBusBars:Num_Loads+NumBusBars-1]
    Gen_aux_Buses = New_Buses[NumBusBars+Num_Loads:Num_Loads+Num_Gens+NumBusBars-1]
    Conn_aux_Buses = New_Buses[Num_Gens+Num_Loads+NumBusBars:Num_Loads+Num_Gens+Num_Lines+NumBusBars-1]

    substation_selected.SubStation_Buses = vcat(BusID,New_Buses)
    substation_selected.Aux_Buses = vcat(Load_aux_Buses,Gen_aux_Buses,Conn_aux_Buses)
    substation_selected.Non_Aux_Buses = Bus_bars

    #Clearing root bus
    System_converted.BusData[BusID,:Pd] = 0.0
    System_converted.BusData[BusID,:Qd] = 0.0

    #Adding additional auxilliary buses
    for i in 1:NumBusBars-1
        push!(System_converted.BusData,[Bus_bars[i+1] 2 0.0 0.0 1.1 0.9 0.0])
    end

    #Adding Load Buses
    for i in 1:Num_Loads
        push!(System_converted.BusData,[
            Load_aux_Buses[i] 3 substation_selected.SubStation_Loads[i,:Pd] substation_selected.SubStation_Loads[i,:Qd] 1.1 0.9 1.0
            ])
    end

    #Adding Gen Buses
    Gen_ind = findall(x->x==BusID,System_converted.Gen_Data[:,1])
    for i in 1:Num_Gens
        push!(System_converted.BusData,[Gen_aux_Buses[i] 1 0.0 0.0 1.1 0.9 1.0])
        System_converted.Gen_Data[Gen_ind[i],:bus] = Gen_aux_Buses[i]
    end


    #Adding Connector Buses
    for i in 1:Num_Lines
        push!(System_converted.BusData,[Conn_aux_Buses[i] 2 0.0 0.0 1.1 0.9 1.0])
    end

    #------------------------------------------------------------#
    #Adding and modifing lines
    #------------------------------------------------------------#

    #Connecting Lines to Connector buses
    LineIDs = substation_selected.SubStation_NonAuxLines[:ID]
    k = 1
    for i in 1:length(System_converted.LineData[:,1])
        if System_converted.LineData[i,:fbus] == BusID
            System_converted.LineData[i,:fbus] = Conn_aux_Buses[k]
            k = k+1
        elseif System_converted.LineData[i,:tbus] == BusID
            System_converted.LineData[i,:tbus] = Conn_aux_Buses[k]
            k = k+1
        end
    end

    #Adding auxilliary Lines
    N_aux_lines = NumBusBars*(Num_Lines+Num_Gens+Num_Loads)
    fbus_array = Int16.(BusID*ones(Num_Lines+Num_Gens+Num_Loads,))
    for i in 1:NumBusBars-1
        fbus_array=vcat(fbus_array,Int16.(Bus_bars[i+1]*ones(Num_Lines+Num_Gens+Num_Loads,)))
    end

    tbus_array = vcat(Conn_aux_Buses,Load_aux_Buses,Gen_aux_Buses)
    for i in 1:NumBusBars-1
        tbus_array=vcat(tbus_array,tbus_array)
    end

    N_Lines = length(System_converted.LineData[:,1])
    substation_selected.SubStation_AuxLines = DataFrame(fbus = Int64[],tbus = Int64[],r =Float64[],x=Float64[],b = Float64[],rate = Float64[],ID = Int64[],is_aux=Float64[])
    for i in 1:N_aux_lines
        push!(System_converted.LineData,[
        #fbus tbus r x b rate ID is_aux
        fbus_array[i] tbus_array[i] 0 10e-4 0 50 N_Lines+i 1.0
        ])

        push!(substation_selected.SubStation_AuxLines,[
        #fbus tbus r x b rate ID is_aux
        fbus_array[i] tbus_array[i] 0 10e-4 0 50 N_Lines+i 1.0
        ])
    end
    push!(System_converted.LineData,[
    #fbus tbus r x b rate ID is_aux
    Bus_bars[1] Bus_bars[2] 0 10e-4 0 50 N_Lines+N_aux_lines+1 2.0
    ])

    push!(substation_selected.SubStation_AuxLines,[
    #fbus tbus r x b rate ID is_aux
    Bus_bars[2] Bus_bars[2] 0 10e-4 0 50 N_Lines+N_aux_lines+1 2.0
    ])

    Lines_at_bus = Dict()
    for i in substation_selected.SubStation_Buses
        ind = vcat(findall(x->x==i,System_converted.LineData[:fbus]),findall(x->x==i,System_converted.LineData[:tbus]))
        Lines_at_bus[i] = [System_converted.LineData[j,:ID] for j in ind]
    end
    substation_selected.Lines_at_bus = Lines_at_bus
    System_converted.SubStations[string("SubStation_",BusID)] = substation_selected
    System_converted.y_bus,System_converted.b_line = Y_Bus(System_converted.LineData,System_converted.N_bus)
    #System = deepcopy(System_converted)
    System.LineData = System_converted.LineData
    System.BusData = System_converted.BusData
    System.N_bus = System_converted.N_bus
    System.SubStations = System_converted.SubStations
    System.Gen_Data = System_converted.Gen_Data
    System.y_bus = System_converted.y_bus
    System.b_line = System_converted.b_line
    System.Line_Constraints = DataFrame(fbus = System.LineData[!,:fbus],tbus = System.LineData[!,:tbus],SL_Rating = System.LineData[!,:rate])
end

function Solve_OBS_aux(System ::System_Struct,timelimit = Inf)
    THETAMAX = 0.6
    M = 10000
    M2 = 2 * THETAMAX
    G = real(System.y_bus);
    B = imag(System.y_bus);
    b = System.b_line;

    N = System.N_bus
    N_Lines = size(System.Line_Constraints,1)
    Nodes_set = System.BusData[:,:bus_i];
    Gen_set = 1:length(System.Gen_Data[:,:bus]);

    Sbase = System.Sbase;

    BusData = System.BusData;
    GenData = System.Gen_Data;
    BranchData = System.Line_Constraints;

    df_lines = System.LineData[:,[:fbus,:tbus]]
    lines_set = Set.(convert.(Array, row for row in eachrow(df_lines)))

    df_aux = System.LineData[findall(x->x==1.0,System.LineData[:is_aux]),[:fbus,:tbus]]
    aux_line_set = Set.(convert.(Array, row for row in eachrow(df_aux)))

    df_non_aux = System.LineData[findall(x->x==0.0,System.LineData[:is_aux]),[:fbus,:tbus]]
    non_aux_line_set = Set.(convert.(Array, row for row in eachrow(df_non_aux)))

    Nodes_at_Lines = Dict()
    for i in System.LineData[:ID]
        Nodes_at_Lines[i] = Set([System.LineData[i,:fbus],System.LineData[i,:tbus]])
    end

    # B_ = -B[1:N,1:N];

    m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,TimeLimit=timelimit))

    # 2.1.Variables
    JuMP.@variable(m, -THETAMAX≤ δ[i in Nodes_set] ≤ THETAMAX)

    JuMP.@variable(m, GenData[g, :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[g, :Pmax][1,1])

    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set])

    JuMP.@variable(m, z[l=System.LineData[:ID];System.LineData[l,:is_aux]==1.0], Bin)

    #Objective Function
    JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g]+GenData[g,:C0][1,1] for g in Gen_set) )

    #Constraints

    #Current law applied on auxilliary buses
    # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
    #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

    # JuMP.@constraint(m, ReferenceAngle,
    #     (δ[1] ==  0.0))
    #Current law applied on all buses
    JuMP.@constraint(m, Nodal_Balance[i = Nodes_set],
        sum(pij[i,j] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

    #Power flows in non auxilliary lines
    JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set ; Set([i,j]) in non_aux_line_set ],
        pij[i,j] == Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i]-δ[j]))

    JuMP.@constraint(m,total_balance, sum(p) == sum(System.BusData[:,:Pd]))
    #Switching constraint
    JuMP.@constraint(m,connectors[i = Nodes_set; BusData[i,:is_aux] == 1.0],
        sum(z[l] for l in vcat(findall(x->x==i,System.LineData[:,:fbus]),findall(x->x==i,System.LineData[:,:tbus])) if System.LineData[l,:is_aux] == 1.0) == 1)

    #Theta Constraints
    JuMP.@constraint(m,theta1[l = System.LineData[:,:ID] ;System.LineData[l,:is_aux] == 1.0],
        δ[System.LineData[l,:fbus]] <= δ[System.LineData[l,:tbus]] + (1-z[l])*M2)

    JuMP.@constraint(m,theta2[l = System.LineData[:,:ID] ;System.LineData[l,:is_aux] == 1.0],
        δ[System.LineData[l,:fbus]] + (1-z[l])*M2 >= δ[System.LineData[l,:tbus]])

    #Capacity Constraints

    JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID]],
        pij[System.LineData[i,:fbus],System.LineData[i,:tbus]] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus]])

    JuMP.@constraint(m,Line_Capacity_aux1_1[i = System.LineData[:,:ID] ;System.LineData[i,:is_aux] == 1.0],
        -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus]])

    JuMP.@constraint(m,Line_Capacity_aux1_2[i = System.LineData[:,:ID] ;System.LineData[i,:is_aux] == 1.0],
         pij[System.LineData[i,:fbus],System.LineData[i,:tbus]] <= System.LineData[i,:rate]*Sbase*z[i])

    JuMP.@constraint(m,Line_Capacity_aux2_1[i = System.LineData[:,:ID] ;System.LineData[i,:is_aux] == 1.0],
        -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus]])

    JuMP.@constraint(m,Line_Capacity_aux2[i = System.LineData[:,:ID] ;System.LineData[i,:is_aux] == 1.0],
        pij[System.LineData[i,:tbus],System.LineData[i,:fbus]] <= System.LineData[i,:rate]*Sbase*z[i])

    JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID] ;System.LineData[i,:is_aux] == 0.0],
        -System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus]] <= System.LineData[i,:rate]*Sbase)

    JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ;System.LineData[i,:is_aux] == 0.0],
        -System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus]] <= System.LineData[i,:rate]*Sbase)

    optimize!(m)

    p_ij = JuMP.value.(pij)
    pg = JuMP.value.(p)
    z = JuMP.value.(z)
    δ_sol = JuMP.value.(δ)
    obj = JuMP.objective_value(m)

    switching_status = zeros(length(System.LineData[System.LineData[:is_aux].==1.0,:ID]),2)
    switching_status[:,1] = System.LineData[System.LineData[:is_aux].==1.0,:ID]
    switching_status[:,2] = [z[i] for i in System.LineData[System.LineData[:is_aux].==1.0,:ID]]
    #========================#
    Pg = JuMP.value.(p)
    Pij = JuMP.value.(pij)
    δ = JuMP.value.(δ)
    V = ones(N)
    Pd = System.BusData[!,3]
    Qd = System.BusData[!,4]
    System.Operating_Cost = JuMP.objective_value(m)
    System.Gen_Data[:Pg] = Pg.data
    Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1])]) : 0 for i in Nodes_set]
    df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
     Pg = Pg , Pd = Pd)
    System.BusData_output = df

    Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
    Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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
    System.LineLoading[:ID] = deepcopy(System.LineData[:ID])
    #========================#
    reduce_System!(System,switching_status)
    removed_lines = sort(collect(setdiff(Set(System.LineLoading[:ID]),System.LineData[:ID])))
    deleterows!(System.LineLoading,removed_lines)
    System.LineLoading[:FromBus] = deepcopy(System.LineData[:fbus])
    System.LineLoading[:ToBus] = deepcopy(System.LineData[:tbus])
    return obj,System,z.data,pg.data
end

function reduce_System!(System ::System_Struct,switching_status)
    aux_lines_ind = switching_status[:,1]
    aux_lines_stat = switching_status[:,2]

    #Removing switched off lines from line data
    for i in 1:length(aux_lines_stat)
        if aux_lines_stat[i] == 0
            deleterows!(System.LineData,findall(x->x==aux_lines_ind[i],System.LineData[:,:ID]))
        end
    end


    #Remove unnecessary buses
    unconnected_buses = []
    for i in 1:length(System.SubStations)

        substation_selected = System.SubStations[string("SubStation_",i)]
        Lines_at_bus = Dict()
        for k in substation_selected.SubStation_Buses
            ind = vcat(findall(x->x==k,System.LineData[:fbus]),findall(x->x==k,System.LineData[:tbus]))
            Lines_at_bus[k] = [System.LineData[j,:ID] for j in ind]
        end

        System.SubStations[string("SubStation_",i)].Lines_at_bus = Lines_at_bus
        for j in System.SubStations[string("SubStation_",i)].SubStation_Buses
            if length(System.SubStations[string("SubStation_",i)].Lines_at_bus[j]) == 0
                push!(unconnected_buses,j)
            end
        end
    end

    for i in 1:length(unconnected_buses)
        deleterows!(System.BusData,findall(x->x==unconnected_buses[i],System.BusData[:,:bus_i]))
        if !isempty(System.BusData_output)
            deleterows!(System.BusData_output,findall(x->x==unconnected_buses[i],System.BusData_output[:,:Bus]))
        end
    end

    Nodes_at_Lines = Dict()
    for i in 1:length(System.LineData[:ID])
        Nodes_at_Lines[System.LineData[i,:ID]] = Set([System.LineData[i,:fbus],System.LineData[i,:tbus]])
    end
    #merging buses
    for substation in System.SubStations
        non_aux_buses = substation[2].Non_Aux_Buses

        for root_bus in non_aux_buses
            bus_list = System.LineData[System.LineData[:fbus].==root_bus,:tbus]
            if length(bus_list) != 0
                merge_buses!(System,root_bus,bus_list)
            end
        end
    end


    #updating System variables
    System.N_bus = length(System.BusData[:,1])
    rename_buses!(System)
    sort!(System.BusData,(:bus_i))
    sort!(System.BusData_output,(:Bus))
end

function reduce_System_coupled!(System ::System_Struct,switching_status,switching_status_coupler)
    aux_lines_ind = switching_status[:,1]
    aux_lines_stat = switching_status[:,2]

    coupler_lines_ind = switching_status_coupler[:,1]
    coupler_lines_stat = switching_status_coupler[:,2]

    #Removing switched off lines from line data
    for i in 1:length(aux_lines_stat)
        if aux_lines_stat[i] == 0
            deleterows!(System.LineData,findall(x->x==aux_lines_ind[i],System.LineData[:,:ID]))
        end
    end

    for i in 1:length(coupler_lines_stat)
        if coupler_lines_stat[i] == 0
            deleterows!(System.LineData,findall(x->x==coupler_lines_ind[i],System.LineData[:,:ID]))
        end
    end


    #Remove unnecessary buses
    unconnected_buses = []
    for i in 1:length(System.SubStations)

        substation_selected = System.SubStations[string("SubStation_",i)]
        Lines_at_bus = Dict()
        for k in substation_selected.SubStation_Buses
            ind = vcat(findall(x->x==k,System.LineData[:fbus]),findall(x->x==k,System.LineData[:tbus]))
            Lines_at_bus[k] = [System.LineData[j,:ID] for j in ind]
        end

        System.SubStations[string("SubStation_",i)].Lines_at_bus = Lines_at_bus
        for j in System.SubStations[string("SubStation_",i)].SubStation_Buses
            if length(System.SubStations[string("SubStation_",i)].Lines_at_bus[j]) == 0
                push!(unconnected_buses,j)
            end
        end
    end

    for i in 1:length(unconnected_buses)
        deleterows!(System.BusData,findall(x->x==unconnected_buses[i],System.BusData[:,:bus_i]))
        if !isempty(System.BusData_output)
            deleterows!(System.BusData_output,findall(x->x==unconnected_buses[i],System.BusData_output[:,:Bus]))
        end
    end

    Nodes_at_Lines = Dict()
    for i in 1:length(System.LineData[:ID])
        Nodes_at_Lines[System.LineData[i,:ID]] = Set([System.LineData[i,:fbus],System.LineData[i,:tbus]])
    end
    #merging buses
    for substation in System.SubStations
        non_aux_buses = substation[2].Non_Aux_Buses
        for root_bus in non_aux_buses
            bus_list = vcat(System.LineData[System.LineData[:fbus].==root_bus,:tbus],System.LineData[System.LineData[:tbus].==root_bus,:fbus])
            for bus in bus_list
                if bus in non_aux_buses
                    filter!(x->x≠ bus,bus_list)
                end
            end
            if length(bus_list) != 0
                merge_buses_coupled!(System,root_bus,bus_list,non_aux_buses)
            end
        end
    end


    #updating System variables
    System.N_bus = length(System.BusData[:,1])
    rename_buses!(System)
    sort!(System.BusData,(:bus_i))
    sort!(System.BusData_output,(:Bus))
end

function merge_buses!(System ::System_Struct,root_bus,bus_list)

    #updating generation data
    for i in bus_list
        System.Gen_Data[findall(x->x==i,System.Gen_Data[:bus]),:bus] = root_bus
    end

    #updating bus data

    rows = []
    for i in bus_list
        push!(rows,findall(x->x==i,System.BusData[:bus_i]))
    end
    System.BusData[findall(x->x==root_bus,System.BusData[:,:bus_i]),:Pd] = sum(System.BusData[i,:Pd] for i in rows)
    System.BusData[findall(x->x==root_bus,System.BusData[:,:bus_i]),:Qd] = sum(System.BusData[i,:Qd] for i in rows)
    System.BusData_output[findall(x->x==root_bus,System.BusData_output[:,:Bus]),:Pd] = sum(System.BusData_output[i,:Pd] for i in rows)
    System.BusData_output[findall(x->x==root_bus,System.BusData_output[:,:Bus]),:Pg] = sum(System.BusData_output[i,:Pg] for i in rows)
    for i in bus_list
        deleterows!(System.BusData,findall(x->x==i,System.BusData[:,:bus_i]))
        deleterows!(System.BusData_output,findall(x->x==i,System.BusData_output[:,:Bus]))
    end

    #updating line data
    deleterows!(System.LineData,findall(x->x==root_bus,System.LineData[:,:fbus]))
    for i in 1:length(System.LineData[:,1])
        if System.LineData[i,:fbus] in bus_list
            System.LineData[i,:fbus] = root_bus
        elseif System.LineData[i,:tbus] in bus_list
            System.LineData[i,:tbus] = root_bus
        end
    end
end

function merge_buses_coupled!(System ::System_Struct,root_bus,bus_list,non_aux_buses)

    #updating generation data
    for i in bus_list
        System.Gen_Data[findall(x->x==i,System.Gen_Data[:bus]),:bus] = root_bus
    end

    #updating bus data

    rows = []
    for i in bus_list
        push!(rows,findall(x->x==i,System.BusData[:bus_i]))
    end
    System.BusData[findall(x->x==root_bus,System.BusData[:,:bus_i]),:Pd] = sum(System.BusData[i,:Pd] for i in rows)
    System.BusData[findall(x->x==root_bus,System.BusData[:,:bus_i]),:Qd] = sum(System.BusData[i,:Qd] for i in rows)
    System.BusData_output[findall(x->x==root_bus,System.BusData_output[:,:Bus]),:Pd] = sum(System.BusData_output[i,:Pd] for i in rows)
    System.BusData_output[findall(x->x==root_bus,System.BusData_output[:,:Bus]),:Pg] = sum(System.BusData_output[i,:Pg] for i in rows)
    for i in bus_list
        deleterows!(System.BusData,findall(x->x==i,System.BusData[:,:bus_i]))
        deleterows!(System.BusData_output,findall(x->x==i,System.BusData_output[:,:Bus]))
    end

    #updating line data
    rows_to_remove = findall(x->x==root_bus,System.LineData[:,:fbus])
    for row in rows_to_remove
        if System.LineData[row,:tbus] in non_aux_buses
            filter!(x->x≠ row,rows_to_remove)
        end
    end
    deleterows!(System.LineData,rows_to_remove)
    for i in 1:length(System.LineData[:,1])
        if System.LineData[i,:fbus] in bus_list
            System.LineData[i,:fbus] = root_bus
        elseif System.LineData[i,:tbus] in bus_list
            System.LineData[i,:tbus] = root_bus
        end
    end
end

function find_bus_substation(bus,System ::System_Struct)


    for substation in System.SubStations

        current_substation_buses = substation[2].SubStation_Buses
        # println(current_substation_buses)
        # substation_root_bus = current_substation_buses[1]
        if bus in current_substation_buses
            return substation[2]
        end

    end
end

function replace_buses!(main_bus,replacement,System ::System_Struct)

    #replace in bus data
    main_ind_bus_data = findall(x->x==main_bus,System.BusData[:bus_i])
    for i in main_ind_bus_data
        System.BusData[i,:bus_i] = replacement
        System.BusData_output[i,:Bus] = replacement
    end

    #replace in line data/loading
    main_ind_fbus = findall(x->x==main_bus,System.LineData[:fbus])
    for i in main_ind_fbus
        System.LineData[i,:fbus] = replacement
        System.LineLoading[i,:FromBus] = replacement
    end

    main_ind_tbus = findall(x->x==main_bus,System.LineData[:tbus])
    for i in main_ind_tbus
        System.LineData[i,:tbus] = replacement
        System.LineLoading[i,:ToBus] = replacement
    end

    #replace in Gen data
    main_ind_gen = findall(x->x==main_bus,System.Gen_Data[:bus])
    for i in main_ind_gen
        System.Gen_Data[i,:bus] = replacement
    end
end

function rename_buses!(System ::System_Struct)


    duplicate_buses_busdata = System.BusData[findall(x->x>length(System.SubStations),System.BusData[:bus_i]),:bus_i]
    duplicate_buses_gendata = System.Gen_Data[findall(x->x>length(System.SubStations),System.Gen_Data[:bus]),:bus]
    duplicate_buses = unique(vcat(duplicate_buses_busdata,duplicate_buses_gendata))

    for bus in duplicate_buses
        corresponding_substation = find_bus_substation(bus,System)
        corresponding_rootbus = corresponding_substation.Non_Aux_Buses[1]
        rootbus_exists = (corresponding_rootbus in System.BusData[:bus_i])
        is_bus_nonaux = (bus in corresponding_substation.Non_Aux_Buses)
        root_nonaux_exists = (corresponding_substation.Non_Aux_Buses[2] in System.BusData[:bus_i])
        if is_bus_nonaux
            if !rootbus_exists
                replace_buses!(bus,corresponding_rootbus,System)
            end
        else
            if rootbus_exists
                if !root_nonaux_exists
                    replace_buses!(bus,corresponding_substation.Non_Aux_Buses[2],System)
                end
            else
                replace_buses!(bus,corresponding_rootbus,System)
            end
        end
    end

end

function find_splitted_buses(System ::System_Struct)

    System_self = deepcopy(System)
    N_lines = length(System_self.LineLoading[:,1])
    filthy_buses_array = []
    for line in 1:N_lines
        if System_self.LineLoading[line,:Utilization] == 0
            ind_1 = findall(x->x==System_self.LineLoading[line,2],System_self.BusData_output[:Bus])
            if !(System_self.BusData_output[ind_1,:Pd] == 0) && (System_self.BusData_output[ind_1,:Pg] == 0)
                push!(filthy_buses_array,find_bus_substation(System_self.LineLoading[line,2],System).SubStationID)
            else
                push!(filthy_buses_array,find_bus_substation(System_self.LineLoading[line,1],System).SubStationID)
            end
        end
    end

    splitted_buses = []
    N_substations = length(System.SubStations)
    N_system = System.N_bus
    for i in N_substations+1:N_system

        splitted_bus = find_bus_substation(System.BusData[i,:bus_i],System).Non_Aux_Buses[1]

        if ! (splitted_bus in filthy_buses_array)
            push!(splitted_buses,splitted_bus)
        end
    end
    return unique(splitted_buses),filthy_buses_array
end

function find_splitted_buses_2(System ::System_Struct)

    splitted_buses = []
    N_substations = length(System.SubStations)
    N_system = System.N_bus
    for i in N_substations+1:N_system
        push!(splitted_buses,find_bus_substation(System.BusData[i,:bus_i],System).Non_Aux_Buses[1])
    end
    return unique(splitted_buses)
end

function Solve_Min_dev_ED(System ::System_Struct,ED_Pg)
    THETAMAX = 0.6
    M = 10000
    M2 = 2 * THETAMAX
    G = real(System.y_bus);
    B = imag(System.y_bus);
    b = System.b_line;

    N = System.N_bus
    N_Lines = size(System.Line_Constraints,1)
    Nodes_set = System.BusData[:,:bus_i];
    Gen_set = System.Gen_Data[:,:bus];

    Sbase = System.Sbase;

    BusData = System.BusData;
    GenData = System.Gen_Data;
    BranchData = System.Line_Constraints;

    df_lines = System.LineData[:,[:fbus,:tbus]]
    lines_set = Set.(convert.(Array, row for row in eachrow(df_lines)))

    df_aux = System.LineData[findall(x->x==1.0,System.LineData[:is_aux]),[:fbus,:tbus]]
    aux_line_set = Set.(convert.(Array, row for row in eachrow(df_aux)))

    df_non_aux = System.LineData[findall(x->x==0.0,System.LineData[:is_aux]),[:fbus,:tbus]]
    non_aux_line_set = Set.(convert.(Array, row for row in eachrow(df_non_aux)))

    Nodes_at_Lines = Dict()
    for i in System.LineData[:ID]
        Nodes_at_Lines[i] = Set([System.LineData[i,:fbus],System.LineData[i,:tbus]])
    end

    aux_lines_at_aux_bus = Dict()
    for i in System.BusData[:bus_i]
        if System.BusData[i,:is_aux] == 1.0
            lines_IDs = vcat(findall(x->x==i,System.LineData[:,:fbus]),findall(x->x==i,System.LineData[:,:tbus]))
            aux_lines = findall(x->x==1.0,System.LineData[lines_IDs,:is_aux])
            aux_lines_at_aux_bus[i] = lines_IDs[aux_lines]
        end
    end



    B_ = -B[1:N,1:N];

    m = JuMP.Model(JuMP.with_optimizer(Ipopt.Optimizer))

    # 2.1.Variables
    JuMP.@variable(m, -0.6≤ δ[i in Nodes_set] ≤ 0.6)

    JuMP.@variable(m, GenData[findall(x->x==g,System.Gen_Data[:,:bus]), :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[findall(x->x==g,System.Gen_Data[:,:bus]), :Pmax][1,1])

    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set])

    JuMP.@variable(m, 0 ≤ z[l=System.LineData[:ID];System.LineData[l,:is_aux]==1.0] ≤ 1) #capacity multiplier
    # JuMP.@variable(m,0 ≤ Bij_[i in Nodes_set,j in Nodes_set;Set([i,j]) in aux_line_set] ≤ 0.1)

    #Objective Function
    JuMP.@objective(m,Min,
        sum(GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:C1][1,1]*p[g]+GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:C0][1,1] for g in Gen_set))

    #Constraints

    #Current law applied on auxilliary buses
    # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
    #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

    JuMP.@constraint(m, ReferenceAngle,
        (δ[1] ==  0.0))
    #Current law applied on all buses
    JuMP.@constraint(m, Nodal_Balance[i = Nodes_set],
        sum(pij[i,j] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

    #Power flows in non auxilliary lines
    JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set ; Set([i,j]) in non_aux_line_set ],
        pij[i,j] == Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i]-δ[j]))

    #Power flows in auxilliary lines
    # JuMP.@constraint(m,Pl_aux[i in Nodes_set,j in Nodes_set ; Set([i,j]) in aux_line_set],
    #     pij[i,j] == Sbase*Bij_[i,j]*(δ[i]-δ[j]))

    JuMP.@constraint(m,total_balance, sum(p) == sum(System.BusData[:,:Pd]))

    #Switching constraints/capacity constraints
    JuMP.@NLconstraint(m,connectors_0[i = Nodes_set; BusData[i,:is_aux] == 1.0],
        0.45≤ (z[aux_lines_at_aux_bus[i][1]]-0.5)^2 + (z[aux_lines_at_aux_bus[i][2]]-0.5)^2 ≤ 0.5)

    JuMP.@NLconstraint(m,connectors_1[i = Nodes_set; BusData[i,:is_aux] == 1.0],
        (z[aux_lines_at_aux_bus[i][1]]-z[aux_lines_at_aux_bus[i][2]])^2 ≥ 0.9)

    JuMP.@constraint(m,connectors_2[i = Nodes_set; BusData[i,:is_aux] == 1.0],
         z[aux_lines_at_aux_bus[i][1]] + z[aux_lines_at_aux_bus[i][2]] == 1)

    #Switching constraints/admittance constraints
    # aux_lines_at_aux_bus[i][1] ----- aux_lines_at_aux_bus[i][2]
    # JuMP.@constraint(m,Bij_symmetry[i in Nodes_set,j in Nodes_set;Set([i,j]) in aux_line_set],Bij_[i,j] == Bij_[j,i])
    #
    # JuMP.@NLconstraint(m,connectors2_1[i = Nodes_set; BusData[i,:is_aux] == 1.0],
    #     (Bij_[collect(Nodes_at_Lines[aux_lines_at_aux_bus[i][1]])[1],collect(Nodes_at_Lines[aux_lines_at_aux_bus[i][1]])[2]]
    #     -Bij_[collect(Nodes_at_Lines[aux_lines_at_aux_bus[i][2]])[1],collect(Nodes_at_Lines[aux_lines_at_aux_bus[i][2]])[2]])^2 ≥ 0.05)
    #
    # JuMP.@constraint(m,connectors2_2[i = Nodes_set; BusData[i,:is_aux] == 1.0],
    #      Bij_[collect(Nodes_at_Lines[aux_lines_at_aux_bus[i][1]])[1],collect(Nodes_at_Lines[aux_lines_at_aux_bus[i][1]])[2]]
    #       + Bij_[collect(Nodes_at_Lines[aux_lines_at_aux_bus[i][2]])[1],collect(Nodes_at_Lines[aux_lines_at_aux_bus[i][2]])[2]] == 0.1)

    #Theta Constraints
    JuMP.@constraint(m,theta1[l = System.LineData[:,:ID] ;System.LineData[l,:is_aux] == 1.0],
        δ[System.LineData[l,:fbus]] <= δ[System.LineData[l,:tbus]] + (1-z[l])*M2)

    JuMP.@constraint(m,theta2[l = System.LineData[:,:ID] ;System.LineData[l,:is_aux] == 1.0],
        δ[System.LineData[l,:fbus]] + (1-z[l])*M2 >= δ[System.LineData[l,:tbus]])

    #Capacity Constraints

    JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID]],
        pij[System.LineData[i,:fbus],System.LineData[i,:tbus]] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus]])

    JuMP.@constraint(m,Line_Capacity_aux1_1[i = System.LineData[:,:ID] ;System.LineData[i,:is_aux] == 1.0],
        -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus]])

    JuMP.@constraint(m,Line_Capacity_aux1_2[i = System.LineData[:,:ID] ;System.LineData[i,:is_aux] == 1.0],
         pij[System.LineData[i,:fbus],System.LineData[i,:tbus]] <= System.LineData[i,:rate]*Sbase*z[i])

    JuMP.@constraint(m,Line_Capacity_aux2_1[i = System.LineData[:,:ID] ;System.LineData[i,:is_aux] == 1.0],
        -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus]])

    JuMP.@constraint(m,Line_Capacity_aux2[i = System.LineData[:,:ID] ;System.LineData[i,:is_aux] == 1.0],
        pij[System.LineData[i,:tbus],System.LineData[i,:fbus]] <= System.LineData[i,:rate]*Sbase*z[i])

    JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID] ;System.LineData[i,:is_aux] == 0.0],
        -System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus]] <= System.LineData[i,:rate]*Sbase)

    JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ;System.LineData[i,:is_aux] == 0.0],
        -System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus]] <= System.LineData[i,:rate]*Sbase)

    JuMP.optimize!(m)

    p_ij = JuMP.value.(pij)
    pg = JuMP.value.(p)
    z = JuMP.value.(z)
    δ_sol = JuMP.value.(δ)
    obj = JuMP.objective_value(m)

    switching_status = zeros(length(System.LineData[System.LineData[:is_aux].==1.0,:ID]),2)
    switching_status[:,1] = System.LineData[System.LineData[:is_aux].==1.0,:ID]
    switching_status[:,2] = round.([z[i] for i in System.LineData[System.LineData[:is_aux].==1.0,:ID]])


    Line_Duals = []
    for i in System.LineData[:,:ID]
        if System.LineData[i,:is_aux] == 1.0
            current_dual_1 = JuMP.dual(Line_Capacity_aux1_1[i])
            current_dual_2 = JuMP.dual(Line_Capacity_aux2_1[i])
            if current_dual_1 != 0
                Line_Duals = append!(Line_Duals,current_dual_1)
            elseif current_dual_2 != 0
                Line_Duals = append!(Line_Duals,current_dual_2)
            else
                Line_Duals = append!(Line_Duals,0)
            end
        end
    end

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
    println("==================================================================================================")
    println("|          System Summary                                                                        |")
    println("==================================================================================================")
    println("Objective Function Value = ", round(System.Operating_Cost,digits=3), " USD/hr")
    println(System.BusData_output)
    println()
    println(System.LineLoading)
    println("=========================== END OF REPORT ========================================")


    #reduce_System!(System ::System_Struct,switching_status)

    return System,z.data,pg.data,obj,Line_Duals
end

function Solve_OBS_SC(System ::System_Struct,method, Delta_P = 0,timelimit=Inf,warm_start = true)
    THETAMAX = 0.6
    M = 10000
    M2 = 2 * THETAMAX
    G = real(System.y_bus);
    B = imag(System.y_bus);
    b = System.b_line;

    N = System.N_bus
    N_Lines = size(System.Line_Constraints,1)
    Nodes_set = System.BusData[:,:bus_i];
    # Gen_set = System.Gen_Data[:,:bus];

    Sbase = System.Sbase;

    BusData = System.BusData;
    GenData = System.Gen_Data;
    BranchData = System.Line_Constraints;

    df_lines = System.LineData[:,[:fbus,:tbus]]
    lines_set = Set.(convert.(Array, row for row in eachrow(df_lines)))

    df_aux = System.LineData[findall(x->x==1.0,System.LineData[:is_aux]),[:fbus,:tbus]]
    aux_line_set = Set.(convert.(Array, row for row in eachrow(df_aux)))

    df_non_aux = System.LineData[findall(x->x==0.0,System.LineData[:is_aux]),[:fbus,:tbus]]
    non_aux_line_set = Set.(convert.(Array, row for row in eachrow(df_non_aux)))

    Nodes_at_Lines = Dict()
    for i in System.LineData[:ID]
        Nodes_at_Lines[i] = Set([System.LineData[i,:fbus],System.LineData[i,:tbus]])
    end

    #B_ = -B[1:N,1:N];

    subs_ = [k for k in keys(System.SubStations)]
    root_buses = [parse(Int,subs_[i][12:length(subs_[i])]) for i in 1:length(subs_)]

    if method == "P-SCOPF"
        a,K_ = Generate_Contingencies(System,method);
        K = 1:K_;
        Gen_set = 1:length(System.Gen_Data[:,:bus])

        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,TimeLimit=timelimit))

        #Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set,k in K] ≤ THETAMAX)

        JuMP.@variable(m, GenData[g ,:Pmin] ≤ p[g in Gen_set] ≤ GenData[g, :Pmax])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set,k in K])

        JuMP.@variable(m, z[l=System.LineData[:ID];System.LineData[l,:is_aux]==1.0], Bin)

        #Objective Function
        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g]+GenData[g,:C0][1,1] for g in Gen_set) )

        #Constraints

        #Current law applied on auxilliary buses
        # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
        #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

        # JuMP.@constraint(m, ReferenceAngle[k in K],
        #     (δ[Nodes_set[1],k] ==  0.0))
        #Current law applied on all buses
        JuMP.@constraint(m, Nodal_Balance[i = Nodes_set,k in K],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

        #Power flows in non auxilliary lines
        JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set,k in K ; Set([i,j]) in non_aux_line_set ],
            pij[i,j,k] == a[i,j,k]*Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i,k]-δ[j,k]))

        #JuMP.@constraint(m,total_balance, sum(p) == sum(System.BusData[:,:Pd]))
        #Switching constraint
        JuMP.@constraint(m,connectors[i = Nodes_set; BusData[i,:is_aux] == 1.0],
            sum(z[l] for l in vcat(findall(x->x==i,System.LineData[:,:fbus]),findall(x->x==i,System.LineData[:,:tbus])) if System.LineData[l,:is_aux] == 1.0) == 1)

        #Theta Constraints
        JuMP.@constraint(m,theta1[l = System.LineData[:,:ID],k in K ;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] <= δ[System.LineData[l,:tbus],k] + (1-z[l])*M2)

        JuMP.@constraint(m,theta2[l = System.LineData[:,:ID] ,k in K;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] + (1-z[l])*M2 >= δ[System.LineData[l,:tbus],k])

        #Capacity Constraints

        JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID],k in K],
            pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
             pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= System.LineData[i,:rate]*Sbase*z[i])

        JuMP.@constraint(m,Line_Capacity_aux2_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= System.LineData[i,:rate]*Sbase*z[i])

        JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 0.0],
            -a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase)

        JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ,k in K ;System.LineData[i,:is_aux] == 0.0],
            -a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase)

            if warm_start
                for aux_line in System.LineData[:ID]
                    if System.LineData[aux_line,:is_aux]==1.0
                        fbus = System.LineData[aux_line,:fbus]
                        if fbus in root_buses
                            JuMP.setvalue(z[aux_line],1)
                        else
                            JuMP.setvalue(z[aux_line],0)
                        end
                    end
                end
            end
        optimize!(m)

        p_ij = JuMP.value.(pij)
        pg = JuMP.value.(p)
        z = JuMP.value.(z)
        δ_sol = JuMP.value.(δ)
        obj = JuMP.objective_value(m)

        switching_status = zeros(length(System.LineData[System.LineData[:is_aux].==1.0,:ID]),2)
        switching_status[:,1] = System.LineData[System.LineData[:is_aux].==1.0,:ID]
        switching_status[:,2] = [z[i] for i in System.LineData[System.LineData[:is_aux].==1.0,:ID]]

        #========================#
        Pg = JuMP.value.(p)
        Pij = JuMP.value.(pij)
        δ = JuMP.value.(δ)
        V = ones(N)
        Pd = System.BusData[!,3]
        Qd = System.BusData[!,4]
        System.Operating_Cost = JuMP.objective_value(m)
        System.Gen_Data[:Pg] = Pg.data
        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1])]) : 0 for i in Nodes_set]
        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
         Pg = Pg , Pd = Pd)
        System.BusData_output = df

        Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
        Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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
        System.LineLoading[:ID] = deepcopy(System.LineData[:ID])
        #========================#
        reduce_System!(System,switching_status)
        removed_lines = sort(collect(setdiff(Set(System.LineLoading[:ID]),System.LineData[:ID])))
        deleterows!(System.LineLoading,removed_lines)
        System.LineLoading[:FromBus] = deepcopy(System.LineData[:fbus])
        System.LineLoading[:ToBus] = deepcopy(System.LineData[:tbus])

        #System.Gen_Data[:Pg] = pg.data
        return obj,System,z.data,pg.data
    elseif method == "OBS-P-1"
        a,K_ = Generate_Contingencies(System,"P-SCOPF");
        K = 1:K_;
        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,TimeLimit=timelimit))
        Gen_set = 1:length(System.Gen_Data[:,:bus])
        #Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set,k in K] ≤ THETAMAX)

        JuMP.@variable(m, GenData[g, :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[g, :Pmax][1,1])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set,k in K])

        JuMP.@variable(m, z[l=System.LineData[:ID],k in K;System.LineData[l,:is_aux]==1.0], Bin)

        #Objective Function
        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g]+GenData[g,:C0][1,1] for g in Gen_set) )

        #Constraints

        #Current law applied on auxilliary buses
        # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
        #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

        # JuMP.@constraint(m, ReferenceAngle[k in K],
        #     (δ[Nodes_set[1],k] ==  0.0))
        #Current law applied on all buses
        JuMP.@constraint(m, Nodal_Balance[i = Nodes_set,k in K],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

        #Power flows in non auxilliary lines
        JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set,k in K ; Set([i,j]) in non_aux_line_set ],
            pij[i,j,k] == a[i,j,k]*Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i,k]-δ[j,k]))

        JuMP.@constraint(m,total_balance, sum(p) == sum(System.BusData[:,:Pd]))
        #Switching constraint
        JuMP.@constraint(m,connectors[i = Nodes_set,k in K; BusData[i,:is_aux] == 1.0],
            sum(z[l,k] for l in vcat(findall(x->x==i,System.LineData[:,:fbus]),findall(x->x==i,System.LineData[:,:tbus])) if System.LineData[l,:is_aux] == 1.0) == 1)

        #Theta Constraints
        JuMP.@constraint(m,theta1[l = System.LineData[:,:ID],k in K ;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] <= δ[System.LineData[l,:tbus],k] + (1-z[l,k])*M2)

        JuMP.@constraint(m,theta2[l = System.LineData[:,:ID] ,k in K;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] + (1-z[l,k])*M2 >= δ[System.LineData[l,:tbus],k])

        #Capacity Constraints

        JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID],k in K],
            pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i,k] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
             pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= System.LineData[i,:rate]*Sbase*z[i,k])

        JuMP.@constraint(m,Line_Capacity_aux2_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i,k] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= System.LineData[i,:rate]*Sbase*z[i,k])

        JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 0.0],
            -a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase)

        JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ,k in K ;System.LineData[i,:is_aux] == 0.0],
            -a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase)

        if warm_start
            for k in K
                for aux_line in System.LineData[:ID]
                    if System.LineData[aux_line,:is_aux]==1.0
                        fbus = System.LineData[aux_line,:fbus]
                        if fbus in root_buses
                            JuMP.setvalue(z[aux_line,k],1)
                        else
                            JuMP.setvalue(z[aux_line,k],0)
                        end
                    end
                end
            end
        end
        optimize!(m)

        p_ij = JuMP.value.(pij)
        pg = JuMP.value.(p)
        z = JuMP.value.(z)
        δ_sol = JuMP.value.(δ)
        obj = JuMP.objective_value(m)

        switching_status = zeros(length(System.LineData[System.LineData[:is_aux].==1.0,:ID]),2)
        switching_status[:,1] = System.LineData[System.LineData[:is_aux].==1.0,:ID]
        switching_status[:,2] = [z[i,1] for i in System.LineData[System.LineData[:is_aux].==1.0,:ID]]

        #========================#
        Pg = JuMP.value.(p)
        Pij = JuMP.value.(pij)
        δ = JuMP.value.(δ)
        V = ones(N)
        Pd = System.BusData[!,3]
        Qd = System.BusData[!,4]
        System.Operating_Cost = JuMP.objective_value(m)
        System.Gen_Data[:Pg] = Pg.data
        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1])]) : 0 for i in Nodes_set]
        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
         Pg = Pg , Pd = Pd)
        System.BusData_output = df

        Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
        Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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
        System.LineLoading[:ID] = deepcopy(System.LineData[:ID])
        #========================#
        reduce_System!(System,switching_status)
        removed_lines = sort(collect(setdiff(Set(System.LineLoading[:ID]),System.LineData[:ID])))
        deleterows!(System.LineLoading,removed_lines)
        System.LineLoading[:FromBus] = deepcopy(System.LineData[:fbus])
        System.LineLoading[:ToBus] = deepcopy(System.LineData[:tbus])

        #System.Gen_Data[:Pg] = pg.data
        return obj,System,z.data,pg.data
    elseif method == "OBS-P-2"
        subs_ = [k for k in keys(System.SubStations)]
        root_buses = [parse(Int,subs_[i][12:length(subs_[i])]) for i in 1:length(subs_)]
        a,K_ = Generate_Contingencies(System,"P-SCOPF");
        K = 1:K_;
        Gen_set = 1:length(System.Gen_Data[:,:bus])
        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,TimeLimit=timelimit))

        #Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set,k in K] ≤ THETAMAX)

        JuMP.@variable(m, GenData[g, :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[g, :Pmax][1,1])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set,k in K])

        JuMP.@variable(m, z[l=System.LineData[:ID],k in K;System.LineData[l,:is_aux]==1.0], Bin)

        #Objective Function
        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g]+GenData[g,:C0][1,1] for g in Gen_set) )

        #Constraints

        #Current law applied on auxilliary buses
        # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
        #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

        # JuMP.@constraint(m, ReferenceAngle[k in K],
        #     (δ[Nodes_set[1],k] ==  0.0))
        #Current law applied on all buses
        JuMP.@constraint(m, Nodal_Balance[i = Nodes_set,k in K],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

        #Power flows in non auxilliary lines
        JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set,k in K ; Set([i,j]) in non_aux_line_set ],
            pij[i,j,k] == a[i,j,k]*Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i,k]-δ[j,k]))

        JuMP.@constraint(m,total_balance, sum(p) == sum(System.BusData[:,:Pd]))
        #Switching constraint
        JuMP.@constraint(m,connectors[i = Nodes_set,k in K; BusData[i,:is_aux] == 1.0],
            sum(z[l,k] for l in vcat(findall(x->x==i,System.LineData[:,:fbus]),findall(x->x==i,System.LineData[:,:tbus])) if System.LineData[l,:is_aux] == 1.0) == 1)

        #Theta Constraints
        JuMP.@constraint(m,precontingency[l=System.LineData[:ID]; (System.LineData[l,:is_aux] == 1.0) && (System.LineData[l,:fbus] in root_buses)],
            z[l,1] == 1)
        JuMP.@constraint(m,precontingency2[l=System.LineData[:ID]; (System.LineData[l,:is_aux] == 1.0) && !(System.LineData[l,:fbus] in root_buses)],
            z[l,1] == 0)

        JuMP.@constraint(m,theta1[l = System.LineData[:,:ID],k in K ;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] <= δ[System.LineData[l,:tbus],k] + (1-z[l,k])*M2)

        JuMP.@constraint(m,theta2[l = System.LineData[:,:ID] ,k in K;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] + (1-z[l,k])*M2 >= δ[System.LineData[l,:tbus],k])

        #Capacity Constraints

        JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID],k in K],
            pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i,k] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
             pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= System.LineData[i,:rate]*Sbase*z[i,k])

        JuMP.@constraint(m,Line_Capacity_aux2_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i,k] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= System.LineData[i,:rate]*Sbase*z[i,k])

        JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 0.0],
            -a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase)

        JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ,k in K ;System.LineData[i,:is_aux] == 0.0],
            -a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase)

            if warm_start
                for k in K[2:length(K)]
                    for aux_line in System.LineData[:ID]
                        if System.LineData[aux_line,:is_aux]==1.0
                            fbus = System.LineData[aux_line,:fbus]
                            if fbus in root_buses
                                JuMP.setvalue(z[aux_line,k],1)
                            else
                                JuMP.setvalue(z[aux_line,k],0)
                            end
                        end
                    end
                end
            end
        optimize!(m)

        p_ij = JuMP.value.(pij)
        pg = JuMP.value.(p)
        z = JuMP.value.(z)
        δ_sol = JuMP.value.(δ)
        obj = JuMP.objective_value(m)

        switching_status = zeros(length(System.LineData[System.LineData[:is_aux].==1.0,:ID]),2)
        switching_status[:,1] = System.LineData[System.LineData[:is_aux].==1.0,:ID]
        switching_status[:,2] = [z[i,1] for i in System.LineData[System.LineData[:is_aux].==1.0,:ID]]

        #========================#
        Pg = JuMP.value.(p)
        Pij = JuMP.value.(pij)
        δ = JuMP.value.(δ)
        V = ones(N)
        Pd = System.BusData[!,3]
        Qd = System.BusData[!,4]
        System.Operating_Cost = JuMP.objective_value(m)
        System.Gen_Data[:Pg] = Pg.data
        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1])]) : 0 for i in Nodes_set]
        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
         Pg = Pg , Pd = Pd)
        System.BusData_output = df

        Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
        Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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
        System.LineLoading[:ID] = deepcopy(System.LineData[:ID])
        #========================#
        reduce_System!(System,switching_status)
        removed_lines = sort(collect(setdiff(Set(System.LineLoading[:ID]),System.LineData[:ID])))
        deleterows!(System.LineLoading,removed_lines)
        System.LineLoading[:FromBus] = deepcopy(System.LineData[:fbus])
        System.LineLoading[:ToBus] = deepcopy(System.LineData[:tbus])

        #System.Gen_Data[:Pg] = pg.data
        return obj,System,z.data,pg.data

    elseif method == "OBS-P-3"
        
        subs_ = [k for k in keys(System.SubStations)]
        root_buses = [parse(Int,subs_[i][12:length(subs_[i])]) for i in 1:length(subs_)]
        a,K_ = Generate_Contingencies(System,"P-SCOPF");
        K = 1:K_;
        Gen_set = 1:length(System.Gen_Data[:,:bus])
        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,TimeLimit=timelimit))

        #Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set,k in K] ≤ THETAMAX)

        JuMP.@variable(m, GenData[g, :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[g, :Pmax][1,1])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set,k in K])

        JuMP.@variable(m, z[l=System.LineData[:ID] ;System.LineData[l,:is_aux]==1.0], Bin)

        JuMP.@variable(m, z_coupler[l=System.LineData[:ID],k in K ;System.LineData[l,:is_aux]==2.0], Bin)
        #Objective Function
        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g]+GenData[g,:C0][1,1] for g in Gen_set) )

        #Constraints

        #Current law applied on auxilliary buses
        # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
        #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

        # JuMP.@constraint(m, ReferenceAngle[k in K],
        #     (δ[Nodes_set[1],k] ==  0.0))
        #Current law applied on all buses
        JuMP.@constraint(m, Nodal_Balance[i = Nodes_set,k in K],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

        #Power flows in non auxilliary lines
        JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set,k in K ; Set([i,j]) in non_aux_line_set ],
            pij[i,j,k] == a[i,j,k]*Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i,k]-δ[j,k]))

        JuMP.@constraint(m,total_balance, sum(p) == sum(System.BusData[:,:Pd]))

        #Switching constraint
        JuMP.@constraint(m,connectors[i = Nodes_set; BusData[i,:is_aux] == 1.0],
            sum(z[l] for l in vcat(findall(x->x==i,System.LineData[:,:fbus]),findall(x->x==i,System.LineData[:,:tbus])) if System.LineData[l,:is_aux] == 1.0) == 1)

        JuMP.@constraint(m,couplers[l = System.LineData[:ID]; System.LineData[l,:is_aux]==2.0],
            z_coupler[l,1] == 1.0)
        #Theta Constraints


        JuMP.@constraint(m,theta1[l = System.LineData[:,:ID],k in K ;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] <= δ[System.LineData[l,:tbus],k] + (1-z[l])*M2)

        JuMP.@constraint(m,theta2[l = System.LineData[:,:ID] ,k in K;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] + (1-z[l])*M2 >= δ[System.LineData[l,:tbus],k])

        JuMP.@constraint(m,theta3[l = System.LineData[:,:ID],k in K ;System.LineData[l,:is_aux] == 2.0],
            δ[System.LineData[l,:fbus],k] <= δ[System.LineData[l,:tbus],k] + (1-z_coupler[l,k])*M2)

        JuMP.@constraint(m,theta4[l = System.LineData[:,:ID] ,k in K;System.LineData[l,:is_aux] == 2.0],
            δ[System.LineData[l,:fbus],k] + (1-z_coupler[l,k])*M2 >= δ[System.LineData[l,:tbus],k])

        #Capacity Constraints

        JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID],k in K],
            pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        #Auxilliary lines constraints
        JuMP.@constraint(m,Line_Capacity_aux1_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
             pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= System.LineData[i,:rate]*Sbase*z[i])

        JuMP.@constraint(m,Line_Capacity_aux2_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= System.LineData[i,:rate]*Sbase*z[i])

        # Coupler constraints
        JuMP.@constraint(m,Line_Capacity_coupler1_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 2.0],
            -System.LineData[i,:rate]*Sbase*z_coupler[i,k] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k])

        JuMP.@constraint(m,Line_Capacity_coupler1_2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 2.0],
             pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= System.LineData[i,:rate]*Sbase*z_coupler[i,k])

        JuMP.@constraint(m,Line_Capacity_coupler2_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 2.0],
            -System.LineData[i,:rate]*Sbase*z_coupler[i,k] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_coupler2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 2.0],
            pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= System.LineData[i,:rate]*Sbase*z_coupler[i,k])

        #Non-auxilliary lines constraints
        JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 0.0],
            -a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase)

        JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ,k in K ;System.LineData[i,:is_aux] == 0.0],
            -a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase)

            if warm_start
                for aux_line in System.LineData[:ID]
                    if System.LineData[aux_line,:is_aux]==1.0
                        fbus = System.LineData[aux_line,:fbus]
                        if fbus in root_buses
                            JuMP.setvalue(z[aux_line],1)
                        else
                            JuMP.setvalue(z[aux_line],0)
                        end
                    end
                end
            end
        optimize!(m)

        p_ij = JuMP.value.(pij)
        pg = JuMP.value.(p)
        z = JuMP.value.(z)
        z_coupler = JuMP.value.(z_coupler)
        δ_sol = JuMP.value.(δ)
        obj = JuMP.objective_value(m)

        switching_status = zeros(length(System.LineData[System.LineData[:is_aux].==1.0,:ID]),2)
        switching_status[:,1] = System.LineData[System.LineData[:is_aux].==1.0,:ID]
        switching_status[:,2] = [z[i] for i in System.LineData[System.LineData[:is_aux].==1.0,:ID]]

        switching_status_coupler = zeros(length(System.LineData[System.LineData[:is_aux].==2.0,:ID]),2)
        switching_status_coupler[:,1] = System.LineData[System.LineData[:is_aux].==2.0,:ID]
        switching_status_coupler[:,2] = [z_coupler[i,1] for i in System.LineData[System.LineData[:is_aux].==2.0,:ID]]

        #========================#
        Pg = JuMP.value.(p)
        Pij = JuMP.value.(pij)
        δ = JuMP.value.(δ)
        V = ones(N)
        Pd = System.BusData[!,3]
        Qd = System.BusData[!,4]
        System.Operating_Cost = JuMP.objective_value(m)
        System.Gen_Data[:Pg] = Pg.data
        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1])]) : 0 for i in Nodes_set]
        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
         Pg = Pg , Pd = Pd)
        System.BusData_output = df

        Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
        Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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
        System.LineLoading[:ID] = deepcopy(System.LineData[:ID])
        #========================#
        reduce_System_coupled!(System,switching_status,switching_status_coupler)
        removed_lines = sort(collect(setdiff(Set(System.LineLoading[:ID]),System.LineData[:ID])))
        deleterows!(System.LineLoading,removed_lines)
        System.LineLoading[:FromBus] = deepcopy(System.LineData[:fbus])
        System.LineLoading[:ToBus] = deepcopy(System.LineData[:tbus])

        #System.Gen_Data[:Pg] = pg.data
        return obj,System,z.data,pg.data, z_coupler.data

    elseif method == "C-SCOPF"
        a_g,a_l,K_ = Generate_Contingencies(System,method);
        K = 1:K_;
        Gen_set = 1:length(System.Gen_Data[:,:bus]);
        #Gen_buses = System.Gen_Data[:,:bus]
        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,TimeLimit=timelimit))

        # 2.1.Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set,k in K] ≤ THETAMAX)

        JuMP.@variable(m, a_g[g,1,k]*GenData[g, :Pmin][1,1] ≤ p[g in Gen_set,k in K] ≤ a_g[g,1,k]*GenData[g, :Pmax][1,1])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set,k in K])

        JuMP.@variable(m, z[l=System.LineData[:ID];System.LineData[l,:is_aux]==1.0], Bin)

        #Objective Function
        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g,1]+GenData[g,:C0][1,1] for g in Gen_set) )

        #Constraints

        #Current law applied on auxilliary buses
        # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
        #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

        # JuMP.@constraint(m, ReferenceAngle[k in K],
        #     (δ[Nodes_set[1],k] ==  0.0))
        #Current law applied on all buses
        JuMP.@constraint(m, Nodal_Balance[i = Nodes_set,k in K],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g,k] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

        #Power flows in non auxilliary lines
        JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set,k in K ; Set([i,j]) in non_aux_line_set ],
            pij[i,j,k] == a_l[i,j,k]*Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i,k]-δ[j,k]))

        #JuMP.@constraint(m,total_balance[k in K], sum(p[g,k] for g in Gen_set) == sum(System.BusData[:,:Pd]))
        #Switching constraint
        JuMP.@constraint(m,connectors[i = Nodes_set; BusData[i,:is_aux] == 1.0],
            sum(z[l] for l in vcat(findall(x->x==i,System.LineData[:,:fbus]),findall(x->x==i,System.LineData[:,:tbus])) if System.LineData[l,:is_aux] == 1.0) == 1)

        #Theta Constraints
        JuMP.@constraint(m,theta1[l = System.LineData[:,:ID],k in K ;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] <= δ[System.LineData[l,:tbus],k] + (1-z[l])*M2)

        JuMP.@constraint(m,theta2[l = System.LineData[:,:ID] ,k in K;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] + (1-z[l])*M2 >= δ[System.LineData[l,:tbus],k])

        JuMP.@constraint(m,Lim_Gen_up[i in Gen_set,k = 2:length(K)],p[i,k]-p[i,1] <= Delta_P[i])
        JuMP.@constraint(m,Lim_Gen_Dw[i in Gen_set,k = 2:length(K)],p[i,1]-p[i,k] <= Delta_P[i])
        #Capacity Constraints

        JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID],k in K],
            pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
             pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= System.LineData[i,:rate]*Sbase*z[i])

        # JuMP.@constraint(m,Line_Capacity_aux2_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
        #     -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])
        #
        # JuMP.@constraint(m,Line_Capacity_aux2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
        #     pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= System.LineData[i,:rate]*Sbase*z[i])

        JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 0.0],
            -a_l[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= a_l[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase)

        JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ,k in K ;System.LineData[i,:is_aux] == 0.0],
            -a_l[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= a_l[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase)

            if warm_start
                for aux_line in System.LineData[:ID]
                    if System.LineData[aux_line,:is_aux]==1.0
                        fbus = System.LineData[aux_line,:fbus]
                        if fbus in root_buses
                            JuMP.setvalue(z[aux_line],1)
                        else
                            JuMP.setvalue(z[aux_line],0)
                        end
                    end
                end
            end

        optimize!(m)

        p_ij = JuMP.value.(pij)
        pg = JuMP.value.(p)
        z = JuMP.value.(z)
        δ_sol = JuMP.value.(δ)
        obj = JuMP.objective_value(m)

        switching_status = zeros(length(System.LineData[System.LineData[:is_aux].==1.0,:ID]),2)
        switching_status[:,1] = System.LineData[System.LineData[:is_aux].==1.0,:ID]
        switching_status[:,2] = [z[i] for i in System.LineData[System.LineData[:is_aux].==1.0,:ID]]
        #========================#
        Pg = JuMP.value.(p)
        Pij = JuMP.value.(pij)
        δ = JuMP.value.(δ)
        V = ones(N)
        Pd = System.BusData[!,3]
        Qd = System.BusData[!,4]
        System.Operating_Cost = JuMP.objective_value(m)
        System.Gen_Data[:Pg] = Pg.data[:,1]
        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1]),1]) : 0 for i in Nodes_set]
        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
         Pg = Pg , Pd = Pd)
        System.BusData_output = df

        Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
        Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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
        System.LineLoading[:ID] = deepcopy(System.LineData[:ID])
        #========================#
        reduce_System!(System ::System_Struct,switching_status)
        removed_lines = sort(collect(setdiff(Set(System.LineLoading[:ID]),System.LineData[:ID])))
        deleterows!(System.LineLoading,removed_lines)
        System.LineLoading[:FromBus] = deepcopy(System.LineData[:fbus])
        System.LineLoading[:ToBus] = deepcopy(System.LineData[:tbus])

        return obj,System,z.data,pg.data
    elseif method == "OBS-C-1"
        a_g,a_l,K_ = Generate_Contingencies(System,"C-SCOPF");
        K = 1:K_;
        Gen_set = 1:length(System.Gen_Data[:,:bus]);
        #Gen_buses = System.Gen_Data[:,:bus]
        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,TimeLimit=timelimit))

        # 2.1.Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set,k in K] ≤ THETAMAX)

        JuMP.@variable(m, a_g[g,1,k]*GenData[g, :Pmin][1,1] ≤ p[g in Gen_set,k in K] ≤ a_g[g,1,k]*GenData[g, :Pmax][1,1])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set,k in K])

        JuMP.@variable(m, z[l=System.LineData[:ID],k in K;System.LineData[l,:is_aux]==1.0], Bin)

        #Objective Function
        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g,1]+GenData[g,:C0][1,1] for g in Gen_set) )

        #Constraints

        #Current law applied on auxilliary buses
        # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
        #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

        # JuMP.@constraint(m, ReferenceAngle[k in K],
        #     (δ[Nodes_set[1],k] ==  0.0))
        #Current law applied on all buses
        JuMP.@constraint(m, Nodal_Balance[i = Nodes_set,k in K],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g,k] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

        #Power flows in non auxilliary lines
        JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set,k in K ; Set([i,j]) in non_aux_line_set ],
            pij[i,j,k] == a_l[i,j,k]*Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i,k]-δ[j,k]))

        #JuMP.@constraint(m,total_balance[k in K], sum(p[g,k] for g in Gen_set) == sum(System.BusData[:,:Pd]))
        #Switching constraint
        JuMP.@constraint(m,connectors[i = Nodes_set,k in K; BusData[i,:is_aux] == 1.0],
            sum(z[l,k] for l in vcat(findall(x->x==i,System.LineData[:,:fbus]),findall(x->x==i,System.LineData[:,:tbus])) if System.LineData[l,:is_aux] == 1.0) == 1)

        #Theta Constraints
        JuMP.@constraint(m,theta1[l = System.LineData[:,:ID],k in K ;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] <= δ[System.LineData[l,:tbus],k] + (1-z[l,k])*M2)

        JuMP.@constraint(m,theta2[l = System.LineData[:,:ID] ,k in K;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] + (1-z[l,k])*M2 >= δ[System.LineData[l,:tbus],k])

        JuMP.@constraint(m,Lim_Gen_up[i in Gen_set,k = 2:length(K)],p[i,k]-p[i,1] <= Delta_P[i])
        JuMP.@constraint(m,Lim_Gen_Dw[i in Gen_set,k = 2:length(K)],p[i,1]-p[i,k] <= Delta_P[i])
        #Capacity Constraints

        JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID],k in K],
            pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i,k] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
             pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= System.LineData[i,:rate]*Sbase*z[i,k])

        JuMP.@constraint(m,Line_Capacity_aux2_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i,k] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= System.LineData[i,:rate]*Sbase*z[i,k])

        JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 0.0],
            -a_l[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= a_l[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase)

        JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ,k in K ;System.LineData[i,:is_aux] == 0.0],
            -a_l[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= a_l[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase)

            if warm_start
                for k in K
                    for aux_line in System.LineData[:ID]
                        if System.LineData[aux_line,:is_aux]==1.0
                            fbus = System.LineData[aux_line,:fbus]
                            if fbus in root_buses
                                JuMP.setvalue(z[aux_line,k],1)
                            else
                                JuMP.setvalue(z[aux_line,k],0)
                            end
                        end
                    end
                end
            end
        optimize!(m)

        p_ij = JuMP.value.(pij)
        pg = JuMP.value.(p)
        z = JuMP.value.(z)
        δ_sol = JuMP.value.(δ)
        obj = JuMP.objective_value(m)

        switching_status = zeros(length(System.LineData[System.LineData[:is_aux].==1.0,:ID]),2)
        switching_status[:,1] = System.LineData[System.LineData[:is_aux].==1.0,:ID]
        switching_status[:,2] = [z[i,1] for i in System.LineData[System.LineData[:is_aux].==1.0,:ID]]
        #println(switching_status)
        #========================#
        Pg = JuMP.value.(p)
        Pij = JuMP.value.(pij)
        δ = JuMP.value.(δ)
        V = ones(N)
        Pd = System.BusData[!,3]
        Qd = System.BusData[!,4]
        System.Operating_Cost = JuMP.objective_value(m)
        System.Gen_Data[:Pg] = Pg.data[:,1]
        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1]),1]) : 0 for i in Nodes_set]
        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
         Pg = Pg , Pd = Pd)
        System.BusData_output = df

        Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
        Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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
        System.LineLoading[:ID] = deepcopy(System.LineData[:ID])
        #========================#
        reduce_System!(System,switching_status)
        removed_lines = sort(collect(setdiff(Set(System.LineLoading[:ID]),System.LineData[:ID])))
        deleterows!(System.LineLoading,removed_lines)
        System.LineLoading[:FromBus] = deepcopy(System.LineData[:fbus])
        System.LineLoading[:ToBus] = deepcopy(System.LineData[:tbus])

        return obj,System,z.data,pg.data
    elseif method == "OBS-C-2"

        subs_ = [k for k in keys(System.SubStations)]
        root_buses = [parse(Int,subs_[i][12:length(subs_[i])]) for i in 1:length(subs_)]
        a_g,a_l,K_ = Generate_Contingencies(System,"C-SCOPF");
        K = 1:K_;
        Gen_set = 1:length(System.Gen_Data[:,:bus]);
        #Gen_buses = System.Gen_Data[:,:bus]
        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,TimeLimit=timelimit))

        # 2.1.Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set,k in K] ≤ THETAMAX)

        JuMP.@variable(m, a_g[g,1,k]*GenData[g, :Pmin][1,1] ≤ p[g in Gen_set,k in K] ≤ a_g[g,1,k]*GenData[g, :Pmax][1,1])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set,k in K])

        JuMP.@variable(m, z[l=System.LineData[:ID],k in K;System.LineData[l,:is_aux]==1.0], Bin)

        #Objective Function
        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g,1]+GenData[g,:C0][1,1] for g in Gen_set) )

        #Constraints

        #Current law applied on auxilliary buses
        # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
        #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

        # JuMP.@constraint(m, ReferenceAngle[k in K],
        #     (δ[Nodes_set[1],k] ==  0.0))
        #Current law applied on all buses
        JuMP.@constraint(m, Nodal_Balance[i = Nodes_set,k in K],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g,k] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

        #Power flows in non auxilliary lines
        JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set,k in K ; Set([i,j]) in non_aux_line_set ],
            pij[i,j,k] == a_l[i,j,k]*Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i,k]-δ[j,k]))

        #JuMP.@constraint(m,total_balance[k in K], sum(p[g,k] for g in Gen_set) == sum(System.BusData[:,:Pd]))
        #Switching constraint
        JuMP.@constraint(m,precontingency[l=System.LineData[:ID]; (System.LineData[l,:is_aux] == 1.0) && (System.LineData[l,:fbus] in root_buses)],
            z[l,1] == 1)
        JuMP.@constraint(m,precontingency2[l=System.LineData[:ID]; (System.LineData[l,:is_aux] == 1.0) && !(System.LineData[l,:fbus] in root_buses)],
            z[l,1] == 0)
        JuMP.@constraint(m,connectors[i = Nodes_set,k in K; BusData[i,:is_aux] == 1.0],
            sum(z[l,k] for l in vcat(findall(x->x==i,System.LineData[:,:fbus]),findall(x->x==i,System.LineData[:,:tbus])) if System.LineData[l,:is_aux] == 1.0) == 1)

        #Theta Constraints
        JuMP.@constraint(m,theta1[l = System.LineData[:,:ID],k in K ;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] <= δ[System.LineData[l,:tbus],k] + (1-z[l,k])*M2)

        JuMP.@constraint(m,theta2[l = System.LineData[:,:ID] ,k in K ;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] + (1-z[l,k])*M2 >= δ[System.LineData[l,:tbus],k])

        JuMP.@constraint(m,Lim_Gen_up[i in Gen_set,k = 2:length(K)],p[i,k]-p[i,1] <= Delta_P[i])
        JuMP.@constraint(m,Lim_Gen_Dw[i in Gen_set,k = 2:length(K)],p[i,1]-p[i,k] <= Delta_P[i])
        #Capacity Constraints

        JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID],k in K],
            pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_1[i = System.LineData[:,:ID],k in K;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i,k] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
             pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= System.LineData[i,:rate]*Sbase*z[i,k])

        JuMP.@constraint(m,Line_Capacity_aux2_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i,k] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= System.LineData[i,:rate]*Sbase*z[i,k])

        JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 0.0],
            -a_l[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= a_l[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase)

        JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ,k in K ;System.LineData[i,:is_aux] == 0.0],
            -a_l[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= a_l[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase)

        optimize!(m)

        if warm_start
            for k in K[2:length(K)]
                for aux_line in System.LineData[:ID]
                    if System.LineData[aux_line,:is_aux]==1.0
                        fbus = System.LineData[aux_line,:fbus]
                        if fbus in root_buses
                            JuMP.setvalue(z[aux_line,k],1)
                        else
                            JuMP.setvalue(z[aux_line,k],0)
                        end
                    end
                end
            end
        end
        p_ij = JuMP.value.(pij)
        pg = JuMP.value.(p)
        z = JuMP.value.(z)
        δ_sol = JuMP.value.(δ)
        obj = JuMP.objective_value(m)

        switching_status = zeros(length(System.LineData[System.LineData[:is_aux].==1.0,:ID]),2)
        switching_status[:,1] = System.LineData[System.LineData[:is_aux].==1.0,:ID]
        switching_status[:,2] = [z[i,1] for i in System.LineData[System.LineData[:is_aux].==1.0,:ID]]
        #println(switching_status)
        #========================#
        Pg = JuMP.value.(p)
        Pij = JuMP.value.(pij)
        δ = JuMP.value.(δ)
        V = ones(N)
        Pd = System.BusData[!,3]
        Qd = System.BusData[!,4]
        System.Operating_Cost = JuMP.objective_value(m)
        System.Gen_Data[:Pg] = Pg.data[:,1]
        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1]),1]) : 0 for i in Nodes_set]
        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
         Pg = Pg , Pd = Pd)
        System.BusData_output = df

        Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
        Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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
        System.LineLoading[:ID] = deepcopy(System.LineData[:ID])
        #========================#
        reduce_System!(System ::System_Struct,switching_status)
        removed_lines = sort(collect(setdiff(Set(System.LineLoading[:ID]),System.LineData[:ID])))
        deleterows!(System.LineLoading,removed_lines)
        System.LineLoading[:FromBus] = deepcopy(System.LineData[:fbus])
        System.LineLoading[:ToBus] = deepcopy(System.LineData[:tbus])

        return obj,System,z.data,pg.data
    end
end

function Solve_OBS_SC_coupled(System ::System_Struct,method, Delta_P = 0,timelimit=Inf,warm_start = true)
    THETAMAX = 0.6
    M = 10000
    M2 = 2 * THETAMAX
    G = real(System.y_bus);
    B = imag(System.y_bus);
    b = System.b_line;

    N = System.N_bus
    N_Lines = size(System.Line_Constraints,1)
    Nodes_set = System.BusData[:,:bus_i];
    # Gen_set = System.Gen_Data[:,:bus];

    Sbase = System.Sbase;

    BusData = System.BusData;
    GenData = System.Gen_Data;
    BranchData = System.Line_Constraints;

    df_lines = System.LineData[:,[:fbus,:tbus]]
    lines_set = Set.(convert.(Array, row for row in eachrow(df_lines)))

    df_aux = System.LineData[findall(x->x==1.0,System.LineData[:is_aux]),[:fbus,:tbus]]
    aux_line_set = Set.(convert.(Array, row for row in eachrow(df_aux)))

    df_non_aux = System.LineData[findall(x->x==0.0,System.LineData[:is_aux]),[:fbus,:tbus]]
    non_aux_line_set = Set.(convert.(Array, row for row in eachrow(df_non_aux)))

    Nodes_at_Lines = Dict()
    for i in System.LineData[:ID]
        Nodes_at_Lines[i] = Set([System.LineData[i,:fbus],System.LineData[i,:tbus]])
    end

    #B_ = -B[1:N,1:N];

    subs_ = [k for k in keys(System.SubStations)]
    root_buses = [parse(Int,subs_[i][12:length(subs_[i])]) for i in 1:length(subs_)]

    if method == "P-SCOPF"
        a,K_ = Generate_Contingencies(System,method);
        K = 1:K_;
        Gen_set = 1:length(System.Gen_Data[:,:bus])

        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,TimeLimit=timelimit))

        #Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set,k in K] ≤ THETAMAX)

        JuMP.@variable(m, GenData[g ,:Pmin] ≤ p[g in Gen_set] ≤ GenData[g, :Pmax])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set,k in K])

        JuMP.@variable(m, z[l=System.LineData[:ID];System.LineData[l,:is_aux]==1.0], Bin)

        #Objective Function
        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g]+GenData[g,:C0][1,1] for g in Gen_set) )

        #Constraints

        #Current law applied on auxilliary buses
        # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
        #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

        # JuMP.@constraint(m, ReferenceAngle[k in K],
        #     (δ[Nodes_set[1],k] ==  0.0))
        #Current law applied on all buses
        JuMP.@constraint(m, Nodal_Balance[i = Nodes_set,k in K],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

        #Power flows in non auxilliary lines
        JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set,k in K ; Set([i,j]) in non_aux_line_set ],
            pij[i,j,k] == a[i,j,k]*Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i,k]-δ[j,k]))

        #JuMP.@constraint(m,total_balance, sum(p) == sum(System.BusData[:,:Pd]))
        #Switching constraint
        JuMP.@constraint(m,connectors[i = Nodes_set; BusData[i,:is_aux] == 1.0],
            sum(z[l] for l in vcat(findall(x->x==i,System.LineData[:,:fbus]),findall(x->x==i,System.LineData[:,:tbus])) if System.LineData[l,:is_aux] == 1.0) == 1)

        #Theta Constraints
        JuMP.@constraint(m,theta1[l = System.LineData[:,:ID],k in K ;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] <= δ[System.LineData[l,:tbus],k] + (1-z[l])*M2)

        JuMP.@constraint(m,theta2[l = System.LineData[:,:ID] ,k in K;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] + (1-z[l])*M2 >= δ[System.LineData[l,:tbus],k])

        #Capacity Constraints

        JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID],k in K],
            pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
             pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= System.LineData[i,:rate]*Sbase*z[i])

        JuMP.@constraint(m,Line_Capacity_aux2_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= System.LineData[i,:rate]*Sbase*z[i])

        JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 0.0],
            -a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase)

        JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ,k in K ;System.LineData[i,:is_aux] == 0.0],
            -a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase)

            if warm_start
                for aux_line in System.LineData[:ID]
                    if System.LineData[aux_line,:is_aux]==1.0
                        fbus = System.LineData[aux_line,:fbus]
                        if fbus in root_buses
                            JuMP.setvalue(z[aux_line],1)
                        else
                            JuMP.setvalue(z[aux_line],0)
                        end
                    end
                end
            end
        optimize!(m)

        p_ij = JuMP.value.(pij)
        pg = JuMP.value.(p)
        z = JuMP.value.(z)
        δ_sol = JuMP.value.(δ)
        obj = JuMP.objective_value(m)

        switching_status = zeros(length(System.LineData[System.LineData[:is_aux].==1.0,:ID]),2)
        switching_status[:,1] = System.LineData[System.LineData[:is_aux].==1.0,:ID]
        switching_status[:,2] = [z[i] for i in System.LineData[System.LineData[:is_aux].==1.0,:ID]]

        #========================#
        Pg = JuMP.value.(p)
        Pij = JuMP.value.(pij)
        δ = JuMP.value.(δ)
        V = ones(N)
        Pd = System.BusData[!,3]
        Qd = System.BusData[!,4]
        System.Operating_Cost = JuMP.objective_value(m)
        System.Gen_Data[:Pg] = Pg.data
        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1])]) : 0 for i in Nodes_set]
        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
         Pg = Pg , Pd = Pd)
        System.BusData_output = df

        Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
        Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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
        System.LineLoading[:ID] = deepcopy(System.LineData[:ID])
        #========================#
        reduce_System!(System,switching_status)
        removed_lines = sort(collect(setdiff(Set(System.LineLoading[:ID]),System.LineData[:ID])))
        deleterows!(System.LineLoading,removed_lines)
        System.LineLoading[:FromBus] = deepcopy(System.LineData[:fbus])
        System.LineLoading[:ToBus] = deepcopy(System.LineData[:tbus])

        #System.Gen_Data[:Pg] = pg.data
        return obj,System,z.data,pg.data
    elseif method == "OBS-P-1"
        a,K_ = Generate_Contingencies(System,"P-SCOPF");
        K = 1:K_;
        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,TimeLimit=timelimit))
        Gen_set = 1:length(System.Gen_Data[:,:bus])
        #Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set,k in K] ≤ THETAMAX)

        JuMP.@variable(m, GenData[g, :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[g, :Pmax][1,1])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set,k in K])

        JuMP.@variable(m, z[l=System.LineData[:ID],k in K;System.LineData[l,:is_aux]==1.0], Bin)

        #Objective Function
        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g]+GenData[g,:C0][1,1] for g in Gen_set) )

        #Constraints

        #Current law applied on auxilliary buses
        # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
        #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

        # JuMP.@constraint(m, ReferenceAngle[k in K],
        #     (δ[Nodes_set[1],k] ==  0.0))
        #Current law applied on all buses
        JuMP.@constraint(m, Nodal_Balance[i = Nodes_set,k in K],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

        #Power flows in non auxilliary lines
        JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set,k in K ; Set([i,j]) in non_aux_line_set ],
            pij[i,j,k] == a[i,j,k]*Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i,k]-δ[j,k]))

        JuMP.@constraint(m,total_balance, sum(p) == sum(System.BusData[:,:Pd]))
        #Switching constraint
        JuMP.@constraint(m,connectors[i = Nodes_set,k in K; BusData[i,:is_aux] == 1.0],
            sum(z[l,k] for l in vcat(findall(x->x==i,System.LineData[:,:fbus]),findall(x->x==i,System.LineData[:,:tbus])) if System.LineData[l,:is_aux] == 1.0) == 1)

        #Theta Constraints
        JuMP.@constraint(m,theta1[l = System.LineData[:,:ID],k in K ;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] <= δ[System.LineData[l,:tbus],k] + (1-z[l,k])*M2)

        JuMP.@constraint(m,theta2[l = System.LineData[:,:ID] ,k in K;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] + (1-z[l,k])*M2 >= δ[System.LineData[l,:tbus],k])

        #Capacity Constraints

        JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID],k in K],
            pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i,k] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
             pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= System.LineData[i,:rate]*Sbase*z[i,k])

        JuMP.@constraint(m,Line_Capacity_aux2_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i,k] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= System.LineData[i,:rate]*Sbase*z[i,k])

        JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 0.0],
            -a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase)

        JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ,k in K ;System.LineData[i,:is_aux] == 0.0],
            -a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase)

        if warm_start
            for k in K
                for aux_line in System.LineData[:ID]
                    if System.LineData[aux_line,:is_aux]==1.0
                        fbus = System.LineData[aux_line,:fbus]
                        if fbus in root_buses
                            JuMP.setvalue(z[aux_line,k],1)
                        else
                            JuMP.setvalue(z[aux_line,k],0)
                        end
                    end
                end
            end
        end
        optimize!(m)

        p_ij = JuMP.value.(pij)
        pg = JuMP.value.(p)
        z = JuMP.value.(z)
        δ_sol = JuMP.value.(δ)
        obj = JuMP.objective_value(m)

        switching_status = zeros(length(System.LineData[System.LineData[:is_aux].==1.0,:ID]),2)
        switching_status[:,1] = System.LineData[System.LineData[:is_aux].==1.0,:ID]
        switching_status[:,2] = [z[i,1] for i in System.LineData[System.LineData[:is_aux].==1.0,:ID]]

        #========================#
        Pg = JuMP.value.(p)
        Pij = JuMP.value.(pij)
        δ = JuMP.value.(δ)
        V = ones(N)
        Pd = System.BusData[!,3]
        Qd = System.BusData[!,4]
        System.Operating_Cost = JuMP.objective_value(m)
        System.Gen_Data[:Pg] = Pg.data
        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1])]) : 0 for i in Nodes_set]
        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
         Pg = Pg , Pd = Pd)
        System.BusData_output = df

        Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
        Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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
        System.LineLoading[:ID] = deepcopy(System.LineData[:ID])
        #========================#
        reduce_System!(System,switching_status)
        removed_lines = sort(collect(setdiff(Set(System.LineLoading[:ID]),System.LineData[:ID])))
        deleterows!(System.LineLoading,removed_lines)
        System.LineLoading[:FromBus] = deepcopy(System.LineData[:fbus])
        System.LineLoading[:ToBus] = deepcopy(System.LineData[:tbus])

        #System.Gen_Data[:Pg] = pg.data
        return obj,System,z.data,pg.data
    elseif method == "OBS-P-2"
        subs_ = [k for k in keys(System.SubStations)]
        root_buses = [parse(Int,subs_[i][12:length(subs_[i])]) for i in 1:length(subs_)]
        a,K_ = Generate_Contingencies(System,"P-SCOPF");
        K = 1:K_;
        Gen_set = 1:length(System.Gen_Data[:,:bus])
        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,TimeLimit=timelimit))

        #Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set,k in K] ≤ THETAMAX)

        JuMP.@variable(m, GenData[g, :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[g, :Pmax][1,1])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set,k in K])

        JuMP.@variable(m, z[l=System.LineData[:ID] ;System.LineData[l,:is_aux]==1.0], Bin)

        JuMP.@variable(m, z_coupler[l=System.LineData[:ID],k in K ;System.LineData[l,:is_aux]==2.0], Bin)
        #Objective Function
        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g]+GenData[g,:C0][1,1] for g in Gen_set) )

        #Constraints

        #Current law applied on auxilliary buses
        # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
        #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

        # JuMP.@constraint(m, ReferenceAngle[k in K],
        #     (δ[Nodes_set[1],k] ==  0.0))
        #Current law applied on all buses
        JuMP.@constraint(m, Nodal_Balance[i = Nodes_set,k in K],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

        #Power flows in non auxilliary lines
        JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set,k in K ; Set([i,j]) in non_aux_line_set ],
            pij[i,j,k] == a[i,j,k]*Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i,k]-δ[j,k]))

        JuMP.@constraint(m,total_balance, sum(p) == sum(System.BusData[:,:Pd]))

        #Switching constraint
        JuMP.@constraint(m,connectors[i = Nodes_set; BusData[i,:is_aux] == 1.0],
            sum(z[l] for l in vcat(findall(x->x==i,System.LineData[:,:fbus]),findall(x->x==i,System.LineData[:,:tbus])) if System.LineData[l,:is_aux] == 1.0) == 1)

        JuMP.@constraint(m,couplers[l = System.LineData[:ID]; System.LineData[l,:is_aux]==2.0],
            z_coupler[l,1] == 1.0)
        #Theta Constraints


        JuMP.@constraint(m,theta1[l = System.LineData[:,:ID],k in K ;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] <= δ[System.LineData[l,:tbus],k] + (1-z[l])*M2)

        JuMP.@constraint(m,theta2[l = System.LineData[:,:ID] ,k in K;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] + (1-z[l])*M2 >= δ[System.LineData[l,:tbus],k])

        JuMP.@constraint(m,theta3[l = System.LineData[:,:ID],k in K ;System.LineData[l,:is_aux] == 2.0],
            δ[System.LineData[l,:fbus],k] <= δ[System.LineData[l,:tbus],k] + (1-z_coupler[l,k])*M2)

        JuMP.@constraint(m,theta4[l = System.LineData[:,:ID] ,k in K;System.LineData[l,:is_aux] == 2.0],
            δ[System.LineData[l,:fbus],k] + (1-z_coupler[l,k])*M2 >= δ[System.LineData[l,:tbus],k])

        #Capacity Constraints

        JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID],k in K],
            pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        #Auxilliary lines constraints
        JuMP.@constraint(m,Line_Capacity_aux1_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
             pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= System.LineData[i,:rate]*Sbase*z[i])

        JuMP.@constraint(m,Line_Capacity_aux2_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= System.LineData[i,:rate]*Sbase*z[i])

        # Coupler constraints
        JuMP.@constraint(m,Line_Capacity_coupler1_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 2.0],
            -System.LineData[i,:rate]*Sbase*z_coupler[i,k] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k])

        JuMP.@constraint(m,Line_Capacity_coupler1_2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 2.0],
             pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= System.LineData[i,:rate]*Sbase*z_coupler[i,k])

        JuMP.@constraint(m,Line_Capacity_coupler2_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 2.0],
            -System.LineData[i,:rate]*Sbase*z_coupler[i,k] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_coupler2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 2.0],
            pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= System.LineData[i,:rate]*Sbase*z_coupler[i,k])

        #Non-auxilliary lines constraints
        JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 0.0],
            -a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase)

        JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ,k in K ;System.LineData[i,:is_aux] == 0.0],
            -a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase)

            if warm_start
                for aux_line in System.LineData[:ID]
                    if System.LineData[aux_line,:is_aux]==1.0
                        fbus = System.LineData[aux_line,:fbus]
                        if fbus in root_buses
                            JuMP.setvalue(z[aux_line],1)
                        else
                            JuMP.setvalue(z[aux_line],0)
                        end
                    end
                end
            end
        optimize!(m)

        p_ij = JuMP.value.(pij)
        pg = JuMP.value.(p)
        z = JuMP.value.(z)
        z_coupler = JuMP.value.(z_coupler)
        δ_sol = JuMP.value.(δ)
        obj = JuMP.objective_value(m)

        switching_status = zeros(length(System.LineData[System.LineData[:is_aux].==1.0,:ID]),2)
        switching_status[:,1] = System.LineData[System.LineData[:is_aux].==1.0,:ID]
        switching_status[:,2] = [z[i] for i in System.LineData[System.LineData[:is_aux].==1.0,:ID]]

        switching_status_coupler = zeros(length(System.LineData[System.LineData[:is_aux].==2.0,:ID]),2)
        switching_status_coupler[:,1] = System.LineData[System.LineData[:is_aux].==2.0,:ID]
        switching_status_coupler[:,2] = [z_coupler[i,1] for i in System.LineData[System.LineData[:is_aux].==2.0,:ID]]

        #========================#
        Pg = JuMP.value.(p)
        Pij = JuMP.value.(pij)
        δ = JuMP.value.(δ)
        V = ones(N)
        Pd = System.BusData[!,3]
        Qd = System.BusData[!,4]
        System.Operating_Cost = JuMP.objective_value(m)
        System.Gen_Data[:Pg] = Pg.data
        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1])]) : 0 for i in Nodes_set]
        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
         Pg = Pg , Pd = Pd)
        System.BusData_output = df

        Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
        Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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
        System.LineLoading[:ID] = deepcopy(System.LineData[:ID])
        #========================#
        reduce_System_coupled!(System,switching_status,switching_status_coupler)
        removed_lines = sort(collect(setdiff(Set(System.LineLoading[:ID]),System.LineData[:ID])))
        deleterows!(System.LineLoading,removed_lines)
        System.LineLoading[:FromBus] = deepcopy(System.LineData[:fbus])
        System.LineLoading[:ToBus] = deepcopy(System.LineData[:tbus])

        #System.Gen_Data[:Pg] = pg.data
        return obj,System,z.data,pg.data, z_coupler.data
    elseif method == "C-SCOPF"
        a_g,a_l,K_ = Generate_Contingencies(System,method);
        K = 1:K_;
        Gen_set = 1:length(System.Gen_Data[:,:bus]);
        #Gen_buses = System.Gen_Data[:,:bus]
        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,TimeLimit=timelimit))

        # 2.1.Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set,k in K] ≤ THETAMAX)

        JuMP.@variable(m, a_g[g,1,k]*GenData[g, :Pmin][1,1] ≤ p[g in Gen_set,k in K] ≤ a_g[g,1,k]*GenData[g, :Pmax][1,1])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set,k in K])

        JuMP.@variable(m, z[l=System.LineData[:ID];System.LineData[l,:is_aux]==1.0], Bin)

        #Objective Function
        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g,1]+GenData[g,:C0][1,1] for g in Gen_set) )

        #Constraints

        #Current law applied on auxilliary buses
        # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
        #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

        # JuMP.@constraint(m, ReferenceAngle[k in K],
        #     (δ[Nodes_set[1],k] ==  0.0))
        #Current law applied on all buses
        JuMP.@constraint(m, Nodal_Balance[i = Nodes_set,k in K],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g,k] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

        #Power flows in non auxilliary lines
        JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set,k in K ; Set([i,j]) in non_aux_line_set ],
            pij[i,j,k] == a_l[i,j,k]*Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i,k]-δ[j,k]))

        #JuMP.@constraint(m,total_balance[k in K], sum(p[g,k] for g in Gen_set) == sum(System.BusData[:,:Pd]))
        #Switching constraint
        JuMP.@constraint(m,connectors[i = Nodes_set; BusData[i,:is_aux] == 1.0],
            sum(z[l] for l in vcat(findall(x->x==i,System.LineData[:,:fbus]),findall(x->x==i,System.LineData[:,:tbus])) if System.LineData[l,:is_aux] == 1.0) == 1)

        #Theta Constraints
        JuMP.@constraint(m,theta1[l = System.LineData[:,:ID],k in K ;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] <= δ[System.LineData[l,:tbus],k] + (1-z[l])*M2)

        JuMP.@constraint(m,theta2[l = System.LineData[:,:ID] ,k in K;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] + (1-z[l])*M2 >= δ[System.LineData[l,:tbus],k])

        JuMP.@constraint(m,Lim_Gen_up[i in Gen_set,k = 2:length(K)],p[i,k]-p[i,1] <= Delta_P[i])
        JuMP.@constraint(m,Lim_Gen_Dw[i in Gen_set,k = 2:length(K)],p[i,1]-p[i,k] <= Delta_P[i])
        #Capacity Constraints

        JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID],k in K],
            pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
             pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= System.LineData[i,:rate]*Sbase*z[i])

        # JuMP.@constraint(m,Line_Capacity_aux2_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
        #     -System.LineData[i,:rate]*Sbase*z[i] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])
        #
        # JuMP.@constraint(m,Line_Capacity_aux2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
        #     pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= System.LineData[i,:rate]*Sbase*z[i])

        JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 0.0],
            -a_l[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= a_l[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase)

        JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ,k in K ;System.LineData[i,:is_aux] == 0.0],
            -a_l[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= a_l[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase)

            if warm_start
                for aux_line in System.LineData[:ID]
                    if System.LineData[aux_line,:is_aux]==1.0
                        fbus = System.LineData[aux_line,:fbus]
                        if fbus in root_buses
                            JuMP.setvalue(z[aux_line],1)
                        else
                            JuMP.setvalue(z[aux_line],0)
                        end
                    end
                end
            end

        optimize!(m)

        p_ij = JuMP.value.(pij)
        pg = JuMP.value.(p)
        z = JuMP.value.(z)
        δ_sol = JuMP.value.(δ)
        obj = JuMP.objective_value(m)

        switching_status = zeros(length(System.LineData[System.LineData[:is_aux].==1.0,:ID]),2)
        switching_status[:,1] = System.LineData[System.LineData[:is_aux].==1.0,:ID]
        switching_status[:,2] = [z[i] for i in System.LineData[System.LineData[:is_aux].==1.0,:ID]]
        #========================#
        Pg = JuMP.value.(p)
        Pij = JuMP.value.(pij)
        δ = JuMP.value.(δ)
        V = ones(N)
        Pd = System.BusData[!,3]
        Qd = System.BusData[!,4]
        System.Operating_Cost = JuMP.objective_value(m)
        System.Gen_Data[:Pg] = Pg.data[:,1]
        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1]),1]) : 0 for i in Nodes_set]
        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
         Pg = Pg , Pd = Pd)
        System.BusData_output = df

        Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
        Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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
        System.LineLoading[:ID] = deepcopy(System.LineData[:ID])
        #========================#
        reduce_System!(System ::System_Struct,switching_status)
        removed_lines = sort(collect(setdiff(Set(System.LineLoading[:ID]),System.LineData[:ID])))
        deleterows!(System.LineLoading,removed_lines)
        System.LineLoading[:FromBus] = deepcopy(System.LineData[:fbus])
        System.LineLoading[:ToBus] = deepcopy(System.LineData[:tbus])

        return obj,System,z.data,pg.data
    elseif method == "OBS-C-1"
        a_g,a_l,K_ = Generate_Contingencies(System,"C-SCOPF");
        K = 1:K_;
        Gen_set = 1:length(System.Gen_Data[:,:bus]);
        #Gen_buses = System.Gen_Data[:,:bus]
        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,TimeLimit=timelimit))

        # 2.1.Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set,k in K] ≤ THETAMAX)

        JuMP.@variable(m, a_g[g,1,k]*GenData[g, :Pmin][1,1] ≤ p[g in Gen_set,k in K] ≤ a_g[g,1,k]*GenData[g, :Pmax][1,1])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set,k in K])

        JuMP.@variable(m, z[l=System.LineData[:ID],k in K;System.LineData[l,:is_aux]==1.0], Bin)

        #Objective Function
        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g,1]+GenData[g,:C0][1,1] for g in Gen_set) )

        #Constraints

        #Current law applied on auxilliary buses
        # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
        #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

        # JuMP.@constraint(m, ReferenceAngle[k in K],
        #     (δ[Nodes_set[1],k] ==  0.0))
        #Current law applied on all buses
        JuMP.@constraint(m, Nodal_Balance[i = Nodes_set,k in K],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g,k] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

        #Power flows in non auxilliary lines
        JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set,k in K ; Set([i,j]) in non_aux_line_set ],
            pij[i,j,k] == a_l[i,j,k]*Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i,k]-δ[j,k]))

        #JuMP.@constraint(m,total_balance[k in K], sum(p[g,k] for g in Gen_set) == sum(System.BusData[:,:Pd]))
        #Switching constraint
        JuMP.@constraint(m,connectors[i = Nodes_set,k in K; BusData[i,:is_aux] == 1.0],
            sum(z[l,k] for l in vcat(findall(x->x==i,System.LineData[:,:fbus]),findall(x->x==i,System.LineData[:,:tbus])) if System.LineData[l,:is_aux] == 1.0) == 1)

        #Theta Constraints
        JuMP.@constraint(m,theta1[l = System.LineData[:,:ID],k in K ;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] <= δ[System.LineData[l,:tbus],k] + (1-z[l,k])*M2)

        JuMP.@constraint(m,theta2[l = System.LineData[:,:ID] ,k in K;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] + (1-z[l,k])*M2 >= δ[System.LineData[l,:tbus],k])

        JuMP.@constraint(m,Lim_Gen_up[i in Gen_set,k = 2:length(K)],p[i,k]-p[i,1] <= Delta_P[i])
        JuMP.@constraint(m,Lim_Gen_Dw[i in Gen_set,k = 2:length(K)],p[i,1]-p[i,k] <= Delta_P[i])
        #Capacity Constraints

        JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID],k in K],
            pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i,k] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
             pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= System.LineData[i,:rate]*Sbase*z[i,k])

        JuMP.@constraint(m,Line_Capacity_aux2_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i,k] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= System.LineData[i,:rate]*Sbase*z[i,k])

        JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 0.0],
            -a_l[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= a_l[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase)

        JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ,k in K ;System.LineData[i,:is_aux] == 0.0],
            -a_l[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= a_l[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase)

            if warm_start
                for k in K
                    for aux_line in System.LineData[:ID]
                        if System.LineData[aux_line,:is_aux]==1.0
                            fbus = System.LineData[aux_line,:fbus]
                            if fbus in root_buses
                                JuMP.setvalue(z[aux_line,k],1)
                            else
                                JuMP.setvalue(z[aux_line,k],0)
                            end
                        end
                    end
                end
            end
        optimize!(m)

        p_ij = JuMP.value.(pij)
        pg = JuMP.value.(p)
        z = JuMP.value.(z)
        δ_sol = JuMP.value.(δ)
        obj = JuMP.objective_value(m)

        switching_status = zeros(length(System.LineData[System.LineData[:is_aux].==1.0,:ID]),2)
        switching_status[:,1] = System.LineData[System.LineData[:is_aux].==1.0,:ID]
        switching_status[:,2] = [z[i,1] for i in System.LineData[System.LineData[:is_aux].==1.0,:ID]]
        #println(switching_status)
        #========================#
        Pg = JuMP.value.(p)
        Pij = JuMP.value.(pij)
        δ = JuMP.value.(δ)
        V = ones(N)
        Pd = System.BusData[!,3]
        Qd = System.BusData[!,4]
        System.Operating_Cost = JuMP.objective_value(m)
        System.Gen_Data[:Pg] = Pg.data[:,1]
        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1]),1]) : 0 for i in Nodes_set]
        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
         Pg = Pg , Pd = Pd)
        System.BusData_output = df

        Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
        Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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
        System.LineLoading[:ID] = deepcopy(System.LineData[:ID])
        #========================#
        reduce_System!(System,switching_status)
        removed_lines = sort(collect(setdiff(Set(System.LineLoading[:ID]),System.LineData[:ID])))
        deleterows!(System.LineLoading,removed_lines)
        System.LineLoading[:FromBus] = deepcopy(System.LineData[:fbus])
        System.LineLoading[:ToBus] = deepcopy(System.LineData[:tbus])

        return obj,System,z.data,pg.data
    elseif method == "OBS-C-2"

        subs_ = [k for k in keys(System.SubStations)]
        root_buses = [parse(Int,subs_[i][12:length(subs_[i])]) for i in 1:length(subs_)]
        a_g,a_l,K_ = Generate_Contingencies(System,"C-SCOPF");
        K = 1:K_;
        Gen_set = 1:length(System.Gen_Data[:,:bus]);
        #Gen_buses = System.Gen_Data[:,:bus]
        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,TimeLimit=timelimit))

        # 2.1.Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set,k in K] ≤ THETAMAX)

        JuMP.@variable(m, a_g[g,1,k]*GenData[g, :Pmin][1,1] ≤ p[g in Gen_set,k in K] ≤ a_g[g,1,k]*GenData[g, :Pmax][1,1])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set,k in K])

        JuMP.@variable(m, z[l=System.LineData[:ID],k in K;System.LineData[l,:is_aux]==1.0], Bin)

        #Objective Function
        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g,1]+GenData[g,:C0][1,1] for g in Gen_set) )

        #Constraints

        #Current law applied on auxilliary buses
        # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
        #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

        # JuMP.@constraint(m, ReferenceAngle[k in K],
        #     (δ[Nodes_set[1],k] ==  0.0))
        #Current law applied on all buses
        JuMP.@constraint(m, Nodal_Balance[i = Nodes_set,k in K],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g,k] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

        #Power flows in non auxilliary lines
        JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set,k in K ; Set([i,j]) in non_aux_line_set ],
            pij[i,j,k] == a_l[i,j,k]*Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i,k]-δ[j,k]))

        #JuMP.@constraint(m,total_balance[k in K], sum(p[g,k] for g in Gen_set) == sum(System.BusData[:,:Pd]))
        #Switching constraint
        JuMP.@constraint(m,precontingency[l=System.LineData[:ID]; (System.LineData[l,:is_aux] == 1.0) && (System.LineData[l,:fbus] in root_buses)],
            z[l,1] == 1)
        JuMP.@constraint(m,precontingency2[l=System.LineData[:ID]; (System.LineData[l,:is_aux] == 1.0) && !(System.LineData[l,:fbus] in root_buses)],
            z[l,1] == 0)
        JuMP.@constraint(m,connectors[i = Nodes_set,k in K; BusData[i,:is_aux] == 1.0],
            sum(z[l,k] for l in vcat(findall(x->x==i,System.LineData[:,:fbus]),findall(x->x==i,System.LineData[:,:tbus])) if System.LineData[l,:is_aux] == 1.0) == 1)

        #Theta Constraints
        JuMP.@constraint(m,theta1[l = System.LineData[:,:ID],k in K ;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] <= δ[System.LineData[l,:tbus],k] + (1-z[l,k])*M2)

        JuMP.@constraint(m,theta2[l = System.LineData[:,:ID] ,k in K ;System.LineData[l,:is_aux] == 1.0],
            δ[System.LineData[l,:fbus],k] + (1-z[l,k])*M2 >= δ[System.LineData[l,:tbus],k])

        JuMP.@constraint(m,Lim_Gen_up[i in Gen_set,k = 2:length(K)],p[i,k]-p[i,1] <= Delta_P[i])
        JuMP.@constraint(m,Lim_Gen_Dw[i in Gen_set,k = 2:length(K)],p[i,1]-p[i,k] <= Delta_P[i])
        #Capacity Constraints

        JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID],k in K],
            pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_1[i = System.LineData[:,:ID],k in K;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i,k] <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k])

        JuMP.@constraint(m,Line_Capacity_aux1_2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
             pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= System.LineData[i,:rate]*Sbase*z[i,k])

        JuMP.@constraint(m,Line_Capacity_aux2_1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            -System.LineData[i,:rate]*Sbase*z[i,k] <= pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_aux2[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 1.0],
            pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= System.LineData[i,:rate]*Sbase*z[i,k])

        JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID],k in K ;System.LineData[i,:is_aux] == 0.0],
            -a_l[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= a_l[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase)

        JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ,k in K ;System.LineData[i,:is_aux] == 0.0],
            -a_l[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= a_l[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase)

        optimize!(m)

        if warm_start
            for k in K[2:length(K)]
                for aux_line in System.LineData[:ID]
                    if System.LineData[aux_line,:is_aux]==1.0
                        fbus = System.LineData[aux_line,:fbus]
                        if fbus in root_buses
                            JuMP.setvalue(z[aux_line,k],1)
                        else
                            JuMP.setvalue(z[aux_line,k],0)
                        end
                    end
                end
            end
        end
        p_ij = JuMP.value.(pij)
        pg = JuMP.value.(p)
        z = JuMP.value.(z)
        δ_sol = JuMP.value.(δ)
        obj = JuMP.objective_value(m)

        switching_status = zeros(length(System.LineData[System.LineData[:is_aux].==1.0,:ID]),2)
        switching_status[:,1] = System.LineData[System.LineData[:is_aux].==1.0,:ID]
        switching_status[:,2] = [z[i,1] for i in System.LineData[System.LineData[:is_aux].==1.0,:ID]]
        #println(switching_status)
        #========================#
        Pg = JuMP.value.(p)
        Pij = JuMP.value.(pij)
        δ = JuMP.value.(δ)
        V = ones(N)
        Pd = System.BusData[!,3]
        Qd = System.BusData[!,4]
        System.Operating_Cost = JuMP.objective_value(m)
        System.Gen_Data[:Pg] = Pg.data[:,1]
        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x -> x==i,System.Gen_Data[!,1]),1]) : 0 for i in Nodes_set]
        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
         Pg = Pg , Pd = Pd)
        System.BusData_output = df

        Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
        Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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
        System.LineLoading[:ID] = deepcopy(System.LineData[:ID])
        #========================#
        reduce_System!(System ::System_Struct,switching_status)
        removed_lines = sort(collect(setdiff(Set(System.LineLoading[:ID]),System.LineData[:ID])))
        deleterows!(System.LineLoading,removed_lines)
        System.LineLoading[:FromBus] = deepcopy(System.LineData[:fbus])
        System.LineLoading[:ToBus] = deepcopy(System.LineData[:tbus])

        return obj,System,z.data,pg.data
    end
end

function Generate_Contingencies(System :: System_Struct,  type ::String,Contingency_Order = 1)

    N_Line = size(System.Line_Constraints,1);
    N_Gen = size(System.Gen_Data,1);
    N_Nodes = System.N_bus;
    contData = System.Line_Constraints;

    if type == "P-SCOPF"

        k = N_Line + 1;
        a = ones(N_Nodes,N_Nodes,k)
        for i in 1:N_Nodes , j in 1:N_Nodes , c in 2:k
            if (contData[c-1,:fbus] == i) && (contData[c-1,:tbus] ==j)
                a[i,j,c] = 0
                a[j,i,c] = 0
            end
        end
            return a,k

    elseif type == "C-SCOPF"

        k = N_Line + N_Gen + 1;
        a_l = ones(N_Nodes,N_Nodes,k)
        a_g = ones(N_Gen,1,k)

        for i in 1:N_Nodes , j in 1:N_Nodes , c in 2:N_Line
            if (contData[c-1,:fbus] == i) && (contData[c-1,:tbus] ==j)
                a_l[i,j,c] = 0
                a_l[j,i,c] = 0
            end
        end

        i = 1;
        for c in N_Line+2:k
            a_g[i,1,c] = 0;
            i = i+1;
        end

        return a_g,a_l,k
    end
end

function solve_ED(System ::System_Struct)
    Gen_set = 1:length(System.Gen_Data[:,:bus]);
    GenData = System.Gen_Data;
    m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,OutputFlag=1))
    JuMP.@variable(m, GenData[g, :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[g, :Pmax][1,1])
    JuMP.@constraint(m,Market_Clearing,sum(p[g] for g in Gen_set) == sum(System.BusData[:Pd]))
    JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g]+GenData[g,:C0][1,1] for g in Gen_set) )
    JuMP.optimize!(m)

    return JuMP.objective_value(m),JuMP.value.(p).data,JuMP.dual(Market_Clearing)
end

function Solve_SCOPF!(System ::System_Struct,method,Delta_P=0)

    THETAMAX = 0.6
    M = 10000
    M2 = 2 * THETAMAX
    G = real(System.y_bus);
    B = imag(System.y_bus);
    b = System.b_line;

    N = System.N_bus
    N_Lines = size(System.Line_Constraints,1)
    Nodes_set = System.BusData[:,:bus_i];
    Gen_set = 1:length(System.Gen_Data[:,:bus]);

    Sbase = System.Sbase;

    BusData = System.BusData;
    GenData = System.Gen_Data;
    BranchData = System.Line_Constraints;

    df_lines = System.LineData[:,[:fbus,:tbus]]
    lines_set = Set.(convert.(Array, row for row in eachrow(df_lines)))

    df_aux = System.LineData[findall(x->x==1.0,System.LineData[:is_aux]),[:fbus,:tbus]]
    aux_line_set = Set.(convert.(Array, row for row in eachrow(df_aux)))

    df_non_aux = System.LineData[findall(x->x==0.0,System.LineData[:is_aux]),[:fbus,:tbus]]
    non_aux_line_set = Set.(convert.(Array, row for row in eachrow(df_non_aux)))

    Nodes_at_Lines = Dict()
    for i in System.LineData[:ID]
        Nodes_at_Lines[i] = Set([System.LineData[i,:fbus],System.LineData[i,:tbus]])
    end

    #B_ = -B[1:N,1:N];

    if method == "P-SCOPF"
        a,K = Generate_Contingencies(System,method);
        K = 1:K;
        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,OutputFlag=1))

        #Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set,k in K] ≤ THETAMAX)

        JuMP.@variable(m, GenData[g, :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[g, :Pmax][1,1])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set,k in K])

        #Objective
        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g]+GenData[g,:C0][1,1] for g in Gen_set) )

        #Constraints:

        # JuMP.@constraint(m, ReferenceAngle[k in K],
        #     (δ[Nodes_set[1],k] ==  0.0))
        #Current law applied on all buses
        JuMP.@constraint(m, Nodal_Balance[i = Nodes_set,k in K],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

        #Power flows in non auxilliary lines
        JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set,k in K;Set([i,j]) in lines_set ],
            pij[i,j,k] == a[i,j,k]*Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i,k]-δ[j,k]))

        #JuMP.@constraint(m,total_balance[k in K], sum(p[:,k]) == sum(System.BusData[:,:Pd]))

        #Capacity Constraints

        JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID],k in K],
            pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])

        JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID],k in K],
            -a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= a[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase)

        JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ,k in K],
            -a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= a[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase)

        optimize!(m)
        p_ij = JuMP.value.(pij)
        pg = JuMP.value.(p)
        δ_sol = JuMP.value.(δ)
        obj = JuMP.objective_value(m)

        #========================#
        Pg = JuMP.value.(p)
        Pij = JuMP.value.(pij)
        δ = JuMP.value.(δ)
        V = ones(N)
        Pd = System.BusData[!,3]
        Qd = System.BusData[!,4]
        System.Operating_Cost = JuMP.objective_value(m)
        System.Gen_Data[:Pg] = Pg.data
        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x->x == i,GenData[:,:bus]),1]) : 0 for i in Nodes_set]
        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
         Pg = Pg , Pd = Pd)
        System.BusData_output = df

        Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
        Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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
        return obj,System


    elseif method == "C-SCOPF"
        a_g,a_l,K_ = Generate_Contingencies(System,method);
        K = 1:K_;
        Gen_set = 1:length(System.Gen_Data[:,:bus]);
        #Gen_buses = System.Gen_Data[:,:bus]
        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,OutputFlag=1))

        # 2.1.Variables
        JuMP.@variable(m, -THETAMAX ≤ δ[i in Nodes_set,k in K] ≤ THETAMAX)

        JuMP.@variable(m, a_g[g,1,k]*GenData[g, :Pmin][1,1] ≤ p[g in Gen_set,k in K] ≤ a_g[g,1,k]*GenData[g, :Pmax][1,1])

        JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set,k in K])

        #Objective Function
        JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g,1]+GenData[g,:C0][1,1] for g in Gen_set) )

        #Constraints

        #Current law applied on auxilliary buses
        # JuMP.@constraint(m, Nodal_Balance[i = BusData[findall(x->x==1.0,System.BusData[:,:is_aux]),:bus_i]],
        #     sum(pij[i,j] for j = Nodes_set if [i,j] in lines_set || [j,i] in lines_set ) == sum(p[g] for g in Gen_set if GenData[findall(x->x==g,System.Gen_Data[:,:bus]),:bus][1,1]==i) - BusData[i,:Pd])

        # JuMP.@constraint(m, ReferenceAngle[k in K],
        #     (δ[Nodes_set[1],k] ==  0.0))
        #Current law applied on all buses
        JuMP.@constraint(m, Nodal_Balance[i = Nodes_set,k in K],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in lines_set ) == sum(p[g,k] for g in Gen_set if GenData[g,:bus][1,1]==i) - BusData[i,:Pd])

        #Power flows in non auxilliary lines
        JuMP.@constraint(m,Pl[i in Nodes_set,j in Nodes_set,k in K ; Set([i,j]) in lines_set ],
            pij[i,j,k] == a_l[i,j,k]*Sbase*(1/(System.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])],:x][1,1]))*(δ[i,k]-δ[j,k]))

        JuMP.@constraint(m,Lim_Gen_up[i in Gen_set,k = 2:length(K)],p[i,k]-p[i,1] <= Delta_P[i])
        JuMP.@constraint(m,Lim_Gen_Dw[i in Gen_set,k = 2:length(K)],p[i,1]-p[i,k] <= Delta_P[i])
        #Capacity Constraints

        JuMP.@constraint(m,lines_consistency[i in System.LineData[:,:ID],k in K],
            pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] == -pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k])


        JuMP.@constraint(m,Line_Capacity_non_aux1[i = System.LineData[:,:ID],k in K ],
            -a_l[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase <= pij[System.LineData[i,:fbus],System.LineData[i,:tbus],k] <= a_l[System.LineData[i,:fbus],System.LineData[i,:tbus],k]*System.LineData[i,:rate]*Sbase)

        JuMP.@constraint(m,Line_Capacity_non_aux2[i = System.LineData[:,:ID] ,k in K ],
            -a_l[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase <=pij[System.LineData[i,:tbus],System.LineData[i,:fbus],k] <= a_l[System.LineData[i,:tbus],System.LineData[i,:fbus],k]*System.LineData[i,:rate]*Sbase)

        optimize!(m)

        p_ij = JuMP.value.(pij)
        pg = JuMP.value.(p)
        δ_sol = JuMP.value.(δ)
        obj = JuMP.objective_value(m)

        Pg = JuMP.value.(p)
        Pij = JuMP.value.(pij)
        δ = JuMP.value.(δ)
        V = ones(N)
        Pd = System.BusData[!,3]
        Qd = System.BusData[!,4]
        System.Operating_Cost = JuMP.objective_value(m)
        System.Gen_Data[:Pg] = Pg.data[:,1]
        Pg = [i in System.Gen_Data[!,1] ? sum(Pg.data[findall(x->x == i,GenData[:,:bus]),1]) : 0 for i in Nodes_set]
        df = DataFrame(Bus = System.BusData[!,1], V = V, δ = δ.data[:,1],
         Pg = Pg , Pd = Pd)
        System.BusData_output = df

        Pij_1 = [ Pij[System.LineData[i,:fbus],System.LineData[i,:tbus],1] for i in 1:size(System.LineData,1) ]
        Pij_2 = [ Pij[System.LineData[i,:tbus],System.LineData[i,:fbus],1] for i in 1:size(System.LineData,1) ]
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


        return obj,System
    end
end

function Verify_SCP_PF(System ::System_Struct)

    B = imag(System.y_bus);
    N = System.N_bus
    N_Lines = size(System.LineData[:,1],1)
    Sbase = System.Sbase;
    BranchData = System.LineData;
    Pg = [ i in System.Gen_Data[:bus] ? sum(System.Gen_Data[findall(x->x==i,System.Gen_Data[:bus]),:Pg]./Sbase) : 0 for i in 1:System.N_bus]
    Pd = System.BusData[:Pd]/Sbase
    Pnet = Pg-Pd
    flag = 0

    #Line removal loop
    for i in 1:N_Lines
        BB = deepcopy(B)
        fbus = BranchData[i,:fbus]
        tbus = BranchData[i,:tbus]
        BB[fbus,tbus] = 0
        BB[tbus,fbus] = 0
        B_ = -BB[2:N,2:N]
        δ = inv(B_)*Pnet[2:N]
        δ = append!([0.0],δ)
        PL = zeros(N_Lines,1)*1.0
        #line flow loop
        for j in 1:N_Lines
            x = BranchData[j,:x]
            fbus_ = BranchData[j,:fbus]
            tbus_ = BranchData[j,:tbus]
            PL[j] = abs((1/x)*(δ[fbus_]-δ[tbus_]))
            println([i,j,PL[j],BranchData[j,:rate]])
            if PL[j] > BranchData[j,:rate]
                flag = 1
            end
        end

        if flag == 1
            println("Failed!!!")
            break
        end

    end
    if flag != 1
        println("Success!")
    end
end

function Verify_SCC_OPF(System ::System_Struct,Δg)
    #OPF solver

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

    BranchData = deepcopy(System.LineData);

    df_lines = BranchData[:,[:fbus,:tbus]]
    lines_set = Set.(convert.(Array, row for row in eachrow(df_lines)))
    lines_rates = Dict()
    for i in 1:length(BranchData[:,1])
        lines_rates[lines_set[i]]=BranchData[i,:rate]
    end
    Nodes_at_Lines = Dict()
    for i in BranchData[:ID]
        Nodes_at_Lines[i] = Set([BranchData[i,:fbus],BranchData[i,:tbus]])
    end
    B_ = -B;

    Pmax = [minimum(vcat(GenData[g,:Pmax],GenData[g,:Pg]+Δg[g])) for g in Gen_set]
    Pmin = [maximum(vcat(GenData[g,:Pmin],GenData[g,:Pg]-Δg[g])) for g in Gen_set]

    flag = 0
    ag,a,K = Generate_Contingencies(System,"C-SCOPF")

    for k in 1:K


        #JuMP.with_optimizer(Gurobi.Optimizer,OutputFlag=0),quiet=true

        m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,OutputFlag=0))


            # 2.1.Variables
            JuMP.@variable(m,-0.6 <= δ[i in Nodes_set] <= 0.6)

            JuMP.@variable(m, Pmin[g]*ag[g,1,k] ≤ p[g in Gen_set] ≤ Pmax[g]*ag[g,1,k])

            JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set ; Set([i,j]) in lines_set])

            # 2.2. Constraints
            JuMP.@constraint(m, ReferenceAngle,
                (δ[Nodes_set[1]] ==  0.0))

            JuMP.@constraint(m,Nodal_balance[i in Nodes_set],
                    sum(pij[i,j] for j = Nodes_set if Set([i,j]) in lines_set) == sum(p[g] for g in Gen_set if System.Gen_Data[g,:bus][1,1] == i) - BusData[i,:Pd])
            JuMP.@constraint(m,pl[i in Nodes_set,j in Nodes_set; Set([i,j]) in lines_set],
                    pij[i,j] == a[i,j,k]*Sbase*(1/(BranchData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])][1],:x]))*(δ[i]-δ[j]))

            JuMP.@constraint(m,pl_rate[i in Nodes_set,j in Nodes_set; Set([i,j]) in lines_set],
                    -Sbase*a[i,j,k]*lines_rates[Set([i,j])] ≤ pij[i,j] ≤ Sbase*a[i,j,k]*lines_rates[Set([i,j])] )

            JuMP.@objective(m,Min,sum(GenData[g,:C1][1,1]*p[g]+GenData[g,:C0][1,1] for g in Gen_set))
            JuMP.optimize!(m)

            solution_status = raw_status(m)
            if solution_status == "Model was proven to be infeasible."
                flag = flag + 1
                println(k)
                println(ag[:,1,k])
            end

        end
        if flag == 0
            println("C-OPF/OBS is verified!!")
        else
            println("Failed!!!")
        end
end

function Comparative_Studies(System ::System_Struct,Splittables)
    System_OBS_C = deepcopy(System)
    System_OBS_CP = deepcopy(System)
    System_OBS_CPP = deepcopy(System)
    System_C_SCOPF = deepcopy(System)

    for i in Splittables
        convert_to_aux!(System_OBS_C,i)
        convert_to_aux!(System_OBS_CP,i)
        convert_to_aux!(System_OBS_CPP,i)
    end

    Delta_P = 0.2*System.Gen_Data[:Pmax]

    obj_DCOPF_C,System_C_SCOPF, = Solve_SCOPF!(System_C_SCOPF,"C-SCOPF",Delta_P)
    obj_OBS_C,System_OBS_C,z0,pg_OBS_C = Solve_OBS_SC(System_OBS_C,"C-SCOPF",Delta_P)
    obj_OBS_CPP,System_OBS_CPP,z1,pg_OBS_CPP = Solve_OBS_SC(System_OBS_CPP,"OBS-C-1",Delta_P)
    obj_OBS_CP,System_OBS_CP,z2,pg_OBS_CP = Solve_OBS_SC(System_OBS_CP,"OBS-C-2",Delta_P)

    return obj_DCOPF_C,obj_OBS_C,obj_OBS_CPP,obj_OBS_CP
end

function find_vulnerable_loads(System ::System_Struct)


    list_of_loads = System.BusData[findall(x->x!=0,System.BusData[:Pd]),:bus_i]

    list_of_unsecure_loads = []

    for i in list_of_loads
        list_of_lines_at_load = vcat(System.LineData[findall(x->x==i,System.LineData[:fbus]),:fbus],System.LineData[findall(x->x==i,System.LineData[:tbus]),:tbus])
        if length(list_of_lines_at_load) <= 1
            push!(list_of_unsecure_loads,i)
        end
    end

    return list_of_unsecure_loads
end

function find_vulnerable_generators(System ::System_Struct)

    list_of_generators = System.Gen_Data[:bus]
    list_of_unsecure_gens = []
    for i in list_of_generators
        list_of_lines_at_gen = vcat(System.LineData[findall(x->x==i,System.LineData[:fbus]),:fbus],System.LineData[findall(x->x==i,system.LineData[:tbus]),:tbus])
        if length(list_of_lines_at_gen) <= 1
            push!(list_of_unsecure_gens,i)
        end
    end
    return list_of_unsecure_gens
end

function find_redundant_lines(System ::System_Struct)

    redundant_lines = []
    for line in System.LineData[:ID]
        buses = Set([System.LineData[line,:fbus],System.LineData[line,:tbus]])
        for line_inner in System.LineData[:ID]
            if line != line_inner
                if Set([System.LineData[line_inner,:fbus],System.LineData[line_inner,:tbus]]) == buses
                    push!(redundant_lines,buses)
                end
            end
        end
    end
    return unique!(redundant_lines)

end
