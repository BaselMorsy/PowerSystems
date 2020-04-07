include("Ybus.jl")
import JuMP, Ipopt,Juniper
using Parameters
using Ipopt, Cbc,Juniper, LinearAlgebra


module System_Definition
    using Parameters
    using DataFrames
    export System_Struct

    @with_kw mutable struct System_Struct
        y_bus
        b_line
        LineData
        BusData
        N_bus
        Sbase

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

    end
end  # module

using Main.System_Definition;

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
        Lines = DataFrame(CSV.read(LineData))
        N_bus = length(unique!(vcat(unique!(Vector(Lines[:fbus])),unique!(Vector(Lines[:tbus])))))
        YBUS,bLine = Y_Bus(Lines,N_bus)
        S_rating_lines = DataFrame(fbus = Lines[!,:fbus],tbus = Lines[!,:tbus]
        ,SL_Rating = Lines[!,:rate])
    end


    if BusData_sheet != nothing
        busData = DataFrame(CSV.read(BusData_sheet))
    end

    if GenData != nothing
        gendata = DataFrame(CSV.read(GenData))
        gen_data = DataFrame(bus = gendata[:bus],C2 = gendata[:c2],
            C1 = gendata[:c1],C0 = gendata[:c0],
            Pmax = gendata[:Pmax],
            Pmin = gendata[:Pmin],
            Qmax = gendata[:Qmax],
            Qmin = gendata[:Qmin]);
    end

    Sys = System_Struct(y_bus = YBUS,b_line = bLine, LineData = Lines,BusData =busData
        ,N_bus=N_bus,Sbase = Sbase)
    append!(Sys.Line_Constraints,S_rating_lines);
    Sys.Gen_Data=gen_data;
    Sys.Pd_total = sum(busData[:Pd]);
    return Sys
end


function Solve_OPF!(System ::System_Struct, method = "ACOPF", print_status = true)

    #OPF solver

    G = real(System.y_bus);
    B = imag(System.y_bus);
    b = System.b_line;

    N = System.N_bus
    N_Lines = size(System.Line_Constraints,1)
    Nodes_set = 1:N;
    Gen_set = System.Gen_Data[!,:bus];

    Sbase = System.Sbase;

    BusData = System.BusData;
    GenData = System.Gen_Data;
    BranchData = System.Line_Constraints;
    S = zeros(N,N)

        for l in 1:N_Lines
            S[Int64(BranchData[l,:fbus]),Int64(BranchData[l,:tbus])] = BranchData[l,:SL_Rating]
            S[Int64(BranchData[l,:tbus]),Int64(BranchData[l,:fbus])] = BranchData[l,:SL_Rating]
        end

    if method == "ACOPF"


            m = JuMP.Model(JuMP.with_optimizer(Ipopt.Optimizer, print_level =0))

            # 2.1.Variables
            JuMP.@variable(m, BusData[i, :Vmin] ≤ v[i in Nodes_set] ≤ BusData[i, :Vmax])
            JuMP.@variable(m, -2*π ≤ δ[i in Nodes_set] ≤ 2*π)

            JuMP.@variable(m, GenData[(GenData[:bus].==g), :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[(GenData[:bus].==g), :Pmax][1,1])
            JuMP.@variable(m, GenData[(GenData[:bus].==g), :Qmin][1,1] ≤ q[g in Gen_set] ≤ GenData[(GenData[:bus].==g), :Qmax][1,1])

            JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set])
            JuMP.@variable(m, qij[i = Nodes_set, j = Nodes_set])

            # 2.2. Constraints
            #JuMP.@constraint(m, ReferenceAngle,
            #    (δ[1] ==  0.0))

            # ACTIVE POWER THROUGH LINE N-M
            JuMP.@NLconstraint(m, p_line[i in Nodes_set, j in Nodes_set],
                 (pij[i,j] ==  Sbase*(v[i]*v[j]*(G[i,j]*cos(δ[i]-δ[j])+B[i,j]*sin(δ[i]-δ[j])) -(v[i]^2)*G[i,j] )))

            # REACTIVE POWER THROUGH LINE N-M
            JuMP.@NLconstraint(m, q_line[i in Nodes_set, j in Nodes_set],
                 (qij[i,j] ==  Sbase*(v[i]*v[j]*(G[i,j]*sin(δ[i]-δ[j]) - B[i,j]*cos(δ[i]-δ[j])) +(v[i]^2)*(B[i,j]-b[i,j]) )));

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

            elseif method =="DCOPF"

                B_ = -B[1:N,1:N];

                m = JuMP.Model(JuMP.with_optimizer(Ipopt.Optimizer, print_level =0))

                # 2.1.Variables
                JuMP.@variable(m, -2*π ≤ δ[i in Nodes_set] ≤ 2*π)

                JuMP.@variable(m, GenData[(GenData[:bus].==g), :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[(GenData[:bus].==g), :Pmax][1,1])

                JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set])

                # 2.2. Constraints
                JuMP.@constraint(m, ReferenceAngle,
                    (δ[1] ==  0.0))

                JuMP.@constraint(m,Nodal_balance[i in Nodes_set],
                        sum(pij[i,j] for j = Nodes_set) == sum(p[g] for g in Gen_set if GenData[(GenData[:bus].==g),:bus][1,1]==i) - BusData[i,:Pd])
                JuMP.@constraint(m,pl[i in Nodes_set,j in Nodes_set],
                        pij[i,j] == Sbase*B_[i,j]*(δ[i]-δ[j]))
                JuMP.@constraint(m,pl_rate[i in Nodes_set,j in Nodes_set],
                        -Sbase*S[i,j] ≤ pij[i,j] ≤ Sbase*S[i,j])

                JuMP.@objective(m,Min,sum(GenData[(GenData[:bus].==g),:C1][1,1]*p[g]+GenData[(GenData[:bus].==g),:C0][1,1] for g in Gen_set))
                JuMP.optimize!(m)

                Pg = JuMP.value.(p)
                Pij = JuMP.value.(pij)
                δ = JuMP.value.(δ)
                V = ones(N)
                Pd = System.BusData[!,3]
                Qd = System.BusData[!,4]
                System.Operating_Cost = JuMP.objective_value(m)

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
end

function Solve_OTS!(System ::System_Struct,method = "ACOTS",use_dc_flag = true , print_status = true)

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
                    (δ[1] ==  0.0))


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
