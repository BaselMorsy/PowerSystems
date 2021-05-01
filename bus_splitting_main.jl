using Parsers
import LinearAlgebra,DataFrames,CSV;
import Ipopt;
using Ipopt;
using LinearAlgebra,DataFrames,CSV;
import JuMP
using JuMP
using Printf
import Gurobi
using Gurobi
using StatsPlots



include("SystemAssembly.jl")
include("Aux_Method.jl")
include("report.jl")
Sbase = 100;

#Load case
show_cases()
system = load_case(49,Sbase)
set_line_capacity_multiplier!(system,0.9)
set_generation_multiplier!(system,5)
system.LineData[!,:b] = vec(zeros(length(system.LineData[!,:b]),1))

println(system.Gen_Data)
system.BusData
show(stdout, MIME("text/latex"),system.BusData)

#Solve DCOPF
System_DCOPF = deepcopy(system)
Pg_DCOPF,_ = Solve_DCOPF!(System_DCOPF)
Pg_DCOPF
minimum_cost,Pg_ED,MCP = solve_ED(System_DCOPF)
System_DCOPF.Gen_Data
println(System_DCOPF.LineLoading)
System_DCOPF.Gen_Data
show(stdout, MIME("text/latex"),System_DCOPF.Gen_Data[:,[1,9]])

SU_DCOPF = sum(System_DCOPF.LineLoading[:Utilization])/6

obj_DCOPF = System_DCOPF.Operating_Cost
#Solve N-1 DCOPF
System_DCOPF_P = deepcopy(system)
System_DCOPF_C = deepcopy(system)

obj_DCOPF_P,System_DCOPF_P, = Solve_SCOPF!(System_DCOPF_P,"P-SCOPF")

show(stdout, MIME("text/latex"),System_DCOPF_P.LineLoading)
show(stdout, MIME("text/latex"),System_DCOPF_P.Gen_Data[:,[1,9]])
println(System_DCOPF_P.Gen_Data)
println(System_DCOPF_P.BusData)
println(System_DCOPF_P.LineLoading)
obj_DCOPF_P

Delta_P = 0.2*System_DCOPF_C.Gen_Data[:Pmax]
obj_DCOPF_C,System_DCOPF_C, = Solve_SCOPF!(System_DCOPF_C,"C-SCOPF",Delta_P)
obj_DCOPF_C

show(stdout, MIME("text/latex"),System_DCOPF_C.LineLoading)
show(stdout, MIME("text/latex"),System_DCOPF_C.Gen_Data[:,[1,9]])
println(System_DCOPF_C.Gen_Data)
println(System_DCOPF_C.BusData)
println(System_DCOPF_C.LineLoading)

SU_DCOPF_P = sum(System_DCOPF_P.LineLoading[:Utilization])/6
SU_DCOPF_C = sum(System_DCOPF_C.LineLoading[:Utilization])/6
#OBS
System_OBS = deepcopy(system)
for i in 1:length(system.SubStations)
    convert_to_aux!(System_OBS,i)
end

obj_OBS,System_OBS,z,pg = Solve_OBS_aux(System_OBS)
obj_OBS
Pg_ED-pg

println(System_OBS.Gen_Data)
println(System_OBS.BusData)
println(System_OBS.LineLoading)

SU_OBS = sum(System_OBS.LineLoading[:Utilization])/6

show(stdout, MIME("text/latex"),System_OBS.LineLoading)
show(stdout, MIME("text/latex"),System_OBS.Gen_Data[:,[1,9]])
#OBS-P
System_OBS_P = deepcopy(system)
for i in 1:length(System_OBS_P.SubStations)
    convert_to_aux!(System_OBS_P,i)
end
obj_OBS_P,System_OBS_P,z,pg = Solve_OBS_SC(System_OBS_P,"P-SCOPF")
obj_OBS_P
show(stdout, MIME("text/latex"),System_OBS_P.LineLoading)
show(stdout, MIME("text/latex"),System_OBS_P.Gen_Data[:,[1,9]])
println(System_OBS_P.Gen_Data)
println(System_OBS_P.LineLoading)
println(System_OBS_P.BusData)

SU_OBS_P = sum(System_OBS_P.LineLoading[:Utilization])/6

#OBS-C
System_OBS_C = deepcopy(system)
for i in 1:length(System_OBS_C.SubStations)
    convert_to_aux!(System_OBS_C,i)
end


Delta_P = 0.2*System_OBS_C.Gen_Data[:Pmax]
obj_OBS_C,System_OBS_C,z,pg = Solve_OBS_SC(System_OBS_C,"C-SCOPF",Delta_P)
obj_OBS_C
show(stdout, MIME("text/latex"),System_OBS_C.LineLoading)
show(stdout, MIME("text/latex"),System_OBS_C.Gen_Data[:,[1,9]])
println(System_OBS_C.Gen_Data)
println(System_OBS_C.BusData)
println(System_OBS_C.LineLoading)
SU_OBS_C = sum(System_OBS_C.LineLoading[:Utilization])/6


System_DCOPF_P.Gen_Data[:Pg]
System_DCOPF_P.Gen_Data[[:Pmax,:Pmin]]
System_DCOPF_C.Gen_Data[1,:Pg]+Delta_P[1]
Delta_P[1]
System_DCOPF_C.Gen_Data[1,:Pg]
System_DCOPF_C.Gen_Data[1,:Pmax]

Pmax = [minimum(vcat(System_DCOPF_C.Gen_Data[g,:Pmax],System_DCOPF_C.Gen_Data[g,:Pg]+Delta_P[g])) for g in 1:5]
Pmin = [maximum(vcat(System_DCOPF_C.Gen_Data[g,:Pmin],System_DCOPF_C.Gen_Data[g,:Pg]-Delta_P[g])) for g in 1:5]
Verify_SCP_PF(System_DCOPF_P)
Verify_SCC_OPF(System_DCOPF_C,Delta_P)
System_OBS_C.BusData
System_OBS_C.Gen_Data
sort!(System_OBS_C.BusData,[:bus_i])
include("Aux_Method.jl")

Verify_SCC_OPF(System_OBS_C,Delta_P)

names = repeat(["No Security","P-Scheme","C-Scheme"], outer = 2)
groups = repeat(["DCOPF", "OBS"],outer = 3)
groupedbar(names,[obj_DCOPF, obj_DCOPF_P, obj_DCOPF_C,obj_OBS ,obj_OBS_P, obj_OBS_C],
    group = groups,ylim = (0,29000),xlabel = "Security Scheme",ylabel = "Cost (\$/h)")
savefig("costs.png")

#================OBS-C-1==================#
#pre and post contingency splitting
System_OBS_C_1 = deepcopy(system)
for i in 1:length(System_OBS_C_1.SubStations)
    convert_to_aux!(System_OBS_C_1,i)
end


Delta_P = 0.2*System_OBS_C_1.Gen_Data[:Pmax]
obj_OBS_C_1,System_OBS_C_1,z,pg = Solve_OBS_SC(System_OBS_C_1,"OBS-C-1",Delta_P)
pg
println(System_OBS_C_1.BusData)
System_OBS_C_1.BusData_output
println(System_OBS_C_1.Gen_Data)
println(System_OBS_C_1.LineLoading)

#only post-contingency splitting
System_OBS_C_2 = deepcopy(system)
for i in 1:length(System_OBS_C_2.SubStations)
    convert_to_aux!(System_OBS_C_2,i)
end

Delta_P = 0.2*System_OBS_C_2.Gen_Data[:Pmax]
obj_OBS_C_2,System_OBS_C_2,z,pg_OBS_C_2 = Solve_OBS_SC(System_OBS_C_2,"OBS-C-2",Delta_P)
System_OBS_C_2.BusData_output
pg
System_OBS_C_2.LineLoading
println(System_OBS_C_2.Gen_Data)
println(System_OBS_C_2.LineLoading)
println(System_OBS_C_2.BusData)

obj_scopf_c,obj_obs_c,obj_obs_cpp,obj_obs_cp = Comparative_Studies(system)

#===================24 Bus System===========================#
show_cases()
system = load_case(20,Sbase)
set_line_capacity_multiplier!(system,1)
set_generation_multiplier!(system,1)
set_load_multiplier!(system,1)
Delta_P = 0.5*system.Gen_Data[:Pmax]

System_DCOPF = deepcopy(system)
Solve_DCOPF!(System_DCOPF)

list_of_loads = system.BusData[findall(x->x!=0,system.BusData[:Pd]),:bus_i]

list_of_unsecure_loads = []

for i in list_of_loads
    list_of_lines_at_load = vcat(system.LineData[findall(x->x==i,system.LineData[:fbus]),:fbus],system.LineData[findall(x->x==i,system.LineData[:tbus]),:tbus])
    if length(list_of_lines_at_load) <= 1
        push!(list_of_unsecure_loads,i)
    end
end
list_of_unsecure_loads

list_of_generators = system.Gen_Data[:bus]
list_of_unsecure_gens = []
for i in list_of_generators
    list_of_lines_at_gen = vcat(system.LineData[findall(x->x==i,system.LineData[:fbus]),:fbus],system.LineData[findall(x->x==i,system.LineData[:tbus]),:tbus])
    if length(list_of_lines_at_gen) <= 1
        push!(list_of_unsecure_gens,i)
    end
end
list_of_unsecure_gens

println(system.Gen_Data)


System_OBS = deepcopy(system)

for i in 1:length(System_OBS.SubStations)
    convert_to_aux!(System_OBS,i)
end


obj_OBS,System_OBS,z,pg = Solve_OBS_aux(System_OBS,300)
splitted_buses,fb = find_splitted_buses(System_OBS)

sum(System_OBS.Gen_Data[:Pg])
sum(System_DCOPF.Gen_Data[:Pg])
obj_OBS
println(System_OBS.LineLoading[:Utilization]-System_DCOPF.LineLoading[:Utilization])

println(System_OBS.LineLoading)
println(System_DCOPF.BusData_output)
println(System_OBS.BusData_output)
println(System_OBS.Gen_Data)
println(System_DCOPF.Gen_Data)

#Corrective Scheme
obj_DCOPF_C,System_DCOPF_C, = Solve_SCOPF!(System_DCOPF,"C-SCOPF",Delta_P)


System_OBS_C = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_C,i)
end
obj_OBS_C,System_OBS_C,z,pg = Solve_OBS_SC(System_OBS_C,"C-SCOPF",Delta_P)

obj_OBS_C

System_OBS_C_1 = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_C_1,i)
end

obj_OBS_C_1,System_OBS_C_1,z,pg = Solve_OBS_SC(System_OBS_C_1,"OBS-C-1",Delta_P,800)
obj_OBS_C_1
System_OBS_C_1.BusData

System_OBS_C_2 = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_C_2,i)
end

obj_OBS_C_2,System_OBS_C_2,z,pg = Solve_OBS_SC(System_OBS_C_2,"OBS-C-2",Delta_P,800)
obj_OBS_C_2
obj_scopf_c,obj_obs_c,obj_obs_cpp,obj_obs_cp = Comparative_Studies(system,splittables)

#Preventive Scheme
System_P_SCOPF = deepcopy(system)
obj_DCOPF_P,System_P_SCOPF, = Solve_SCOPF!(System_P_SCOPF,"P-SCOPF")

System_P_SCOPF.LineLoading
println(System_P_SCOPF.BusData_output)
println(system.LineData)

obj_DCOPF_P

System_OBS_P = deepcopy(system)

for i in splitted_buses
    convert_to_aux!(System_OBS_P,i)
end
obj_OBS_P, System_OBS_P,_,_ = Solve_OBS_SC(System_OBS_P,"P-SCOPF",0,300)

obj_OBS_P

find_splitted_buses(System_OBS_P)

Verify_SCP_PF(System_P_SCOPF)

sum(System_P_SCOPF.Gen_Data[:Pg])
sum(System_OBS_P.Gen_Data[:Pg])
sum(System_OBS.Gen_Data[:Pg])

println(System_OBS_P.BusData_output)
print(System_OBS_P.LineData)

System_OBS_P_1 = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_P_1,i)
end


obj_OBS_P_1, System_OBS_P_1,_,_ = Solve_OBS_SC(System_OBS_P_1,"OBS-P-1",0,800)

obj_OBS_P_1
find_splitted_buses(System_OBS_P_1)

println(System_OBS_P_1.BusData_output)
println(System_OBS_P_1.LineLoading)
println(System_OBS_P_1.Gen_Data)


System_OBS_P_2 = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_P_2,i)
end
obj_OBS_P_2, System_OBS_P_2,_,_ = Solve_OBS_SC(System_OBS_P_2,"OBS-P-2",0,800)
obj_OBS_P_2

println(System_OBS_P_2.BusData)
println(System_OBS_P_2.LineLoading)
println(System_OBS_P_2.Gen_Data)
#=========================================================#

show_cases()
system = load_case(20,Sbase)
set_line_capacity_multiplier!(system,1)
set_generation_multiplier!(system,1.5)
set_load_multiplier!(system,1)

#Solve DCOPF
System_DCOPF = deepcopy(system)
Pg_DCOPF,_ = Solve_DCOPF!(System_DCOPF)
Pg_DCOPF
minimum_cost,Pg_ED,MCP = solve_ED(System_DCOPF)


SU_DCOPF = sum(System_DCOPF.LineLoading[:Utilization])/6

obj_DCOPF = System_DCOPF.Operating_Cost
#Solve N-1 DCOPF
System_DCOPF_P = deepcopy(system)
System_DCOPF_C = deepcopy(system)

obj_DCOPF_P,System_DCOPF_P, = Solve_SCOPF!(System_DCOPF_P,"P-SCOPF")

obj_DCOPF_P

Delta_P = 0.2*System_DCOPF_C.Gen_Data[:Pmax]
obj_DCOPF_C,System_DCOPF_C, = Solve_SCOPF!(System_DCOPF_C,"C-SCOPF",Delta_P)
obj_DCOPF_C



SU_DCOPF_P = sum(System_DCOPF_P.LineLoading[:Utilization])/6
SU_DCOPF_C = sum(System_DCOPF_C.LineLoading[:Utilization])/6

#OBS
System_OBS = deepcopy(system)
for i in 1:length(system.SubStations)
    convert_to_aux!(System_OBS,i)
end

obj_OBS,System_OBS,z,pg = Solve_OBS_aux(System_OBS)
obj_OBS
Pg_ED-pg
splitted_buses,fb = find_splitted_buses(System_OBS)
println(System_OBS.Gen_Data)
println(System_OBS.BusData)
println(System_OBS.LineData)
println(System_OBS.LineLoading)

SU_OBS = sum(System_OBS.LineLoading[:Utilization])/6

#OBS-P

System_OBS_P_2 = deepcopy(system)
for i in splitted_buses
    convert_to_aux_coupled!(System_OBS_P_2,i)
end
obj_OBS_P_2, System_OBS_P_2,_,_ = Solve_OBS_SC_coupled(System_OBS_P_2,"OBS-P-2",0,200)
obj_OBS_P_2


println(System_OBS_P_2.BusData_output)
println(System_OBS_P_2.LineLoading)
println(System_OBS_P_2.Gen_Data)

#OBS-P

System_OBS_P_2 = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_P_2,i)
end
obj_OBS_P_2, System_OBS_P_2,_,_ = Solve_OBS_SC(System_OBS_P_2,"OBS-P-2",0,200)
obj_OBS_P_2

println(System_OBS_P_2.BusData)
println(System_OBS_P_2.LineLoading)
println(System_OBS_P_2.Gen_Data)


#===================30 Bus System===========================#
show_cases()
system_30 = load_case(34,Sbase)
set_line_capacity_multiplier!(system_30,25)
set_generation_multiplier!(system_30,25)
system_30.LineData[!,:b] = vec(zeros(length(system_30.LineData[!,:b]),1))
System_DCOPF_30 = deepcopy(system_30)
Pg_DCOPF_30,_ = Solve_DCOPF!(System_DCOPF_30)

obj_scopf_c30,obj_obs_c30,obj_obs_cpp30,obj_obs_cp30 = Comparative_Studies(system_30)

include("Aux_Method.jl")
System_C_SCOPF = deepcopy(system)
System_C_SCOPF.Gen_Data
Delta_P = 0.2*System_C_SCOPF.Gen_Data[:Pmax]
obj_DCOPF_C,System_C_SCOPF, = Solve_SCOPF!(System_C_SCOPF,"C-SCOPF",Delta_P)


df_lines = system_30.LineData[:,[:fbus,:tbus]]
lines_set = Set.(convert.(Array, row for row in eachrow(df_lines)))
Nodes_at_Lines = Dict()
for i in system_30.LineData[:ID]
    Nodes_at_Lines[i] = Set([system_30.LineData[i,:fbus],system_30.LineData[i,:tbus]])
end
for (id,v) in Nodes_at_Lines
    println(id)
    println(v)
end
Nodes_set = system_30.BusData[:,:bus_i];
for i in Nodes_set
    for j in Nodes_set
        if Set([i,j]) in lines_set
            println([i,j])
            println(system_30.LineData[[id for (id,v) in Nodes_at_Lines if v == Set([i,j])][1],:x])
        end
    end
end

##########################-14 Bus System-############################

#Instantiate System
system_dir = string(pwd(),"\\Systems\\14_Bus");
branch_dir = string(system_dir,"\\branch.csv");
bus_dir = string(system_dir,"\\bus.csv");
gen_dir = string(system_dir,"\\gen.csv");

System_14Bus_ACOPF = SystemAssembly(branch_dir,Sbase,bus_dir,gen_dir)
set_load_multiplier!(System_14Bus_ACOPF,1)
set_line_capacity_multiplier!(System_14Bus_ACOPF,1.5)
Solve_OPF!(System_14Bus_ACOPF,"DCOPF")


t_14 = @elapsed Sys_new = OBS_optimizer_unsorted!(deepcopy(System_14Bus_ACOPF),
    "DCOPF",1,Inf,"dsp_s",0.8)

t_14 = @elapsed Sys_new = OBS_optimizer_unsorted!(deepcopy(System_14Bus_ACOPF),
    "DCOPF",1,Inf,"dsp_s_a",0.1)

t_14_2 = @elapsed Sys_new_2 = OBS_optimizer_unsorted!(deepcopy(System_14Bus_ACOPF),
    "DCOPF",1,Inf,"dsp",1)

Y = System_14Bus_ACOPF.Y_bus_inc

abs(det(Y))


det(Sys_new.y_bus)
det(System_14Bus_ACOPF.y_bus)


round.(System_14Bus_ACOPF.Bus_Duals,digits = 3)
round.(System_14Bus_ACOPF.Line_Duals,digits = 3)

frequency = freqtable(round.(System_14Bus_ACOPF.Line_Duals,digits = 3))
frequency2= freqtable(round.(System_14Bus_ACOPF.Bus_Duals,digits = 3))


round.(Sys_new.Bus_Duals,digits = 3)
round.(Sys_new.Line_Duals,digits = 3)

frequency_new = freqtable(round.(Sys_new.Line_Duals,digits = 3))
frequency2_new= freqtable(round.(Sys_new.Bus_Duals,digits = 3))

##########################-5 Bus System-############################

system_dir_5 = string(pwd(),"\\Systems\\5_Bus_splitting");
branch_dir_5 = string(system_dir_5,"\\branch5.csv");
bus_dir_5 = string(system_dir_5,"\\bus5.csv");
gen_dir_5 = string(system_dir_5,"\\gen5.csv");

System_5Bus_ACOPF = SystemAssembly(branch_dir_5,Sbase,bus_dir_5,gen_dir_5)


set_load_multiplier!(System_5Bus_ACOPF,0.8)
Solve_OPF!(System_5Bus_ACOPF,"ACOPF")

Sys_new_5 = OBS_optimizer_unsorted!(deepcopy(System_5Bus_ACOPF),
    "ACOPF",1,Inf,"dsp_s_a",1)
Sys_new_5.Line_Constraints

det(Sys_new_5.y_bus)
det(System_5Bus_ACOPF.y_bus)

round.(System_5Bus_ACOPF.Bus_Duals,digits = 3)
round.(System_5Bus_ACOPF.Line_Duals,digits = 3)

frequency_5 = freqtable(round.(System_5Bus_ACOPF.Line_Duals,digits = 3))
frequency2_5= freqtable(round.(System_5Bus_ACOPF.Bus_Duals,digits = 3))


round.(Sys_new_5.Bus_Duals,digits = 3)
round.(Sys_new_5.Line_Duals,digits = 3)

frequency_new_5 = freqtable(round.(Sys_new_5.Line_Duals,digits = 3))
frequency2_new_5 = freqtable(round.(Sys_new_5.Bus_Duals,digits = 3))

##########################-118 Bus System-############################

system_dir_118 = string(pwd(),"\\Systems\\118_Bus");
branch_dir_118 = string(system_dir_118,"\\branch118.csv");
bus_dir_118 = string(system_dir_118,"\\bus118.csv");
gen_dir_118 = string(system_dir_118,"\\gen118.csv");

System_118Bus_ACOPF = SystemAssembly(branch_dir_118,Sbase,bus_dir_118,gen_dir_118)

set_load_multiplier!(System_118Bus_ACOPF,4.5)
set_generation_multiplier!(System_118Bus_ACOPF,2.5)
set_line_capacity_multiplier!(System_118Bus_ACOPF,0.7)
Solve_OPF!(System_118Bus_ACOPF,"DCOPF")

t = @elapsed Sys_new_118 = OBS_optimizer_unsorted!(deepcopy(System_118Bus_ACOPF),
    "DCOPF",0,Inf,"dsp_s_a",0.9)

t2 = @elapsed Sys_new_118_2 = OBS_optimizer_unsorted!(deepcopy(System_118Bus_ACOPF),
    "DCOPF",1,Inf,"dsp_s",0.02)


System_118Bus_ACOPF.Line_Constraints[1:5,:]
size(Sys_new_118_2.A_ex)
vals,vecs = eigen(System_118Bus_ACOPF.y_bus)
vals_abs = abs.(vals)
maximum(vals_abs)

abs(det(System_118Bus_ACOPF.Y_bus_inc))
abs(det(Sys_new_118_2.Y_bus_inc))
det(transpose(Sys_new_118_2.A_ex)*Sys_new_118_2.Y_pr*Sys_new_118_2.A_ex)
prod(vals)
vals_s,vecs_s = eigen(Sys_new_118_2.Y_bus_inc)
vals_abs_s = abs.(vals_s)
maximum(vals_abs_s)
minimum(vals_abs)
minimum(vals_abs_s)

abs(sum(vals))
abs(sum(vals_s))

size(vals_abs)
size(vals_abs_s)
size(Sys_new_118_2.Y_bus_inc)
size(System_118Bus_ACOPF.y_bus)

round.(System_118Bus_ACOPF.Bus_Duals,digits = 3)
round.(System_118Bus_ACOPF.Line_Duals,digits = 3)

frequency_118 = freqtable(round.(System_118Bus_ACOPF.Line_Duals,digits = 3))
frequency2_118 = freqtable(round.(System_118Bus_ACOPF.Bus_Duals,digits = 3))


round.(Sys_new_118.Bus_Duals,digits = 3)
round.(Sys_new_118.Line_Duals,digits = 3)

frequency_new_118 = freqtable(round.(Sys_new_118.Line_Duals,digits = 3))
frequency2_new_118 = freqtable(round.(Sys_new_118.Bus_Duals,digits = 3))
