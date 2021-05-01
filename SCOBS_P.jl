import LinearAlgebra,DataFrames,CSV;
import Ipopt;
using Ipopt;
using LinearAlgebra,DataFrames,CSV;
import JuMP
using JuMP
using Printf
import Gurobi
using Gurobi

include("SystemAssembly.jl")
include("Aux_Method.jl")
include("report.jl")
Sbase = 100;

#24-bus
show_cases()
system = load_case(19,Sbase)


set_line_capacity_multiplier!(system,1)
set_generation_multiplier!(system,1.5)
set_load_multiplier!(system,1)


print_latex_data(system)
System_DCOPF = deepcopy(system)

Solve_DCOPF!(System_DCOPF)

find_vulnerable_loads(system)
find_vulnerable_generators(system)
find_redundant_lines(System_DCOPF)
println(system.BusData)

System_DCOPF.Operating_Cost

System_OBS = deepcopy(system)

for i in 1:length(System_OBS.SubStations)
    convert_to_aux!(System_OBS,i)
end




obj_OBS,System_OBS,z,pg = Solve_OBS_aux(System_OBS,300)
obj_OBS
splitted_buses,fb = find_splitted_buses(System_OBS)
splitted_buses = find_splitted_buses_2(System_OBS)
println(System_OBS.BusData_output)
println(System_DCOPF.LineData)

#Preventive Scheme
System_P_SCOPF = deepcopy(system)
obj_DCOPF_P,System_P_SCOPF, = Solve_SCOPF!(System_P_SCOPF,"P-SCOPF")


println(System_P_SCOPF.BusData_output)
println(System_P_SCOPF.LineLoading)

obj_DCOPF_P


System_OBS_P_2 = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_P_2,i)
end
obj_OBS_P_2, System_OBS_P_2,_,_ = Solve_OBS_SC(System_OBS_P_2,"OBS-P-2",0,300)
obj_OBS_P_2

println(System_OBS_P_2.BusData_output)
println(System_OBS_P_2.LineLoading)
println(System_OBS_P_2.Gen_Data)


System_OBS_P_3 = deepcopy(system)
for i in splitted_buses
    convert_to_aux_coupled!(System_OBS_P_3,i)
end
obj_OBS_P_3, System_OBS_P_3,_,_ = Solve_OBS_SC_coupled(System_OBS_P_3,"OBS-P-2",0,300)
obj_OBS_P_3


println(System_OBS_P_3.BusData_output)
println(System_OBS_P_3.LineLoading)
println(System_OBS_P_3.Gen_Data)


#=========================================#
#39-Bus
show_cases()
system = load_case(40,Sbase)


set_line_capacity_multiplier!(system,1.2)
set_generation_multiplier!(system,1)
set_load_multiplier!(system,1)
Delta_P = 0.8*system.Gen_Data[:Pmax]


print_latex_data(system)
System_DCOPF = deepcopy(system)

Solve_DCOPF!(System_DCOPF)

find_vulnerable_loads(system)
find_vulnerable_generators(system)
find_redundant_lines(System_DCOPF)
println(system.BusData)

System_DCOPF.Operating_Cost

System_OBS = deepcopy(system)

for i in 1:length(System_OBS.SubStations)
    convert_to_aux!(System_OBS,i)
end




obj_OBS,System_OBS,z,pg = Solve_OBS_aux(System_OBS,300)
obj_OBS
splitted_buses,fb = find_splitted_buses(System_OBS)
splitted_buses = find_splitted_buses_2(System_OBS)
println(System_OBS.BusData_output)
println(System_DCOPF.LineData)

#Preventive Scheme
System_P_SCOPF = deepcopy(system)
obj_DCOPF_P,System_P_SCOPF, = Solve_SCOPF!(System_P_SCOPF,"P-SCOPF")


println(System_P_SCOPF.BusData_output)
println(System_P_SCOPF.LineLoading)

obj_DCOPF_P


System_OBS_P_2 = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_P_2,i)
end
obj_OBS_P_2, System_OBS_P_2,_,_ = Solve_OBS_SC(System_OBS_P_2,"OBS-P-2",0,300)
obj_OBS_P_2

println(System_OBS_P_2.BusData_output)
println(System_OBS_P_2.LineLoading)
println(System_OBS_P_2.Gen_Data)


System_OBS_P_3 = deepcopy(system)
for i in splitted_buses
    convert_to_aux_coupled!(System_OBS_P_3,i)
end
obj_OBS_P_3, System_OBS_P_3,_,_ = Solve_OBS_SC_coupled(System_OBS_P_3,"OBS-P-2",0,600)
obj_OBS_P_3


println(System_OBS_P_3.BusData_output)
println(System_OBS_P_3.LineLoading)
println(System_OBS_P_3.Gen_Data)

println(System_OBS_P_3.BusData)
println(System_OBS_P_3.LineData)


#=========================================#
#5-Bus
show_cases()
system = load_case(51,Sbase)


set_line_capacity_multiplier!(system,0.9)
set_generation_multiplier!(system,5)
set_load_multiplier!(system,1)


print_latex_data(system)
System_DCOPF = deepcopy(system)

Solve_DCOPF!(System_DCOPF)

find_vulnerable_loads(system)
find_vulnerable_generators(system)
find_redundant_lines(System_DCOPF)
println(system.BusData)

System_DCOPF.Operating_Cost

System_OBS = deepcopy(system)

for i in 1:length(System_OBS.SubStations)
    convert_to_aux!(System_OBS,i)
end




obj_OBS,System_OBS,z,pg = Solve_OBS_aux(System_OBS,300)
obj_OBS
splitted_buses,fb = find_splitted_buses(System_OBS)
splitted_buses = find_splitted_buses_2(System_OBS)
println(System_OBS.BusData_output)
println(System_DCOPF.LineData)

#Preventive Scheme
System_P_SCOPF = deepcopy(system)
obj_DCOPF_P,System_P_SCOPF, = Solve_SCOPF!(System_P_SCOPF,"P-SCOPF")


println(System_P_SCOPF.BusData_output)
println(System_P_SCOPF.LineLoading)

obj_DCOPF_P


System_OBS_P_2 = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_P_2,i)
end
obj_OBS_P_2, System_OBS_P_2,_,_ = Solve_OBS_SC(System_OBS_P_2,"OBS-P-2",0,300)
obj_OBS_P_2

println(System_OBS_P_2.BusData_output)
println(System_OBS_P_2.LineLoading)
println(System_OBS_P_2.Gen_Data)


System_OBS_P_3 = deepcopy(system)
for i in splitted_buses
    convert_to_aux_coupled!(System_OBS_P_3,i)
end
obj_OBS_P_3, System_OBS_P_3,_,_,z_coupler= Solve_OBS_SC_coupled(System_OBS_P_3,"OBS-P-2",0,600)
obj_OBS_P_3

length(z_coupler)
println(z_coupler)
for i in 1:length(z_coupler)/2
    if z_coupler[(17,i)] == 0 && z_coupler[(28,i)] == 0
        println("Possible Error, i = "*string(i))
    end
    println([i,z_coupler[(17,i)],z_coupler[(28,i)]])
end

a,k = Generate_Contingencies(System_OBS_P_3,"P-SCOPF")

i = 29
line_ind = Set([findall(x->x==0,a[:,:,i])[1][1],findall(x->x==0,a[:,:,i])[1][2]])


println(System_OBS_P_3.LineData)
println(System_OBS_P_3.BusData_output)
println(System_OBS_P_3.LineLoading)
println(System_OBS_P_3.Gen_Data)|
System_OBS_P_3.
println(System_OBS_P_3.BusData)
println(System_OBS_P_3.LineData)
