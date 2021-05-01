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


show_cases()
system = load_case(19,Sbase)


set_line_capacity_multiplier!(system,1)
set_generation_multiplier!(system,1.5)
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

System_OBS_P = deepcopy(system)

for i in splitted_buses
    convert_to_aux!(System_OBS_P,i)
end
obj_OBS_P, System_OBS_P,_,_ = Solve_OBS_SC(System_OBS_P,"P-SCOPF",0,600,false)

obj_OBS_P

find_splitted_buses_2(System_OBS_P)

Verify_SCP_PF(System_P_SCOPF)

find_bus_substation(65,System_OBS_P).SubStationID

sum(System_P_SCOPF.Gen_Data[:Pg])
sum(System_OBS_P.Gen_Data[:Pg])
sum(System_OBS.Gen_Data[:Pg])

println(System_OBS_P.BusData_output)
print(System_OBS_P.LineLoading)

System_OBS_P_1 = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_P_1,i)
end


obj_OBS_P_1, System_OBS_P_1,_,_ = Solve_OBS_SC(System_OBS_P_1,"OBS-P-1",0,600)

obj_OBS_P_1

find_splitted_buses_2(System_OBS_P_1)


findall(x->x==72,System_OBS_P_1.BusData_output[:Bus])[1,1]
find_bus_substation(72,System_OBS_P_1).SubStationID


println(System_OBS_P_1.BusData_output)
println(System_OBS_P_1.LineLoading)
println(System_OBS_P_1.Gen_Data)

findall(x->x==13,System_OBS_P_1.Gen_Data[:bus])

System_OBS_P_2 = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_P_2,i)
end
obj_OBS_P_2, System_OBS_P_2,_,_ = Solve_OBS_SC(System_OBS_P_2,"OBS-P-2",0,300)
obj_OBS_P_2

println(System_OBS_P_2.BusData_output)
println(System_OBS_P_2.LineLoading)
println(System_OBS_P_2.Gen_Data)


System_OBS_P_2 = deepcopy(system)
for i in splitted_buses
    convert_to_aux_coupled!(System_OBS_P_2,i)
end
obj_OBS_P_2, System_OBS_P_2,_,_ = Solve_OBS_SC_coupled(System_OBS_P_2,"OBS-P-2",0,300)
obj_OBS_P_2


println(System_OBS_P_2.BusData_output)
println(System_OBS_P_2.LineLoading)
println(System_OBS_P_2.Gen_Data)

#Corrective Scheme
obj_DCOPF_C,System_DCOPF_C, = Solve_SCOPF!(System_DCOPF,"C-SCOPF",Delta_P)


System_OBS_C = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_C,i)
end
obj_OBS_C,System_OBS_C,z,pg = Solve_OBS_SC(System_OBS_C,"C-SCOPF",Delta_P,300)

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


#=========================================#
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

System_OBS_P = deepcopy(system)

for i in splitted_buses
    convert_to_aux!(System_OBS_P,i)
end
obj_OBS_P, System_OBS_P,_,_ = Solve_OBS_SC(System_OBS_P,"P-SCOPF",0,600,false)

obj_OBS_P

find_splitted_buses_2(System_OBS_P)

Verify_SCP_PF(System_P_SCOPF)

find_bus_substation(65,System_OBS_P).SubStationID

sum(System_P_SCOPF.Gen_Data[:Pg])
sum(System_OBS_P.Gen_Data[:Pg])
sum(System_OBS.Gen_Data[:Pg])

println(System_OBS_P.BusData_output)
print(System_OBS_P.LineLoading)

System_OBS_P_1 = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_P_1,i)
end


obj_OBS_P_1, System_OBS_P_1,_,_ = Solve_OBS_SC(System_OBS_P_1,"OBS-P-1",0,600)

obj_OBS_P_1

find_splitted_buses_2(System_OBS_P_1)


findall(x->x==72,System_OBS_P_1.BusData_output[:Bus])[1,1]
find_bus_substation(72,System_OBS_P_1).SubStationID


println(System_OBS_P_1.BusData_output)
println(System_OBS_P_1.LineLoading)
println(System_OBS_P_1.Gen_Data)

findall(x->x==13,System_OBS_P_1.Gen_Data[:bus])

System_OBS_P_2 = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_P_2,i)
end
obj_OBS_P_2, System_OBS_P_2,_,_ = Solve_OBS_SC(System_OBS_P_2,"OBS-P-2",0,300)
obj_OBS_P_2

println(System_OBS_P_2.BusData_output)
println(System_OBS_P_2.LineLoading)
println(System_OBS_P_2.Gen_Data)


System_OBS_P_2 = deepcopy(system)
for i in splitted_buses
    convert_to_aux_coupled!(System_OBS_P_2,i)
end
obj_OBS_P_2, System_OBS_P_2,_,_ = Solve_OBS_SC_coupled(System_OBS_P_2,"OBS-P-2",0,600)
obj_OBS_P_2


println(System_OBS_P_2.BusData_output)
println(System_OBS_P_2.LineLoading)
println(System_OBS_P_2.Gen_Data)

println(System_OBS_P_2.BusData)
println(System_OBS_P_2.LineData)

#Corrective Scheme
obj_DCOPF_C,System_DCOPF_C, = Solve_SCOPF!(System_DCOPF,"C-SCOPF",Delta_P)


System_OBS_C = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_C,i)
end
obj_OBS_C,System_OBS_C,z,pg = Solve_OBS_SC(System_OBS_C,"C-SCOPF",Delta_P,300)

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

#=========================================#
show_cases()
system = load_case(51,Sbase)


set_line_capacity_multiplier!(system,0.9)
set_generation_multiplier!(system,5)
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

System_OBS_P = deepcopy(system)

for i in splitted_buses
    convert_to_aux!(System_OBS_P,i)
end
obj_OBS_P, System_OBS_P,_,_ = Solve_OBS_SC(System_OBS_P,"P-SCOPF",0,600,false)

obj_OBS_P

find_splitted_buses_2(System_OBS_P)

Verify_SCP_PF(System_P_SCOPF)

find_bus_substation(65,System_OBS_P).SubStationID

sum(System_P_SCOPF.Gen_Data[:Pg])
sum(System_OBS_P.Gen_Data[:Pg])
sum(System_OBS.Gen_Data[:Pg])

println(System_OBS_P.BusData_output)
print(System_OBS_P.LineLoading)

System_OBS_P_1 = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_P_1,i)
end


obj_OBS_P_1, System_OBS_P_1,_,_ = Solve_OBS_SC(System_OBS_P_1,"OBS-P-1",0,600)

obj_OBS_P_1

find_splitted_buses_2(System_OBS_P_1)


findall(x->x==72,System_OBS_P_1.BusData_output[:Bus])[1,1]
find_bus_substation(72,System_OBS_P_1).SubStationID


println(System_OBS_P_1.BusData_output)
println(System_OBS_P_1.LineLoading)
println(System_OBS_P_1.Gen_Data)

findall(x->x==13,System_OBS_P_1.Gen_Data[:bus])

System_OBS_P_2 = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_P_2,i)
end
obj_OBS_P_2, System_OBS_P_2,_,_ = Solve_OBS_SC(System_OBS_P_2,"OBS-P-2",0,300)
obj_OBS_P_2

println(System_OBS_P_2.BusData_output)
println(System_OBS_P_2.LineLoading)
println(System_OBS_P_2.Gen_Data)


System_OBS_P_2 = deepcopy(system)
for i in splitted_buses
    convert_to_aux_coupled!(System_OBS_P_2,i)
end
obj_OBS_P_2, System_OBS_P_2,_,_ = Solve_OBS_SC_coupled(System_OBS_P_2,"OBS-P-2",0,600)
obj_OBS_P_2



println(System_OBS_P_2.BusData_output)
println(System_OBS_P_2.LineLoading)
println(System_OBS_P_2.Gen_Data)

println(System_OBS_P_2.BusData)
println(System_OBS_P_2.LineData)

#Corrective Scheme
obj_DCOPF_C,System_DCOPF_C, = Solve_SCOPF!(System_DCOPF,"C-SCOPF",Delta_P)


System_OBS_C = deepcopy(system)
for i in splitted_buses
    convert_to_aux!(System_OBS_C,i)
end
obj_OBS_C,System_OBS_C,z,pg = Solve_OBS_SC(System_OBS_C,"C-SCOPF",Delta_P,300)

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
