using LinearAlgebra,DataFrames,DataFramesMeta,CSV;
using Printf
import JuMP, Ipopt, Cbc
using JuMP, Ipopt, Cbc


include("SystemAssembly.jl")
include("report.jl")
Sbase = 100;

#Instantiate System

system_dir = string(pwd(),"\\Systems\\5_Bus_splitting");
branch_dir = string(system_dir,"\\branch5.csv");
bus_dir = string(system_dir,"\\bus5.csv");
gen_dir = string(system_dir,"\\gen5.csv");

branch_dir_5S = string(system_dir,"\\branch6.csv");
branch_dir_5S_2 = string(system_dir,"\\branch6_2.csv");
branch_dir_5S_3 = string(system_dir,"\\branch6_3.csv");
bus_dir_5S = string(system_dir,"\\bus6.csv");
gen_dir_5S = string(system_dir,"\\gen6.csv");

System_5Bus_ACOPF = SystemAssembly(branch_dir,Sbase,bus_dir,gen_dir)
System_6Bus_ACOPF = SystemAssembly(branch_dir_5S,Sbase,bus_dir_5S,gen_dir_5S)
System_6Bus_ACOPF_2 = SystemAssembly(branch_dir_5S_2,Sbase,bus_dir_5S,gen_dir_5S)
System_6Bus_ACOPF_3 = SystemAssembly(branch_dir_5S_3,Sbase,bus_dir_5S,gen_dir_5S)

get_network_matrices!(System_5Bus_ACOPF)
System_5Bus_ACOPF.Y_bus_inc
System_5Bus_ACOPF.Y_bus_inc == System_5Bus_ACOPF.y_bus
System_5Bus_ACOPF.I
System_5Bus_ACOPF.A


create_extended_incidence_matrix!(System_5Bus_ACOPF)
System_5Bus_ACOPF.A_ex
System_5Bus_ACOPF.z
System_5Bus_ACOPF.ρ
System_5Bus_ACOPF.Δ
System_5Bus_ACOPF.D

get_connection_status(System_5Bus_ACOPF)


system_perturber(System_5Bus_ACOPF,6)

Solve_OPF!(System_5Bus_ACOPF,"ACOPF")
Solve_OPF!(System_6Bus_ACOPF,"DCOPF")
Solve_OPF!(System_6Bus_ACOPF_2,"ACOPF")
Solve_OPF!(System_6Bus_ACOPF_3,"ACOPF")

System_5Bus_ACOPF.Operating_Cost
System_6Bus_ACOPF.Operating_Cost
System_6Bus_ACOPF_2.Operating_Cost
System_6Bus_ACOPF_3.Operating_Cost
System_5Bus_ACOPF.
Solve_OTS!(System_5Bus_ACOPF,"DCOTS",false)
Solve_OTS!(System_6Bus_ACOPF,"DCOTS",false)
Solve_OTS!(System_6Bus_ACOPF_2,"DCOTS",false)
Solve_OTS!(System_6Bus_ACOPF_3,"DCOTS",false)

ybus5 = System_5Bus_ACOPF.y_bus
ybus6 = System_6Bus_ACOPF.y_bus
ybus6_2 = System_6Bus_ACOPF_2.y_bus
ybus6_3 = System_6Bus_ACOPF_3.y_bus

A = Matrix(1I, 6, 5)
A[5,5] = 0
A[6,5] = 1
A * ybus5 * transpose(A)

det_5 = det(ybus5)
det_6 = det(ybus6)
det_6_2 = det(ybus6_2)
det_6_3 = det(ybus6_3)

abs(det_5)
abs(det_6)
abs(det_6_2)
abs(det_6_3)

eig_values_5,eig_vec_5 = eigen(ybus5)
eig_abs_5 = abs.(eig_values_5)

eig_values_6,eig_vec_6 = eigen(ybus6)
eig_abs_6 = abs.(eig_values_6)

eig_values_6_2,eig_vec_6_2 = eigen(ybus6_2)
eig_abs_6_2 = abs.(eig_values_6_2)

eig_values_6_3,eig_vec_6_3 = eigen(ybus6_3)
eig_abs_6_3 = abs.(eig_values_6_3)

sum(eig_abs_5)
sum(eig_abs_6)
sum(eig_abs_6_2)
sum(eig_abs_6_3)

system_dir_3 = string(pwd(),"\\Systems\\3_Bus");
branch_dir_3 = string(system_dir_3,"\\branch3.csv");
bus_dir_3 = string(system_dir_3,"\\bus3.csv");
gen_dir_3 = string(system_dir_3,"\\gen3.csv");

System_3Bus_ACOPF = SystemAssembly(branch_dir_3,Sbase,bus_dir_3,gen_dir_3)

Solve_OPF!(System_3Bus_ACOPF,"ACOPF")
remove_line!(System_3Bus_ACOPF,[0,1,1])

Solve_OTS!(System_3Bus_ACOPF,"ACOTS",false)

System_5Bus_ACOPF.LineData


System_5Bus_ACOPF.LineData[:,1]
findfirst(x->x==5,System_5Bus_ACOPF.LineData[:,1])
size(System_5Bus_ACOPF.LineData)
System_5Bus_ACOPF.N_bus
LineData = System_5Bus_ACOPF.LineData
delete!(LineData[:],1)

a = [1, 2,3]

b = [4,5,6]

c = vec(zeros(5,1))
count(x->x == 3,a)
df = DataFrame(A = a,B = b)
b=[]
