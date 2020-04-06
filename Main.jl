using LinearAlgebra,DataFrames,DataFramesMeta,CSV;
import JuMP, Ipopt, Cbc
import Juniper
using Juniper,JuMP, Ipopt, Cbc

include("SystemAssembly.jl")
include("report.jl")
Sbase = 1000;

#Instantiate System
System_ACOPF = SystemAssembly("branch.csv",Sbase,"bus.csv","gen.csv")
System_ACOPF_2 = SystemAssembly("branch.csv",Sbase,"bus.csv","gen.csv")


line_status=Solve_OTS!(System_ACOPF,0,"ACOTS")
Solve_OPF!(System_ACOPF_2,"ACOPF")
System_ACOPF_2.BusData_output


System_14_ACOPF = SystemAssembly("branch14.csv",1000,"bus14.csv","gen14.csv")
System_14_ACOPF_M = SystemAssembly("branch14.csv",1000,"bus14.csv","gen14.csv")

line_status_14 = Solve_OTS!(System_14_ACOPF,1,"ACOTS")
Solve_OPF!(System_14_ACOPF_M,"ACOPF")
System_14_ACOPF_M.LineLoading
System_14_ACOPF_M.LineLoading[:,1]

for i in Nodes_set, j in Nodes_set if ([i,j] in lines_set_1 || [i,j] in lines_set_2)
    println(JuMP.start_value(v))
    println(JuMP.start_value(δ))
    println(JuMP.start_value(p))
    println(JuMP.start_value(q))
    println(JuMP.start_value(pij))
    println(JuMP.start_value(qij))
    println(JuMP.start_value(a))
