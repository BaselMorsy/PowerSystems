using LinearAlgebra,DataFrames,DataFramesMeta,CSV;
import JuMP, Ipopt, Cbc
import Juniper
using Juniper,JuMP, Ipopt, Cbc

include("SystemAssembly.jl")
include("report.jl")
Sbase = 100;

#Instantiate System

system_dir = string(pwd(),"\\Systems\\5_Bus");
branch_dir = string(system_dir,"\\branch5.csv");
bus_dir = string(system_dir,"\\bus5.csv");
gen_dir = string(system_dir,"\\gen5.csv");

System_5Bus_ACOPF = SystemAssembly(branch_dir,Sbase,bus_dir,gen_dir)
System_5Bus_DCOPF = SystemAssembly(branch_dir,Sbase,bus_dir,gen_dir)
System_5Bus_ACOTS_noDC = SystemAssembly(branch_dir,Sbase,bus_dir,gen_dir)
System_5Bus_ACOTS_DC = SystemAssembly(branch_dir,Sbase,bus_dir,gen_dir)
System_5Bus_DCOTS = SystemAssembly(branch_dir,Sbase,bus_dir,gen_dir)

Solve_OPF!(System_5Bus_ACOPF,"ACOPF")
Solve_OPF!(System_5Bus_DCOPF,"DCOPF")
Solve_OTS!(System_5Bus_ACOTS_noDC,"ACOTS",false)
Solve_OTS!(System_5Bus_ACOTS_DC,"ACOTS",true)
_,_,rem_info = Solve_OTS!(System_5Bus_DCOTS,"DCOTS")
remove_line!(System_5Bus_DCOTS,rem_info)
