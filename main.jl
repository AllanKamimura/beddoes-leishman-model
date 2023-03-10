import Dierckx
import PyCall
import Printf
import MAT
import LinearAlgebra
import BenchmarkTools
import DataFrames
import Statistics
import CSV


PyCall.pushfirst!(PyCall.pyimport("sys")."path", "");

PyCall.py"""
from load_data.load_dataframe import load_dataframe

df = load_dataframe()
"""

print("aaaa")
df = DataFrames.DataFrame();

include("./load_data/load_dataframe.jl");
pd_to_df!(PyCall.py"df", df);

println("Load BL-params")
println(df)

include("./load_data.jl");
include("./BL_dxdt.jl");
include("./BL_RKF45.jl");
include("./BL_model.jl");

Plots.default(fmt = :png)
BL_model(10022, df, true)
