import PyCall
import Printf
import DataFrames

PyCall.pushfirst!(PyCall.pyimport("sys")."path", "");

PyCall.py"""
from load_data.load_dataframe import load_dataframe

print("bbbb")
df = load_dataframe()
"""

print("aaaa")
df = DataFrames.DataFrame();


include("./load_data/load_dataframe.jl");
pd_to_df!(PyCall.py"df", df);

println("Load BL-params")
println(df)

include("./load_data.jl");
include("./BL_model.jl");
include("./BL_RKF45.jl");
include("./BL_model.jl");

#Plots.default(fmt = :png)