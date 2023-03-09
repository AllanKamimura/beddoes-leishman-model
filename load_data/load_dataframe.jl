function pd_to_df!(df_pd, df)
    for col in df_pd.columns
        df[!, col] = getproperty(df_pd, col).values
    end
    df
end