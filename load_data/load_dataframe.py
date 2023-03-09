import pandas as pd
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

def load_dataframe(excel_url = "https://docs.google.com/spreadsheets/d/1XoKzawITVjYBsEfcukOobY4Z7rUq-bkbidEYGKQFOx4/edit#gid=0"):
    excel_url = excel_url.replace("/edit#gid=", "/export?format=xlsx&gid=")

    df = pd.read_excel(excel_url, 
                    usecols="A:AQ",
                    nrows = 8)

    return df