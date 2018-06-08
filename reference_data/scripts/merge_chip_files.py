import pandas as pd
from functools import reduce


file_list = [line.strip() for line in open("file_list.txt", 'r')]
file_list = [file + "_cleaned.txt" for file in file_list]

df_names = ["df_" + str(i) for i in range(len(file_list))]

for i in range(len(df_names)):
    df_names[i] = pd.read_table(file_list[i])

for i in range(len(df_names)):
    df_names[i] = df_names[i].round(decimals=4)


df_final = reduce(lambda left, right: pd.merge(left, right, how='outer', on='!Sample_description'), df_names)
df_final = df_final.fillna(value=-1)
df_final.to_csv('test-merged.txt', sep="\t")