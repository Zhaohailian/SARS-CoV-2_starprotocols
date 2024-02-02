import os,sys
import pandas as pd

df = pd.read_csv(sys.argv[1], sep='\t', header=None)
sorted_df = df.sort_values(by=9)  ## 按照pvalue排序

#sorted_df.to_csv('sorted_file.txt', sep='\t', header=None, index=False)

print(sorted_df)

min_value = sorted_df.iloc[:5, 0].str[1:].astype(int).min()  ##输出pvalue排名前5中energy最小的结构
#print(min_value)

def extract_energy_range(file_path, start_index):
 with open(file_path, 'r') as file:
  lines = file.readlines()
  energy_indices = [i for i, line in enumerate(lines) if 'ENERGY' in line]
  if start_index < len(energy_indices):
   start_energy_index = energy_indices[start_index - 1]
   end_energy_index = energy_indices[start_index] if start_index < len(energy_indices) else len(lines)
   result_lines = lines[start_energy_index:end_energy_index]
   result = ''.join(result_lines)
   print(result)
  else:
   print("Invalid start index.")

extract_energy_range(sys.argv[2], min_value)


