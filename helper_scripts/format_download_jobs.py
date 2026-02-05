
import re

MC_name = 'DQ_data'
# MC_name = 'HF'
# MC_name = 'genpurp'

jobs_file = 'jobs.txt'
out_file = 'jobs_formatted.csv'

# Read all lines from the jobs file
with open('input_data/' + MC_name + '/' + jobs_file, 'r') as f:
    lines = f.readlines()

    # Join every 7 lines into a single line
    pattern = r"\d{10}\n\t"
    lines = re.split(pattern, ''.join(lines))

    lines = [line.replace('\n', '').replace('⬇️', '').strip() for line in lines] # Remove newlines and extra spaces
    lines = [line.split('\t') for line in lines]

with open('input_data/' + MC_name + '/' + out_file, 'w') as f:
    for line in lines:
        if line[-1] == 'Done':
            f.write(','.join(line) + '\n')