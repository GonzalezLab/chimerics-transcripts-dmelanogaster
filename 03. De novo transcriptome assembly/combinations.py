import sys
import itertools

file = sys.argv[1]
output = sys.argv[2]

filenames = open(file,'r').read().splitlines()

file_contents = []
for filename in filenames:
  with open(filename, 'r') as file:
    contents = file.read().splitlines()
    file_contents.append(contents)

comb=list(itertools.product(*file_contents))

with open(output, 'w') as file:
    for element in comb:
        line = ' '.join(element) + '\n'
        file.write(line)
