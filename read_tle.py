from tletools import TLE

filename = "gp.txt"

with open(filename) as f:
    content = f.readlines()

num_lines = len(content)

tles = []

i = 0
while i < num_lines:

  temp = ["Name"]
  temp.append(content[i:i+1][0])
  temp.append(content[i+1:i+2][0])
  tles.append(TLE.from_lines(*temp))
  i += 2

print(tles[1].argp)
