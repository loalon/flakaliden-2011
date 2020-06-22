
file="/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/data/ko2Species.tsv"
outfile = open('/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/data/ko2speciesMap.tsv', 'w', newline='\n')

myDict = {}

def invert_dict(d): 
    inverse = dict() 
    for key in d: 
        # Go through the list that is saved in the dict:
        for item in d[key]:
            # Check if in the inverted dict the key exists
            if item not in inverse: 
                # If not create a new list
                inverse[item] = [key] 
            else: 
                inverse[item].append(key) 
    return inverse
    
with open(file, 'r', encoding='utf-8') as inFile:
  for line in inFile:
    splits=line.rstrip().split('\t')
    #gene = splits[0]
    species = splits[9]
    kos = splits[1].split(',')
    myDict[species]=kos


#invMap = {v: k for k, v in myDict.items()}
invMap = invert_dict(myDict)

for k,value in invMap.items():
  outfile.write(k+'\t')
  
  counter=1
  for gen in value:
    outfile.write(gen)
    if counter < len(value):
      outfile.write('|')
      counter+=1
  outfile.write('\n')
			
outfile.close()
print(invMap)
