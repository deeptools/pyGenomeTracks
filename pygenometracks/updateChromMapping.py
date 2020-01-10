from os import listdir
import os.path

my_dir = 'ChromosomeMappings'
files = [os.path.join(my_dir, f) for f in listdir(my_dir)
         if os.path.isfile(os.path.join(my_dir, f))]

bigDic = {}
for file in files:
    with open(file, 'r') as f:
        for line in f.readlines():
            try:
                a, b = line.strip().split("\t")
            except ValueError:
                pass
            else:
                if a != b:
                    if a not in bigDic:
                        bigDic[a] = []
                    if b not in bigDic[a] and a != b:
                        bigDic[a].append(b)

print("There are {} entries in the dictionary".format(len(bigDic)))
numberOfCor = [len(bigDic[a]) for a in bigDic]
print("The occurence of each is:")
print({x: numberOfCor.count(x) for x in set(numberOfCor)})
threshold = 1
large = [x for x in bigDic if len(bigDic[x]) > threshold]
print("The ones with more than one correspondance are:")
for k in large:
    print('{}:{}'.format(k, bigDic[k]))

with open('pygenometracks/_chromMapping.py', 'w') as f:
    f.write('__chromMapping__ = ' + str(bigDic))
