# File for use with testing of the scripts
import os 
print("Path at terminal when executing this file")
print(os.getcwd() + "\n")

print("This file full path (following symlinks)")
full_path = os.path.realpath(__file__)
print(full_path + "\n")
print("This file directory only")
print(os.path.dirname(full_path))

'''
with open("/Users/gonzalgo/Downloads/grland_04583_16033_003_160324_ALTTBB_HH_03.ann", "r") as f:
        lines = f.readlines()
        print(type(lines))
        a = []
        print(type(a))
        print("start: " + lines[61])
        print(lines[62])
        print(lines[63])
        print(lines[64])
        x = 0
        
        for line in lines:
            #print(line + str(x))
            
            
            if line.startswith("GRD Latitude Lines"):
                l = f.readline()
                print(l)
                print(x)
                
            
            x += 1
        
'''
            