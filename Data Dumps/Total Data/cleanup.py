import os

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path += '/Timestamps/'
print(dir_path)

list = os.listdir(dir_path) # dir is your directory path
number_files = len(list)
print(number_files)

for i in range(1,number_files+1):
    if i%100 != 0 and i<int(number_files/100)*100 and i>100:
        os.remove(dir_path + '%d.zip' %i)
    else:
        print("%d.zip not deleted." %i)
