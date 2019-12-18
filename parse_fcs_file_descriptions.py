import re
import pandas
array = []

flag = False
with open("merged-file.txt") as open_file:
    for line in open_file:
        line=line.strip()
        if line.startswith("Huseyin"):
            flag = True
            array2 = []
        if flag and line:
            array2.append(line)
            if line.startswith("micro"):
                flag = False
                array.append(array2)
array
columns=["filename","date","strain","plasmid","backbone","iptg"]
df=pandas.DataFrame([[row[i] for i in [0, 0, 2, 4, 4, 5]] for row in array],columns=columns)

df.filename=df.filename.replace({'.mqd':'.fcs'}, regex=True)
df.date=df.date.replace({'\.[0-9]{4}.mqd':''}, regex=True)
df.date=df.date.replace({'Huseyin':''}, regex=True)
df.backbone=df.backbone.replace({'::.*':''}, regex=True)
df.plasmid=df.plasmid.replace({'.*::':''}, regex=True)


df.strain=df.strain.replace('CC118λPir','CC118Lpir')
df.strain=df.strain.replace('DH5α','DH5alpha')
df.strain=df.strain.replace('KT 2440','KT2440')

df.plasmid=df.plasmid.str.replace("-","_")
df.plasmid=df.plasmid.str.capitalize()
df.plasmid=df.plasmid.replace('Phif_pone','Phif_p1')
df.iptg=df.iptg.replace(0,1)

df.to_csv("file_description.csv",index=False)
