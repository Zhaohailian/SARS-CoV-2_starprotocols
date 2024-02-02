import os,sys
fh1=open(sys.argv[1],'r')
fh2=open(sys.argv[2],'r')

notoverlap={}
for line in fh1:
 s=line.strip().split('\t')
 if s[6] not in notoverlap:
  notoverlap[s[6]]=1

for line in fh2:
 s=line.strip().split('\t')
 info=s[0].split('_')
 readid='_'.join([info[0],info[1],info[2]])
 if readid in notoverlap:
  print(line.strip())

fh1.close()
fh2.close()

