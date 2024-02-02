import os,sys
fh=open(sys.argv[1],'r')
length=int(sys.argv[2])

for line in fh:
 s=line.strip().split('\t')
 if s[7]=="VirusHeadGenomeTail":
  if s[5]=="+":
   
   if s[6]=="Plus":
    strand="+"
   elif s[6]=="Minus":
    strand="-"
   outline='\t'.join([s[3],s[4],str(int(s[4])+length),s[8],"255",strand])
   print(outline)

  elif s[5]=="-":
   
   if s[6]=="Plus":
    strand="+"
   elif s[6]=="Minus":
    strand="-"
   outline='\t'.join([s[3],str(int(s[4])-length),s[4],s[8],"255",strand])
   print(outline)
 
 elif s[7]=="GenomeHeadVirusTail":
  if s[5]=="+":
   if s[6]=="Plus":
    strand="+"
   elif s[6]=="Minus":
    strand="-"
   outline='\t'.join([s[3],str(int(s[4])-length),s[4],s[8],"255",strand])
   print(outline)
  elif s[5]=="-":
   if s[6]=="Plus":
    strand="+"
   elif s[6]=="Minus":
    strand="-"
   outline='\t'.join([s[3],s[4],str(int(s[4])+length),s[8],"255",strand])
   print(outline)

fh.close()
