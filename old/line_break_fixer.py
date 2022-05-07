with open ("output.txt", "r") as myfile: 
    s = myfile.read().replace('\n','') 
s_new = re.sub(r'\dT\d+-\d', r'\n\g<0>',s)
with open('output_fixed.txt','a') as myfile: 
    myfile.write(s_new) 
