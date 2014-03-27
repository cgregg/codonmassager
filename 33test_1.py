#

input = open('ini5.txt') #'input' is file handle
GBoutputtest = open('outputtest.txt','w')
a=1

for line in input:
    #print a
    if a % 2 == 1:
        wordlist = line.split() #list of strings
        #print wordlist
        
        #break #break out of for loop (opposite is continue where you stayin the in for loop)
        print wordlist[-2:]
        jwordlist = wordlist[-2:]
        jwordlist = " ".join(jwordlist)
        print jwordlist
        #break
        GBoutputtest.write(jwordlist+'\n')
    a = a + 1

        
GBoutputtest.close()
