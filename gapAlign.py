
#This program was adapted from my work on lab 6.1.

#This function starts everything going and calls globalAlign.
def main():
  #This gives the five inputs of globalAlign. Change them here. Note that for the function to work, s1 must be at least as long as s2.
  s1 = 'CTCGAGTCTAGAGCATCTCGAGTCTAGAGCATCTCGAGTCTAGAGCATCTCGAGTCTAGAGCATCTCGAGTCTAGA'
  s2 = 'GCATGCATGCATGCAT'
  #gapStart and gapContinue are the penalties assigned to starting and continuing a gap.
  gapStart = -1
  gapContinue = -0.1
  match = 1
  mismatch = -1
  print('Computing global alignment of strings')
  align1, align2 = gapAlign(s1, s2, match, mismatch, gapStart, gapContinue)
  print(align1)
  print(align2)
  
  return # done with main() function

#gapAlign aligns two sequences, preferring to continue gaps than open them
def gapAlign(s1, s2, match, mismatch, gapStart, gapContinue):
  #This bit builds a table to compute the alignment with
  centralTable = initializeTable(s1, s2)
  #This builds the backtrack table to compute the alignment with
  centralBacktrack = initializeTable(s1[1:], s2[1:])
  #This builds the table to track gaps on the first string
  upperTable = initializeTable(s1, s2)
  #This builds the table to map decisions relating to gaps on the first string
  upperBacktrack = initializeTable(s1[1:], s2[1:])
  #This builds the table to track gaps on the second string
  lowerTable = initializeTable(s1, s2)
  #This builds the table to map decisions relating to gaps on the second string
  lowerBacktrack = initializeTable(s1[1:], s2[1:])
  size1 = len(s1)
  size2 = len(s2)
  i = 0
  j = 0
  #tableBuilder fills the 0,0 of the blank tables
  upperTable, centralTable, lowerTable, upperBacktrack, centralBacktrack, lowerBacktrack = tableBuilder(0, 0, upperTable, centralTable, lowerTable, upperBacktrack, centralBacktrack, lowerBacktrack, match, mismatch, gapStart, gapContinue, s1, s2)
  #This section ensures that the edges of  the lower and upper tables are not favored, to prevent out-of-range errors
  for k in range(size1+1):
    upperTable[k][0] = -100
  for l in range(size2+1):
    lowerTable[0][l] = -100
  #This section completely fills the blank tables
  while j < size2-1:
    j = j + 1
    tempj = j
    tempi = i
    while tempj >= 0:
      upperTable, centralTable, lowerTable, upperBacktrack, centralBacktrack, lowerBacktrack = tableBuilder(tempi, tempj, upperTable, centralTable, lowerTable, upperBacktrack, centralBacktrack, lowerBacktrack, match, mismatch, gapStart, gapContinue, s1, s2)
      tempi = tempi + 1
      tempj = tempj - 1
  while i < size1-size2:
    i = i + 1
    tempi = i
    tempj = j
    while tempj >= 0:
      upperTable, centralTable, lowerTable, upperBacktrack, centralBacktrack, lowerBacktrack = tableBuilder(tempi, tempj, upperTable, centralTable, lowerTable, upperBacktrack, centralBacktrack, lowerBacktrack, match, mismatch, gapStart, gapContinue, s1, s2)
      tempi = tempi + 1
      tempj = tempj - 1
  while i < size1:
    i = i + 1
    tempi = i
    tempj = j
    while tempi < size1:
      upperTable, centralTable, lowerTable, upperBacktrack, centralBacktrack, lowerBacktrack = tableBuilder(tempi, tempj, upperTable, centralTable, lowerTable, upperBacktrack, centralBacktrack, lowerBacktrack, match, mismatch, gapStart, gapContinue, s1, s2)
      tempi = tempi + 1
      tempj = tempj - 1
  print('table built')
  #This section prints the tables, and is not necessary unless debugging is occuring
  #print('backtrack tables:')
  #print('upper backtrack:')
  #for i in range(len(upperBacktrack)):
    #print(upperBacktrack[i])
  #print('central backtrack:')
  #for i in range(len(centralBacktrack)):
    #print(centralBacktrack[i])
  #print('lower backtrack:')
  #for i in range(len(lowerBacktrack)):
    #print(lowerBacktrack[i])
  #print('reference tables:')
  #print('upper table:')
  #for i in range(len(upperTable)):
    #print(upperTable[i])
  #print('central table:')
  #for i in range(len(centralTable)):
    #print(centralTable[i])
  #print('lower table:')
  #for i in range(len(lowerTable)):
    #print(lowerTable[i])
        
  #These are the alignments we construct
  alignment1 = ''
  alignment2 = ''
  s1Position = len(s1)
  s2Position = len(s2)
  backtrackLevel = 0
  #This section constructs the alignment
  while s1Position > 0 and s2Position > 0:
    #print('s1Position, s2Position: ', s1Position, s2Position)
    if backtrackLevel == 0:
      if centralBacktrack[s1Position-1][s2Position-1] == 'm':
        alignment1 = s1[s1Position-1] + alignment1
        alignment2 = s2[s2Position-1] + alignment2
        s1Position = s1Position - 1
        s2Position = s2Position - 1
        #print('central m')
      elif centralBacktrack[s1Position-1][s2Position-1] == 'u':
        backtrackLevel = 1
        alignment1 = '-' + alignment1
        alignment2 = s2[s2Position-1] + alignment2
        #s1Position = s1Position - 1
        #print('central u')
      elif centralBacktrack[s1Position-1][s2Position-1] == 'd':
        backtrackLevel = -1
        alignment2 = '-' + alignment2
        alignment1 = s1[s1Position-1] + alignment1
        #s2Position = s2Position - 1
       #print('central d')
      else:
        print('error: central backtrack value of ', centralBacktrack[s1Position][s2Position])
    elif backtrackLevel == 1:
      if upperBacktrack[s1Position-1][s2Position-1] == 'd':
        backtrackLevel = 0
        s2Position = s2Position - 1
        #print('upper d')
      elif upperBacktrack[s1Position-1][s2Position-1] == 'w':
        alignment1 = '-' + alignment1
        alignment2 = s2[s2Position-1] + alignment2
        s2Position = s2Position - 1
        #print('upper w')
      else:
        print('error: upper backtrack value of ', upperBacktrack[s1Position][s2Position])
    elif backtrackLevel == -1:
      if lowerBacktrack[s1Position-1][s2Position-1] == 'u':
        backtrackLevel = 0
        s1Position = s1Position - 1
        #print('lower u')
      elif lowerBacktrack[s1Position-1][s2Position-1] == 'n':
        alignment2 = '-' + alignment2
        alignment1 = s1[s1Position-1] + alignment1
        s1Position = s1Position - 1
        #print('lower n')
      else:
        print('error: lower backtrack value of ', lowerBacktrack[s1Position][s2Position])
      #print(s1Position, s2Position)
      #print('edge reached with alignments: ', alignment1, alignment2)
  #This section ensures both strings are fully incorporated in the alignment
  while s1Position > 0:
    #print('completing s1')
    s1Position = s1Position - 1
    #print('s1Position: ', s1Position)
    alignment1 = s1[s1Position] + alignment1
    #print(alignment1)
  while s2Position > 0:
    #print('completing s2')
    s2Position = s2Position - 1
    #print('s2Position: ', s2Position)
    alignment2 = s2[s2Position] + alignment2
    #print(alignment2)
    
  #This makes sure both strings are the same length
  while len(alignment1) != len(alignment2):
    if len(alignment1) > len(alignment2):
      alignment2 = '-' + alignment2
    else:
      alignment1 = '-' + alignment1
  print('final score = ', centralTable[len(s1)][len(s2)])
  return alignment1, alignment2
  
#This function returns the match score if the two given characters 
#match and the mismatch otherwise.
def subscore(one, two, match, mismatch):
  if one == two:
    return match
  else:
    return mismatch

#This function takes coordinates and computes values for all six tables at those coordinates, referencing previous coordinates
def tableBuilder(i, j, upperTable, centralTable, lowerTable, upperBacktrack, centralBacktrack, lowerBacktrack, match, mismatch, gapStart, gapContinue, s1, s2):
  #print(i, j)
  lowerTable[i+1][j+1] = max((lowerTable[i][j+1] + gapContinue), (centralTable[i][j+1] + gapStart))
  upperTable[i+1][j+1] = max((upperTable[i+1][j] + gapContinue), (centralTable[i+1][j] + gapStart))
  centralTable[i+1][j+1] = max(lowerTable[i+1][j+1], upperTable[i+1][j+1], (centralTable[i][j] + subscore(s1[i], s2[j], match, mismatch)))
  if (lowerTable[i][j+1] + gapContinue) >= (centralTable[i][j+1] + gapStart):
    lowerBacktrack[i][j] = 'n'
  else:
    lowerBacktrack[i][j] = 'u'
  if (upperTable[i+1][j] + gapContinue) >= (centralTable[i+1][j] + gapStart):
    upperBacktrack[i][j] = 'w'
  else:
    upperBacktrack[i][j] = 'd'
  if (centralTable[i][j] + subscore(s1[i], s2[j], match, mismatch)) >= lowerTable[i+1][j+1] and (centralTable[i][j] + subscore(s1[i], s2[j], match, mismatch)) >= upperTable[i+1][j+1]:
    centralBacktrack[i][j] = 'm'
  elif upperTable[i+1][j+1] >= lowerTable[i+1][j+1]:
    centralBacktrack[i][j] = 'u'
  else:
    centralBacktrack[i][j] = 'd'
  return upperTable, centralTable, lowerTable, upperBacktrack, centralBacktrack, lowerBacktrack
  
  #This builds blank tables for tableBuilder to fill
def initializeTable(nrows, ncols):
  table = []
  for i in range(len(nrows)+1):
    addition = []
    for j in range(len(ncols)+1):
      addition = addition + [0]
    table = table + [addition]
  return table
  
main()
