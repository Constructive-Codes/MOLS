from pathlib import Path

for file_path in Path('.').glob('*.txt'):
    with open(str(file_path),'r',encoding='utf-8') as file:
        mols = []
        square = []
        firstflag = True
        for line in file:
            if line[0]!='#':
                row = [int(x) for x in line.split()]
                if firstflag:
                    [n,m] = row
                    firstflag = False
                elif len(row)>0:
                    square.append(row)
                    if len(square)==n:
                        mols.append(square)
                        square = []


    for a in range(0,m):
        for i in range(n):
            col = {}
            row = {}
            for j in range(n):
                if int(mols[a][i][j])!=mols[a][i][j] or mols[a][i][j]<0 or mols[a][i][j]>=n :
                    print("Violation: " + str(file_path))
                if mols[a][i][j] in row:
                    print("Violation: " + str(file_path))
                row[mols[a][i][j]] = 1
                if mols[a][j][i] in col:
                    print("Violation: " + str(file_path))
                col[mols[a][j][i]] = 1
                    
    for a in range(0,m):
        for b in range(a+1,m):
            pairs = {}
            for i in range(n):
                for j in range(n):
                    pair = (mols[a][i][j],mols[b][i][j])
                    if pair in pairs:
                        print("Violation: " + str(file_path))
                    pairs[pair] = 1
            if len(pairs)!=n*n:
                print("Violation: " + str(file_path))
    print(str(file_path) + " completed " + str(n) + " " + str(m))
    
            
        

    
