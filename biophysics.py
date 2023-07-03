import math
import numpy as np

rna = open("RNA.asa", "r")
protein = open("PROTEIN.asa", "r")

rnatext=rna.read()
rna.close()
proteintext=protein.read()
protein.close()

rs = rnatext.split("\n")
ps= proteintext.split("\n")
j=0
aa_list= ["ALA", "ARG", "ASN", "ASP","CYS", "GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","GLU"]
base_list = ["A","U","G","C"]
k = {'A': {"ALA":0, "ARG":0, "ASN":0, "ASP":0,"CYS":0, "GLN":0,"GLY":0,"HIS":0,"ILE":0,"LEU":0,"LYS":0,"MET":0,"PHE":0,"PRO":0,"SER":0,"THR":0,"TRP":0,"TYR":0,"VAL":0,"GLU":0},
     'U':{"ALA":0, "ARG":0, "ASN":0, "ASP":0,"CYS":0, "GLN":0,"GLY":0,"HIS":0,"ILE":0,"LEU":0,"LYS":0,"MET":0,"PHE":0,"PRO":0,"SER":0,"THR":0,"TRP":0,"TYR":0,"VAL":0,"GLU":0},
     'C':{"ALA":0, "ARG":0, "ASN":0, "ASP":0,"CYS":0, "GLN":0,"GLY":0,"HIS":0,"ILE":0,"LEU":0,"LYS":0,"MET":0,"PHE":0,"PRO":0,"SER":0,"THR":0,"TRP":0,"TYR":0,"VAL":0,"GLU":0},
     'G':{"ALA":0, "ARG":0, "ASN":0, "ASP":0,"CYS":0, "GLN":0,"GLY":0,"HIS":0,"ILE":0,"LEU":0,"LYS":0,"MET":0,"PHE":0,"PRO":0,"SER":0,"THR":0,"TRP":0,"TYR":0,"VAL":0,"GLU":0}} 
# row: 1. ala 2. arg 3. asm 4. asp 5. cys 6. gln 7. gly 8.his 9.ile 10.leu 11.lys 12.met 13.phe 14.pro 15.ser 16.thr 17.trp 18.tyr 19.val 20.glu
# column: 1.A 2.U 3.C 4.G

na = 0 
nu=0 
nc=0 
ng =0

for i in range(len(rs)):
    for j in range(len(ps)):
        rsa = rs[i].split()
        psa = ps[j].split()
        if(rsa[9] and psa[9]):
            a=[float(rsa[6]), float(rsa[7]), float(rsa[8])]
            b=[float(psa[6]),float(psa[7]), float(psa[8])]
            # threshold = pow((pow((float(psa[6])-float(rsa[6])),2) + pow((float(psa[7])-float(rsa[7])),2) + pow((float(psa[8])-float(rsa[8])),2)), 0.5)
            if math.dist(a,b)<=4:
                match rsa[3]:
                    case "A":
                        na=na+1
                        k['A'][psa[3]]=k['A'][psa[3]]+1
                    case "U":
                        nu=nu+1
                        k['U'][psa[3]]=k['U'][psa[3]]+1
                    case "G":
                        ng=ng+1
                        k['G'][psa[3]]=k['G'][psa[3]]+1
                    case "C":
                        nc=nc+1
                        k['C'][psa[3]]=k['C'][psa[3]]+1                      
        j=j+1
propensity_matrix = {}

# NOW WE FIND NSQ,NIPQ,SUMS, ETC
# all possible interactions sumNIpq
# sumNIpq = 0
# for rna_base in k.keys():
#     for aa in k[rna_base].keys():
#         sumNIpq += k[rna_base][aa]

# for nsq 
for rna_base in base_list:
    freqs = list(k[rna_base].values())
    total_aa = sum(freqs)
     # number of times base interact any aa
    # print("number of times ", base, " interact with any amino acid " , nsq)
    # Nsp finding for all the amino acids
    # for aa in k['A'].keys():
    #     Nsp = 0
    #     for rna_base in k.keys():
    #         freqs = list(k[rna_base].values())
    #         nsq = sum(freqs)
    #         if aa in k[rna_base]:
    #             Nsp += k[rna_base][aa]
        # print("Total number of interactions of", aa, "with any RNA molecule (Nsp):", total_interactions)
        # propensities = [ NIpq/sumNIpq / ((Nsp / sumNIpq) * (nsq / sumNIpq)) for NIpq in freqs if Nsp!=0]

    propensities = [f/total_aa for f in freqs]
    propensity_matrix[rna_base] = dict(zip(k[rna_base].keys(), propensities))

# print("the propensity matrix is \n", propensity_matrix)
# Nlpq = 15
# Nsp = 

# def calculate_propensity(Nlpq, Nsp, sumpNsp, Nsq, sumqNsq):
#     propensity = (Nlpq / (Nsp / sumpNsp)) * (Nsq / sumqNsq)
#     return propensity

# propensity = calculate_propensity(Nlpq, Nsp, sumpNsp, Nsq, sumqNsq)
# print(propensity)



# Define the gas constant and temperature
RT = 0.59

# Define the interaction matrix
P = np.array([[propensity_matrix[key][aa] for aa in aa_list] for key in base_list])
# Calculate deltaG^statpq
deltaG = -RT * np.log(P)
 
# Print the result in tabular form
print("   ", "   ".join(aa_list))
for i, b in enumerate(base_list):
    row = [f"{deltaG[i,j]:.2f}" for j in range(len(aa_list))]
    print(b, "   ".join(row))