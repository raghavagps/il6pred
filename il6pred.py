##############################################################################
# IL6Pred is developed for predicting, desigining and scanning interleukin-6  #
# inducing peptides. It is developed by Prof G. P. S. Raghava's group. Please #
# cite: IL6Pred
# ############################################################################
import sys
import subprocess
import pkg_resources
required = {'joblib','pandas','numpy','tqdm'}
installed = {pkg.key for pkg in pkg_resources.working_set}
python = sys.executable
subprocess.check_call([python, '-m', 'pip', 'install', *required], stdout=subprocess.DEVNULL)
import argparse
import os
import sys
import numpy as np
import pandas as pd
import math
import itertools
import re
import glob
import time
import warnings
from time import sleep
from tqdm import tqdm
from sklearn.externals import joblib
from collections import Counter
from collections import Sequence
from numpy.core.umath_tests import inner1d
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description='Please provide following arguments')
## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence in FASTA format or single sequence per line in single letter code")
parser.add_argument("-o", "--output",type=str, help="Output: File for saving results by default outfile.csv")
parser.add_argument("-j", "--job",type=int, choices = [1,2,3], help="Job Type: 1:predict, 2:design and 3:scan, by default 1 (predict)")
parser.add_argument("-t","--threshold", type=float, help="Threshold: Value between 0 to 1 by default 0.5")
parser.add_argument("-w","--winleng", type=int, choices =range(5, 30), help="Window Length: 5 to 30 (scan mode only), by default 10")
parser.add_argument("-d","--display", type=int, choices = [1,2], help="Display: 1:Interleukin-6 inducing peptide, 2: All peptides, by default 1")
args = parser.parse_args()

# Function for generating all possible mutants
def seq_mutants(aa_seq):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    aa_seq = aa_seq.upper()
    mut_seq=[]
    for pos in range(0,len(aa_seq)):
        for aa in std:
            mut_seq += [aa_seq[:pos] + aa + aa_seq[pos+1:]]
    return mut_seq
# Function for generating pattern of a given length
def seq_pattern(aa_seq,win_len):
    aa_seq == aa_seq.upper
    seq_pat=[]
    for i1 in range(0, (len(aa_seq) + 1 - win_len)):
        i2 = i1 + int(win_len)
        seq_pat += [aa_seq[i1:i2]]
    return seq_pat

def feature_gen(file):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    sf = ['hydrogen_bonds','Grantham_gap3','Polarity_2','N_perc','Charge_3','NVWV_3','L','Hydrophobicity_1','A','SS_3']
    def aac_comp(file):
        filename, file_extension = os.path.splitext(file)
        f = open(filename+".aac", 'w')
        sys.stdout = f
        df = pd.read_csv(file, header = None)
        zz = df.iloc[:,0]
        print("A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,")
        for j in zz:
            for i in std:
                count = 0
                for k in j:
                    temp1 = k
                    if temp1 == i:
                        count += 1
                    composition = (count/len(j))*100
                print("%.2f"%composition, end = ",")
            print("")
        f.truncate()
    
    def bond(file) :
        tota = []
        hy = []
        Si = []
        Du = []
        b1 = []
        b2 = []
        b3 = []
        b4 = []
        bb = pd.DataFrame()
        filename, file_extension = os.path.splitext(file)
        df12 = pd.read_csv(file, header = None)
        df = pd.DataFrame(df12[0].str.upper())
        bonds=pd.read_csv("Data/bonds.csv")
        for i in range(0,len(df)) :
            tot = 0
            h = 0
            S = 0
            D = 0
            tota.append([i])
            hy.append([i])
            Si.append([i])
            Du.append([i])
            for j in range(0,len(df[0][i])) :
                temp = df[0][i][j]
                for k in range(0,len(bonds)) :
                    if bonds.iloc[:,0][k] == temp :
                        tot = tot + bonds.iloc[:,1][k]
                        h = h + bonds.iloc[:,2][k]
                        S = S + bonds.iloc[:,3][k]
                        D = D + bonds.iloc[:,4][k]
            tota[i].append(tot)
            hy[i].append(h)
            Si[i].append(S)
            Du[i].append(D)
        for m in range(0,len(df)) :
            b1.append(tota[m][1])
            b2.append(hy[m][1])
            b3.append(Si[m][1])
            b4.append(Du[m][1])
        
        bb["total_bonds"] = b1 
        bb["hydrogen_bonds"] = b2
        bb["single_bond"] = b3
        bb["double_bond"] = b4 
        
        #bb.to_csv(filename+".bonds_info")
        return bb
    
    ###########################atom###############
    def atc(file) :
        import pandas as pd
        atom=pd.read_csv("Data/atom.csv",header=None)
        filename, file_extension = os.path.splitext(file)
        at=pd.DataFrame()
        i = 0
        C_atom = []
        H_atom = []
        N_atom = []    
        O_atom = []
        S_atom = []
     
        while i < len(atom):
            C_atom.append(atom.iloc[i,1].count("C"))
            H_atom.append(atom.iloc[i,1].count("H"))
            N_atom.append(atom.iloc[i,1].count("N"))
            O_atom.append(atom.iloc[i,1].count("O"))
            S_atom.append(atom.iloc[i,1].count("S"))
            i += 1
        atom["C_atom"]=C_atom  
        atom["O_atom"]=O_atom  
        atom["H_atom"]=H_atom  
        atom["N_atom"]=N_atom 
        atom["S_atom"]=S_atom 
    ##############read file ##########	
        df12 = pd.read_csv(file, header = None)
        test1 = pd.DataFrame(df12[0].str.upper())
        #test1 = pd.read_csv(file,header=None)
        dd = []
        for i in range(0, len(test1)):
            dd.append(test1[0][i].upper())
        test = pd.DataFrame(dd)
        count_C = 0
        count_H = 0
        count_N = 0
        count_O = 0
        count_S = 0
        count = 0
        i1 = 0
        j = 0
        k = 0
        C_ct = []
        H_ct = []
        N_ct = []
        O_ct = []
        S_ct = []
        while i1 < len(test) :
            while j < len(test[0][i1]) :
                while k < len(atom) :
                    if test.iloc[i1,0][j]==atom.iloc[k,0].replace(" ","") :
                        count_C = count_C + atom.iloc[k,2]
                        count_H = count_H + atom.iloc[k,3]
                        count_N = count_N + atom.iloc[k,4]
                        count_O = count_O + atom.iloc[k,5]
                        count_S = count_S + atom.iloc[k,6]
                    #count = count_C + count_H + count_S + count_N + count_O
                    k += 1
                k = 0
                j += 1
            C_ct.append(count_C)
            H_ct.append(count_H)
            N_ct.append(count_N)
            O_ct.append(count_O)
            S_ct.append(count_S)
            count_C = 0
            count_H = 0
            count_N = 0
            count_O = 0
            count_S = 0
            j = 0
            i1 += 1
        test["C_count"]=C_ct  
        test["H_count"]=H_ct  
        test["N_count"]=N_ct 
        test["O_count"]=O_ct  
        test["S_count"]=S_ct     
    
        ct_total = []
        m = 0
        while m < len(test) :
            ct_total.append(test.iloc[m,1] + test.iloc[m,2] + test.iloc[m,3] + test.iloc[m,4] + test.iloc[m,5])
            m += 1
        test["count"]=ct_total    
    ##########final output#####
        final = pd.DataFrame()
        n = 0
        p = 0
        C_p = []
        H_p = []
        N_p = []
        O_p = []
        S_p = []
        while n < len(test):
            C_p.append((test.iloc[n,1]/test.iloc[n,6])*100)
            H_p.append((test.iloc[n,2]/test.iloc[n,6])*100)
            N_p.append((test.iloc[n,3]/test.iloc[n,6])*100)
            O_p.append((test.iloc[n,4]/test.iloc[n,6])*100)
            S_p.append((test.iloc[n,5]/test.iloc[n,6])*100)
            n += 1
        final["C_perc"] = C_p
        final["H_perc"] = H_p
        final["N_perc"] = N_p
        final["O_perc"] = O_p
        final["S_perc"] = S_p
    
        #(final.round(2)).to_csv(filename+".atc")
        return final.round(2)
    
    ########################################atom_bond#################
    def atom_bond(file) :
        filename, file_extension = os.path.splitext(file)
        file_atc=atc(file)
        file_bonds=bond(file)
        df1=file_atc.iloc[:,0:6]
        df2=file_bonds.iloc[:,0:5]
        df3 = pd.concat([df1,df2],axis=1)
        df3.to_csv(filename+".atom_bond", index=None)	
    	
    ###################################soc#######################################
    def soc(file,gap):
        ff = []
        filename, file_extension = os.path.splitext(file)
        df = pd.read_csv(file, header = None)
        df2 = pd.DataFrame(df[0].str.upper())
        for i in range(0,len(df2)):
            ff.append(len(df2[0][i]))
        if min(ff) < gap:
            print("Error: All sequences' length should be higher than :", gap)
            return 0
        mat1 = pd.read_csv("Data/Schneider-Wrede.csv", index_col = 'Name')
        mat2 = pd.read_csv("Data/Grantham.csv", index_col = 'Name')
        h1 = []
        h2 = []
        for n in range(1, gap+1):
            h1.append('Schneider_gap' + str(n))
        for n in range(1, gap + 1):
            h2.append('Grantham_gap' + str(n))
        s1 = []
        s2 = []
        for i in range(0,len(df2)):
            for n in range(1, gap+1):
                sum = 0
                sum1 =0
                sum2 =0
                sum3 =0
                for j in range(0,(len(df2[0][i])-n)):
                    sum = sum + (mat1[df2[0][i][j]][df2[0][i][j+n]])**2
                    sum1 = sum/(len(df2[0][i])-n)
                    sum2 = sum2 + (mat2[df2[0][i][j]][df2[0][i][j+n]])**2
                    sum3 = sum2/(len(df2[0][i])-n)
                s1.append(sum1)
                s2.append(sum3)
        zz = np.array(s1).reshape(len(df2),gap)
        zz2 = np.array(s2).reshape(len(df2),gap)
        zz3 = round(pd.concat([pd.DataFrame(zz, columns = h1),pd.DataFrame(zz2,columns = h2)], axis = 1),4)
        zz3.to_csv(filename+".soc", index = None, encoding = 'utf-8')
    #######################################CeTD#######################
    def ctd(file):
        attr=pd.read_csv("Data/aa_attr_group.csv", sep = "\t")
        filename, file_extension = os.path.splitext(file)
        df1 = pd.read_csv(file, header = None)
        df = pd.DataFrame(df1[0].str.upper())
        n = 0
        stt1 = []
        m = 1
        for i in range(0,len(attr)) :
            st =[]
            stt1.append([])
            for j in range(0,len(df)) :
                stt1[i].append([])
                for k in range(0,len(df[0][j])) :        
                    while m < 4 :
                        while n < len(attr.iloc[i,m]) :
                            if df[0][j][k] == attr.iloc[i,m][n] :
                                st.append(m)
                                stt1[i][j].append(m)
                            n += 2
                        n = 0  
                        m += 1
                    m = 1 
    #####################Composition######################     
        f = open(filename+".comp", 'w')
        sys.stdout = f
        std = [1,2,3]
        print("1,2,3,")
        for p in range (0,len(df)) :
            for ii in range(0,len(stt1)) :
                #for jj in stt1[ii][p]:
                for pp in std :
                    count = 0
                    for kk in stt1[ii][p] :
                        temp1 = kk
                        if temp1 == pp :
                            count += 1
                        composition = (count/len(stt1[ii][p]))*100    
                    print("%.2f"%composition, end = ",")
                print("")  
        f.truncate() 
    
    #################################Transition#############
        tt = []
        tr=[]
        kk =0
        for ii in range(0,len(stt1)) :
            tt = []
            tr.append([])
            for p in range (0,len(df)) :
                tr[ii].append([])
                while kk < len(stt1[ii][p]) :
                    if kk+1 <len(stt1[ii][p]):
                    #if  stt1[ii][p][kk] < stt1[ii][p][kk+1] or stt1[ii][p][kk] > stt1[ii][p][kk+1]: # condition for adjacent values
                        tt.append(stt1[ii][p][kk])
                        tt.append(stt1[ii][p][kk+1])
                        tr[ii][p].append(stt1[ii][p][kk])
                        tr[ii][p].append(stt1[ii][p][kk+1])
                        
                    kk += 1
                kk = 0    
            
        pp = 0
        xx = []
        xxx = []
        for mm in range(0,len(tr)) :
            xx = []
            xxx.append([])
            for nn in range(0,len(tr[mm])):
                xxx[mm].append([])
                while pp < len(tr[mm][nn]) :
                    xx .append(tr[mm][nn][pp:pp+2])
                    xxx[mm][nn].append(tr[mm][nn][pp:pp+2])
                    pp+=2
                pp = 0
            
        f1 = open(filename+'.trans', 'w')
        sys.stdout = f1
        std1 = [[1,1],[1,2],[1,3],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3]]
        print("1->1,1->2,1->3,2->1,2->2,2->3,3->1,3->2,3->3,")
        for rr in range(0,len(df)) :
            for qq in range(0,len(xxx)):
                for tt in std1 :
                    count = 0
                    for ss in xxx[qq][rr] :
                        temp2 = ss
                        if temp2 == tt :
                            count += 1
                    print(count, end = ",")  
                print("")
        f1.truncate()       
       
        #################################Distribution#############
        c_11 = []
        c_22 = []
        c_33 = []
        zz = []
        print("0% 25% 50% 75% 100%")
        for x in range(0,len(stt1)) :
            #c_11.append([])
            c_22.append([])
            #c_33.append([])
            yy_c_1 = []
            yy_c_2 = []
            yy_c_3 = []
            ccc = []
        
            k = 0
            j = 0
            for y in range(0,len(stt1[x])):
                #c_11[x].append([])
                c_22[x].append([])
                for i in range(1,4) :
                    cc = []
                    c1 = [index for index,value in enumerate(stt1[x][y]) if value == i]
                    c_22[x][y].append(c1)
        cc = []
        for ss in range(0,len(df)):
            for uu in range(0,len(c_22)):
                for mm in range(0,3):
                    for ee in range(0,101,25):
                        k = (ee*(len(c_22[uu][ss][mm])))/100
                        cc.append(math.floor(k))
        f2 = open(filename+'.dist', 'w')
        sys.stdout = f2
        print("0% 25% 50% 75% 100%")
        for i in range (0,len(cc),5):
            print(*cc[i:i+5])
        f2.truncate()           
        head = []
        header1 = ['Hydrophobicity','NVWV','Polarity','Polarizability','Charge','SS','SA']
        for i in header1:
            for j in range(1,4):
                head.append(i+'_'+str(j))
        df11 = pd.read_csv(filename+".comp")
        df_1 = df11.iloc[:,:-1]
        zz = pd.DataFrame()
        for i in range(0,len(df_1),7):
            zz = zz.append(pd.DataFrame(pd.concat([df_1.loc[i],df_1.loc[i+1],df_1.loc[i+2],df_1.loc[i+3],df_1.loc[i+4],df_1.loc[i+5],df_1.loc[i+6]],axis=0)).transpose()).reset_index(drop=True)
        zz.columns = head
        zz.to_csv(filename+".ctd_comp", index=None, encoding='utf-8')
        head2 = []
        header2 = ['1->1','1->2','1->3','2->1','2->2','2->3','3->1','3->2','3->3']
        for i in header2:
            for j in range(1,8):
                head2.append(i+'_'+str(j))
        df12 = pd.read_csv(filename+".trans")
        df_2 = df12.iloc[:,:-1]
        ss = pd.DataFrame()
        for i in range(0,len(df_2),7):
            ss = ss.append(pd.DataFrame(pd.concat([df_2.loc[i],df_2.loc[i+1],df_2.loc[i+2],df_2.loc[i+3],df_2.loc[i+4],df_2.loc[i+5],df_2.loc[i+6]],axis=0)).transpose()).reset_index(drop=True)
        ss.columns = head2
        ss.to_csv(filename+".ctd_trans", index=None, encoding='utf-8')
        head3 = []
        header3 = ['0%','25%','50%','75%','100%']
        header4 = ['Hydrophobicity','NVWV','Polarity','Polarizability','Charge','SS','SA']
        for j in range(1,4):
            for k in header4:
                for i in header3:
                    head3.append(i+'_'+str(j)+'_'+k)
        df_3 = pd.read_csv(filename+".dist", sep=" ")
        rr = pd.DataFrame()
        for i in range(0,len(df_3),21):
            rr = rr.append(pd.DataFrame(pd.concat([df_3.loc[i],df_3.loc[i+1],df_3.loc[i+2],df_3.loc[i+3],df_3.loc[i+4],df_3.loc[i+5],df_3.loc[i+6],df_3.loc[i+7],df_3.loc[i+8],df_3.loc[i+9],df_3.loc[i+10],df_3.loc[i+11],df_3.loc[i+12],df_3.loc[i+13],df_3.loc[i+14],df_3.loc[i+15],df_3.loc[i+16],df_3.loc[i+17],df_3.loc[i+18],df_3.loc[i+19],df_3.loc[i+20]],axis=0)).transpose()).reset_index(drop=True)
        rr.columns = head3
        rr.to_csv(filename+".ctd_dist", index=None, encoding='utf-8')    
    ##############################################################################
    filename, file_extension = os.path.splitext(file)
    aac_comp(file)
    atom_bond(file)
    ctd(file)
    soc(file,3)
    df1 = pd.read_csv(filename+".aac")
    df11 = pd.read_csv(filename+".atom_bond")
    df41 = pd.read_csv(filename+".soc")
    df60 = pd.read_csv(filename+".ctd_comp")
    df61 = pd.read_csv(filename+".ctd_trans")
    df62 = pd.read_csv(filename+".ctd_dist")
    df = pd.concat([df1.iloc[:,:-1],df11.iloc[:,:]],axis=1)
    df1 = pd.concat([df41.iloc[:,:],df60.iloc[:,:-1],df61.iloc[:,:-1],df62.iloc[:,:]],axis=1)        
    dff = pd.read_csv(file)
    df2 = pd.concat([df,df1],axis=1)		
    xx = len(dff)+1
    df3 = df2.iloc[:xx,:]
    df4 = df3[sf]
    df4.columns = ['BTC_H','SCO3_G3','CeTD_PO2','ATC_N','CeTD_CH3','CeTD_VW3','AAC_L','CeTD_HB1','AAC_A','CeTD_SS3']
    ############################################		
    file_path = file
    while not os.path.exists(file_path):
        time.sleep(10)
    
    if os.path.isfile(file_path):
    #Verifies CSV file was created, then deletes unneeded files.
        for CleanUp in glob.glob(filename+'*'):
            #print(CleanUp)
            if not CleanUp.endswith(file):
                os.remove(CleanUp)	
    				
    filelist=glob.glob("Binary*")
    for file_1 in filelist:
        os.remove(file_1)
    filelist1=glob.glob("AAIndex*")
    for file1 in filelist1:
        os.remove(file1)
    return df4

def adjusted_classes(y_scores, t):
    return [1 if y >= t else -1 for y in y_scores]

def Perform_testing(clf,name,X,t):
    Y_pred = clf.predict(X)
    Y_scores=[]
    Y_scores=clf.predict_proba(X)[:,-1]
    Y_pred = adjusted_classes(Y_scores,t)
    return Y_pred,Y_scores

print('##############################################################################')
print('# This program IL6Pred is developed for predicting, desigining and scanning  #')
print('# interleukin-6 inducing peptides, developed by Prof G. P. S. Raghava group. #')
print('# ############################################################################')

# Parameter initialization or assigning variable for command level arguments

Sequence= args.input        # Input variable 
 
# Output file 
 
if args.output == None:
    result_filename= "outfile.csv" 
else:
    result_filename = args.output
         
# Threshold 
if args.threshold == None:
        Threshold = 0.11
else:
        Threshold= float(args.threshold)
# Job Type 
if args.job == None:
        Job = int(1)
else:
        Job = int(args.job)
# Window Length 
if args.winleng == None:
        Win_len = int(10)
else:
        Win_len = int(args.winleng)

# Display
if args.display == None:
        dplay = int(1)
else:
        dplay = int(args.display)


print('Summary of Parameters:')
print('Input File: ',Sequence,'; Threshold: ', Threshold,'; Job Type: ',Job)
print('Output File: ',result_filename,'; Window Length: ',Win_len,'; Display: ',dplay)



    
#------------------ Read input file ---------------------
def load_model(path):
    clf = joblib.load(path)
    return clf

f=open(Sequence,"r")
len1 = f.read().count('>')
f.close()

f=open(Sequence,"r")
seqs=[]
seqid=[]
str1='';
header=0
line=0
if len1 >= 1: # read fasta file
    for l in f:
        if l.startswith('>'):
            if header != 0:
                seqs += [str1]
                str1 = '';
            header = 1
            line+=1
            seqid += [l.rstrip()]
        else:
            str1 += l.rstrip()
    seqs += [str1]
else: # read single line file    
    for l in f:
        if len(l) >= 8 and len(l)<= 25:
            seqs+=[l.strip()]
            seqid += ['>Seq_' + str(line)]
        line+=1
f.close()

fout= open(result_filename,"w+")

i1 = 0
#======================= Prediction Module start from here =====================
if Job == 1:
    fout.write('# Sequence_ID, Sequence, Prediction, Score\n')
    clf=load_model('RF_model')
    for Sequence in tqdm(seqs):
        header = seqid[i1]
        if len(Sequence) >= 8: 
            if len(Sequence) >= 26:
                Sequence = Sequence[0:25]
            ss = []
            ss.append(Sequence)
            pd.DataFrame(ss).to_csv(Sequence,index=None,header=False)
            X = feature_gen(Sequence)
            Y_pred,Y_score=Perform_testing(clf,'RF',X,Threshold)
            Y_score = np.round(Y_score,2)
            flag=""
            if Y_pred[0]==1:
                flag='IL-6 inducer'
            else:
                flag='IL-6 non-inducer'
            if dplay == 1:
                if Y_pred[0]==1:
                    fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
            else:
                fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
        os.remove(Sequence)
        i1 = i1 +1
#===================== Design Model Start from Here ======================
elif Job == 2:
    fout.write('# Sequence_ID, Sequence, Prediction, Score\n')
    print('==== Designing Peptides: Processing sequences please wait ...')
    i1 = 0
    for Sequence1 in seqs:
        pat_seq=[]
        header = seqid[i1]
        print('SeqID : ',header,'# under process ...')
        if len(Sequence1) >= 8: 
            fout.write("#Main Sequence (SequenceId, Sequence): %s,%s\n" % (header,Sequence1))
            fout.write('# Sequence_ID, pattern, Prediction, Score\n')
            if len(Sequence1) >= 26:
                Sequence1 = Sequence1[0:25]
            pat_seq = seq_mutants(Sequence1)
            for Sequence in tqdm(pat_seq):
                clf=load_model('RF_model')
                ss = []
                ss.append(Sequence)
                pd.DataFrame(ss).to_csv(Sequence,index=None,header=False)
                X = feature_gen(Sequence)
                Y_pred,Y_score=Perform_testing(clf,'RF',X,Threshold)
                Y_score = np.round(Y_score,2)
                flag=""
                if Y_pred[0]==1:
                    flag='IL-6 inducer'
                else:
                    flag='IL-6 non-inducer'
                if dplay == 1:
                    if Y_pred[0]==1:
                        fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                else:
                    fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                os.remove(Sequence)
        i1 = i1 + 1
#=============== Scan Model start from here ==================

else:
    print('==== Scanning Peptides: Processing sequences please wait ...')
    i1 = 0
    for Sequence1 in seqs:
        pat_seq=[]
        header = seqid[i1]
        print('SeqID : ',header,'# under process ...')
        if len(Sequence1) >= Win_len: 
            fout.write("#Main Sequence (SequenceId, Sequence): %s,%s\n" % (header,Sequence1))
            fout.write('# Sequence_ID, pattern, Prediction, Score\n')
            pat_seq = seq_pattern(Sequence1,Win_len)
            for Sequence in tqdm(pat_seq):
                clf=load_model('RF_model')
                ss = []
                ss.append(Sequence)
                pd.DataFrame(ss).to_csv(Sequence,index=None,header=False)
                X = feature_gen(Sequence)
                Y_pred,Y_score=Perform_testing(clf,'RF',X,Threshold)
                Y_score = np.round(Y_score,2)
                flag=""
                if Y_pred[0]==1:
                    flag='IL-6 inducer'
                else:
                    flag='IL-6 non-inducer'
                if dplay == 1:
                    if Y_pred[0]==1:
                        fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                else:
                    fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                os.remove(Sequence)
        i1 = i1 +1
fout.close()

print('\n======= Thanks for using IL6Pred. Your results are stored in file :',result_filename,' =====\n\n')
print('Please cite: IL6Pred\n')
