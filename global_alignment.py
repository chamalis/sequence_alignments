from matrices import *

def do_global_alignment( seqs, exchangeMatrix, penalty, print_Matrix):
    """ do pairwise global alignment using DP """
    # get the two sequences by calling the function Sequence
    seq1 = seqs[0].Sequence
    seq2 = seqs[1].Sequence
    #put "-" infront of the 2 sequences
    seq1 = "-" + seq1 #comment for later use : the first sequence will be put on the rows 
    seq2 = "-" + seq2 #comment for later use : the second sequence will be put as columns
    #the length of the two sequences is being put in the variables len1 and len2
    len1, len2 = len(seq1), len(seq2)
    #create the score_matrix, the matrix that all the calculated scores will be scored
    #create the trace_matrix that will be used for the traceback path
    score_matrix, trace_matrix = [], []
    #fill in the score_matrix, trace_matrix with len1 rows while appending it with len2 elements: len2 columns
    for i in range(len1): 
         #the score_matrix, trace_matrix will have as many as columns as the 1st sequence
         score_matrix.append([0]*len2)
         trace_matrix.append([0]*len2)
    #initialise the [0][0] element of the matrix with 0
    score_matrix[0][0] = 0
    #for every column beginning from the second column until the last column of the matrix
    for j in range(1,len2):
       #calculate the first row 
        score_matrix[0][j]= score_matrix[0][j-1]-penalty
    #for every row beggining from the second row until the last row
    for i in range(1,len1):
       # calculate the first column of the matrix
        score_matrix[i][0]= score_matrix[i-1][0]-penalty
        
   ###############################Calculate DP table######################################### 
    #for every row in the matrix beginning from the second row, as the score_matrix[0][0]=0  
    for i in range(1,len1):
        #for every column beginning from the second one
        for j in range(1,len2):
            #calculate the match
            match = score_matrix[i-1][j-1]+ exchangeMatrix[ord(seq1[i])-ord('A')][ord(seq2[j])-ord('A')]         
            #calculate the deletion
            delete = score_matrix[i-1][j] - penalty
            #calculate the insertion
            insert = score_matrix[i][j-1] - penalty
            #get the maximum value of the matching, deletion and the insertion
            score_matrix[i][j] = max(match, delete, insert)
    # show the matrix in an appearable way
    if print_Matrix == True:
        for k in range(len(score_matrix)):
            print score_matrix[k]
  
    return score_matrix

###############################Traceback################################################################################################################
#function for showing the traceback of the global alignment with input the 2 sequences, the score_matrix from the previous dunction and the penalty value
def give_global_traceback(seqs,exchangeMatrix, score_matrix, penalty):
    
    #get the seq1,seq2; put the "-"as first element in the 2 sequences; 
    seq1 = seqs[0].Sequence
    seq2 = seqs[1].Sequence
    seq1 = "-" + seq1
    seq2 = "-" + seq2
    #find the length of the two sequences
    len1, len2 = len(seq1), len(seq2)
    #initialise two strings : aseq1 and aseq2
    aseq1 =""
    aseq2 =""
    #initialise the string line_division
    line_division = ""
    # initialise i and j
    i = len1 - 1
    j = len2 - 1
    final_score = score_matrix[i][j]
    #loop over the score_matrix as long as the i and j is not 0
    while (i!=0 or j!=0):
        #check for equal residues between the 2 given sequences, if so, print a vertical line between them 
        if seq1[i] == seq2[j]:
            line_division = "|" + line_division
        # if they are not similar, print an empty " "
        else:
            line_division = " " + line_division   
        #get variable score as the element of th earray that is in the position [i][j]
        score = score_matrix[i][j]
        #get the diagonal score
        score_diag = score_matrix[i-1][j-1]
        #get the score that comes from the exactly above element of the array
        score_up = score_matrix[i-1][j]
        #get the score from the left element of the array 
        score_left = score_matrix[i][j-1]
        #check wether the score value to see where the current value came from
        if score == (score_left - penalty):
            #if the current value came from the calculation of the left score - penalty then fill in the aseq1 with a gap and the aseq2 with the value of the seq2[j]
            aseq1 = '-' + aseq1
            aseq2 = seq2[j] + aseq2
            # reduce j by 1, so that to go to the previous column 
            j -= 1
        # check wether the score comes from the calculation of the score_up decreased by the penalty
        elif score == (score_up - penalty):
            #if score equal to the value of the above calculation then calculate aseq1
            aseq1 = seq1[i] + aseq1
            #since the score comes from the value up, means that there will be a gap in the sequence 2, therefore add "-" to aseq2
            aseq2 = "-" + aseq2
            #REduce i by one to move to the previous row
            i -= 1
        #check wether the score comes from the previous diagonal element after a specific calculation is been done including th eexchangeMatrix in it.
        elif score == score_diag + exchangeMatrix[ord(seq1[i])-ord('A')][ord(seq2[j])-ord('A')]:
            # If yes, then put the the value of seq1[1] to aseq1 and do the same for aseq2
            aseq1 = seq1[i] + aseq1
            aseq2 = seq2[j] + aseq2
            #Reduce i and j by one so that to move to the previous diagonal position in the matrix, in the next loop
            i -= 1
            j -= 1
        #if nothing of the above if and elif is not true then print "ERROR"
        else:
            print "ERROR"
    #print the aligned sequence 1 
    print aseq1
    #print the vertical line if there are similar residues
    print line_division
    #print the aligned sequence 2
    print aseq2
    #print the score 
    print "score =",final_score
    #the function returns all the following elements everytime we cal her from a different point of the program 
    return aseq1, aseq2, score, seq1, seq2

########################################################## semiglobal alignment##################################
