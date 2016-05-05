'''
collection of functions for rule-based approach.

Created on Jan 8, 2014

@author: HELE
'''

blurbp1=15 #for difference of breakpoint loci between c1 and c2. 10pb
blurbp2=100 #for difference of breakpoint loci between two continuous lines. 150bp
blurbp3=15 #for novel insertion detection
delimeter="\t"
htmlStartDiv='<div>'
htmlEndDiv='</div>\n'

class SVtypes:
    error=0
    tandem=1
    deletion=2
    inversion=4
    interspersedDuplication=tandem+deletion
    translocation=8
    novelInsertion=16
    interchromosomal=32

def getTypeFromSting(classification):
    if classification=="DUP":
        return SVtypes.tandem
    elif classification=="DEL":
        return SVtypes.deletion
    elif classification=="INV":
        return SVtypes.inversion
    elif classification=="INTDUP":
        return SVtypes.interspersedDuplication
    elif classification=="TRX":
        return SVtypes.translocation
    elif classification=="INS":
        return SVtypes.novelInsertion
    elif classification=="INTRX":
        return SVtypes.interchromosomal
    else:
        return SVtypes.error
    
# def getResultString(result):
    # string=""
    # bpstr=""
    # r=result[0]
    # if r==-1: # deleted by merging two lines
        # return
    # if r==-2:
        # return
        # string="deleted for translocation"
    # if r==SVtypes.tandem:
        # string="DUP"
    # if r==SVtypes.deletion:
        # string="DEL"
    # if r==SVtypes.inversion:
        # string="INV"
    # if r==SVtypes.interspersedDuplication:
        # string="INTDUP"
    # if r==SVtypes.translocation:
        # string="TRX"
    # if r==SVtypes.novelInsertion:
        # string="INS"
    # if r==SVtypes.error:
        # string="Unknown"
    # if r==SVtypes.interchromosomal:
        # string="INTRX"
# #    for breakpoint in result[1]:
# #        bpstr+=delimeter+breakpoint
    # return string#+bpstr
    
def getResultType(result):
    string=""
    r=result[0]
    if r==SVtypes.tandem:
        string="DUP"
    if r==SVtypes.deletion:
        string="DEL"
    if r==SVtypes.inversion:
        string="INV"
    if r==SVtypes.interspersedDuplication:
        string="INTDUP"
    if r==SVtypes.translocation:
        string="TRX"
    if r==SVtypes.novelInsertion:
        string="INS"
    if r==SVtypes.error:
        string="Unknown"
    if r==SVtypes.interchromosomal:
        string="INTRX"
    return string

def printResultWithLine(result,content):
    string=""
    bpstr=""
    r=result[0]
    if r==-1:
        return
    if r==SVtypes.tandem:
        string="DUP"
    if r==SVtypes.deletion:
        string="DEL"
    if r==SVtypes.inversion:
        string="INV"
    if r==SVtypes.interspersedDuplication:
        string="INTDUP"
    if r==SVtypes.translocation:
        string="TRX"
    if r==SVtypes.novelInsertion:
        string="INS"
    if r==SVtypes.error:
        string="Unknown"
    if r==SVtypes.interchromosomal:
        string="INTRX"
    for breakpoint in result[1]:
        bpstr+=delimeter+breakpoint
    print string+bpstr
#     print string+":"+content[l].split("\t")[0].split(":")[1]

def detect (prevSV,prevResult,sv):
    if prevResult:
        prevResult=prevResult[0]
    result=-1
    #anchor_chrom,c1_anchor,c1_anchor_dir,realign_chrom,c1_realign,c1_realign_dir = sv
    anchor_chrom,c1_anchor,c1_anchor_dir = sv['chr1'],sv['pos1'],sv['dir1']
    realign_chrom,c1_realign,c1_realign_dir = sv['chr2'],sv['pos2'],sv['dir2']

    if realign_chrom == anchor_chrom:
        if c1_realign < c1_anchor:
            anchor_chrom,c1_anchor,c1_anchor_dir = sv['chr2'],sv['pos2'],sv['dir2']
            realign_chrom,c1_realign,c1_realign_dir = sv['chr1'],sv['pos1'],sv['dir1']

#    line=line.strip().split(delimeter) 
#    c1_realign=line[0].split(':')[1]
#    c1_realign_dir=line[1]
#    c1_realign_consensus=line[2]
#    c1_anchor=line[3].split(':')[1]
#    c1_anchor_dir=line[4]
#    c1_anchor_consensus=line[5]
#    c1_long_support=line[6]
#    c1_long_support_bases=line[7]
#    c1_short_support=line[8]
#    c1_short_support_bases=line[9]
#    c1_short_support_max_len=line[10]
#    c1_avg_realign_mapq=line[11]
#    c2_realign=line[12].split(':')[1]
#    c2_realign_dir=line[13]
#    c2_realign_consensus=line[14]
#    c2_anchor=line[15].split(':')[1]
#    c2_anchor_dir=line[16]
#    c2_anchor_consensus=line[17]
#    c2_long_support=line[18]
#    c2_long_support_bases=line[19]
#    c2_short_support=line[20]
#    c2_short_support_bases=line[21]
#    c2_short_support_max_len=line[22]
#    c2_avg_realign_mapq=line[23]
#    otherInfo=line[24]
    
    if realign_chrom != anchor_chrom:
        result=SVtypes.interchromosomal
#    elif c1_realign_dir!=c1_anchor_dir and c2_realign_dir!=c2_anchor_dir:
    elif c1_realign_dir!=c1_anchor_dir:
#         not inversion
        
        #if +anchor_loci>-realign_loci and -anchor_loci<+realign_loci, is tandem
        #if +anchor_loci<-realign_loci and -anchor_loci>+realign_loci, is deletion
        if int(c1_anchor)<int(c1_realign) and c1_anchor_dir=="-" and c1_realign_dir=="+":
            result=SVtypes.tandem
        elif int(c1_anchor)<int(c1_realign) and c1_anchor_dir=="+" and c1_realign_dir=="-":
            result=SVtypes.deletion
#        if (int(c1_anchor)>int(c1_realign) and c1_anchor_dir=="+") and (int(c2_anchor)<int(c2_realign) and c2_anchor_dir=="-"):
#            result=SVtypes.tandem
#        elif (int(c2_anchor)>int(c2_realign) and c2_anchor_dir=="+") and (int(c1_anchor)<int(c1_realign) and c1_anchor_dir=="-"):
#            result=SVtypes.tandem
#        elif (int(c2_anchor)>int(c2_realign) and c2_anchor_dir=="-") and (int(c1_anchor)<int(c1_realign) and c1_anchor_dir=="+"):
#            result=SVtypes.deletion
#        elif (int(c1_anchor)>int(c1_realign) and c1_anchor_dir=="-") and (int(c2_anchor)<int(c2_realign) and c2_anchor_dir=="+"):
#            result=SVtypes.deletion
        elif abs(int(c1_anchor)-int(c1_realign))<blurbp3:
            result=SVtypes.novelInsertion
        else:
            result=SVtypes.error
        
        #if two lines, interspersed duplication = del + tandem
        if prevSV:
            panchor_chrom,pc1_anchor,pc1_anchor_dir = prevSV['chr1'],prevSV['pos1'],prevSV['dir1']
            prealign_chrom,pc1_realign,pc1_realign_dir = prevSV['chr2'],prevSV['pos2'],prevSV['dir2']
            if prealign_chrom == anchor_chrom:
                if pc1_realign < pc1_anchor:
                    anchor_chrom,pc1_anchor,pc1_anchor_dir = prevSV['chr2'],prevSV['pos2'],prevSV['dir2']
                    prealign_chrom,pc1_realign,pc1_realign_dir = prevSV['chr1'],prevSV['pos1'],prevSV['dir1']
            distance1 = pc1_realign-int(c1_realign)
            distance2 = pc1_anchor-int(c1_anchor)
            if abs(distance1)<blurbp2 or abs(distance2)<blurbp2:#if this two events are nearby.
                if result+prevResult==SVtypes.interspersedDuplication:
                    result=SVtypes.interspersedDuplication
        
#    elif c1_realign_dir==c1_anchor_dir and c2_realign_dir==c2_anchor_dir:
    elif c1_realign_dir==c1_anchor_dir and c1_anchor_dir==c1_realign_dir:
#         inversion
        result=SVtypes.inversion
    else:
        result=SVtypes.error

#    startPoint=min(c1_realign,c1_anchor,c2_realign,c2_anchor)
#    endPoint=max(c1_realign,c1_anchor,c2_realign,c2_anchor)
    startPoint=min(c1_realign,c1_anchor,c1_anchor,c1_realign)
    endPoint=max(c1_realign,c1_anchor,c1_anchor,c1_realign)
#    if "nserted sequence" in otherInfo and abs(int(startPoint)-int(endPoint))<5:
#        result=SVtypes.novelInsertion
    return [result,[startPoint,endPoint]]

def detectTransloc(idx,sv_info,tolerance):
    #find the mobile part. r:[interInserType,[p1,p2,p3,p4]]
    if idx-1 < 0:
        return []    
    sv1 = sv_info[idx-1]
    sv2 = sv_info[idx]
    p1 = int(sv1['pos1'])
    p2 = int(sv2['pos1'])
    p3 = int(sv1['pos2'])
    p4 = int(sv2['pos2'])
    mobilePart = [p1,p2] if abs(p2-p1) > tolerance else [p3,p4]
    # flip coords if loc1 > loc2
    mobilePart = mobilePart if mobilePart[0]<mobilePart[1] else [mobilePart[1],mobilePart[0]]
    # try to find if the mobile part is deleted. If so, it is translocation.
    translocs = []
    for i,sv in enumerate(sv_info):
        if sv['classification']==getResultType([SVtypes.deletion]):
            p1 = int(sv['pos1'])
            p2 = int(sv['pos2'])
            if abs(mobilePart[0]-p1)<tolerance and abs(mobilePart[1]-p2)<tolerance:
                translocs = [idx-1,idx,i]
                break
    return translocs
    
def realignLoci(line):
    return int(line.strip().split(delimeter)[0].split(":")[1])
def anchorLoci(line):
    return int(line.strip().split(delimeter)[3].split(":")[1])
def wrapDIV(str):
    return '<p>'+str+'</p>\n'

def wrapColor(str,kind):
    htmlRightWrapper='<div class="alert alert-success">'
    htmlWrongWrapper='<div class="alert alert-danger">'
    htmlOtherWrapper='<div class="alert alert-warning">'
    htmlOther2Wrapper='<div class="alert alert-info">'
    return str
    if kind=="right":
        return htmlRightWrapper+str+htmlEndDiv
    if kind=="wrong":
        return htmlWrongWrapper+str+htmlEndDiv
    if kind=="other":
        return htmlOtherWrapper+str+htmlEndDiv
    if kind=="other2":
        return htmlOther2Wrapper+str+htmlEndDiv
    
def writeComapreResultToHTML(myResult,stand,compareResult,recal,precision):
    myResult=myResult.strip().split("\n")
    stand=stand.strip().split("\n")
    checkroll=[]
    for x in compareResult:
        for y in x:
            checkroll.append(y)
    
    htmlHead=open("head.html").read();
    htmlFoot=open("foot.html").read();
    htmlBody=""
    if not recal==0 and not precision==0:
        htmlBody+="Recal:"+str(recal)+"; Precision:"+str(precision)
    htmlBody+="<tr><td>Detect result</td><td>Source in breakpoint file</td><td>line No. in breakpoint file</td></tr>";

    
    i=0;
    unusedbk=0;
    lastStandPoint=0;
    for cr in compareResult:
        htmlCol1="";
        htmlCol2="";
        htmlclass="";
        if cr==[]: #my result has extra line
            htmlCol1=wrapColor(wrapDIV(myResult[i]), "other2")
            htmlCol2=wrapColor(wrapDIV("None source in breakpoint file"), "other2")
            htmlclass="other2"
            htmlBody+="<tr class='"+htmlclass+"'><td>"+htmlCol1+"</td>"+"<td>"+htmlCol2+"</td><td>"+"-</td></tr>"
        elif cr[0]==-2:
            htmlCol2=wrapColor(wrapDIV(stand[lastStandPoint+1]), "other");
            htmlCol1=wrapColor(wrapDIV("Deletion is part of translocation"), "other");
            htmlclass="other"
            htmlBody+="<tr class='"+htmlclass+"'><td>"+htmlCol1+"</td><td>"+htmlCol2+"</td><td>"+str(lastStandPoint+1)+"</td></tr>"
            checkroll.append(lastStandPoint+1)
        else:
            jump=cr[0]-lastStandPoint
#             print jump
## This is not very good, may miss breakpoints in standard.
            if jump>1: #standard has extra line, which may be  missed by Socrates.
                for x in range(1,jump):
                    if lastStandPoint+x not in checkroll and (lastStandPoint+x) < len(stand):
                        checkroll.append(lastStandPoint+x)
                        htmlCol2=wrapColor(wrapDIV(stand[lastStandPoint+x]), "other");
                        htmlCol1=wrapColor(wrapDIV("Cannot be detected"), "other");
                        htmlclass="other"
                        htmlBody+="<tr class='"+htmlclass+"'><td>"+htmlCol1+"</td><td>"+htmlCol2+"</td><td>"+str(lastStandPoint+x)+"</td></tr>"
                        unusedbk+=1
#                         print "unused dk:",lastStandPoint+x
                lastStandPoint+=jump
                
            if len(cr)>0:
                col2Str=""
                col1Str=""
                col1Str=wrapDIV(myResult[i])
                for x in cr:
                    col2Str+=wrapDIV(stand[x])
                htmlCol1=wrapColor(col1Str, "right");
                htmlCol2=wrapColor(col2Str, "right");
                lastStandPoint=min(cr)
                htmlclass="right"
                htmlBody+="<tr class='"+htmlclass+"'><td>"+htmlCol1+"</td><td>"+htmlCol2+"</td><td>"+str(cr)+"</td></tr>"
        i+=1
    html=htmlHead+htmlBody+htmlFoot
    open("index.html",'w').write(html)
#     print "unused bk amount:", unusedbk
