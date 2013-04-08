#!/usr/bin/env python
# Find and replace likely copyright headers in .cpp/.h files.
# Customise the directories to traverse below.

import os

#excludedir = ["..\\Lib"]

oldheaders=set()
nonheaders=set()

def update_source(filename, copyright):
    global oldheaders
    global nonheaders
    
    utfstr = chr(0xef)+chr(0xbb)+chr(0xbf)
    fdata = file(filename,"r+").read()
    isUTF = False
    if (fdata.startswith(utfstr)):
        isUTF = True
        fdata = fdata[3:]
    
    i=0
    line_comment=False
    block_comment=False
    prev_slash=False # or prev_star
    for c in fdata:
        i+=1
        if block_comment:
            # use prev_slash as prev_star instead
            if prev_slash and c=='/':
                block_comment=False
                prev_slash=False
                break # good behaviour to stop at end of first block comment?
                # (or better to break after an empty line?)
            elif c=='*':
                prev_slash=True
            else:
                prev_slash=False
        elif c in [' ','\t','\f','\v']:
            pass
        elif c in ['\r','\n']:
            line_comment=False
        elif c=='/':
            if prev_slash:
                line_comment=True
                prev_slash=False
            else:
                prev_slash=True
        elif c=='*' and prev_slash:
            prev_slash=False
            block_comment=True
        else:
            i -= 1  # go back (keep this character)
            break   # end of header
    header=fdata[0:i]
    
    if header in oldheaders:
        fdata = fdata[len(header):]
    elif header in nonheaders:
        pass
    else:
        print header
        r=raw_input("Remove above header? (y/N): ")
        if r[0]=='y' or r[0]=='Y':
            oldheaders.add(header)
            fdata = fdata[len(header):]
        else:
            nonheaders.add(header)
    
    if not (fdata.startswith(copyright)):
        print "updating "+filename
        fdata = copyright + fdata
        if (isUTF):
            file(filename,"w").write(utfstr+fdata)
        else:
            file(filename,"w").write(fdata)

def recursive_traversal(dir, copyright):
    global excludedir
    fns = os.listdir(dir)
    #print "listing "+dir
    for fn in fns:
        fullfn = os.path.join(dir,fn)
        #if (fullfn in excludedir):
            #continue
        if (os.path.isdir(fullfn)):
            recursive_traversal(fullfn, copyright)
        else:
            if (fullfn.endswith(".h") or fullfn.endswith(".cpp")):
                update_source(fullfn, copyright)


cright = file("licence-template.txt","r+").read()
recursive_traversal("include", cright)
recursive_traversal("model", cright)
exit()
