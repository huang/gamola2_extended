#!/usr/bin/python

import sys

lines = open(sys.argv[1])
startp = 0
endp = 0
id_ = 1
for line in lines:
    #newline = line.replace(".", "_")
    tokens = line.split("\t")
    #name = tokens[0]
    startp = endp+1
    endp = startp + int(tokens[1])-4
    newline = "    %i  %i  %i  [+1" % (id_, startp, endp)
    print newline.rstrip()
    endp += 103
    id_ += 1
