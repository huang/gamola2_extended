#!/usr/bin/python

import sys

lines = open(sys.argv[1])
for line in lines:
    if line.startswith(">"):
        #newline = line.replace(".", "_")
        newline = ">%s" % (line.split(" ")[1])
        print newline.rstrip()
    else:
        print line.rstrip()
