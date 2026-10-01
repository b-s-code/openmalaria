#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This file is part of OpenMalaria.

Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
Copyright (C) 2020-2026 University of Basel
Copyright (C) 2025-2026 The Kids Research Institute Australia

OpenMalaria is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.
 
This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
"""

import os
import re

#excludedir = ["..\\Lib"]

oldheaders=set()
nonheaders=set()

def update_source(filename, copyright):
    global oldheaders
    global nonheaders
    
    utfstr = chr(0xef)+chr(0xbb)+chr(0xbf)
    fdata = open(filename,"r+").read()
    isUTF = False
    if (fdata.startswith(utfstr)):
        isUTF = True
        fdata = fdata[3:]
    
    i=0
    line_comment=False
    block_comment=False
    prev_slash=False # or prev_star
    want_block=True # false once we have a block comment or set of line comments
    for c in fdata:
        i+=1
        if block_comment:
            # use prev_slash as prev_star instead
            if prev_slash and c=='/':
                block_comment=False
                prev_slash=False
                want_block=False
                # or: break to stop now (don't eat any new-lines)
            elif c=='*':
                prev_slash=True
            else:
                prev_slash=False
        elif c in [' ','\t','\f','\v']:
            pass
        elif c in ['\r','\n']:
            line_comment=False # end of line
        elif want_block and c=='/':
            if prev_slash:
                line_comment=True
                prev_slash=False
            else:
                prev_slash=True
        elif want_block and c=='*' and prev_slash:
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
        print(header)
        r=input("Remove above header? (y/N): ")
        if r[0]=='y' or r[0]=='Y':
            oldheaders.add(header)
            fdata = fdata[len(header):]
        else:
            nonheaders.add(header)
    
    if not (fdata.startswith(copyright)):
        print(("updating "+filename))
        fdata = copyright + fdata
        if (isUTF):
            open(filename,"w").write(utfstr+fdata)
        else:
            open(filename,"w").write(fdata)

def recursive_traversal(dir, suffixes, update, copyright):
    """Call update(file, copyright) on every file under dir whose name ends with one of suffixes."""
    fns = os.listdir(dir)
    #print "listing "+dir
    for fn in fns:
        fullfn = os.path.join(dir,fn)
        if (os.path.isdir(fullfn)):
            if not (fn.startswith('.') or fn.startswith('build')):
                # Skip hidden dirs, build dirs.
                recursive_traversal(fullfn, suffixes, update, copyright)
        elif fn.endswith(suffixes):
            update(fullfn, copyright)

def line_comment(c_header, prefix):
    """Convert the licence template, which uses c-style comments, to line comments starting with prefix (e.g. '#' or 'REM')."""
    lines = []
    for line in c_header.splitlines():
        if line.startswith(' */'):
            break
        elif line.startswith('/*'):
            line = line[2:]
        elif line.startswith(' *'):
            line = line[2:]
        lines.append(prefix + line)
    return '\n'.join(lines) + '\n'

def update_line_comment_header(filename, copyright):
    """Replaces the copyright header in a cmake/python/bash/bat file i.e. a file which uses line comments, not c-style comments."""
    # files without the full licence block (e.g. old short notices) are left alone
    fdata = open(filename, newline='').read()  # newline='' keeps CRLF line endings as-is
    if '\r\n' in fdata:
        copyright = copyright.replace('\n', '\r\n')
    full_block = re.compile(r'^(#|REM) This file is part of OpenMalaria\.\r?\n.*?Boston, MA 02110-1301, USA\.\r?\n', re.DOTALL | re.MULTILINE)
    new = full_block.sub(lambda m: copyright, fdata, count=1)
    if new != fdata:
        print("updating "+filename)
        open(filename, "w", newline='').write(new)

cright = open("util/licence-template.txt","r+").read()
recursive_traversal("model", (".h", ".cpp"), update_source, cright)
recursive_traversal("unittest", (".h", ".cpp"), update_source, cright)
recursive_traversal(".", ("CMakeLists.txt", ".sh"), update_line_comment_header, line_comment(cright, "#"))
recursive_traversal(".", (".bat",), update_line_comment_header, line_comment(cright, "REM"))

print('Remember to update the text for --version in model/CommandLine.cpp!')

exit()
