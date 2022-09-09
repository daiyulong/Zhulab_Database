#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# @Time   : 2022/9/7
# @Author : Dai Yulong

# Add inchikey and smiles information to msp files.
import os
import sys
from blockGenerator import blockGenerator

def parse_note_to_dict(file):
    f = open(file,"r",encoding="utf-8").readlines()[:]
    header = f[0]
    data_dict = {}
    for line in f[1:]:
        line_list = line.strip().split("\t")
        name = line_list[0]
        inchikey = line_list[1]
        smiles = line_list[2]
        if name not in data_dict :
            data_dict[name] = [inchikey,smiles]
        else:
            print("1")
            sys.exit(1)
    return data_dict

if __name__ == "__main__":
    path = r"E:\projectResearch\03Zhulab\output\CombineInchi.Smile-neg.txt"
    data = parse_note_to_dict(path)
    print(data)
    msp_path = r"E:\projectResearch\03Zhulab\output\Zhulab-neg_all_meta.msp"
    output = r"E:\projectResearch\03Zhulab\output\result-neg.msp"
    w = open(output,"w",encoding="utf-8")
    for block in blockGenerator(msp_path):
        name = block[0].strip().split("NAME: ")[1]
        if name in data:
            inchikey = data[name][0]
            smiles = data[name][1]
            block.insert(4,inchikey + "\n")
            block.insert(5,smiles + "\n")
            w.write("".join(block) +"\n")
        else:
            print(name)


