from blockGenerator import blockGenerator

if __name__ == "__main__":
    path = r"E:\projectResearch\03Zhulab\output\Zhulab-pos_all_meta.msp"

    for block in blockGenerator(path):
        print(block[0].strip().split("NAME: ")[1])


