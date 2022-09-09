def blockGenerator(msp_file_path=""):
    f = open(msp_file_path, "r", encoding="utf-8")
    block = []

    def block_check(block):
        if "NAME" in block[0]:
            return True, None
        else:
            return False, ["NAME: None\n"] + block

    while True:
        line = f.readline()
        if line:
            if line != " \n":
                block.append(line)
                continue
            else:
                if len(block) > 1:
                    x, y = block_check(block)
                    if x:
                        yield block
                        block = []
                        continue
                    else:
                        yield y
                        block = []
                else:
                    continue

        else:
            break