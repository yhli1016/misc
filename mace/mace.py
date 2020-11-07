import re


def include(filename, nl0=None, nl1=None):
    with open(filename, "r") as infile:
        content = infile.readlines()
    nl_start = nl0 if nl0 is not None else 1
    nl_end   = nl1 if nl1 is not None else len(content)
    longline = "".join(content[(nl_start-1):nl_end])
    return longline


def expand_line(macro, line):
    flag = False
    match_result = re.search(r"<[a-zA-Z0-9_]+>", line)
    if match_result is not None:
        mac_name = match_result.group()[1:-1]
        if mac_name in macro.keys():
            mac_text = str(macro[mac_name]).lstrip("\n").rstrip("\n")
            line = re.sub(r"<%s>" % mac_name, mac_text, line)
            flag = True
    return line, flag


def main(macro, input, output):
    with open(input, "r") as in_file:
        content_raw = in_file.readlines()
    content_expanded = []
    for line in content_raw:
        while True:
            line, flag = expand_line(macro, line)
            if flag == False:
                break
        content_expanded.append(line)
    with open(output, "w") as out_file:
        for line in content_expanded:
            out_file.write(line)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,  action="store", required=True)
    parser.add_argument("-o", "--output", type=str, action="store", required=True)
    parser.add_argument("-m", "--macro", type=str, action="store", nargs="*")
    args = parser.parse_args()
    macro = dict()
    for argv in args.macro:
        s = argv.split("=")
        if len(s) == 1:
            macro[s[0]] = ""
        else:
            macro[s[0]] = s[1]
    main(macro, args.input, args.output)
