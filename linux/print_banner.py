#! /usr/bin/env python


def print_banner_line(banner, width=80):
    """
    Print a banner like '# --------------- FOO ------------------' to stdout.

    The number 4 in this piece of code counts for the two spaces wrapping the
    central text and an '# ' ahead of the banner.

    :param banner: the central text in the banner
    :param width: total width of the banner
    :return: None
    """
    num_marks_total = width - len(banner) - 4
    num_marks_left = num_marks_total // 2 
    num_marks_right = num_marks_total - num_marks_left
    banner_with_marks = "# "
    mark = "-"
    for i in range(num_marks_left):
        banner_with_marks += mark
    banner_with_marks += " %s " % banner
    for i in range(num_marks_right):
        banner_with_marks += mark
    print(banner_with_marks)


def print_banner_block(banner, width=80):
    """
    Print a banner like 
    #-----------------------------------
    #               FOO 
    #-----------------------------------
    to stdout.


    :param banner: the central text in the banner
    :param width: total width of the banner
    :return: None
    """
    num_marks_total = width - len(banner) - 1
    num_marks_left = num_marks_total // 2 
    banner_with_marks = "#"
    mark = " "
    for i in range(num_marks_left):
        banner_with_marks += mark
    banner_with_marks += banner
    banner_border = "#"
    for i in range(width-1):
        banner_border += "-"
    print(banner_border)
    print(banner_with_marks)
    print(banner_border)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("banner", type=str,  action="store")
    parser.add_argument("--width", "-w",   type=int, action="store", default="80")
    parser.add_argument("--type", "-t", type=str, action="store", default="line")
    args = parser.parse_args()
    if args.type == "line":
        print_banner_line(args.banner, args.width)
    elif args.type == "block":
        print_banner_block(args.banner, args.width)
    else:
        raise ValueError("ERROR: Unknown banner type '%s'" % args.flavor)
