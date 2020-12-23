#! /usr/bin/env python


def print_banner_line(banner=None, width=80, mark="-"):
    """
    Print a banner like '#--------------- FOO ---------------' to stdout.

    The magic number 3 accounts for a '#' ahead of the banner and two spaces
    wrapping the banner text.

    :param banner: string, central text in the banner
    :param width: integer, total width of the banner
    :param mark: character, marker
    :return: None
    """
    num_marks_total = width - len(banner) - 3
    num_marks_left = num_marks_total // 2
    num_marks_right = num_marks_total - num_marks_left
    banner_with_marks = "#" + "".join([mark for _ in range(num_marks_left)])
    banner_with_marks += " %s " % banner
    banner_with_marks += "".join([mark for _ in range(num_marks_right)])
    print(banner_with_marks)


def print_banner_block(banner, width=80, mark="-"):
    """
    Print a banner like
    #----------------------------------#
    #               FOO                #
    #----------------------------------#
    to stdout.

    :param banner: string, central text in the banner
    :param width: integer, total width of the banner
    :param mark: character, marker
    :return: None
    """
    num_spaces_total = width - len(banner) - 2
    num_spaces_left = num_spaces_total // 2
    num_spaces_right = num_spaces_total - num_spaces_left
    banner_with_spaces = "#" + "".join([" " for _ in range(num_spaces_left)])
    banner_with_spaces += banner
    banner_with_spaces += "".join([" " for _ in range(num_spaces_right)]) + "#"
    border = "#" + "".join([mark for _ in range(width-2)]) + "#"
    print(border)
    print(banner_with_spaces)
    print(border)


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
